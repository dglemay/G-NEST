--#########################################################################
-- Copyright (C) 2009-2012 William F. Martin and Danielle G. Lemay
--
-- This program is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the
-- Free Software Foundation;
--
-- This program is distributed in the hope that it will be useful, but
-- WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
-- See the GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License along
-- with this program; if not, write to the Free Software Foundation, Inc.,
-- 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
--#########################################################################
SET dynamic_library_path TO '/local/dglemay:$libdir';

---------------------------------------------------------------------------
--  Create a new schema with the same name as the project.
--  If the schema already exists, a number suffix is added to make it unique.
--  Creation of the schema without first testing for the availability of 
--  the name is done to avoid a race condition.
--
--  RETURN the name of the schema actually created.
---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION create_schema(
    p_project_name           varchar)
    RETURNS varchar AS $$
DECLARE
  v_i                        integer;
  v_schema                   varchar;
BEGIN
  v_i := 0;

  WHILE (v_schema IS NULL)  LOOP
    v_schema := p_project_name ||
      CASE WHEN v_i = 0 THEN '' ELSE '_' || v_i::varchar END;

    BEGIN
      EXECUTE 'CREATE SCHEMA ' || v_schema;
    EXCEPTION
      WHEN  duplicate_schema  THEN  v_schema := NULL;
    END;

    v_i := v_i + 1;
  END LOOP;

  EXECUTE 'SET search_path TO ' || v_schema || ',public';

  RETURN v_schema;

END;
$$ LANGUAGE plpgsql;


-------------------------------------------------------------------------
-- Create a directory owned by the 'postgres' user but readable by anyone.
-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION create_tmp_dir() RETURNS varchar AS $$
  use File::Temp('tempdir');
  my $dir = tempdir();
  my $cmd = "chmod a+rx $dir";
  system($cmd)==0  or  elog(ERROR, "Can't execute command: $cmd");
  $dir;
$$ LANGUAGE plperlu;


-------------------------------------------------------------------------
-- Create a directory owned by the 'postgres' user but readable by anyone.
-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION create_readable_subdir(
    p_parent_path            varchar,
    p_subdir_path            varchar)
    RETURNS void AS $$

  my $cmd = "cd $_[0]; mkdir -m 0755 $_[1]";
  system($cmd)==0  or  elog(ERROR, "Can't execute command: $cmd");

$$ LANGUAGE plperlu;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION purge_tmp_dir(varchar) RETURNS varchar AS $$
  my $cmd = "rm -fr $_[0]";
  system($cmd)==0  or  elog(ERROR, "Can't execute command: $cmd");
$$ LANGUAGE plperlu;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION blk_stats_query( 
    p_blk_size               integer, 
    p_max_win_size           integer, 
    p_mode                   varchar,  -- 'TNS', 'ANC', or 'pvalue' 
    p_chromosome             varchar,  -- null for all 
    p_chrom_length           integer)
    RETURNS varchar AS $$

  my ($blk_size, $max_win_size, $mode, $chromosome, $chrom_length)
      = @_;
  $mode = lc($mode);

  my $fnx = $mode eq 'pvalue' ? 'min' : 'max';

  my $query = <<QUERY_1;
WITH
  pass_1 AS (
    SELECT
      nh_parm AS win_size,
      (chr_start_pos / $blk_size) AS blk_index,
      $mode AS value
    FROM nh_stats JOIN gene_nh USING(nh_id)
    WHERE (NOT by_gene_counts)  AND chromosome = '$chromosome'
    ),
  pass_2 AS (
    SELECT
      win_size,
      blk_index,
      $fnx(value) AS value
    FROM pass_1
    GROUP BY win_size, blk_index
    )
SELECT
  win_size,
  blk_index + add_blk_index  AS blk_index,
  $fnx(value) AS value
FROM
  pass_2,
  generate_series(0, $max_win_size/$blk_size) AS add_blk_index
WHERE add_blk_index*$blk_size < win_size  AND
      (blk_index + add_blk_index) * $blk_size < $chrom_length
GROUP BY win_size, (blk_index + add_blk_index);

QUERY_1

  $query;

$$ LANGUAGE plperl;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gnest_fill_nh_by_win_size(
    p_min_win_size           integer,
    p_max_win_size           integer,
    p_win_size_increment     integer)
    RETURNS void AS 

  'gnest_fill_nh_by_win_size.so'

LANGUAGE C;
    

-------------------------------------------------------------------------
--CREATE OR REPLACE FUNCTION gnest_dot_product(  
--    p_gene_1_expr_ranks      integer[],
--    p_gene_2_expr_ranks      integer[])
--    RETURNS integer AS $$
--
--  SELECT sum( $1[i] * $2[i] )::integer
--  FROM generate_series(1,array_length($1,1)) AS i;
--
--$$ LANGUAGE SQL;

-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gnest_dot_product(  
    p_gene_1_expr_ranks      integer[],
    p_gene_2_expr_ranks      integer[])
    RETURNS integer AS 

  'gnest_dot_product.so'

LANGUAGE C;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gen_nh_by_gc(
    p_min_gene_count         integer,
    p_max_gene_count         integer)
    RETURNS SETOF gb_nh_info AS $$
DECLARE
  v_chrom                    RECORD;  -- chromosome_info
  v_out                      gb_nh_info;
BEGIN
  FOR v_chrom IN (SELECT * FROM chromosome_info) LOOP
    v_out.chromosome := v_chrom.chromosome;
    FOR v_start_gene_pos_index IN 1..(v_chrom.num_genes-p_min_gene_count+1)
    LOOP
      v_out.start_gene_pos_index := v_start_gene_pos_index;
      FOR v_num_genes IN p_min_gene_count ..
          LEAST(p_max_gene_count, v_chrom.num_genes-v_start_gene_pos_index)
      LOOP
        v_out.num_genes := v_num_genes;
        RETURN NEXT v_out;
      END LOOP;
    END LOOP;
  END LOOP;
END;
$$ LANGUAGE plpgsql;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_nh_stats_grid(
    p_blk_size               integer,
    p_mode                   varchar,  -- 'TNS', 'ANC', or 'pvalue'
    p_chromosome             varchar)
    RETURNS SETOF nh_stats_grid_cell AS $$
DECLARE
  v_max_win_size             integer;
  v_chrom_length             integer;
BEGIN
  SELECT max(nh_parm) INTO v_max_win_size
  FROM nh_stats WHERE NOT by_gene_counts;

  SELECT chrom_length INTO v_chrom_length
  FROM chromosome_info WHERE chromosome = p_chromosome;

  EXECUTE 'CREATE TEMP TABLE blk_stats AS ' || blk_stats_query(p_blk_size,
      v_max_win_size, p_mode, p_chromosome, v_chrom_length);
  CREATE INDEX blk_stat_key ON blk_stats(win_size, blk_index);

  ---------------------------------------------------------------------------
  --  This query is takes the sparse matrix of values to graph and fills in
  --  The holes with null values.
  ---------------------------------------------------------------------------
  RETURN QUERY
    WITH
      win AS (
        SELECT DISTINCT nh_parm AS win_size
        FROM nh_stats
        WHERE NOT by_gene_counts
        ORDER BY win_size
        ),
      blks(blk_index) AS (
        SELECT generate_series(1, (SELECT max(blk_index) FROM blk_stats))
        )
    SELECT win_size, blk_index, value
    FROM blks CROSS JOIN win
         LEFT OUTER JOIN blk_stats USING(win_size, blk_index)
    ORDER BY win_size, blk_index;

  DROP TABLE blk_stats;

END;
$$ LANGUAGE plpgsql;


--------------------------------------------------------------------------
--  Given a cluster of genes that overlap each other (including transitively),
--  Return the highest expressed non-overlapping gene(s).
--
--  If A overlaps B overlaps C (but A doesn't overlap C), then if A has the
--  highest expression, both A and C should be returned.
--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION select_genes_from_cluster(
    p_genes                  varchar[])
    RETURNS SETOF varchar
    LANGUAGE plpgsql AS $$
DECLARE
  v_gene_name                varchar;
  v_start_pos                integer;
  v_end_pos                  integer;
BEGIN
  TRUNCATE tmp_gene_clust;

  INSERT INTO tmp_gene_clust(gene_name, start_pos, end_pos, max_expr)
    WITH tmp_ids(gene_name) AS (SELECT unnest(p_genes))
    SELECT
      gene_name,
      start_pos,
      end_pos,
      max(expr) AS max_expr
    FROM tmp_ids
    JOIN raw_gene_info g USING(gene_name)
    LEFT OUTER JOIN raw_expr_data e USING(gene_name)
    GROUP BY gene_name, start_pos, end_pos
    ORDER BY max_expr DESC;

  LOOP
    SELECT gene_name, start_pos, end_pos
    INTO v_gene_name, v_start_pos, v_end_pos
    FROM tmp_gene_clust
    ORDER BY max_expr DESC LIMIT 1;

    EXIT WHEN NOT FOUND;

    RETURN NEXT v_gene_name;

    DELETE FROM tmp_gene_clust
    WHERE GREATEST(start_pos, v_start_pos) <= LEAST(v_end_pos, end_pos);
  END LOOP;

END;
$$;


--------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_non_overlapping_genes()
    RETURNS SETOF varchar
    LANGUAGE plpgsql AS $$
DECLARE
  v_genes                    varchar[];
  v_clust_end_pos            integer;
  v_chromosome               varchar;
  v_rec                      RECORD;
BEGIN
  CREATE TEMP TABLE tmp_gene_clust(
    gene_name                varchar,
    start_pos                integer,
    end_pos                  integer,
    max_expr                 real
  );

  FOR v_chromosome IN SELECT chromosome FROM raw_chrom_lengths  LOOP
    v_genes := NULL;
    FOR v_rec IN
        SELECT * FROM raw_gene_info
          WHERE chromosome = v_chromosome
          ORDER BY start_pos, end_pos DESC
    LOOP
      ----------------------------------------------------------------
      --  Past the current cluster?  Handle it and start a new one.
      ----------------------------------------------------------------
      IF v_rec.start_pos > v_clust_end_pos THEN
        IF array_length(v_genes,1) = 1 THEN
          RETURN NEXT v_genes[1];
        ELSE
          RETURN QUERY SELECT select_genes_from_cluster(v_genes);
        END IF;

        v_genes := NULL;
        v_clust_end_pos := 0;
      END IF;

      v_genes := v_genes || v_rec.gene_name;
      v_clust_end_pos := GREATEST(v_clust_end_pos, v_rec.end_pos);
    END LOOP;  -- end gene

    IF array_length(v_genes,1) = 1 THEN
      RETURN NEXT v_genes[1];
    ELSE
      RETURN QUERY SELECT select_genes_from_cluster(v_genes);
    END IF;

  END LOOP; -- end chromosome

  DROP TABLE tmp_gene_clust;
END;
$$;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION export_db(
    p_schema                 varchar,
    p_include_corr_matrix    boolean,
    p_output_file            varchar)
    RETURNS void AS $$

  my ($schema, $include_corr_matrix, $output_file) = @_;
  
  my @tables = qw(chromosome_info samples gene_info expr_data
      gene_nh gene_to_nh nh_stats nh_synteny syntenic_blocks);

  push @tables, 'corr_matrix'  if ($include_corr_matrix);

  my $t_list = join(' ', map("-t $schema.$_", @tables) );

  my $sed_cmds = '/^CREATE TABLE/,/^);/ p; /^COPY /,/^\\\./ p;';

  my $cmd = "pg_dump $t_list gnest | sed -n '{$sed_cmds}' > $output_file";
  my $rv = system($cmd);
  $rv==0  or  elog(ERROR, "Error($rv): system('$cmd')");

  $cmd = "chmod a+r $output_file";
  my $rv = system($cmd);
  $rv==0  or  elog(ERROR, "Error($rv): system('$cmd')");

$$ LANGUAGE plperlu;


-------------------------------------------------------------------------
--  Format a query to return the statistic values for each non-empty
--  cell in the graph.
-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION raw_grid_query(
    p_chromosome             varchar,
    p_chrom_length           integer,
    p_chr_first_gene_id      integer,
    p_min_nh_parm            integer,
    p_by_gene_counts         boolean,
    p_x_axis_bp              boolean)
    RETURNS varchar AS $$

  my ($chromosome, $chrom_length, $chr_first_gene_id, 
      $min_nh_parm, $p_by_gene_counts, $p_x_axis_bp) = @_;

  my $by_gene_counts = (lc(substr($p_by_gene_counts,0,1)) eq 't');
  my $x_axis_bp = (lc(substr($p_x_axis_bp,0,1)) eq 't');


  #---------------------------------------------------------------------
  #  This is to compute the x-coordinates range of neighborhood.
  #
  #  If the x-axis is in base pairs (units displayed as MB):
  #    If the y-axis is sliding windows, then to determine the x-coordinate
  #    span of a sliding window, use the full width based upon the win_size
  #    centered on the center of the genes span.
  #    However, the values must be truncated at the length of the chromosome.
  #
  #    IF the y-axis is gene counts, then use the span of the genes in
  #    the neighborhood.
  #
  #  If the x-axis is gene indices, use the gene_id range, shifted so that
  #  the chromosome starts with zero.
  #---------------------------------------------------------------------
  my $x_value = $x_axis_bp ?
    ($by_gene_counts ?
      'generate_series( ' .
      '(lower_bound::real/100000.0)::integer * 100000,' .
      '(upper_bound::real/100000.0)::integer * 100000,' .
      ' 100000) '
      :

    'generate_series( ' .
    '((lower_bound+upper_bound-nh_parm)::real/200000.0)::integer * 100000,' .
    'LEAST(' .
    '((lower_bound+upper_bound+nh_parm)::real/200000.0)::integer * 100000,' .
    " $chrom_length - ($chrom_length % 100000)), " .
    ' 100000) '
    )

    :  # (X-axis as gene indices)
 'generate_series(start_gene_pos_index + 1, start_gene_pos_index + num_genes)'
    ;

  my $not_by_gc = $by_gene_counts ? '' : 'NOT';

  my $query = <<RAW_GRID_QUERY;
    WITH 
      pass_1 AS (
        SELECT
          $x_value AS x,
          nh_parm AS y,
          tns, anc, pvalue
        FROM nh_stats JOIN gene_nh USING(nh_id)
        WHERE ($not_by_gc by_gene_counts) AND chromosome = '$chromosome' AND
              anc > 0
        )
    SELECT
      x,
      y,
      min(pvalue) AS pvalue,
      max(tns)    AS tns,
      max(anc)    AS anc
    FROM pass_1
    GROUP BY x, y
RAW_GRID_QUERY

  $query;

$$ LANGUAGE plperl;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION build_grid(
    p_chromosome             varchar,
    p_chr_num_genes          integer,
    p_chr_first_gene_id      integer,
    p_chr_length             integer,
    p_min_nh_parm            integer,
    p_max_nh_parm            integer,
    p_by_gene_counts         boolean,
    p_x_axis_bp              boolean)
    RETURNS void AS $$
DECLARE
  v_min_x                    integer;
  v_max_x                    integer;
  v_step_x                   integer;
  v_min_y                    integer;
  v_max_y                    integer;
  v_step_y                   integer;
BEGIN

  EXECUTE 'CREATE TEMP TABLE tmp_raw_grid AS ' || raw_grid_query(
      p_chromosome, p_chr_length, p_chr_first_gene_id, 
      p_min_nh_parm, p_by_gene_counts, p_x_axis_bp);
  CREATE INDEX raw_grid_index ON tmp_raw_grid(x,y);

  ------------------------------------------------------------------------
  --  Make the grid just large enough to encompass all non-empty cells
  --  in the sparse grid.  This in effect truncates the portion of the
  --  chromosome beyond the right-most neighborhood (for base-pair graphs).
  ------------------------------------------------------------------------
  SELECT max(x), max(y)  INTO v_max_x, v_max_y  FROM tmp_raw_grid;

  IF p_x_axis_bp THEN
    v_min_x := 0;
    v_step_x := 100000;
  ELSE
    v_min_x := 1;
    v_step_x := 1;
  END IF;

  v_step_y := CASE WHEN p_by_gene_counts THEN 1 ELSE p_min_nh_parm END;

  ---------------------------------------------------------------------------
  --  This query takes the sparse matrix of values to graph and fills in
  --  the holes with null values.
  ---------------------------------------------------------------------------
  INSERT INTO tmp_grid_values(x, y, pvalue, tns, anc)
    WITH
      x_values AS (SELECT generate_series(v_min_x, v_max_x, v_step_x) AS x),
      y_values AS (
        SELECT generate_series(p_min_nh_parm, v_max_y, v_step_y) AS y )
    SELECT x, y, pvalue, tns, anc
    FROM x_values CROSS JOIN y_values
         LEFT OUTER JOIN tmp_raw_grid USING(x,y)
    ORDER BY x,y;

  DROP TABLE tmp_raw_grid;

END;
$$ LANGUAGE plpgsql;



---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION generate_legends(
    output_dir               varchar,
    graphics_file_type       varchar)
    RETURNS void AS $$

  for (stat_type in c('pvalue', 'anc_or_tns')) {

    fname = paste(output_dir, '/', stat_type, '_', 'legend', '.',
                  graphics_file_type, sep='')
    do.call(graphics_file_type, list(fname))
  
    my.colors <- if (stat_type == 'pvalue')
         { terrain.colors(24)[c(1,6:24)] } else { rev(heat.colors(12)) }
  
    ColorLevels <- seq(0,1, length=length(my.colors))
  
    layout(matrix(c(1,2)), widths=c(1,2), heights=c(3,1), respect=TRUE)
  
    # Color Scale
    image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),
          nrow=1), col=my.colors, xlab="",ylab="", xaxt="n")
    dev.off()

    cmd = sprintf('chmod a+r %s', fname);
    rv = system(cmd)
    if (rv != 0)  stop(sprintf("Can't execute command: %s", cmd))
  }

$$ LANGUAGE plr;


---------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION generate_graphs(
    main_title               varchar,
    chromosome               varchar,
    chrom_length             integer,
    use_tns                  boolean,
    by_gene_counts           boolean,
    x_axis_bp                boolean,
    output_dir               varchar,
    outfile_fmt_str          varchar,
    graphics_file_type       varchar)
    RETURNS void AS $$

  if (!file.exists(output_dir))   dir.create(output_dir, mode="0755")

  stats_types = if (use_tns) c('pvalue','anc','tns') else c('pvalue','anc')

  #-----------------------------------------------------------------------
  #  Read the grid values from the database.
  #-----------------------------------------------------------------------
  grid_df = pg.spi.exec('SELECT * FROM tmp_grid_values ORDER BY y,x')

# rs = dbSendQuery(con, 'SELECT * FROM tmp_grid_values ORDER BY y,x')
# grid_df = fetch(rs, n=-1)


  #-----------------------------------------------------------------------
  #  Set up x-axis details:
  #    Using base pairs:
  #     - cells are 100000 bp (0.1 MB) wide
  #     - ticks are every 10 MB
  #     - tick labels must be converted from bp to MB.
  #    Using gene indices:
  #     - cells are 1 gene wide
  #     - ticks are every 100 genes 
  #-----------------------------------------------------------------------
  x_values = sort(unique(grid_df$x))
  if (x_axis_bp) {            #  x-axis details by base pairs
    x_ticks  = 10000000 * (0:(max(x_values) %/% 10000000))
    x_tick_labels = x_ticks/1000000
    x_axis_label = sprintf("Location on chromosome %s (MB)", chromosome)
  } else {                    #  x-axis details by gene index
    x_ticks  = seq(0, max(x_values), 100)
    x_tick_labels = x_ticks
    x_axis_label = sprintf("Gene index on chromosome %s", chromosome)
  }

  y_values = sort(unique(grid_df$y))
  if (by_gene_counts) {
    y_ticks = sort(unique(grid_df$y))
    y_tick_labels = y_ticks
    y_axis_label = "Number genes in neighborhood"
    
  } else {
    y_ticks = seq(0, max(grid_df$y), 1000000)
    y_tick_labels = y_ticks / 1000000
    y_axis_label = "Window size (MB) of gene neighborhood"
  }

  grid_dim = c(length(x_values), length(y_values))

  for (stat_type in stats_types) {
    my.colors <- if (stat_type == 'pvalue')
       { terrain.colors(24)[c(1,6:24)] } else { rev(heat.colors(12)) }

    fname = paste(output_dir, '/', 
      sprintf(outfile_fmt_str, stat_type), '.', graphics_file_type, sep='')
    m = grid_df[,stat_type]
    dim(m) = grid_dim

    do.call(graphics_file_type, list(fname))
    image(x_values, y_values, z=m, xlab=x_axis_label, ylab=y_axis_label, 
          zlim=c(0,1), col=my.colors, xaxt="n", yaxt="n")
    axis(1, at=x_ticks, labels=x_tick_labels)
    axis(2, at=y_ticks, labels=y_tick_labels)
    title(main=sprintf("%s : chr %s\n(%s)", main_title, chromosome, stat_type))
    dev.off()
  }

  cmd = sprintf('chmod a+r %s/*', output_dir);
  rv = system(cmd)
  if (rv != 0)  stop(sprintf("Can't execute command: %s", cmd))

$$ LANGUAGE plr;
