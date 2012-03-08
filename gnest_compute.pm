##########################################################################
# Copyright (C) 2009-2012 William F. Martin and Danielle G. Lemay
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation;
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
package gnest_compute;
use strict;
use DBI;
use DBD::Pg;


#--------------------------------------------------------------------------
sub pg_version_supports_array_agg_order_by {
  my ($dbh) = @_;

  my $sth = $dbh->prepare(
    "SELECT regexp_replace(version(), 'PostgreSQL +([0-9]+).*', '\\1')");
  $sth->execute();
  my ($pg_version) = $sth->fetchrow_array();

  $pg_version >= 9;
}


#--------------------------------------------------------------------------
sub create_legacy_gen_expr_arrays {
  my ($dbh) = @_;
 
  my $query = <<'LEGACY_GEN_EXPR_ARRAYS';

CREATE OR REPLACE FUNCTION gen_expr_arrays()
    RETURNS SETOF expr_array AS $$
DECLARE
  v_num_samples              integer;
  v_out                      expr_array;
  v_gene_id                  integer;
  v_expr_rank                integer;
BEGIN
  SELECT count(*) INTO v_num_samples FROM samples;
  v_out.gene_id = 0;

  FOR v_gene_id, v_expr_rank IN (
      SELECT gene_id, expr_rank
      FROM expr_data AS e  JOIN gene_info USING(gene_id)
      WHERE not silent
      ORDER BY gene_id, sample_id)
  LOOP
    IF v_gene_id <> v_out.gene_id THEN
      IF v_out.gene_id <> 0 THEN
        RETURN NEXT v_out;
      END IF;
      v_out.gene_id = v_gene_id;
      v_out.expr_ranks_array = array[]::integer[];
    END IF;

    v_out.expr_ranks_array := array_append(v_out.expr_ranks_array,
      2 * v_expr_rank - v_num_samples - 1);
  END LOOP;

  RETURN NEXT v_out;
END;
$$ LANGUAGE plpgsql;

LEGACY_GEN_EXPR_ARRAYS

  $dbh->do($query)  or  die "Error: Creating function gen_expr_arrays()\n";
}


#--------------------------------------------------------------------------
sub correlation {
  my ($dbh) = @_;

  my $query;

  #-----------------------------------------------------------------------
  #  This function computes the Spearman correlation in expression values
  #  (across samples) between pairs of genes.
  #
  #  To facilitate computation of pairwise correlations, the expression ranks
  #  are put into arrays (so that the computation only involves combinations
  #  of two rows).
  #
  #  Only fill in half the matrix (where gene_1 < gene_2).
  #
  #  Spearman's correlation is just Pearson's correlation computed using ranks.
  #  In this particular problem, several factors of the formula for computing
  #  Pearson's correlation are constant for all gene pairs. Those factors are:
  #    
  #    * standard deviation of ranks for both genes
  #    * N (number of samples)
  #    * mean of ranks
  #
  #  The remaining variable factors are (rank - mean_rank) for each of the
  #  genes.  This term is multiplied by 2 so then it will always be an integer.
  #  The variable terms end up being:
  #       (2 * rank - num_samples - 1)  
  #
  #  (thus keeping the computation in whole integers up until the last step,
  #   when dividing by the 'divisor').
  #  The divisor  =  stddev^2 * 4 * (N-1)
  #    where the 4 is to compensate for the two variable terms being doubled.
  #
  #  All of this provides optimization.
  #-----------------------------------------------------------------------
  my $sth = $dbh->prepare('SELECT count(*) FROM samples');
  $sth->execute();
  my ($num_samples) = $sth->fetchrow_array();

  $sth = $dbh->prepare(
    "SELECT stddev_samp(r) FROM generate_series(1,$num_samples) AS r");
  $sth->execute();
  my $sd = $sth->fetchrow_array();
  my $divisor = $sd * $sd * 4.0 * ($num_samples - 1.0);

  my $expr_arrays_subquery = <<GEN_EXPR_ARRAYS;
      SELECT
        gene_id,
        array_agg(2*expr_rank - $num_samples - 1  ORDER BY sample_id)
            AS expr_ranks_array
      FROM expr_data  JOIN gene_info USING(gene_id)
      WHERE not silent
      GROUP BY gene_id
GEN_EXPR_ARRAYS

  if (not pg_version_supports_array_agg_order_by($dbh)) {
    create_legacy_gen_expr_arrays($dbh);
    $expr_arrays_subquery = 'SELECT * FROM gen_expr_arrays()';
  }
 
  my $query = <<CORR_MATRIX;
INSERT INTO corr_matrix
  WITH 
    expr_arrays AS ( $expr_arrays_subquery )
  SELECT
    t1.gene_id AS gene_id_1,
    t2.gene_id AS gene_id_2,
    gnest_dot_product(t1.expr_ranks_array, t2.expr_ranks_array)::real
        / $divisor
      AS sp_corr
  FROM
    expr_arrays AS t1,
    expr_arrays AS t2
  WHERE
    t1.gene_id < t2.gene_id
CORR_MATRIX
  $dbh->do($query)  or  die "Error in computing correlations\n";
  #gnh_db::error_when_table_empty($dbh, 'corr_matrix');
}


#----------------------------------------------------------------------------
#  All possible gene neighborhoods are created in advance of computing
#  statistics for them.
#  A gene neighborhood is defined as a group of consecutively located genes
#  on a chromosome.
#
#  They can be created on the basis of base pair window sizes (with varying
#  numbers of genes) and/or a range of gene counts (varying window bp sizes).
#----------------------------------------------------------------------------
sub create_nh {
  my ($dbh, $min_win_size, $max_win_size, $min_gene_count, $max_gene_count) 
     = @_;

  my $query = <<QUERY_CREATE;
    CREATE TEMP TABLE tmp_nh_by_win_size(
      chromosome                 varchar,
      win_size                   integer,
      start_gene_pos_index       integer,
      num_genes                  integer
    );
QUERY_CREATE
  $dbh->do($query)  or  die "Error: CREATE tmp_nh_by_win_size\n";

  $query = "SELECT gnest_fill_nh_by_win_size($min_win_size, " .
           "    $max_win_size, $min_win_size)";
  $dbh->do($query)  or  die "Error gnest_fill_nh_by_win_size\n";

  $query = <<TMP_NH_BY_GC;
CREATE TEMP TABLE tmp_nh_by_gc AS
  SELECT * FROM gen_nh_by_gc($min_gene_count, $max_gene_count);
TMP_NH_BY_GC
    $dbh->do($query)  or  die "Error creating tmp_nh_by_gc\n";
    #gnh_db::error_when_table_empty($dbh, 'tmp_nh_by_gc');


  #--------------------------------------------------------------------------
  #  Generate the unique list of neighborhoods, regardless of whether
  #  they're determined by window size of num genes.
  #
  #  The 'nh_id' column uniquely identifies a particular set of genes,
  #  regardless of how the neighborhood was determined.
  #--------------------------------------------------------------------------

  $query = <<QUERY_4;
INSERT INTO gene_nh(chromosome, start_gene_pos_index, num_genes,
        chr_start_pos, nh_bp_size, genes_list, lower_bound, upper_bound)
  WITH 
    all_nh AS (
      SELECT chromosome, start_gene_pos_index, num_genes
      FROM tmp_nh_by_gc
      UNION
      SELECT chromosome, start_gene_pos_index, num_genes
      FROM tmp_nh_by_win_size
      ),
    pass_1 AS (
      SELECT DISTINCT
        all_nh.chromosome,
        start_gene_pos_index,
        all_nh.num_genes,
        first_gene_id + start_gene_pos_index
            AS start_gene_id,
        first_gene_id + start_gene_pos_index + all_nh.num_genes - 1
            AS end_gene_id
      FROM all_nh JOIN chromosome_info USING(chromosome)
      )
  SELECT
    pass_1.chromosome,
    start_gene_pos_index,
    num_genes,
    g_start.chr_start_pos,
    g_end.chr_start_pos - g_start.chr_start_pos + 1 
        AS nh_bp_size,
    (WITH g AS (SELECT generate_series(start_gene_id, end_gene_id) AS gene_id
     ORDER BY gene_id)   SELECT array_agg(gene_id) FROM g)
        AS genes_list,
    g_start.lower_bound,
    g_end.upper_bound
  FROM pass_1
    JOIN gene_info AS g_start ON (g_start.gene_id = start_gene_id)
    JOIN gene_info AS g_end   ON (g_end.gene_id   = end_gene_id)
  ORDER BY pass_1.chromosome, start_gene_pos_index, num_genes;
QUERY_4
  $dbh->do($query)  or  die "Error inserting into gene_nh\n";
  #gnh_db::error_when_table_empty($dbh, 'gene_nh');

  $query = <<GENE_TO_NH;
INSERT INTO gene_to_nh(nh_id, gene_id)
  SELECT nh_id, unnest(genes_list) AS gene_id
  FROM gene_nh
GENE_TO_NH
  $dbh->do($query)  or  die "Error: insert into gene_to_nh\n";

  #---------------------------------------------------------------------------
  #  Create a table of the neighborhoods by window sizes, with the nh_id added.
  #---------------------------------------------------------------------------
  $query = <<QUERY_5;
INSERT INTO nh_by_win_size(nh_id, win_size)
  SELECT nh_id, win_size
  FROM gene_nh JOIN tmp_nh_by_win_size 
      USING(chromosome, start_gene_pos_index, num_genes);
QUERY_5
  $dbh->do($query)  or  die "Error: INSERT INTO nh_by_win_size\n";

}


#--------------------------------------------------------------------------
#
#--------------------------------------------------------------------------
sub compute_nh_synteny_scores {
  my ($dbh) = @_;
  my $query;

  #--------------------------------------------------------------------------
  #  The number of genomes is used as a divisor in computing syntenic score.
  #--------------------------------------------------------------------------
  my $sth = $dbh->prepare(
      'SELECT count(DISTINCT target_taxon_id) FROM syntenic_blocks');
  $sth->execute();
  my $num_genomes = $sth->fetchrow_array();

  if ($num_genomes > 0) {
    $query = <<INSERT_NH_SYNTENY;
      INSERT INTO nh_synteny(nh_id, target_taxon_id)
        SELECT
          nh_id,
          target_taxon_id
        FROM 
          gene_nh AS g,
          syntenic_blocks AS s
        WHERE
          g.chromosome = s.chromosome  AND
          chr_start_pos BETWEEN start_pos and end_pos  AND
          (chr_start_pos + nh_bp_size - 1) BETWEEN start_pos and end_pos;
INSERT_NH_SYNTENY
    $dbh->do($query)  or  die "Error: INSERT INTO nh_synteny\n";

    $query = <<SYN_COUNTS;
      CREATE TEMP TABLE nh_syn_counts AS 
        SELECT nh_id, count(*) as synteny_count
        FROM nh_synteny GROUP BY nh_id;
SYN_COUNTS
    $dbh->do($query)  or  die "Error: CREATE TABLE nh_syn_counts\n";

    $query = <<SYN_UPDATE_1;
      UPDATE gene_nh AS g
      SET synteny_score = c.synteny_count::real / ${num_genomes}::real
      FROM nh_syn_counts AS c
      WHERE c.nh_id = g.nh_id;
SYN_UPDATE_1
    $dbh->do($query)  or  die "Error: UPDATE gene_nh (#1)\n";

    $query = <<SYN_UPDATE_2;
UPDATE gene_nh AS g
SET synteny_score = 0.0
WHERE synteny_score IS NULL;
SYN_UPDATE_2
    $dbh->do($query)  or  die "Error: UPDATE gene_nh (#2)\n";
  }

}


#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
sub compute_nh_stats {
  my ($dbh, $schema, $num_perms, $min_gc, $max_gc) = @_;

  my $path = $ENV{GNEST_BIN} ?  $ENV{GNEST_BIN} . '/' : '';

  my $cmd = "${path}gnest_anc --num_perms $num_perms " .
            "--min_gc $min_gc --max_gc $max_gc " .
      ($schema ? " --schema $schema" : '');
  my $rv = system($cmd);
  $rv==0  or  die "Error($rv) in gnest_anc\n";
}


#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
sub compute_tns {
  my ($dbh) = @_;

  my $query = <<UPDATE_NH_STATS;
    UPDATE nh_stats
    SET tns = CASE WHEN pvalue >= 0.05 THEN 0
                   ELSE synteny_score * (1.0 - pvalue) * anc END
    FROM gene_nh
    WHERE gene_nh.nh_id = nh_stats.nh_id;
UPDATE_NH_STATS
  $dbh->do($query)  or  die "Error: UPDATE nh_stats (synteny_score)\n";

}


#---------------------------------------------------------------------------
sub compute_nh_win_rank {
  my ($dbh, $use_tns) = @_;

  my $tns_clause = $use_tns ? 'tns DESC,' : '';
  
  my $query = <<CREATE_NH_RANK;
    CREATE TEMP TABLE tmp_nh_win_rank AS
      SELECT
        nh_id,
        nh_parm,
        rank() OVER (PARTITION BY nh_id
                     ORDER BY $tns_clause pvalue, nh_parm)
          AS nh_win_rank
      FROM nh_stats
      WHERE NOT by_gene_counts
CREATE_NH_RANK
  $dbh->do($query)  or  die "Error: CREATE TABLE tmp_nh_win_rank\n";

  $dbh->do('CREATE INDEX tmp_rank_index ON tmp_nh_win_rank(nh_id, nh_parm)')
      or  die "Error creating index on tmp_nh_win_rank\n";

  $query = <<UPDATE_NH_STATS_2;
    UPDATE nh_stats AS ns
    SET nh_win_rank = r.nh_win_rank
    FROM tmp_nh_win_rank AS r
    WHERE not ns.by_gene_counts            AND
          ns.nh_id          = r.nh_id      AND
          ns.nh_parm        = r.nh_parm;
UPDATE_NH_STATS_2
  $dbh->do($query)  or  die "Error: UPDATE nh_stats (nh_win_rank)\n";

  $query = <<UPDATE_NH_STATS_3;
    UPDATE nh_stats AS ns
    SET nh_win_rank = 1
    WHERE by_gene_counts
UPDATE_NH_STATS_3
  $dbh->do($query)  or  die "Error: UPDATE nh_stats #2 (nh_win_rank)\n";

  $dbh->do('DROP TABLE tmp_nh_win_rank') or die "Error dropping tmp_nh_ramk\n";
}

1;
