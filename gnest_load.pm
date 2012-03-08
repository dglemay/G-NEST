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
package gnest_load;
use strict;
use DBI;
use DBD::Pg;
use gnest_table_loader;


our $chromosome_pattern = qr/"?(chr)?[12]?[\dXY]"?/i;
my $real_pattern = qr/[-+]?\d+(\.\d+)([Ee][-+]\d+)?/;
my $samples_info_loaded = 0;
my $genes_filtered = 0;


#------------------------------------------------------------------------
sub validate_options {
  my ($h_opts) = @_;

  die "You cannot filter on both MAS5 calls and expression threshold\n"
      if ($h_opts->{filter_on_mas5} and $h_opts->{filter_min_expr});

}


#------------------------------------------------------------------------
sub remove_quotes {
  local($_) = shift;

  if (/^".*"$/) {
     s/\A"(.*)"\z/\1/;
     s/\\"//g;
  }
  $_;
}


#------------------------------------------------------------------------
sub normalize_chromosome_name {
  my ($chrom) = @_;

  $chrom =~ remove_quotes($chrom);
  $chrom =~ s/^chr//i;

  uc($chrom);
}


#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
sub load_genes_and_chromosomes {
  my ($dbh, $galaxy_mode, $taxon_id, $chrom_length_file, $genes_info_file)
      = @_;

  my $query = <<CREATE_CHROM_LEN_QUERY;
    CREATE TEMP TABLE raw_chrom_lengths(
      chromosome                 varchar,
      chrom_length               integer
    )
CREATE_CHROM_LEN_QUERY
  $dbh->do($query);

  #------------------------------------------------------------------------
  #  If a chromosome lengths file was provided, then load it.
  #------------------------------------------------------------------------
  if ($chrom_length_file) {

    gnest_table_loader::load_table($dbh, {
      file_name => $chrom_length_file,
      file_tag => $galaxy_mode ? 'Chromosomes file' : $chrom_length_file,
      table_name => 'raw_chrom_lengths',
      columns => [
        { col_name => 'chromosome',   syntax => '\S+', might_be_quoted => 1 },
        { col_name => 'chrom_length', syntax => '\d+' }
        ],
      row_preprocessor => sub {
        my $row = shift;
        return undef
            unless ($row->{chromosome} =~ /\A$chromosome_pattern\z/);
        $row->{chromosome} = normalize_chromosome_name($row->{chromosome});
        $row;
        }
      } );
  }

  #------------------------------------------------------------------------
  #  Otherwise use the chromosome info for this taxon from the library.
  #  The scenario of the user not providing a file and it not being in the
  #  library is fatal.
  #------------------------------------------------------------------------
  else {
    $query = <<FROM_CHR_LIB;
      INSERT INTO raw_chrom_lengths(chromosome, chrom_length)
        SELECT chromosome, chrom_length
        FROM chrom_lengths_lib
        WHERE project_taxon_id = ?
FROM_CHR_LIB
    my $sth = $dbh->prepare($query);
    $sth->execute($taxon_id);

    if ($sth->rows() == 0) {
       die "You need to supply chromosomes lengths for taxon=$taxon_id";
    }
  }

  $query = <<CREATE_GENE_INFO;
    CREATE TEMP TABLE raw_gene_info(
      gene_name      varchar,
      chromosome     varchar,
      strand         integer,
      start_pos      integer,
      end_pos        integer,
      silent         boolean  DEFAULT false
    )
CREATE_GENE_INFO
  $dbh->do($query);

  gnest_table_loader::load_table($dbh, {
    file_name => $genes_info_file,
    file_tag => $galaxy_mode ? 'Genes info file' : $genes_info_file,
    table_name => 'raw_gene_info',
    columns => [
      { col_name => 'gene_name',  syntax => '\S+',  might_be_quoted => 1 },
      { col_name => 'chromosome', syntax => '\S+',  might_be_quoted => 1 },
      { col_name => 'strand',     syntax => '[-+]?1?',
        might_be_quoted => 1 },
      { col_name => 'start_pos',  syntax => '\d+',     },
      { col_name => 'end_pos',    syntax => '\d+',     },
      ],
    row_preprocessor => sub {
      my $row = shift;
      return undef 
          unless ($row->{chromosome} =~ /\A$chromosome_pattern\z/);
      $row->{strand} .= '1'  unless substr($row->{strand}, -1) eq '1';
      $row->{chromosome} = normalize_chromosome_name($row->{chromosome});
      $row;
    }
    } );
}


#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
sub load_synteny {
  my ($dbh, $galaxy_mode, $project_taxon_id, @synteny_inputs) = @_;

  foreach (0 .. scalar(@synteny_inputs)/2 - 1) {
    my $target_taxon_id = $synteny_inputs[2*$_];
    my $synteny_file = $synteny_inputs[2*$_+1];

    gnest_table_loader::load_table($dbh, {
      file_name => $synteny_file,
      file_tag => $galaxy_mode ? "Synteny file (taxon=$target_taxon_id)"
          : $synteny_file,
      table_name => 'syntenic_blocks',
      columns => [
        { col_name => 'chromosome', syntax => '\S+' },
        { col_name => 'start_pos',  syntax => '\d+', },
        { col_name => 'end_pos',    syntax => '\d+', },
        ],
      table_cols => [ qw(target_taxon_id chromosome start_pos end_pos) ],
      row_preprocessor => sub {
        my $row = shift;
        $row->{chromosome} = normalize_chromosome_name($row->{chromosome});
        $row->{target_taxon_id} = $target_taxon_id;

        $row;
        }
      } );
  }

  my $query = <<INSERT_SYNTENIC;
    INSERT INTO syntenic_blocks(target_taxon_id, chromosome,
                                start_pos, end_pos)
      WITH already_loaded AS (
        SELECT DISTINCT target_taxon_id FROM syntenic_blocks)
      SELECT 
        target_taxon_id,
        chromosome,
        start_pos,
        end_pos
      FROM syntenic_blocks_lib AS sbl
      WHERE project_taxon_id = $project_taxon_id  AND
            NOT EXISTS (
              SELECT NULL 
              FROM already_loaded AS al
              WHERE al.target_taxon_id = sbl.target_taxon_id
              )
INSERT_SYNTENIC
  $dbh->do($query);
}


#-------------------------------------------------------------------------
#  Optionally, load sample information
#-------------------------------------------------------------------------
sub load_samples {
  my ($dbh, $galaxy_mode, $samples_file) = @_;

  $samples_info_loaded = 1;

  my $query = <<CREATE_SAMPLES;
    CREATE TEMP TABLE raw_samples_info(
      sample_name                varchar,
      bio_state                  varchar,
      replicate                  integer
    )
CREATE_SAMPLES
  $dbh->do($query);

  gnest_table_loader::load_table($dbh, {
    file_name => $samples_file,
    table_name => 'raw_samples_info',
    columns => [
      { col_name => 'sample_name', syntax => '\S+',           },
      { col_name => 'bio_state',   syntax => '[^\t]+',        },
      { col_name => 'replicate',   syntax => qr/(rep)?\d+/i,  },
      ],
    row_preprocessor => sub {
      my $row = shift;
      $row->{sample_name} =~ s/.CEL$//i;
      $row->{replicate} =~ s/^rep//i;

      $row;
      }
    } );

  $query = <<INSERT_TMP_SAMPLES;
    INSERT INTO samples(sample_name, bio_state, replicate)
      SELECT sample_name, bio_state, replicate
      FROM raw_samples_info;
INSERT_TMP_SAMPLES
  $dbh->do($query);

  $dbh->do('DROP TABLE raw_samples_info');
}


#-------------------------------------------------------------------------
#  Parse the header for the expression file or the MAS5 file.
#-------------------------------------------------------------------------
sub parse_expr_header {
  my ($h_tab, $header) = @_;

  $h_tab->{first_data_processed} = 0;
  $h_tab->{sample_names} = [ split /\t/, $header ];
  foreach (@{$h_tab->{sample_names}}) {
    $_ = remove_quotes($_);
    s/\.CEL$//i;
  }

  $h_tab->{repeated_col_pattern} = $h_tab->{table_name} eq 'mas5_calls' ?
          '[APM]' : $real_pattern;

  $h_tab->{row_pattern} = '"?\S+"?' . "(\t$h_tab->{repeated_col_pattern})+";
}


#-------------------------------------------------------------------------
#  Parse a data row for the expression file or the MAS5 file.
#-------------------------------------------------------------------------
sub parse_expr_row {
  my ($h_tab, $gene_name, @expr_cols) = @_;
  my %expr_values;

  unless ($h_tab->{first_data_processed}) {
    shift @{$h_tab->{sample_names}}
        if (scalar(@{$h_tab->{sample_names}}) == scalar(@expr_cols) + 1);

    die "$h_tab->{file_tag}: Number of data columns doesn't match header\n"
        unless (scalar(@{$h_tab->{sample_names}}) == scalar(@expr_cols));
  
    $h_tab->{row_pattern} = '"?\S+\"?' .
            "(\t$h_tab->{repeated_col_pattern})" .
            '{' . scalar(@{$h_tab->{sample_names}}) . '}';

    $h_tab->{first_data_processed} = 1;
  }

  @expr_values{@{$h_tab->{sample_names}}} = @expr_cols;
  $gene_name = remove_quotes($gene_name);

  my @out_rows = map(
    { gene_name => $gene_name, sample_name => $_, expr => $expr_values{$_} },
    keys %expr_values);

  if ($h_tab->{table_name} eq 'mas5_calls') {
    $_->{is_present} = $_->{expr} eq 'P' ? 't' : 'f'   foreach (@out_rows);
  }

  map join("\t", @{ $_ }{ @{$h_tab->{table_cols}} }), @out_rows;
}


#-------------------------------------------------------------------------
#
#-------------------------------------------------------------------------
sub load_expression_data {
  my ($dbh, $galaxy_mode, $expr_data_file) = @_;

  my $query = <<CREATE_EXPR;
CREATE TABLE raw_expr_data(
  gene_name      varchar,
  sample_name    varchar,
  expr           real
);
CREATE_EXPR
  $dbh->do($query);

  gnest_table_loader::load_table($dbh, {
    file_name => $expr_data_file,
    file_tag => $galaxy_mode ? 'Expression file' : $expr_data_file,
    table_name => 'raw_expr_data',
    table_cols => [ qw(gene_name sample_name expr) ],
    process_header => \&parse_expr_header,
    custom_data_handler => \&parse_expr_row,
  } );


  #----------------------------------------------------------------------
  #  If samples information wasn't loaded, just use the names from the
  #  expression data, simply copying the sample_name as bio_state.
  #  In this case, the replicate column is meaningless.
  #----------------------------------------------------------------------
  unless ($samples_info_loaded) {
    $query = <<INSERT_SAMPLES;
      INSERT INTO samples(sample_name, bio_state, replicate)
        SELECT DISTINCT sample_name, sample_name, ''
        FROM raw_expr_data;
INSERT_SAMPLES
    $dbh->do($query);
  }
}


#-----------------------------------------------------------------------------
#  Only use genes that are present in all replicates of at least one bio_state.
#-----------------------------------------------------------------------------
sub filter_by_mas5_data {
  my ($dbh, $galaxy_mode, $mas5_calls_file) = @_;

  $genes_filtered = 1;

  my $query = <<CREATE_MAS5;
    CREATE TEMP TABLE mas5_calls(
      gene_name      varchar,
      sample_name    varchar,
      is_present     boolean
    )
CREATE_MAS5
  $dbh->do($query);

  gnest_table_loader::load_table($dbh, {
    file_name => $mas5_calls_file,
    file_tag => $galaxy_mode ? 'MAS5 calls file' : $mas5_calls_file,
    table_name => 'mas5_calls',
    table_cols => [ qw(gene_name sample_name is_present) ],
    process_header => \&parse_expr_header,
    custom_data_handler => \&parse_expr_row,
  } );


  #----------------------------------------------------------------------
  #  Keep genes that have a bio_state having all its replicates present.
  #----------------------------------------------------------------------
  $query = <<EXPRESSED_INSERT_1;
    CREATE TEMP TABLE expressed_genes AS 
      SELECT DISTINCT gene_name
      FROM mas5_calls JOIN samples USING(sample_name)
      GROUP BY gene_name, bio_state
      HAVING bool_and(is_present);
EXPRESSED_INSERT_1
  $dbh->do($query);

  $dbh->do('DROP TABLE mas5_calls');
}


#-----------------------------------------------------------------------------
#  Only use genes that meet the expression threshold in all replicates of
#  at least one bio_state.
#-----------------------------------------------------------------------------
sub filter_by_min_expr {
  my ($dbh, $threshold) = @_;

  $genes_filtered = 1;

  my $query = <<EXPRESSED_INSERT_2;
    CREATE TEMP TABLE expressed_genes AS 
      SELECT DISTINCT gene_id
      FROM raw_expr_data JOIN samples USING(sample_name)
      GROUP BY gene_id, bio_state;
      HAVING bool_and(expr > $threshold)
EXPRESSED_INSERT_2
  $dbh->do($query);
}


#-------------------------------------------------------------------------
#  The selection of which genes are to be kept (must be expressed) 
#  for analysis can't be done until after expression data is loaded.
#  Where multiple genes overlap, only the one with the highest maximum
#  expression (across all samples) is kept.
#  Those issues are handled here.
#-------------------------------------------------------------------------
sub finish_load {
  my ($dbh) = shift;

  my $query;

  if ($genes_filtered) {
    $dbh->do('UPDATE raw_gene_info SET silent = True')
        or die "Error:  raw_gene_info.silent := True\n";

    $query = <<UNSILENCE;
      UPDATE raw_gene_info AS r
      SET silent = False
      FROM expressed_genes AS e
      WHERE r.gene_name = e.gene_name
UNSILENCE
    $dbh->do($query)  or  die "Error: UNSILENCE\n";

    $dbh->do('DROP TABLE expressed_genes');
  }

  #---------------------------------------------------------------------
  #  If there are any chromosomes with fewer than 2 genes, omit them.
  #---------------------------------------------------------------------
  $query = <<OMIT_CHROMS;
    CREATE TEMP TABLE chroms_to_delete AS
      SELECT chromosome
      FROM raw_gene_info
      GROUP BY chromosome
      HAVING count(*) < 2
OMIT_CHROMS
  $dbh->do($query)  or  die "Error: chroms_to_delete\n";

  if ($dbh->rows() > 0) {
    $query = <<DELETE_CHROMS;
      DELETE FROM chromosome_info AS c
      USING chroms_to_delete AS d
      WHERE c.chromosome = d.chromosome
DELETE_CHROMS
    $dbh->do($query)  or  die "Error: delete_chroms\n";

    $query = <<DELETE_GENES_2;
      DELETE FROM raw_gene_info AS g
      USING chroms_to_delete AS d
      WHERE g.chromosome = d.chromosome
DELETE_GENES_2
    $dbh->do($query);
  }
  $dbh->do('DROP TABLE chroms_to_delete')  or
      die "Error: DROP chroms_to_delete\n";


  #-------------------------------------------------------------------------
  #  Cluster overlapping genes.
  #  Then put select genes in the final 'gene_info' table.
  #  For overlapping genes, use the gene with the highest max expression.
  #
  #  'gene_id' integer key values are automatically populated.
  #-------------------------------------------------------------------------
  $dbh->do('CREATE TEMP TABLE gene_cluster_map AS ' .
           'SELECT * FROM cluster_overlapping_genes()' );
  $dbh->do('CREATE INDEX cluster_index ON gene_cluster_map(gene_name)');

  my $query = <<SELECT_GENES;
    INSERT INTO gene_info(chromosome, gene_name, chr_start_pos, gene_pos_index,
                          lower_bound, upper_bound, silent)
      WITH
        pass_1 AS (
          SELECT
            gene_name,
            max(expr) AS max_expr
          FROM raw_expr_data
          JOIN raw_gene_info USING(gene_name)
          GROUP BY gene_name
          ),
        pass_2 AS (
          SELECT
            cluster_start_pos,
            gene_name,
            row_number() OVER (PARTITION BY chromosome, cluster_start_pos
                         ORDER BY max_expr DESC)
              AS expr_rank
          FROM pass_1 JOIN gene_cluster_map USING(gene_name)
          )
        SELECT
          chromosome,
          gene_name,
          (CASE WHEN strand > 0 THEN start_pos ELSE end_pos END)
              AS chr_start_pos,
          rank() OVER (PARTITION BY chromosome ORDER BY cluster_start_pos)
              AS gene_pos_index,
          start_pos AS lower_bound,
          end_pos   AS upper_bound,
          silent
        FROM pass_2 JOIN raw_gene_info USING(gene_name)
        WHERE expr_rank = 1
        ORDER BY chromosome, cluster_start_pos
SELECT_GENES
  $dbh->do($query)  or  die "Error:  (INSERT INTO gene_info)\n";

  $dbh->do('DROP TABLE gene_cluster_map');

  #----------------------------------------------------------------------
  $query = <<INSERT_CHROMS;
INSERT INTO chromosome_info(chromosome, num_genes, first_gene_id, chrom_length)
  WITH
    counts AS (SELECT chromosome, count(*) AS num_genes 
               FROM gene_info GROUP BY chromosome),
    gene_1 AS (SELECT chromosome, gene_id AS first_gene_id
               FROM gene_info WHERE gene_pos_index = 1)
  SELECT chromosome, num_genes, first_gene_id, chrom_length
  FROM counts JOIN gene_1 USING(chromosome)
       JOIN raw_chrom_lengths USING(chromosome)
INSERT_CHROMS
  $dbh->do($query);


  #----------------------------------------------------------------------
  #  The expression data is filtered to include only the genes of interest
  #  and the samples of interest by joining them to those restricted tables.
  #  Also, the output table uses the internal gene_id instead of the
  #  input gene_name for a key.
  #
  #  Gene expression ranks are computed here, because they will be needed
  #  in computation of Spearman's correlation later on.
  #----------------------------------------------------------------------
  my $query = <<INSERT_EXPR;
    INSERT INTO expr_data(gene_id, sample_id, expr, expr_rank)
      SELECT
        gene_id,
        sample_id,
        expr,
        rank() OVER (PARTITION BY gene_id ORDER BY expr DESC)
          AS expr_rank
      FROM raw_expr_data
           JOIN gene_info USING(gene_name)
           JOIN samples USING(sample_name)
      ORDER BY gene_id, sample_id
INSERT_EXPR
  $dbh->do($query);

  $dbh->do('DROP TABLE raw_chrom_lengths, raw_gene_info, raw_expr_data');

}

1;
