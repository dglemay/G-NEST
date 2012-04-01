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
package gnest_reports;
use strict;
use DBI;
use DBD::Pg;
use File::Temp('tempdir', 'tmpnam');


#------------------------------------------------------------------------
#  Generic report handler.
#------------------------------------------------------------------------
sub report {
  my ($dbh, $label, $dir, $fname, $query) = @_;

  my $cmd = sprintf("COPY (%s) TO '%s' CSV HEADER DELIMITER '\t';",
                    $query, "$dir/$fname");
  $dbh->do($cmd)  or  die "Error in $label report\n";
}


#------------------------------------------------------------------------
sub dump_corr_matrix {
  my ($dbh, $tmp_dir, $project) = @_;

  my $query = <<CORR_MATRIX_QUERY;
    SELECT
      g1.gene_name  AS gene_name_1,
      g2.gene_name  AS gene_name_2,
      CASE WHEN sp_corr = 0 THEN '0' ELSE sp_corr::numeric(3,2) END
          AS correlation
    FROM corr_matrix
    JOIN gene_info AS g1 ON(g1.gene_id = gene_id_1)
    JOIN gene_info AS g2 ON(g2.gene_id = gene_id_2)
    ORDER BY gene_name_1, gene_name_2
CORR_MATRIX_QUERY

  report($dbh, 'Correlation matrix ', $tmp_dir, "${project}_corr_matrix.csv",
         $query);
}


#------------------------------------------------------------------------
sub list_silent_genes {
  my ($dbh, $tmp_dir, $project) = @_;

  report($dbh, 'silent genes', $tmp_dir, "${project}_silent_genes.txt",
         'SELECT gene_name FROM gene_info WHERE silent ORDER BY gene_name');
}


#------------------------------------------------------------------------
#  Report on the gene neighborhoods.
#------------------------------------------------------------------------
sub nh_reports {
  my ($dbh, $tmp_dir, $project, @target_taxons) = @_;

  my $synteny_dependent_cols = '';
  my $genome_bool_cols = '';
  my $genome_from_list = '';

  if (scalar(@target_taxons)) {
    my $case_fmt =
        "  CASE WHEN s_%s.target_taxon_id IS NULL THEN 'N' ELSE 'Y' END\n" .
        "  AS taxon_%d";

    $genome_bool_cols =
      join(",\n", '', map(sprintf($case_fmt, $_, $_), @target_taxons));
    my $join_fmt = "  LEFT OUTER JOIN nh_synteny AS s_%d \n" .
                   "  ON (g.nh_id = s_%d.nh_id and s_%d.target_taxon_id = %d)";
    $genome_from_list = 
      join("\n", map(sprintf($join_fmt, $_, $_, $_, $_), @target_taxons));

    $synteny_dependent_cols = join(',', 'tns::numeric(4,2)',
                         'synteny_score::numeric(4,2)',
                         map("taxon_$_", @target_taxons), '');
  }

  my $synteny_query = <<SYNTENY_QUERY;
    CREATE TEMP TABLE tmp_nh_report_info AS
      WITH  pass_1 AS (
        SELECT
          nh_id, 
          array_to_string(array_agg(gene_name),',') AS genes
        FROM gene_to_nh JOIN gene_info USING(gene_id)
        GROUP BY nh_id
        )
      SELECT
        g.nh_id,
        genes
        $genome_bool_cols
      FROM
        pass_1 AS g
    $genome_from_list
SYNTENY_QUERY
  $dbh->do($synteny_query) 
        or die "Error: CREATE TEMP TABLE tmp_nh_report_info\n";

  $dbh->do('CREATE INDEX rep_nh_id ON tmp_nh_report_info(nh_id)');


  my $query_fmt_str = <<QUERY_FMT_STR;
      WITH pass_1 AS (
          SELECT nh_id, nh_parm, anc, pvalue, tns
          FROM nh_stats 
          WHERE nh_win_rank = 1 AND (%s by_gene_counts)
          )
      SELECT
        g.nh_id,
        chromosome,
        start_gene_pos_index,
        num_genes,
        chr_start_pos,
        nh_bp_size,
        lower_bound,
        upper_bound,
        %s
        anc::numeric(4,2),
        pvalue::numeric(5,3),
        tns::numeric(4,2),
        synteny_score::numeric(4,2),
        $synteny_dependent_cols 
        genes
      FROM pass_1 JOIN gene_nh AS g USING(nh_id)
                  JOIN tmp_nh_report_info USING(nh_id)
      ORDER BY g.nh_id
QUERY_FMT_STR

  report($dbh, 'Neighborhood by gene counts',
    $tmp_dir, $project . '_nh_report_by_gc.csv',
    sprintf($query_fmt_str, '', ''));

  report($dbh, 'Neighborhood by window size',
    $tmp_dir, $project . '_nh_report_by_win.csv',
    sprintf($query_fmt_str, 'NOT', 'nh_parm AS win_size,'));

  if (scalar(@target_taxons)) {
    $dbh->do('DROP TABLE tmp_nh_report_info')  or 
        die "Error: DROP tmp_nh_report_info\n";
  }
}


#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
sub report_genes_best_score {
  my ($dbh, $tmp_dir, $project, $use_tns) = @_;

  my $score_type = $use_tns ? 'tns' : 'pvalue';
  my $max_or_min = $use_tns ? 'max' : 'min';

  my $query_fmt_str = <<TOP_QUERY_FMT_STR;
    WITH pass_1 AS (
      SELECT gene_id, $max_or_min($score_type) AS $score_type 
      FROM nh_stats JOIN gene_to_nh USING(nh_id)
      WHERE %s by_gene_counts
      GROUP BY gene_id
      )
    SELECT chromosome, gene_id, gene_name, ${score_type}::numeric(5,3)
    FROM pass_1 JOIN gene_info USING(gene_id)
    ORDER BY gene_id
TOP_QUERY_FMT_STR

  my $fname_fmt_str = "${project}_best_${score_type}_by_%s.csv";

  my $fname_gc = sprintf($fname_fmt_str, 'gc');
  report($dbh, 'Gene top score by gene counts',
    $tmp_dir, sprintf($fname_fmt_str, 'gc'),
    sprintf($query_fmt_str, '') );

  report($dbh, 'Gene top score by window size',
    $tmp_dir, sprintf($fname_fmt_str, 'win'),
    sprintf($query_fmt_str, 'NOT') );
}


#-------------------------------------------------------------------------
#
#-------------------------------------------------------------------------
sub report_genes_tau {
  my ($dbh, $tmp_dir, $project) = @_;

  my $query = <<QUERY;
WITH pass_1 AS (   -- average the replicates
  SELECT 
    gene_id,
    bio_state,
    avg(expr)::numeric(8,2) AS avg_expr
  FROM expr_data JOIN samples USING(sample_id) 
  GROUP BY bio_state, gene_id)
SELECT   -- compute statistic per gene
  gene_name,
  gene_id,
  (count(*) - sum(avg_expr)/max(avg_expr))/(count(*)-1)  AS tau
FROM pass_1 JOIN gene_info USING(gene_id)
GROUP BY gene_id, gene_name
QUERY

  report($dbh, 'Tau', $tmp_dir, "${project}_tau.csv", $query);
}


#-------------------------------------------------------------------------
#  Dump out most of the tables (definitions and data), allowing the user
#  to recreate the database for his/her own queries.
#
#  The output syntax is specific to postgres.
#
#  The 'corr_matrix' table is omitted because it's huge, and maybe less useful.
#-------------------------------------------------------------------------
sub export_db {
  my ($dbh, $tmp_dir, $project, $schema, $include_corr_matrix) = @_;

  $dbh->do(sprintf("SELECT export_db('%s', %s, '%s/%s_db_export.sql')",
           $schema,
           $include_corr_matrix ? 'True' : 'False',
           $tmp_dir,
           $project))
    or die "Error: calling export_db\n";
}


#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
sub custom_tracks {
  my ($dbh, $use_tns, $project, $track_names, $win_fname, $gc_fname) = @_;

  my $description = $use_tns ? 
      'Total Neighborhood Score (1=good,0=bad)' :
      'Average Neighborhood Correlation (ANC)';

  my $stat_col = $use_tns ? 'tns' : 'anc';

  $track_names ||= $project . '_tracks';

  my $header = "track type=bedGraph name=\"$track_names\" " .
               "description=\"$description\" " .
               'visibility=full color=0,205,0 altColor=255,140,0 ' .
               'autoScale=off gridDefault=on graphType=points ' .
               'viewLimits=0:1 yLineMark=0.8 yLineOnOff=off ' .
               'windowingFunction=maximum smoothingWindow=off';

  my @cols;
  foreach my $not_by_gc ('', 'NOT') {
    my $fname = $not_by_gc ? $win_fname : $gc_fname;
#       sprintf("${project}_track_by_%s.csv", $not_by_gc ? 'win' : 'gc');
    open(OUTF, ">$fname")  or  die "Can't create file $fname\n";
    print OUTF $header, "\n";

    my $query = <<QUERY;
      SELECT
        'chr' || chromosome  AS chromosome,
        lower_bound,
        upper_bound,
        ${stat_col}::numeric(4,2)
      FROM gene_nh AS g
           JOIN nh_stats USING(nh_id)
      WHERE  nh_win_rank = 1 AND ($not_by_gc by_gene_counts)
      ORDER BY g.nh_id
QUERY
    my $sth = $dbh->prepare($query);
    $sth->execute();
    while (@cols = $sth->fetchrow_array()) {
      print OUTF join("\t", @cols), "\n";
    }
    close OUTF;
  }
}


1;
