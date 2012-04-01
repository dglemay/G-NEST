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
package gnest_graphs;
use strict;
use DBI;
use DBD::Pg;
use gnest_db;

#--------------------------------------------------------------------------
#
#--------------------------------------------------------------------------
sub make_graphs {
  my ($dbh, $project, $graphs_dir, $graphs_title, $graphics_format,
      $use_tns,
      $min_win_size, $max_win_size,
      $min_gene_count, $max_gene_count) = @_;

  my $create_legends_sth = $dbh->prepare('SELECT generate_legends(?,?)');
  $create_legends_sth->execute($graphs_dir, $graphics_format);

  my $query = <<CREATE_TMP_GRID;
    CREATE TEMP TABLE tmp_grid_values(
      x          integer,
      y          integer,
      pvalue     real,
      tns        real,
      anc        real
    );
CREATE_TMP_GRID
  $dbh->do($query)  or  die "Error creating tmp_grid_values\n";

  my $build_grid_sth = $dbh->prepare('SELECT build_grid(?,?,?,?,?,?,?,?)');

  my $graph_sth = $dbh->prepare('SELECT generate_graphs(?,?,?,?,?,?,?,?,?)');

  my $chrom_sth = $dbh->prepare('SELECT * FROM chromosome_info');
  $chrom_sth->execute();

  my $subdirs_created = 0;
  while (my $h_chrom = $chrom_sth->fetchrow_hashref()) {

    foreach my $by_gene_counts (0, 1) {  # boolean
      foreach my $x_axis_bp (0, 1) {     # boolean

        my $subdir = sprintf("%s_%s_%s", 
                $project,
                $by_gene_counts ? 'gc' : 'win',
                $x_axis_bp ? 'bp' : 'gi');

        gnest_db::create_readable_subdir($graphs_dir, $subdir)
          unless $subdirs_created;

        $dbh->do('TRUNCATE TABLE tmp_grid_values');

        #-----------------------------------------------------------------
        #  Build a grid of all cells in the rectangular graph.
        #  Cell data includes pvalue, anc, and tns, for 3 different graphs
        #  to be produced from each grid.
        #-----------------------------------------------------------------
        $build_grid_sth->execute(
            $h_chrom->{chromosome},
            $h_chrom->{num_genes},
            $h_chrom->{first_gene_id},
            $h_chrom->{chrom_length},
            $by_gene_counts ? $min_gene_count : $min_win_size,
            $by_gene_counts ? $max_gene_count : $max_win_size,
            $by_gene_counts ? 't' : 'f',
            $x_axis_bp      ? 't' : 'f'
            )   or die "Error: calling build_grid()\n";

        #-----------------------------------------------------------------
        #  Invoke an R function to produce image files.
        #-----------------------------------------------------------------
        $graph_sth->execute(
            $graphs_title,
            $h_chrom->{chromosome},
            $h_chrom->{chrom_length},
            $use_tns ? 't' : 'f',
            $by_gene_counts ? 't' : 'f',
            $x_axis_bp      ? 't' : 'f',
            "$graphs_dir/$subdir",
            sprintf("%s_chr%s_%%s_%s_%s",
                $project,
                $h_chrom->{chromosome},
                $by_gene_counts ? 'gc' : 'win',
                $x_axis_bp ? 'bp' : 'gi',
                ),
            $graphics_format,
            )  or  die "Error: calling graph_nh_stat\n";
      }
    }
    $subdirs_created = 1;
  }

  $dbh->do('DROP TABLE tmp_grid_values')  or
      die "Error dropping tmp_grid_values\n";
}

1;
