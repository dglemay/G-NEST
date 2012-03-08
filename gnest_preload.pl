#!/usr/bin/perl
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
use strict;
BEGIN { unshift @INC, $ENV{GNEST_LIB} if ($ENV{GNEST_LIB}); }

use Getopt::Long;
use gnest_db;
use gnest_table_loader;
use gnest_load;


#-----------------------------------------------------------------------
sub preload_chromosome_lengths {
  my ($dbh, $project_taxon_id, $chrom_length_file) = @_;

   gnest_table_loader::load_table($dbh, {
      file_name => $chrom_length_file,
      file_tag => $chrom_length_file,
      table_name => 'chrom_lengths_lib',
      columns => [
        { col_name => 'chromosome',   syntax => '\S+', might_be_quoted => 1 },
        { col_name => 'chrom_length', syntax => '\d+' }
        ],
      table_cols => [ qw(project_taxon_id chromosome chrom_length) ],
      row_preprocessor => sub {
        my $row = shift;
        $row->{project_taxon_id} = $project_taxon_id;
        return undef
           unless ($row->{chromosome} =~ /\A$gnest_load::chromosome_pattern\z/);
        $row->{chromosome} =
          gnest_load::normalize_chromosome_name($row->{chromosome});
        $row;
        }
      } );
}


#-----------------------------------------------------------------------
sub add_to_syntenic_lib {
  my ($dbh, $project_taxon_id, @synteny_args) = @_;

  my $num_pairs = scalar(@synteny_args);

  while (@synteny_args) {
    my $target_taxon_id = shift @synteny_args;
    my $synteny_file    = shift @synteny_args;

    gnest_table_loader::load_table($dbh, {
      file_name => $synteny_file,
      file_tag => $synteny_file,
      table_name => 'syntenic_blocks_lib',
      columns => [
        { col_name => 'chromosome', syntax => '\S+' },
        { col_name => 'start_pos',  syntax => '\d+', },
        { col_name => 'end_pos',    syntax => '\d+', },
        ],
      table_cols =>
        [ qw(project_taxon_id target_taxon_id chromosome start_pos end_pos) ],
      row_preprocessor => sub {
        my $row = shift;
        $row->{chromosome} =
            gnest_load::normalize_chromosome_name($row->{chromosome});
        $row->{project_taxon_id} = $project_taxon_id;
        $row->{target_taxon_id}  = $target_taxon_id;
        $row;
        }
      } );
  }
}


{ ####  MAIN  ###############################################

  my $dbh = gnest_db::db_connect();

  my %opts;
  GetOptions(\%opts, 'project_taxon_id=i', 'chromosomes=s');

  preload_chromosome_lengths($dbh, $opts{project_taxon_id}, $opts{chromosomes});

  add_to_syntenic_lib($dbh, $opts{project_taxon_id}, @ARGV);
}

