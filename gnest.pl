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
use File::Temp('tmpnam', 'tempdir');
use Cwd('abs_path');

use gnest_db;
use gnest_load;
use gnest_compute;
use gnest_reports;
use gnest_graphs;

my (%opts, @synteny_inputs, $rv, $start_time, $cmd);


#-------------------------------------------------------------------------
sub progress {
  return unless $opts{progress};

  my $tag = shift;

  use integer;
  my $elapsed = time() - $start_time;
  my $minutes = $elapsed / 60;
  my $seconds = $elapsed % 60;
  printf "%d:%02d => $tag\n", $minutes, $seconds;
}


#-------------------------------------------------------------------------
sub check_file_exists {
  local($_) = shift;
  die "File $_ doesn't exist\n"  unless -f;
}


#-------------------------------------------------------------------------
#  Populates global variables  %opts and @synteny_inputs
#-------------------------------------------------------------------------
sub parse_parameters {

  my @option_specs = (
    {
      tag => 'project_taxon_id', suffix => '=i',
      required => 1,
    },
    
    {
      tag => 'project', suffix => '=s',
      validate => sub {
        local($_) = shift; 
        die "Project name must be start with a letter and otherwise be " .
            "alphanumeric or underscore characters"
            unless (/\A[a-z]\w+\z/i);
      },
    },
   
    {
      tag => 'tracks_name', suffix => '=s',
    },

    {
      tag => 'chromosomes', suffix => '=s',
    },

    {
      tag => 'genes', suffix => '=s',
      validate => \&check_file_exists,
      required => 1,
    },

    {
      tag => 'samples', suffix => '=s',
    },

    {
      tag => 'expr_data', suffix => '=s',
      validate => \&check_file_exists,
      required => 1,
    },

    {
      tag => 'filter_on_mas5', suffix => '=s',
      validate => \&check_file_exists,
    },

    {
      tag => 'filter_min_expr', suffix => '=f',
    },

    {
      tag => 'min_win_size', suffix => '=i',
      default => 100000,
    },

    {
      tag => 'max_win_size', suffix => '=i',
      default => 10000000,
    },

    {
      tag => 'min_gene_count', suffix => '=i',
      default => 2,
    },

    {
      tag => 'max_gene_count', suffix => '=i',
      default => 10,
    },

    {
      tag => 'num_permutations', suffix => '=i',
      default => 1000,
    },

    {
      tag => 'graphs_title', suffix => '=s',
    },

    {
      tag => 'graphics_format', suffix => '=s',
      default => 'pdf',
    },

    map ({ tag=>$_, suffix=>'!' },
     qw(corr_matrix allow_overlapping_genes export_db keep_project no_synteny
        progress galaxy)),

    { tag => 'reports_output',    suffix => '=s', },
    { tag => 'graphs_output',     suffix => '=s', },
    { tag => 'tracks_gc_output',  suffix => '=s', },
    { tag => 'tracks_win_output', suffix => '=s', },
  );

  GetOptions(\%opts, map("$_->{tag}$_->{suffix}", @option_specs) );

  foreach my $h_opt_spec (@option_specs) {
    my $val = $opts{$h_opt_spec->{tag}};

    if ($val) {
      &{$h_opt_spec->{validate}}($val)  if ($h_opt_spec->{validate});
    }

    else {
      die "Missing required parameter: $h_opt_spec->{tag}\n"
          if ($h_opt_spec->{required});

      $opts{$h_opt_spec->{tag}} ||= $h_opt_spec->{default}
          if ($h_opt_spec->{default});
    }
  }

  die "Options filter_on_mas5 and filter_min_expr are mutually exclusive\n"
      if ($opts{filter_on_mas5} and $opts{filter_min_expr});

  if ($opts{galaxy}) {

    die "Error: reports_output, graphs_output, tracks_gc_output, and " .
        "tracks_win_output all required in galaxy mode\n"
        unless ($opts{reports_output} and $opts{graphs_output} and
                $opts{tracks_gc_output} and $opts{tracks_win_output});
  }
  else {
    die "Error: reports_output, graphs_output, tracks_gc_output, and " .
        "tracks_win_output are only applicable in galaxy mode\n"
        if ($opts{reports_output} or $opts{graphs_output} or
            $opts{tracks_gc_output} or $opts{tracks_win_output});
  }

  if ($opts{max_win_size} % $opts{min_win_size}) {
    die "Error:  max_win_size must be a multiple of min_win_size\n";
  }

  $opts{project} ||= sprintf("pr_%d", $opts{project_taxon_id});

  return  if ($opts{no_synteny});

  @synteny_inputs = @ARGV;
}


#-------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------
sub make_zip_file {
  my ($in_dir, $out_file) = @_;

  my $out_file = abs_path($out_file);
  my $tmp_out = ($out_file =~ /\.zip$/) ? $out_file : tmpnam() . '.zip';

  $cmd = "cd $in_dir; zip -qr $tmp_out .";
  $rv = system($cmd);
  $rv==0  or  die "Error($rv): $cmd\n";

  unless ($tmp_out eq $out_file) {
    $cmd = "cp $tmp_out $out_file";
    $rv = system($cmd);
    $rv==0  or  die "Error($rv): $cmd\n";
  }
}


{ ####  MAIN  ###############################################
  my @target_taxons = ();

  $start_time = time();
  parse_parameters();

  if ($opts{progress}) {
    $start_time = time();
    my $t = localtime($start_time);
    print "START: $t\n";
  }

  my $dbh = gnest_db::db_connect();
  my $schema = gnest_db::db_project_init($opts{project}, $opts{keep_project});
  progress("Connect to database, schema=$schema");

  gnest_load::load_genes_and_chromosomes($dbh, $opts{galaxy},
      $opts{project_taxon_id}, $opts{chromosomes}, $opts{genes});
  progress("Load genes and chromosomes");

  unless ($opts{no_synteny}) {
    gnest_load::load_synteny($dbh, $opts{galaxy}, $opts{project_taxon_id},
                             @synteny_inputs);
    progress("Load synteny");

    my $sth =
        $dbh->prepare('SELECT DISTINCT target_taxon_id FROM syntenic_blocks');
    $sth->execute();
    my $a_vals = $sth->fetchall_arrayref();
    @target_taxons = map $_->[0], @$a_vals;
  }

  if (exists($opts{samples})) {
    gnest_load::load_samples($dbh, $opts{galaxy}, $opts{samples});
    progress("Load samples info");
  }

  gnest_load::load_expression_data($dbh, $opts{galaxy}, $opts{expr_data});
  progress("Load expression data");

  if ($opts{filter_on_mas5}) {
    gnest_load::filter_by_mas5_data($dbh, $opts{galaxy},
                                    $opts{filter_on_mas5});
    progress("Filter on present mas5 calls");
  }

  if ($opts{filter_min_expr}) {
    gnest_load::filter_by_min_expr($dbh, $opts{filter_min_expr});
    progress("Filter on min expression");
  }

  gnest_load::finish_load($dbh, $opts{allow_overlapping_genes});
  progress("Finish load");

  gnest_compute::correlation($dbh);
  progress("Compute correlations");

  gnest_compute::create_nh($dbh,
      @opts{qw(min_win_size max_win_size min_gene_count max_gene_count)});
  progress("Create neighborhoods");

  gnest_compute::compute_nh_stats($dbh, $schema, $opts{num_permutations},
      $opts{min_gene_count}, $opts{max_gene_count});
  progress("Compute neighborhood statistics");

  if (scalar(@target_taxons) > 0) {
    gnest_compute::compute_nh_synteny_scores($dbh);
    progress("Compute nh synteny scores");

    gnest_compute::compute_tns($dbh);
    progress("Compute tns");
  }

  gnest_compute::compute_nh_win_rank($dbh, scalar(@target_taxons) > 0);
  progress('compute_nh_win_rank');
 

#--------------------------------------------------------------------------
#  All individual report files (except tracks files) are created in a
#  temporary directory writable by the 'postgres' user.
#  At the end, they are zipped up into output files (created by linux user).
#
#  Galaxy insists on specifying the output file names (ending in ".dat"),
#  even if the extension doesn't make sense for the file type.
#  In that case, files must be created with names (extensions) that make
#  sense for the file type and later be renamed to Galaxy's names.
#--------------------------------------------------------------------------
  my $reports_dir = gnest_db::get_tmp_dir();

  gnest_reports::dump_corr_matrix($dbh, $reports_dir, $opts{project})
      if ($opts{corr_matrix});

  gnest_reports::list_silent_genes($dbh, $reports_dir, $opts{project});

  gnest_reports::nh_reports($dbh, $reports_dir, $opts{project},
                            @target_taxons);

  gnest_reports::report_genes_best_score($dbh, $reports_dir, $opts{project},
                                           scalar(@target_taxons)>0 );

  gnest_reports::export_db($dbh, $reports_dir, $opts{project}, $schema,
                           $opts{corr_matrix})
      if ($opts{export_db});

  make_zip_file($reports_dir, 
      $opts{reports_output} || "$opts{project}_reports.zip");
  progress('Reports');


#--------------------------------------------------------------------------
#  Each custom tracks file (2 of them, by gene counts and by window size)
#  is kept separate, but compressed with gzip.
#
#  For galaxy, the gzipped output file must be renamed to match the name
#  passed by galaxy.
#--------------------------------------------------------------------------
  my $scratch_dir = tempdir(CLEANUP=>1);

  my $tracks_file_fmt = "$opts{project}_track_by_%s.csv";
  my $tmp_tracks_win_output = sprintf($tracks_file_fmt, 'win');
  my $tmp_tracks_gc_output  = sprintf($tracks_file_fmt, 'gc');

  gnest_reports::custom_tracks($dbh,
      scalar(@target_taxons)>0, $opts{project}, '',
      "$scratch_dir/$tmp_tracks_win_output",
      "$scratch_dir/$tmp_tracks_gc_output");

  $cmd = "cd $scratch_dir; gzip -c $tmp_tracks_win_output > " .
      Cwd::abs_path($opts{tracks_win_output} ? $opts{tracks_win_output} :
                    "$tmp_tracks_win_output.gz");
  $rv = system($cmd);
  $rv==0  or  die "Error($rv): $cmd\n";

  $cmd = "cd $scratch_dir; gzip -c $tmp_tracks_gc_output > " .
      Cwd::abs_path($opts{tracks_gc_output} ? $opts{tracks_gc_output} :
                    "$tmp_tracks_gc_output.gz");
  $rv = system($cmd);
  $rv==0  or  die "Error($rv): $cmd\n";
  progress('Custom tracks');

#--------------------------------------------------------------------------
#  Create the graphs in a temp directory owned by 'postgres' user
#  (required because stored procedures execute as 'postgres').
#  Then zip the directory contents into a file created by logged in user.
#--------------------------------------------------------------------------
  my $graphs_dir = gnest_db::get_tmp_dir();

  gnest_graphs::make_graphs($dbh, $opts{project}, $graphs_dir,
      $opts{graphs_title} || "$opts{project}: Gene Neighborhoods",
      $opts{graphics_format},
      scalar(@target_taxons)>0,
      $opts{min_win_size}, $opts{max_win_size},
      $opts{min_gene_count}, $opts{max_gene_count});

  make_zip_file($graphs_dir, 
      $opts{graphs_output} || "$opts{project}_graphs.zip");
  progress('Graphs');

}
