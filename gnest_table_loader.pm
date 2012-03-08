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
package gnest_table_loader;
#-------------------------------------------------------------------------
#  This module encapsulates issues of uploading data to database tables.
#  It used "COPY" instead of "INSERT" for efficency.
#-------------------------------------------------------------------------
use strict;
use DBI;
use DBD::Pg;
use File::Spec;
use Carp::Assert;

#-------------------------------------------------------------------------
#  Members of the table information hash :
#   file tag (for use in error messages, if other than the file name)
#   file_name
#   table_name
#   table_cols
#   copy_data_verbatum (boolean)
#   row_preprocessor (call-back function gets row hash and returns it after
#                     modifying it or returns undef to skip a row)
#
#   process_header (call-back function given the header line assigns the 
#                   regexp for subsequent data lines)
#   custom_data_handler (call-back function gets raw columns and returns
#                        list of row hash)
#
#   columns;  an array of hashes with these fields:
#     input_tag (for use in error messages)
#     col_name
#     syntax  (regexp, if null, then '\S+')
#     might_be_quoted  (boolean)
#-------------------------------------------------------------------------


#------------------------------------------------------------------------
sub remove_quotes {
  local($_) = shift;

  if (/\A"(.*)"\z/) {
    $_ = $1;
    s/\\"/"/g;
  }
  $_;
}


#------------------------------------------------------------------------
#--  Return the regexp pattern for the subseqent data rows.
#------------------------------------------------------------------------
sub process_header {
  my ($h_tab, $header) = @_;

  if (exists $h_tab->{process_header}) {
    $h_tab->{copy_data_verbatum} = 0;
    &{$h_tab->{process_header}}($h_tab, $header);
    return;
  }

  $h_tab->{copy_data_verbatum} = not exists $h_tab->{row_preprocessor};
  my @col_patterns = ();
  $h_tab->{input_cols} = [];
  foreach my $col (@{$h_tab->{columns}}) {
     push @{$h_tab->{input_cols}}, $col->{col_name};

     my $pattern = $col->{syntax} || '\S+';
     if ($col->{might_be_quoted}) {
       $pattern = '"?' . $pattern . '"?';
       $h_tab->{copy_data_verbatum} = 0;
     }
     push @col_patterns, $pattern;
  }

  $h_tab->{table_cols} ||= $h_tab->{input_cols};

  $h_tab->{row_pattern} = join("\t", @col_patterns);
}


#------------------------------------------------------------------------
sub handle_syntax_error {
  my ($h_tab, $line_num, $line) = @_;
  
  my @cols = split /\t/, $line;
  foreach my $i (0 .. scalar(@{$h_tab->{columns}})-1) {
    my $col_spec = $h_tab->{columns}[$i];
    if ($col_spec->{syntax} && $cols[$i] !~ /\A$col_spec->{syntax}\z/) {
      die "Syntax error in file (\"$cols[$i]\" doesn't match pattern " .
        "\"$col_spec->{syntax}\" " .
        "($h_tab->{file_tag} line #$line_num, " . 
        "column '$col_spec->{col_name}')\n";
    }
  }
  
  die "Syntax error in file ($h_tab->{file_tag} line #$line_num)\n";
}


#------------------------------------------------------------------------
sub start_copy_command {
  my ($dbh, $h_tab) = @_;

  my $copy_query = "COPY $h_tab->{table_name}(" .
      join(',', @{$h_tab->{table_cols}}) . ") FROM STDIN";

  my $sth = $dbh->prepare($copy_query);
  $sth->execute()  || die $sth->errstr . " table=$h_tab->{table_name}";
}


#------------------------------------------------------------------------
sub prepare_copy_row {
  my ($h_tab, $line) = @_;

  return ( $line )  if ($h_tab->{copy_data_verbatum});

  my @cols_data = split /\t/, $line;

  return &{$h_tab->{custom_data_handler}}($h_tab, @cols_data)
      if (exists $h_tab->{custom_data_handler});

  foreach (0..$#cols_data) {
    $cols_data[$_] = remove_quotes($cols_data[$_])
        if ($h_tab->{columns}[$_]{might_be_quoted});
  }
 
  my $row = {};
  @{$row}{ @{$h_tab->{input_cols}}} = @cols_data;

  $row = &{$h_tab->{row_preprocessor}}($row) 
      if (exists $h_tab->{row_preprocessor});

  $row ? ( join("\t", @{$row}{ @{$h_tab->{table_cols}}}) ) : ();
}


#------------------------------------------------------------------------
sub load_table {
  my ($dbh, $h_tab) = @_;

  $h_tab->{file_tag} ||= $h_tab->{file_name};

  my $abs_path = File::Spec->rel2abs($h_tab->{file_name});
  die "File $abs_path doesn't exist\n"  unless -f $abs_path;
  open(TABLE_INPUT, $abs_path)  or  die "Can't open $abs_path\n";

  my $header = <TABLE_INPUT>;
  $header =~ s/[\r\n]*$//;
  process_header($h_tab, $header);

  start_copy_command($dbh, $h_tab);
  
  my $line_num = 2;
  my $line;
  my @out_lines;
  while ($line = <TABLE_INPUT>) {
    $line =~ s/[\r\n]+$//;

    if ($line !~ /\A$h_tab->{row_pattern}\z/) {
      $dbh->func('pg_putcopyend');
      handle_syntax_error($h_tab, $line_num, $line);
    }

    @out_lines = prepare_copy_row($h_tab, $line);

    foreach (@out_lines) {
      $dbh->func($_ . "\n", 'putline')  or 
          die $dbh->errstr . " ($h_tab->{file_tag} line #$line_num)";
    }
    ++$line_num;
  }

  $dbh->func('pg_putcopyend')  or  die $dbh->errstr;
  close(TABLE_INPUT);
}

1;
