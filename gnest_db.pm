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
package gnest_db;
use strict;
use DBI;
use DBD::Pg;

my ($schema, $dbh, $keep_project, @tmp_dirs);


#-------------------------------------------------------------------------
sub db_connect {
  my $dbname = $ENV{GNEST_DB} || 'gnest';
  my $db_user = $ENV{GNEST_DB_USER} || $ENV{USER};
  my $db_host = $ENV{GNEST_DB_HOST} || 'localhost';
  $dbh = DBI->connect("DBI:Pg:database=$dbname;host=$db_host", $db_user)
            or die "Can't open gnest database";
}


#-------------------------------------------------------------------------
#  Call this after connecting to the 'gnest' database.
#
#  Create a new schema named based upon the project name and configure the 
#  session to choose the new schema by default over the public schema.
#  If no project name is provided, use the public schema instead.
#
#  Create tables that need to be segregated data for the project.
#-------------------------------------------------------------------------
sub db_project_init {
  my ($project);
  ($project, $keep_project) = @_;

  $dbh->do('SET SESSION client_min_messages = WARNING')  
      or die "Error setting client_min_messages\n";

  my $so_dir = $ENV{GNEST_LIB} || $ENV{PWD};
  $dbh->do(sprintf('SET dynamic_library_path TO \'%s:$libdir\'', $so_dir))
      or die "Error setting dynamic_library_path\n";

  if ($project) {
    my $sth = $dbh->prepare('SELECT create_schema(?)');
    $sth->execute($project);
    ($schema) = $sth->fetchrow_array();
  }
  
  $dbh->do("SELECT create_tables()")  or die "Can't create database tables\n";

  $schema;
}


#-------------------------------------------------------------------------
#  Return the path to a temp directory writable by the 'postgres' user.
#-------------------------------------------------------------------------
sub get_tmp_dir {
  my $sth = $dbh->prepare('SELECT create_tmp_dir()');
  $sth->execute();
  push @tmp_dirs, $sth->fetchrow_array();

  $tmp_dirs[-1];
}


#-------------------------------------------------------------------------
sub create_readable_subdir {
  my $sth = $dbh->prepare('SELECT create_readable_subdir(?,?)');
  $sth->execute(@_);
}


#-------------------------------------------------------------------------
#  This is useful for debugging.
#-------------------------------------------------------------------------
sub get_row_count {
  my ($table_name) = @_;

  my $sth = $dbh->prepare("SELECT count(*) FROM $table_name");
  $sth->execute();
  my ($count) = $sth->fetchrow_array();

  $count;
}


#-------------------------------------------------------------------------
#  When finished (successfully or not), delete all the project's data.
#  However, if the 'keep_project' option was selected (and successful), 
#  then preserve the database.
#-------------------------------------------------------------------------
END {
  $dbh->do("SELECT purge_tmp_dir('$_')")  foreach (@tmp_dirs);

  $dbh->do("DROP SCHEMA $schema CASCADE")
      if ($schema and ($? or !$keep_project));

  $dbh->disconnect();
}

1;

