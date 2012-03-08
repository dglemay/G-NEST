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
SET SESSION client_min_messages = WARNING;
CREATE EXTENSION plr;
CREATE LANGUAGE plperl;
CREATE LANGUAGE plperlu;


CREATE TABLE chrom_lengths_lib(
  project_taxon_id           integer,
  chromosome                 varchar,
  chrom_length               integer
);

CREATE TABLE syntenic_blocks_lib(
  project_taxon_id           integer,
  target_taxon_id            integer,
  chromosome                 varchar,
  start_pos                  integer,
  end_pos                    integer
);

CREATE TYPE win_nh_info AS (
  win_size                   integer,
  start_gene_pos_index       integer,
  num_genes                  integer
);

CREATE TYPE nh_stats_grid_cell AS (
  win_size                   integer,
  blk_index                  integer,
  value                      real
);

CREATE TYPE gb_nh_info AS (
  chromosome                 varchar,
  start_gene_pos_index       integer,
  num_genes                  integer
);

CREATE TYPE expr_array AS (
  gene_id                    integer,
  expr_ranks_array           integer[]
);

CREATE TYPE gene_cluster_map AS (
  chromosome                 varchar,
  cluster_start_pos          integer,
  gene_name                  varchar
);



----------------------------------------------------------------------

CREATE OR REPLACE FUNCTION create_tables() RETURNS void AS $$


CREATE TABLE gene_info(
  gene_id                    SERIAL PRIMARY KEY,
  gene_name                  varchar,
  chromosome                 varchar,
  chr_start_pos              integer,
  gene_pos_index             integer,
  lower_bound                integer,
  upper_bound                integer,
  silent                     boolean
);

CREATE TABLE chromosome_info(
  chromosome                 varchar PRIMARY KEY,
  num_genes                  integer,
  first_gene_id              integer,
  chrom_length               integer
);

CREATE TABLE samples(
  sample_id                  SERIAL PRIMARY KEY,
  sample_name                varchar,
  bio_state                  varchar,
  replicate                  varchar
);

CREATE TABLE syntenic_blocks(
  target_taxon_id            integer,
  chromosome                 varchar,
  start_pos                  integer,
  end_pos                    integer
);

CREATE TABLE expr_data(
  gene_id                    integer,
  sample_id                  integer,
  expr                       real,
  expr_rank                  integer
);

CREATE TABLE corr_matrix(
  gene_id_1                  integer,
  gene_id_2                  integer,
  sp_corr                    real
);

CREATE TABLE gene_nh(
  nh_id                      SERIAL PRIMARY KEY,
  chromosome                 varchar,
  start_gene_pos_index       integer,
  num_genes                  integer,
  chr_start_pos              integer,
  nh_bp_size                 integer,
  lower_bound                integer,
  upper_bound                integer,
  genes_list                 integer[],
  synteny_score              real
);

CREATE TABLE gene_to_nh(
  nh_id                      integer,
  gene_id                    integer
);
CREATE INDEX gn_nh_index ON gene_to_nh(nh_id);
CREATE INDEX gn_gene_index ON gene_to_nh(gene_id);

CREATE TABLE nh_synteny(
  nh_id                      integer,
  target_taxon_id            integer
);

CREATE TABLE nh_by_win_size(
  nh_id                      integer,
  win_size                   integer
);

CREATE TABLE nh_stats(
  nh_id                      integer,
  by_gene_counts             boolean,
  nh_parm                    integer,  -- win_size or num_genes
  anc                        real,
  pvalue                     real,
  tns                        real,
  nh_win_rank                integer,

  UNIQUE(nh_id, by_gene_counts, nh_parm)
);

$$ LANGUAGE SQL;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION blk_stats_query( 
    p_blk_size               integer, 
    p_max_win_size           integer, 
    p_mode                   varchar,  -- 'TNS', 'ANC', or 'pvalue' 
    p_chromosome             varchar,  -- null for all 
    p_chrom_length           integer)
    RETURNS varchar AS $$
BEGIN
  RAISE EXCEPTION 'You must load the function "blk_stats_query"';
END;
$$ LANGUAGE plpgsql;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gnh_fill_nh_by_win_size(
  p_min_win_size           integer,
  p_max_win_size           integer,
  p_win_size_increment     integer)
  RETURNS void AS $$
BEGIN
  RAISE EXCEPTION
      'You must load library for C function "gnh_fill_nh_by_win_size"';
END;
$$ LANGUAGE plpgsql;
    

-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION dot_product(  
    p_gene_1_expr_ranks      integer[],
    p_gene_2_expr_ranks      integer[])
    RETURNS integer AS $$
BEGIN
  RAISE EXCEPTION 'You must load function "dot_product"';
END;
$$ LANGUAGE plpgsql;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION make_expr_arrays()
    RETURNS SETOF expr_array AS $$
BEGIN
  RAISE EXCEPTION 'You must load function "make_expr_arrays"';
END;
$$ LANGUAGE plpgsql;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION gen_nh_by_gc(
    p_min_gene_count         integer,
    p_max_gene_count         integer)
    RETURNS SETOF gb_nh_info AS $$
BEGIN
  RAISE EXCEPTION 'You must load function "gen_nh_by_gc"';
END;
$$ LANGUAGE plpgsql;


-------------------------------------------------------------------------
CREATE OR REPLACE FUNCTION get_nh_stats_grid(
    p_blk_size               integer,
    p_mode                   varchar,  -- 'TNS', 'ANC', or 'pvalue'
    p_chromosome             varchar)
    RETURNS SETOF nh_stats_grid_cell AS $$
BEGIN
  RAISE EXCEPTION 'You must load function "get_nh_stats_grid"';
END;
$$ LANGUAGE plpgsql;

