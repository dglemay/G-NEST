/*#########################################################################
*  Copyright (C) 2009-2012 William F. Martin and Danielle G. Lemay
* 
*  This program is free software; you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by the
*  Free Software Foundation;
* 
*  This program is distributed in the hope that it will be useful, but
*  WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*  See the GNU General Public License for more details.
* 
*  You should have received a copy of the GNU General Public License along
*  with this program; if not, write to the Free Software Foundation, Inc.,
*  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*-#########################################################################

CREATE OR REPLACE FUNCTION  gnest_fill_nh_by_win_size(
    p_min_win_size           integer,
    p_max_win_size           integer,
    p_win_size_increment     integer)
    RETURNS void AS $$

*/

#include "postgres.h"
#include <stdlib.h>
#include <sys/stat.h>
#include "executor/spi.h"
#include "utils/builtins.h"
#include "utils/array.h"
#include "fmgr.h"
#include "catalog/pg_type.h"
#include "funcapi.h"
#include "utils/lsyscache.h"
#include "utils/fmgroids.h"

#include "c.h"


static SPIPlanPtr insert_plan;


//-------------------------------------------------------------------------
void get_chromosomes(
    char ***chromosomes,
    int *num_chromosomes)
{
   bool is_null;

   if (SPI_exec("SELECT chromosome FROM chromosome_info ", 1000)
                != SPI_OK_SELECT)
     elog(ERROR, "Error in reading chromosome_info");

   *num_chromosomes = SPI_processed;
   *chromosomes = (char **)palloc(*num_chromosomes * sizeof(char *));

   int i;
   for (i=0; i<*num_chromosomes; ++i) {
     char *val = SPI_getvalue(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1);
     (*chromosomes)[i] = (char *)palloc(strlen(val)+1);
     strcpy((*chromosomes)[i], val);
   }
}


//-------------------------------------------------------------------------
void get_start_positions(
    char *chromosome,
    int **chr_start_positions,
    int *chr_num_genes)
{
  char query[120];
  bool is_null;

  sprintf(query, "SELECT chr_start_pos FROM gene_info WHERE chromosome = '%s' "
                 "ORDER BY chr_start_pos", chromosome);

  if (SPI_exec(query, 500000) != SPI_OK_SELECT)
    elog(ERROR, "Error in reading gene_info");

  *chr_num_genes = SPI_processed;
  *chr_start_positions = (int *)palloc(*chr_num_genes * sizeof(int *));
  int i;
  for (i=0; i<*chr_num_genes; ++i) {
    (*chr_start_positions)[i] = DatumGetInt32( SPI_getbinval(
        SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &is_null));
  }
}


//------------------------------------------------------------------
//  First find the neighborhood for each combination gene
//  (start position on the chromosome) and window size.
//------------------------------------------------------------------
void build_grid(int *chr_start_pos_list, int chr_num_genes,
        int min_win_size, int max_win_size, int win_size_increment,
        int *num_win_sizes, int **grid_ptr)
{
  int s_pos, end_pos, ws_ind, nh_length, grid_bytes;
  int *grid;

  *num_win_sizes = 1 + (max_win_size-min_win_size)/win_size_increment;

  grid_bytes = *num_win_sizes * chr_num_genes * sizeof(int);
  grid = *grid_ptr = (int *)palloc(grid_bytes);
  memset(grid, 0, grid_bytes);

  // Iterate through start positions in the array of genes
  for (s_pos=0; s_pos < chr_num_genes-1; ++s_pos) {

    // Iterate through ranges of genes from largest to smallest
    for (end_pos=chr_num_genes-1; end_pos>s_pos; --end_pos) {
      nh_length = chr_start_pos_list[end_pos] - chr_start_pos_list[s_pos] + 1;

      //----------------------------------------------------------------
      //  Iterate through the list of window sizes:
      //  For each combination of window size and starting gene,
      //  save the largest number of genes that fit in that window.
      //----------------------------------------------------------------
      for (ws_ind=0; ws_ind<*num_win_sizes; ++ws_ind) {
        if (nh_length <= min_win_size + ws_ind*win_size_increment &&
            grid[ ws_ind*chr_num_genes + s_pos] == 0) {

//elog(WARNING, "into grid: ws_ind=%d, chr_num_genes=%d, s_pos=%d, end_pos=%d",
//     ws_ind, chr_num_genes, s_pos, end_pos);

          grid[ ws_ind*chr_num_genes + s_pos] = end_pos - s_pos + 1;
        }
      }
    }
  }
}


//---------------------------------------------------------------------
void save_row(
    Datum *chromosome_datum,
    int win_size,
    int start_gene_pos_index,
    int num_genes)
{
  static Datum param_vals[4];

  if (!insert_plan) {
    Oid params[] = { VARCHAROID, INT4OID, INT4OID, INT4OID };
    insert_plan = SPI_prepare("INSERT INTO tmp_nh_by_win_size("
        "chromosome, win_size, start_gene_pos_index, num_genes) "
        "VALUES($1,$2,$3,$4)", 4, params);
    if (!insert_plan)
      elog(ERROR, "Error preparing INSERT INTO tmp_nh_by_win_size: %d",
           SPI_result);
  }

  param_vals[0] = *chromosome_datum;
  param_vals[1] = Int32GetDatum(win_size);
  param_vals[2] = Int32GetDatum(start_gene_pos_index);
  param_vals[3] = Int32GetDatum(num_genes);

  int ret_val = SPI_execute_plan(insert_plan, param_vals, NULL, 0, 1L);
  if (ret_val != SPI_OK_INSERT) {
    elog(ERROR, "Error inserting into tmp_nh_by_win_size: %d", 
           SPI_result);
  }
}


//---------------------------------------------------------------------
#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif


PG_FUNCTION_INFO_V1(gnest_fill_nh_by_win_size);
Datum gnest_fill_nh_by_win_size(PG_FUNCTION_ARGS) {
  Datum chrom_datum;
  VarChar *vc;
  char **chromosomes;
  int num_chromosomes;
  char *chromosome;
  int *chr_start_positions;
  int num_genes, chr_num_genes;
  int *grid, num_win_sizes;

  insert_plan = 0;
  
  int min_win_size       = PG_GETARG_INT32(0);
  int max_win_size       = PG_GETARG_INT32(1);
  int win_size_increment = PG_GETARG_INT32(2);

  SPI_connect();

  get_chromosomes(&chromosomes, &num_chromosomes);

  int c;
  for (c=0; c<num_chromosomes; ++c) {
    chromosome = chromosomes[c];
    vc = palloc(30);
    strcpy(VARDATA(vc), chromosome);
    SET_VARSIZE(vc, strlen(chromosome) + VARHDRSZ);
    chrom_datum = PointerGetDatum(vc);

    get_start_positions(chromosome, &chr_start_positions, &chr_num_genes);

//elog(WARNING, "chromosome=%s\n", chromosome);
    build_grid(chr_start_positions, chr_num_genes, 
               min_win_size, max_win_size, win_size_increment,
               &num_win_sizes, &grid);

    int s_pos, ws_ind;
    for (s_pos=0; s_pos < chr_num_genes-1; ++s_pos) {
      for (ws_ind=0; ws_ind<num_win_sizes; ++ws_ind) {
        num_genes = grid[chr_num_genes*ws_ind + s_pos];
        if (num_genes >= 2) {
          save_row(&chrom_datum, min_win_size+ws_ind*win_size_increment,
                   s_pos, num_genes);
        }
      }
    }
    pfree(grid);
    pfree(chr_start_positions);
  }

  SPI_finish();
}
