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

 INPUT parameters: 
  1) Database name.
  2) Number of random permutations to be done.
  3) Chromosome list subset (optional; one token with comma separation).

 INPUT TABLES: chromosome_info, gene_info, corr_matrix, gene_nh,
               nh_by_win_size

 OUTPUT:
  nh_stats(nh_id, by_gene_counts, nh_parm, anc, pvalue)
  where:
    by_gene_counts (boolean) indicates if the computation is based upon
                   constant number of gene counts (as opposed to window size)
    nh_parm        (integer) is either gene count or window size, depending
                   upon the value of 'by_gene_counts'

 PREREQUISITES:
  Before running thie program, the following must be done in the database:
  1) Correlation matrix built which includes at least all genes of interest
     (table is 'corr_matrix').
  2) Table 'gene_nh' contains info about all gene neighborhoods to be analyzed.
  3) If analysis is to be done by sliding bp window size, then 'nh_by_win_size'
     needs to be populated.

  The 'corr_matrix' table is abstractly a 2-dim matrix[i,j] with both indices
  being the 'gene_id' values (which start from one) in the 'gene_info' table.
  The 'gene_id' values in all the neighborhood-related tables are one-based
  indices into the matrix.  Given that 'C' used zero based indexing, in this
  program one is subtracted from the gene_ids to be used as indices.

 ASSUMPTIONS about input data:
  1) Gene ids are in order of position on chromosomes (within a chromosome).

  2) All genes in the 'gene_info' table are in the correlation matrix.

  3) Gene neighborhoods IDs are consecutive and in order of
     start_gene_pos_index, num_genes  (within a chromosome).

 ALGORITHM:

  1) Fetch the number of genes in the 'gene_info' table.

  2) Read the correlation matrix into a big array, representing an NxN matrix.
     Verify that the matrix size and indices matches what
     would be expected for the number of genes (from #1).

  3) Determine the next chromosome to process.  For each chromosome:

    A) Retrieve the information about gene neighborhoods from gene_nh table.
       The algorithm computing ANC values is ignorant of whether these
       gene neighborhoods are based upon base pair window sizes or
       upon number of genes.

    B) Iterate through random permutations of the gene ids.
       For each iteration:
       i)  Generate random gene id mapping (except first is identity mapping)
       ii) For each neighborhood:
           a) Compute total of correlation values
           b) Compute and save the ANC.
       (Note that these values can be used to calculate pvalues for either
        a given window size or a given number of genes.)

    C) Compute pvalues (how the real ANC compares to values with randomly
       shuffled genes).

  GENE ID MAPPING to zero-based indices:
    When reading the correlation matrix and gene neighborhood information 
    from the database, substract one from the gene_id values so they can
    be used as zero-based indices into the matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "libpq-fe.h"  // postgres calls

#define GENE_ID_TYPE unsigned short

static PGconn *conn;
static PGresult *res;
static ExecStatusType status;

static short *corr_matrix;
static GENE_ID_TYPE *genes_list = 0;
static int N;         // total number of genes

static int num_allocs = 0;
static void *alloc_list[20];

typedef struct {
  char chromosome[10];
  int  genes_list_offset;
  int  num_genes;
}  CHROM_INFO;

typedef struct {
  int nh_id;
  int start_gene_pos_index;
  int num_genes;
  short anc;
  short *rand_anc;
  float p_value;
} NH_INFO;

typedef struct {
  int      start_gene_pos_index;
  int      num_nh;
  NH_INFO *nh;
} NH_GROUP;


//-----------------------------------------------------------------------
void parse_opts(int argc, char **argv,
    char *schema, int *num_perms,
    int *min_gc, int *max_gc)
{
  char *env;
  int i, c, option_index;
  struct option long_options[] = {
    {"schema",         required_argument,  0, 's'},
    {"num_perms",      required_argument,  0, 'n'},
    {"min_gc",         required_argument,  0, 'm'},
    {"max_gc",         required_argument,  0, 'x'}
  };

  *schema = 0;
  *num_perms = 0;

  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1)
  {
    switch (option_index) {
      case 0:  strcpy(schema, optarg);           break;
      case 1:  *num_perms = atoi(optarg);        break;
      case 2:  *min_gc    = atoi(optarg);        break;
      case 3:  *max_gc    = atoi(optarg);        break;
    }
  }

  if (!*num_perms > 0) {
    printf("Usage:\n  %s  <OPTIONS>  [ chromosome ] ... \n"
      "  where <OPTIONS> are:\n"
      "   [ --schema <?> ] \n"
      "     --num_perms <?> \n"
      "     --min_gc    <?> \n"
      "     --max_gc    <?> \n"
      "\n"
      "   schema         is the default database schema for the session\n"
      "   num_perms      is the number of random gene permutations desired\n"
      "   min_gc         is the minimum gene count in neighborhoods \n"
      "   max_gc         is the minimum gene count in neighborhoods \n",
      argv[0]
      );
    exit(1);
  }

}


//-----------------------------------------------------------------------
static void close_pg() {
  PQfinish(conn);
}


//-----------------------------------------------------------------------
static PGconn *my_connect()
{
  char conn_str[100];
  char *ptr = conn_str;
  PGconn *t_conn;

  char *db_name = getenv("GNEST_DB");
  if (!db_name)  db_name = "gnest";

  char *db_user = getenv("GNEST_DB_USER");
  if (!db_user)  db_user = getenv("USER");

  char *db_host = getenv("GNEST_DB_HOST");
  if (!db_host)  db_host = "localhost";

  ptr += sprintf(ptr, "dbname=%s host=%s user=%s", db_name, db_host, db_user);
  t_conn = PQconnectdb(conn_str);
  if (PQstatus(t_conn) != CONNECTION_OK) {
    fprintf(stderr, "Connection to database failed: %s", PQerrorMessage(conn));
    exit(1);
  }
  return t_conn;
}


//-----------------------------------------------------------------------
//  Memory allocated here is freed after processing each chromosome
//  (rather than having the rest of the application keep track of memory).
//-----------------------------------------------------------------------
static void * my_malloc(int num_bytes) {
  alloc_list[num_allocs] = malloc(num_bytes);
  return (alloc_list[num_allocs++]);
}


//-----------------------------------------------------------------------
static void my_free() {
  int i;

  for (i=0; i<num_allocs; ++i)
    free(alloc_list[i]);

  num_allocs = 0;
}


//-----------------------------------------------------------------------
static int get_number_of_genes() {
  int max_gene_id, num_rows;

  res = PQexec(conn, "SELECT max(gene_id), count(*) FROM gene_info");
  max_gene_id = (int)atoi(PQgetvalue(res, 0, 0));
  num_rows = (int)atoi(PQgetvalue(res, 0, 1));
  PQclear(res);

  if (num_rows > 65535) {   // must change #define GENE_ID_TYPE 
    printf("Number of genes (%d) exceeds capacity of program capacity\n"
           "Program must be recompiled after changing index type\n", num_rows);
    exit(1);
  }

  if (max_gene_id != num_rows) {
    printf("Error in gene_info data: max_gene_id (%d) != count (%d)\n",
           max_gene_id, num_rows);
    exit(1);
  }

  return num_rows;
}


//-----------------------------------------------------------------------
//  corr_matrix  is indexed by gene_ids (minus one to make them zero based).
//-----------------------------------------------------------------------
inline static int corr_index(int i, int j) {
  return N*i + j;
}


//-----------------------------------------------------------------------
//  Read the correlation data into an array (representing a matrix).
//  Get the array count first to allow for an exact memory allocation.
//  Initialize to zeroes; the database stores a sparse matrix such that
//  entries are missing for silent genes (and thus assumed to be zero).
//
//  Compute the number of genes from the number of rows and
//  query the maximum gene_id in the table.
//-----------------------------------------------------------------------
static void read_corr_matrix() {
  int gene_id_1, gene_id_2, num_rows, num_bytes;
  int alloc_len;
  float corr_value;
  char *buf;

  res = PQexec(conn,
      "COPY corr_matrix(gene_id_1, gene_id_2, sp_corr) TO STDIN");
  status = PQresultStatus(res);
  if (status != PGRES_COPY_OUT) {
    printf("Error (COPY corr_matrix): %s\n", PQresStatus(status));
    exit(1);
  }

  alloc_len = N*N;

  corr_matrix = (short *)malloc(alloc_len * sizeof(short));
  memset(corr_matrix, 0, alloc_len * sizeof(short));

  num_rows = 0;
  while ((num_bytes = PQgetCopyData(conn, &buf, 0)) > 0) {
    sscanf(buf, "%d %d %f", &gene_id_1, &gene_id_2, &corr_value);
    corr_matrix[ corr_index(gene_id_1-1, gene_id_2-1) ] =
      corr_matrix[ corr_index(gene_id_2-1, gene_id_1-1) ] =
        (short)rintf((float)32767.0 * corr_value);
    PQfreemem(buf);
    ++num_rows;
  }

  if (num_bytes == -2) {
    printf("COPY ERROR: %s\n", PQerrorMessage(conn));
    exit(1);
  }

  PQclear(res);
}


//-----------------------------------------------------------------------
//  Return NULL if there are no more chromosomes to process.
//-----------------------------------------------------------------------
static CHROM_INFO * get_next_chromosome()
{
  static int first_time = 1;
  static CHROM_INFO *chroms;
  static int num_chroms;
  static int curr_chrom_ind = 0;

  int i;

  if (first_time) {
    res = PQexec(conn, 
        "SELECT chromosome, num_genes, first_gene_id FROM chromosome_info");
    status = PQresultStatus(res);
    if (status != PGRES_TUPLES_OK) {
      printf("Error (read gene_nh): %s\n", PQresStatus(status));
      exit(1);
    }
    num_chroms = PQntuples(res);
    chroms = (CHROM_INFO *)malloc(num_chroms * sizeof(CHROM_INFO));

    for (i=0; i<num_chroms; ++i) {
      strcpy(chroms[i].chromosome,  (char *)PQgetvalue(res, i, 0));
      chroms[i].num_genes         = (int)atoi(PQgetvalue(res, i, 1));
      chroms[i].genes_list_offset =((int)atoi(PQgetvalue(res, i, 2))) - 1;
    }

    PQclear(res);
  }
  
  return( (curr_chrom_ind < num_chroms) ? &chroms[curr_chrom_ind++] : NULL);
}


//--------------------------------------------------------------------------
static void init_genes_list()
{
  int i;

  if (!genes_list)
    genes_list = malloc(N * sizeof(GENE_ID_TYPE));
    
  for (i=0; i<N; ++i) 
    genes_list[i] = (GENE_ID_TYPE) i;
}


//--------------------------------------------------------------------------
static void randomize_genes_list() {
  GENE_ID_TYPE i, j, tmp;

  for (i=N-1; i>=1; --i) {
    j = rand() % (i+1);
    tmp = genes_list[i];
    genes_list[i] = genes_list[j];
    genes_list[j] = tmp;
  }
}
  

//--------------------------------------------------------------------
//  Gene neighborhoods are fetched in the order for which their ANC
//  values will be computed (at least within each chromosome).
//--------------------------------------------------------------------
static void get_gene_nh_for_chrom(CHROM_INFO *chrom_info, int num_perms,
                           NH_INFO **nh_list, int *num_nh, int *min_nh_id,
                           NH_GROUP **nh_groups, int *num_groups)
{
  int i;
  char query[200];
  NH_GROUP *nh_group = 0;
  NH_INFO *nh;
  short *anc_buf;

  sprintf(query, "SELECT nh_id, start_gene_pos_index, num_genes "
          "FROM gene_nh "
          "WHERE chromosome = '%s' "
          "ORDER BY start_gene_pos_index, num_genes", 
          chrom_info->chromosome);
  res = PQexec(conn, query);
  status = PQresultStatus(res);
  if (status != PGRES_TUPLES_OK) {
    printf("Error (read gene_nh): %s\n", PQresStatus(status));
    exit(1);
  }
  *num_nh = PQntuples(res);

  //-----------------------------------------------------------------------
  //  Buffers are allocated here and used in smaller pieces.
  //  Neighborhoods are grouped by starting gene to facilitate processing.
  //-----------------------------------------------------------------------
  anc_buf =    (short *)   my_malloc(sizeof(short) * num_perms * (*num_nh) );
  *nh_list =   (NH_INFO *) my_malloc(sizeof(NH_INFO) * (*num_nh));
  *nh_groups = (NH_GROUP *)my_malloc(sizeof(NH_GROUP) * chrom_info->num_genes);
  *num_groups = 0;
  nh_group = *nh_groups;

  for (i=0, nh=*nh_list;  i<*num_nh;  ++i, ++nh) {

    nh->rand_anc = &anc_buf[i*num_perms];

    nh->nh_id                 = (int)atoi(PQgetvalue(res, i, 0));
    //nh->start_gene_pos_index  =((int)atoi(PQgetvalue(res, i, 1))) - 1;
    nh->start_gene_pos_index  =((int)atoi(PQgetvalue(res, i, 1)));
    nh->num_genes             = (int)atoi(PQgetvalue(res, i, 2));

    if (i==0)   *min_nh_id = nh->nh_id;

    if (i==0 || nh_group->start_gene_pos_index != nh->start_gene_pos_index) {
      nh_group = &(*nh_groups)[(*num_groups)++];
      nh_group->nh = nh;
      nh_group->start_gene_pos_index = nh->start_gene_pos_index;
      nh_group->num_nh = 0;
    }
    ++nh_group->num_nh;
  }

  PQclear(res);
}


//-----------------------------------------------------------------------
//  Correlation totals (used to compute ANC) are computed in stair-step fashion 
//  with the totals of the smaller gene neighborhoods being part of the
//  totals for the larger gene neighborhoods (sharing the same start gene).
//  This avoids re-adding the same pairwise gene correlations.
//-----------------------------------------------------------------------
static void compute_anc_values_for_nh_group(CHROM_INFO *chrom_info,
                                  NH_GROUP *nh_group, int perm_ind)
{
  int i, j;
  int max_num_genes = nh_group->nh[nh_group->num_nh-1].num_genes;
  int min_gene_pos_index =
        chrom_info->genes_list_offset + nh_group->start_gene_pos_index;
  int max_gene_pos_index = min_gene_pos_index + max_num_genes - 1;
  int corr_total;
  int nh_ind=0;
  int nh_num_genes = nh_group->nh[nh_ind].num_genes;
  int last_gene_in_nh =
      min_gene_pos_index + nh_group->nh[nh_ind].num_genes - 1;
  short anc;

  corr_total = 0;
  for (i=min_gene_pos_index+1; i<=max_gene_pos_index; ++i) {
    for (j=min_gene_pos_index; j<i; ++j) {
      corr_total += corr_matrix[ corr_index(genes_list[i], genes_list[j]) ];
    }

    if (i == min_gene_pos_index + nh_num_genes - 1) {  // last gene in nh
      anc = (short) rintf(
          (float)corr_total/(float)(nh_num_genes * (nh_num_genes-1)/2) );
      if (perm_ind==-1)
        nh_group->nh[nh_ind].anc = anc;
      else
        nh_group->nh[nh_ind].rand_anc[perm_ind] = anc;

      if (nh_ind < nh_group->num_nh)
        nh_num_genes = nh_group->nh[++nh_ind].num_genes;
    }
  }
}


//-----------------------------------------------------------------------
//  Neighborhood groups share the same start gene.
//-----------------------------------------------------------------------
static void compute_anc_values(CHROM_INFO *chrom_info,
        NH_GROUP *nh_groups, int num_nh_groups, int num_perms)
{
  int p, grp;

  for (p=-1; p<num_perms; ++p) {
    if (p != -1)   randomize_genes_list();   // -1 for real ANC (not random)

    for (grp=0; grp<num_nh_groups; ++grp) {
      compute_anc_values_for_nh_group(chrom_info, &nh_groups[grp], p);
    }
  }
}


//------------------------------------------------------------------------
static int compare_shorts(const void *v_a, const void *v_b) {
  short *a = (short *)v_a;
  short *b = (short *)v_b;

  return(*a==*b ? 0: (*a<*b ? -1 : 1));
}


//------------------------------------------------------------------------
//  Return the number of elements in sorted array cmpScores that are greater
//  than score.  This is basically a binary search except we don't require
//  an exact match, just the index of the smallest member of cmpScores that
//  is greater than score.
//------------------------------------------------------------------------
static int find_greaters(short score, short *cmpScores, int cmpScoresCount)
{
  int leftI = 0, rightI = cmpScoresCount-1;
  while (rightI > leftI+1) {
    int midI = (rightI + leftI) / 2;
    if (cmpScores[midI] <= score)
        leftI = midI;
    else
        rightI = midI;
  }

  // Just in case score is completely outside of the range of cmpScores:
  if (cmpScores[leftI] > score)
      rightI = leftI;
  if (cmpScores[rightI] <= score)
    rightI++;

  return (cmpScoresCount - rightI);
}


//-----------------------------------------------------------------------
//  "COPY" is used instead of "INSERT" for optimization.
//-----------------------------------------------------------------------
static void save_pvalues(int by_gene_counts, int parm,
                         NH_INFO **sel_nh, int num_sel_nh, int num_perms,
                         short *anc_values)
{
  char *fmt_str = by_gene_counts ?  "%d\tt\t%d\t%f\t%f\n" :
                                    "%d\tf\t%d\t%f\t%f\n";
  int i, rv, anc_buf_len = num_sel_nh*num_perms;
  float f_anc, pvalue;
  NH_INFO *nh;
  char outbuf[50];

  qsort(anc_values, anc_buf_len, sizeof(short), compare_shorts);

  res = PQexec(conn,
      "COPY nh_stats(nh_id, by_gene_counts, nh_parm, anc, pvalue) FROM STDIN");
  status = PQresultStatus(res);
  if (status != PGRES_COPY_IN) {
    printf("Error (COPY nh_stats): %s\n", PQresStatus(status));
    exit(1);
  }
  PQclear(res);

  for (i=0; i<num_sel_nh; ++i) {
    nh = sel_nh[i];
    pvalue = (float)find_greaters(nh->anc, anc_values, anc_buf_len) /
             (float)anc_buf_len;
    f_anc = (float)nh->anc/32767.0;   // convert back from short to float

    sprintf(outbuf, fmt_str, nh->nh_id, parm, f_anc, pvalue);

    rv = PQputCopyData(conn, outbuf, strlen(outbuf));
    if (rv == -1) {
      printf("COPY error on record:%s\n", outbuf);
      exit(1);
    }
  }
  PQputCopyEnd(conn, 0);
}


//-----------------------------------------------------------------------
//  If rows were populated into the 'nh_by_win_size' table, then compute
//  pvalues for the neighborhoods with varying window sizes.
//-----------------------------------------------------------------------
static void compute_pvalues_by_win_size(CHROM_INFO *chrom_info,
        NH_INFO *nh_list, int num_nh, int min_nh_id, int num_perms)
{
  char query[200];
  int i, nh_id, w, rv, num_rows, win_size, num_win_sizes, *win_sizes;
  short *anc_values;
  int num_sel_nh;
  NH_INFO **sel_nh;

  //----------------------------------------------------------------------
  //  Get window sizes in descending size order, used in subsequent queries.
  //----------------------------------------------------------------------
  sprintf(query, "SELECT DISTINCT win_size FROM nh_by_win_size "
                 "ORDER BY win_size DESC");
  res = PQexec(conn, query);
  num_win_sizes = PQntuples(res);

  if (num_win_sizes == 0) 
    return;

  win_sizes = (int *)my_malloc(num_win_sizes * sizeof(int));
  for (i=0; i<num_win_sizes; ++i) {
    win_sizes[i] = atoi(PQgetvalue(res, i, 0));
  }
  PQclear(res);

  //----------------------------------------------------------------------
  //  For each window size, get the gene neighborhoods
  //----------------------------------------------------------------------
  for (w=0; w<num_win_sizes; ++w) {
    win_size = win_sizes[w];
    sprintf(query,
            "SELECT nh_id FROM nh_by_win_size JOIN gene_nh USING(nh_id) "
            "WHERE chromosome = '%s' AND win_size = %d ORDER BY nh_id",
            chrom_info->chromosome, win_size);
    res = PQexec(conn, query);
    num_sel_nh = PQntuples(res);
    
    //----------------------------------------------------------------------
    //  First iteration: allocate space for ANC values from all permutations
    //  of all neighborhoods of this window size.
    //  Also, allocate an array for pointers to neighborhoods of this size.
    //
    //  The first iteration has the largest space requirements, so these
    //  buffers can be safely reused in subsequent iterations.
    //----------------------------------------------------------------------
    if (w==0) {
      anc_values = (short *)malloc(num_sel_nh * num_perms * sizeof(short));
      sel_nh = (NH_INFO **)malloc(num_sel_nh * sizeof(NH_INFO *));
    }

    //----------------------------------------------------------------------
    //  For all neighborhoods of the current size, assemble:
    //    1) An array of pointers to these neighborhoods.
    //    2) All computed random ANC values combined into an array.
    //----------------------------------------------------------------------
    for (i=0; i<num_sel_nh; ++i) {
      nh_id = atoi(PQgetvalue(res, i, 0));
      sel_nh[i] = &nh_list[nh_id - min_nh_id];
      memcpy(&anc_values[i*num_perms], sel_nh[i]->rand_anc,
             num_perms*sizeof(short));
    }
    PQclear(res);

    save_pvalues(0, win_size, sel_nh, num_sel_nh, num_perms, anc_values);
  }

  free(anc_values);   // free these buffers immediately since they are huge
  free(sel_nh);
}


//-----------------------------------------------------------------------
int compare_nh_gc(const void *v_a, const void *v_b) {
  NH_INFO *a = *(NH_INFO **)v_a;
  NH_INFO *b = *(NH_INFO **)v_b;

  return( a->num_genes == b->num_genes ? 0 :
          (a->num_genes < b->num_genes ? -1 : 1) );
}


//-----------------------------------------------------------------------
//  Compute ANC pvalues according to number of genes.
//-----------------------------------------------------------------------
static void compute_pvalues_by_gene_count(NH_INFO *nh_list, int num_nh,
        int min_nh_id, int num_perms,
        int min_gene_count, int max_gene_count)
{
  char query[200];
  int i, nh_id, rv, num_rows, ng, ng_offset, num_sorted_nh;
  //int min_gene_count, max_gene_count;
  short *anc_values;
  NH_INFO **sorted_nh = (NH_INFO **)malloc(num_nh * sizeof(NH_INFO *));
  NH_INFO *nh;

/*
  res = PQexec(conn,
      "SELECT min_gene_count, max_gene_count FROM nh_by_gc_bounds");
  if (PQntuples(res) == 0) 
    return;

  min_gene_count = atoi(PQgetvalue(res, 0, 0));
  max_gene_count = atoi(PQgetvalue(res, 0, 1));
*/

  //---------------------------------------------------------------------
  //  Assemble all the gene neighborhoods with gene counts in range
  //  into an array.  Sort the array in order of increasing gene counts.
  //---------------------------------------------------------------------
  num_sorted_nh = 0;
  for (i=0, nh=nh_list; i<num_nh;  ++i, ++nh)
    if (min_gene_count <= nh->num_genes && nh->num_genes <= max_gene_count)
      sorted_nh[num_sorted_nh++] = nh;

  qsort(sorted_nh, num_sorted_nh, sizeof(NH_INFO *), compare_nh_gc);

  anc_values = (short *)malloc(num_sorted_nh * num_perms * sizeof(short));

  //---------------------------------------------------------------------
  //  Iterate through the gene neighborhoods.  Every time the gene count
  //  increases, compute and save the pvalues for the prior gene count.
  //---------------------------------------------------------------------
  ng = min_gene_count;
  ng_offset = 0;
  for (i=0; i<num_sorted_nh; ++i) {
    nh = sorted_nh[i];

    if (nh->num_genes != ng) {
      save_pvalues(1, ng, &sorted_nh[ng_offset], i-ng_offset,
          num_perms, anc_values);
      ++ng;
      ng_offset = i;
    }

    memcpy(&anc_values[(i-ng_offset)*num_perms], nh->rand_anc,
           num_perms*sizeof(short));

  }
  save_pvalues(1, ng, &sorted_nh[ng_offset], num_sorted_nh-ng_offset,
               num_perms, anc_values);

  free(anc_values);   // free these buffers immediately since they are huge
  free(sorted_nh);
}


//-----------------------------------------------------------------------
int main(int argc, char **argv)
{
  int num_perms;
  CHROM_INFO *chrom_info;
  NH_INFO *nh_list;
  NH_GROUP *nh_groups;
  int num_nh_groups, num_nh, min_nh_id;
  int min_gc, max_gc;
  char schema[30], cmd[100];

  parse_opts(argc, argv, schema, &num_perms, &min_gc, &max_gc);

  atexit(close_pg);

  conn = my_connect();

  if (schema[0]) {
    sprintf(cmd, "SET search_path TO %s,public", schema);

    res = PQexec(conn, cmd);
    status = PQresultStatus(res);
    if (status != PGRES_COMMAND_OK) {
      printf("Error (SET search_path): %s\n", PQresStatus(status));
      exit(1);
    }
  }

  N = get_number_of_genes();

  read_corr_matrix();

  while (chrom_info = get_next_chromosome()) {

    init_genes_list(N);
    
    get_gene_nh_for_chrom(chrom_info, num_perms,
                          &nh_list, &num_nh, &min_nh_id,
                          &nh_groups, &num_nh_groups);

    compute_anc_values(chrom_info, nh_groups, num_nh_groups, num_perms);


    compute_pvalues_by_win_size(chrom_info, nh_list, num_nh, min_nh_id,
                                  num_perms);

    compute_pvalues_by_gene_count(nh_list, num_nh, min_nh_id, num_perms,
                                  min_gc, max_gc);

    my_free();
  }

}
