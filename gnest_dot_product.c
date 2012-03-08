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

FUNCTION gnest_dot_product(
    p_gene_1_expr_ranks      integer[],
    p_gene_2_expr_ranks      integer[],
    RETURNS integer
*/

#include "postgres.h"
#include "utils/builtins.h"
#include "utils/array.h"
#include "catalog/pg_type.h"


//---------------------------------------------------------------------
#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

PG_FUNCTION_INFO_V1(gnest_dot_product);
Datum gnest_dot_product(PG_FUNCTION_ARGS) {

  ArrayType *v;

  v = PG_GETARG_ARRAYTYPE_P(0);
  int *gene_1_expr_ranks = (int *) ARR_DATA_PTR(v);
  int num_samples = ARR_DIMS(v)[0];
  
  v = PG_GETARG_ARRAYTYPE_P(1);
  int *gene_2_expr_ranks = (int *) ARR_DATA_PTR(v);

  int cum = 0;
  for (int i=0; i<num_samples; ++i)
    cum += gene_1_expr_ranks[i] * gene_2_expr_ranks[i];

  PG_RETURN_INT32(cum);
}
