
perl  \
   gnest.pl \
    --project_taxon_id    10090                           \
    --chromosomes         chromosome_info.csv             \
    --genes               genes.csv                       \
    --samples             samples.csv                     \
    --expr_data           expr_data.csv                   \
    --filter_on_mas5      mas5_calls.txt                  \
    --keep_project --progress  --export_db --corr_matrix  \
    9598  mouseVchimp_min10K_max100K_mark2.bed  \
    9913  mouseVcow_min10K_max100K_mark2.bed    \
    9615  mouseVdog_min10K_max100K_mark2.bed    \
    9606  mouseVhuman_min10K_max100K_mark2.bed  \
    10116 mouseVrat_min10K_max100K_mark2.bed 

#  USE DEFAULTS FOR:
#    --filter_min_expr
#    --min_win_size
#    --max_win_size
#    --min_gene_count
#    --max_gene_count
#    --num_permutations
