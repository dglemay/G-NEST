
perl \
   gnest.pl \
    --project_taxon_id    10090                           \
    --chromosomes         mouse_chroms.csv                \
    --genes               genes.csv                       \
    --samples             samples.csv                     \
    --expr_data           expr_data.csv                   \
    --filter_on_mas5      mas5_calls.txt                  \
    --keep_project --progress  --export_db --corr_matrix  \
    9598  mouse/mouseVchimp.bed  \
    9913  mouse/mouseVcow.bed    \
    9615  mouse/mouseVdog.bed    \
    9606  mouse/mouseVhuman.bed  \
    10116 mouse/mouseVrat.bed

#  USE DEFAULTS FOR:
#    --filter_min_expr
#    --min_win_size
#    --max_win_size
#    --min_gene_count
#    --max_gene_count
#    --num_permutations
