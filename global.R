# library(gdata)
# library(tools)



run_experiment <- function(cell, eln, program, centering_fun){
    tga_experiment_path <- glue::glue("{out_file_root}{program}/TGA/EdU/Experiments/RawData/{cell}/{eln}/3D/")
    results <- run_TGA_analysis(tga_experiment_path = tga_experiment_path,
                            centering_fun = tolower(centering_fun)) # method_init
    return(results)
}

load_seurat_obj <- function(path)

validate_seurat_obj <- function(seurat_obj){
}

get_metadata_cols <- function(seurat_obj){
}

create_meta_UMAP_ggplot <- function(seurat_obj){
}

create_gene_expression_UMAP_ggplot <- function(seurat_obj){
}
