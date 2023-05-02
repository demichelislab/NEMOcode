#' hg38tohg19 chain file
#'
#' Chain file to lift positions from hg38 to hg19. It can be used with rtracklayer::liftover function
"chainhg38_hg19"


#' contigs_hg38 chain file
#'
#' Chromosome names in hg38 contig format
"contigs_hg38"


#' panel_design chain file
#'
#' A data.frame with the NEMO panel design, annotating regions subsets, methylation states, genes, and positions
"panel_design"


#' refDist chain file
#'
#' Reference distributions for cell types of interest: ad:Adeno, ne:NE, adm:PBMC, cf:cfDNA
"refDist"


#' control_mat_cfDNA matrix of cfDNA samples
#'
#' Reference samples for tumor content estimation (from Fox-Fisher et al. eLife2021)
"control_mat"


#' control_mat_PBMC matrix of PBMC samples
#'
#' Reference samples for tumor content estimation in tissue (from the BLUEPRINT consortium)
"control_mat_pbmc"
