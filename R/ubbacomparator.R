
measure_batch_effect <- function(
    sce.object,
    batch_label,
    n_dimred=ncol( SingleCellExperiment::reducedDim( sce.object, "PCA" ) ),


) {
    ## kBET analysis ----------------------------------------------------------
    kBET::kBET( SingleCellExperiment::reducedDim( sce.object, "PCA" ), batch=SingleCellExperiment::colData( sce.object )

    ## Silhouette width -------------------------------------------------------
    kBET::batch_sil()

    ## Principcal Componenent Regression --------------------------------------
    kBET::pcRegression()
}

source( "/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/R/ubbanormaliser.R" )
source( "/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/R/ubbadimred.R" )
source( "/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/R/ubbacluster.R" )

pbmc8k <- DropletUtils::read10xCounts( "/Users/natem/codeSandbox/tcr/0_data/pbmc8k/GRCh38", type="sparse", col.names=TRUE )
pbmc4k <- DropletUtils::read10xCounts( "/Users/natem/codeSandbox/tcr/0_data/pbmc4k/GRCh38", type="sparse", col.names=TRUE )
rownames( pbmc8k ) <- SummarizedExperiment::rowData( pbmc8k )$Symbol
rownames( pbmc4k ) <- SummarizedExperiment::rowData( pbmc4k )$Symbol
SummarizedExperiment::colData( pbmc8k )$Batch <- "8k"
SummarizedExperiment::colData( pbmc4k )$Batch <- "4k"
colnames( pbmc8k ) <- paste( colnames( pbmc8k ), SummarizedExperiment::colData( pbmc8k )$Batch, sep="_" )
colnames( pbmc4k ) <- paste( colnames( pbmc4k ), SummarizedExperiment::colData( pbmc4k )$Batch, sep="_" )

combined <- SummarizedExperiment::cbind( pbmc8k, pbmc4k )

comb.norm <- normalise_filter_counts( combined, filter_hvg=FALSE, has_spike=FALSE )
comb.dimred <- dimensionality_reduction( comb.norm$data, limit_hvg=TRUE, num_hvg=2000 )
comb.cluster <- basic_cluster( comb.dimred$data, resolution=1.3, n_dimreds=comb.dimred$info$actual_num_pcs )
length( colnames( comb.dimred$data ) )
colnames( comb.dimred$data )[duplicated( colnames( comb.dimred$data ) )]

