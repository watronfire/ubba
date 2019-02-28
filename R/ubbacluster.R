library( dplyr )
#' @param sce.object The SingleCellExperiment object
#' @param n_dimreds How many PCs to use for further dimensionality reduction. Defauls to all available.
#' @param resolution Value of the resolution parameter for clustering, use a value above 1.0 if you want to obtain a larger number of communities, or vice-versa.
#' @param test Denotes which test to use determine differentially expressed genes. Available options are found in ?Seurat::FindAllMarkers
#' @param features Features to include in featuresplot. Will plot expression across reduced dimensional space.
#' @param plot_reduction Which dimensionality reduction to use for feature plots. Available options are PCA, TSNE, and UMAP.
#' @param verbose Whether to add plots, print progress, etc.
basic_cluster <- function(
    sce.object,
    n_dimreds=ncol( SingleCellExperiment::reducedDim( sce.object, "PCA" ) ),
    resolution=1.0,
    test="wilcox",
    features=c( "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A" ),
    plot_reduction="UMAP",
    verbose=TRUE
) {
    set.seed( 42 )
    cluster_plots <- list()

    ## Clustering -----------------------------------------------------------------------
    se.object <- Seurat::as.Seurat( sce.object )
    se.object <- Seurat::FindNeighbors( object=se.object, reduction="PCA", dims=1:n_dimreds, do.plot=FALSE )
    se.object <- Seurat::FindClusters( object=se.object, resolution=resolution, random.seed=42, verbose=verbose )

    if( verbose ) {
        cluster_plots$basic_clusters <- Seurat::DimPlot( object=se.object, pt.size=1, label=TRUE, reduction=plot_reduction ) +
            ggplot2::ggtitle( label=paste( "Clusters at resolution:", resolution ) )
        cluster_plots$features <- Seurat::FeaturePlot( object=se.object, features=features, reduction=plot_reduction ) +
            ggplot2::ggtitle( label="Common PMBC feature expression" )
    }

    # Assign clusters from Seurat object to SingleCellExperiment object
    SummarizedExperiment::colData( sce.object )$cluster <- se.object@active.ident

    ## Differential Expression ----------------------------------------------------------
    Seurat::Idents( se.object ) <- se.object$RNA_snn_res.1.3
    markers <- Seurat::FindAllMarkers( object=se.object, only.pos=TRUE, verbose=verbose )

    # Could perhaps rename clusters based on most differentially expressed gene. Who knows?

    if( verbose ) {
        top_de_features <- markers %>% group_by( cluster ) %>% top_n( n = 1, wt = avg_logFC )
        cluster_plots$de_feature_plot <- Seurat::FeaturePlot( object=se.object, features=top_de_features$gene, reduction=plot_reduction )

        for( i in seq( 0.1, 2.0, 0.1 ) ) {
            se.object <- Seurat::FindClusters( object=se.object, resolution=i, random.seed=42, verbose=TRUE )
        }

        p1 <- Seurat::DimPlot( se.object, reduction="UMAP", label=TRUE, pt.size=1 ) + ggplot2::theme( legend.position="none" )

        p2 <- se.object@meta.data %>%
            select( starts_with( "RNA_snn_res." ) ) %>%
            clustree::clustree( prefix="RNA_snn_res.", return="plot" ) +
            ggplot2::theme( legend.position="none" )

        cluster_plots$cluster_stability <- scater::multiplot( p1, p2, cols=2 )
    }

    lst( data=sce.object,
         cluster_plots,
         seurat.obj=se.object,
         info=lst( clusters=se.object@active.ident,
                   markers=markers ) )
}

#source( "/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/R/ubbanormaliser.R" )
#source( "/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/R/ubbadimred.R" )

#pbmc8k <- DropletUtils::read10xCounts( "/Users/natem/codeSandbox/tcr/0_data/pbmc8k/GRCh38", type="sparse", col.names=TRUE )
#pbmc8k <- DropletUtils::read10xCounts( "/Users/natem/codeSandbox/tcr/0_data/pbmc8k/GRCh38", type="sparse", col.names=TRUE )
#rownames( pbmc8k ) <- SummarizedExperiment::rowData( pbmc8k )$Symbol
#sce.norm <- normalise_filter_counts( pbmc8k, filter_hvg=FALSE )
#sce.dimred <- dimensionality_reduction( sce.norm$data, limit_hvg=TRUE, num_hvg=2000, n_dimred="auto", verbose=TRUE )
#se.out2 <- basic_cluster( sce.dimred$data, n_dimreds=sce.dimred$info$actual_num_pcs, resolution=0.8, verbose=TRUE )



