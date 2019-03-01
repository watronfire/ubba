#' @param sce.object The SingleCellExperiment object
#' @param on_hvg Whether to perform dimensionality reduction on highly variable features or all features
#' @param has_spike Does this contain spike-ins, for which the feature names are preseded by ERCC
#' @param limit_hvg Whether to limit number of highly variable features for dimensionality reduction
#' @param num_hvg Number of highly variable features to use, if limiting
#' @param n_dimred How many PCs to use for further dimensionality reduction. Defauls to auto which performs Horn's parallel analysis to choose the number of PCs.
#' @param verbose Whether to add plots
dimensionality_reduction <- function(
    sce.object,
    on_hvg=TRUE,
    has_spike=any( grepl( "^ERCC-", rownames( sce.object ) ) ),
    limit_hvg=FALSE,
    num_hvg=1000,
    n_dimred="auto",
    verbose=TRUE
) {
    set.seed( 42 )
    dimred_plots <- list()

    ## Principle Component Analysis -----------------------------------------------------
    if( on_hvg ) {
        var.fit <- scran::trendVar( sce.object, method="loess", use.spikes=has_spike )
        var.out <- scran::decomposeVar( sce.object, var.fit )
        if( limit_hvg ) {
            top.hvg <- head( order( var.out$bio, decreasing=TRUE ), num_hvg )
        } else {
            top.hvg <- order( var.out$bio, decreasing=TRUE )
        }
    } else {
        top.hvg <- rownames( sce.object$data )
    }

    # Pretty standard, eh?
    sce.object <- scater::runPCA( sce.object, ncomponents=50, feature_set=top.hvg, set.seed=42, method="irlba" )

    if( verbose ) {
        dimred_plots$pca <- scater::plotPCA( sce.object, add_ticks=FALSE ) + ggplot2::ggtitle( label="PCA Normalized Gene Counts" )
    }

    ## Calculate Dimensionality ---------------------------------------------------------
    # SLOW
    if( n_dimred=="auto" ) {
        if( verbose ) {
            print( "Determining dimensionality of dataset using Horn’s parallel analysis (Horn 1965). May take a good long while." )
        }
        npcs <- scran::parallelPCA( sce.object, value="n", subset.row=top.hvg, approximate=TRUE )
        actual_num_pcs <- as.integer( npcs )
    } else {
        actual_num_pcs <- as.integer( n_dimred )
    }

    if( verbose ){
        pcs_labels <- colnames( SingleCellExperiment::reducedDim( sce.object, "PCA" ) )
        if( n_dimred=="auto" ) {
            dimred_plots$dimensionality <- ggplot2::ggplot() +
                ggplot2::geom_point( ggplot2::aes( x=1:length( pcs_labels ), y=attr( npcs, "percentVar" )[1:length( pcs_labels )] ) ) +
                ggplot2::geom_vline( xintercept=actual_num_pcs+0.5 ) +
                ggplot2::labs( title="Variance Explained per PC", x="PC #", y="Variance Explained" )
        } else {
            dimred_plots$dimensionality <- ggplot2::ggplot() +
                ggplot2::geom_point( ggplot2::aes( x=1:length( pcs_labels ), y=attr( SingleCellExperiment::reducedDim( sce.object, "PCE" ), "percentVar" )[1:length( pcs_labels )] ) ) +
                ggplot2::geom_vline( xintercept=actual_num_pcs+0.5 ) +
                ggplot2::labs( title="Variance Explained per PC", x="PC #", y="Variance Explained" )
        }
        print( glue::glue( "Assuming dimensionality of dataset to be: {actual_num_pcs}" ) )
    }

    ## t-Distributed Stochastic Neighbor Embedding --------------------------------------
    sce.object <- scater::runTSNE( sce.object, use_dimred="PCA", n_dimred=actual_num_pcs, rand_seed=42 )

    if( verbose ) {
        dimred_plots$tsne <- scater::plotTSNE( sce.object, add_ticks=FALSE ) + ggplot2::ggtitle( label="tSNE Normalized Gene Counts" )
    }

    ## Uniform Manifold Approximation and Projection ------------------------------------
    sce.object <- scater::runUMAP( sce.object, use_dimred="PCA", n_dimred=actual_num_pcs, set.seed=42 )
    if( verbose ) {
        dimred_plots$umap <- scater::plotUMAP( sce.object, add_ticks=FALSE ) + ggplot2::ggtitle( label="UMAP Normalized Gene Counts" )
    }

    lst( data=sce.object, dimred_plots, info=lst( actual_num_pcs ) )

}
