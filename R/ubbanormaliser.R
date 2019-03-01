#' @param sce.object The SingleCellExperiment object
#' @param filter_cells Whether the cells have to be filtered
#' @param filter_features Whether the features have to be filtered
#' @param filter_hvg Whether to filter on highly variable features
#' @param normalisation How to normalise
#' @param has_spike Does this contain spike-ins, for which the feature names are preseded by ERCC
#' @param verbose Whether to add plots
#' @param nmads Number of median deviations for filtering outlier cells
#' @param min_expression Minimum library size for a cell
#' @param min_ave_expression Minimal average expression of a feature
#' @param min_cell Minimum cells a gene must be expressed in.
#' @param hvg_fdr FDR feature filtering cutoff
#' @param hvg_bio Biological feature filtering cutoff
#' @param min_variable_fraction Minimal number of variable features to retain
normalise_filter_counts <- function(
    sce.object,
    filter_cells=TRUE,
    filter_features=TRUE,
    filter_hvg=TRUE,
    normalisation="scran_size_factors",
    has_spike=any( grepl( "^ERCC-", rownames( sce.object ) ) ),
    verbose=TRUE,
    nmads=3,
    min_expression=200,
    min_ave_expression=0.005,
    min_cell=0,
    hvg_fdr=0.05,
    hvg_bio=0.5,
    min_variable_fraction=0.15
) {
    if( verbose ) {
        requireNamespace( "grDevices" )
        requireNamespace( "ggplot2" )
        requireNamespace( "graphics" )
        requireNamespace( "KernSmooth" )
    }

    num_genes <- nrow( sce.object )
    countsObj <- SingleCellExperiment::counts( sce.object )
    normalisation_plots <- list()
    normalisation_steps <- tibble::tibble()

    ## Create data object ---------------------------------------------------------------

    # Remove cells which contain less than minimum express and add to metadata. Throw error if no cells
    # remain
    SummarizedExperiment::colData( sce.object )$library.size <- Matrix::colSums( SingleCellExperiment::counts( sce.object ) )
    sce.object <- sce.object[,SummarizedExperiment::colData( sce.object )$library.size > min_expression]
    if( ncol( sce.object ) == 0 ) {
        stop( "No cell contains gene counts!" )
    }

    # Remove genes which are not expressed in minimum number of cells. Throw errors etc.
    keepFeatures <- Matrix::rowSums( SingleCellExperiment::counts( sce.object ) > 0 ) > min_cell
    sce.object <- sce.object[keepFeatures,]

    # Determines percentage reads which are for mitocondrial genes. Add them to metadata
    mito.genes <- grepl( pattern="^(mt|MT|Mt)-", x=SummarizedExperiment::rowData( sce.object )$Symbol )
    percent.mito <- Matrix::colSums( SingleCellExperiment::counts( sce.object )[mito.genes,] ) / Matrix::colSums( SingleCellExperiment::counts( sce.object ) )
    SummarizedExperiment::colData( sce.object )$percent.mito <- percent.mito
    has_mito <- any( mito.genes )

    # Check if spike-ins are present with specified prefix
    spike <- grepl( pattern="^ERCC-", x=SummarizedExperiment::rowData( sce.object )$Symbol )
    has_spike <- any( spike )

    # Specify what has been learned about dataset so far.
    feature_controls <- list()
    if( has_mito ) feature_controls$Mt <- mito.genes
    if( has_spike ) feature_controls$ERCC <- spike

    # Calculate QC metrics based on what has been determined so far.
    sce.object <- scater::calculateQCMetrics( sce.object, feature_controls=feature_controls, compact=TRUE )

    if( has_spike ) SingleCellExperiment::isSpike( sce.object, "ERCC" ) <- spike

    if( verbose ) {
        normalisation_steps <- tibble::tribble(
            ~type, ~nfeatures, ~ncells,
            "original", nrow( sce.object ), ncol( sce.object )
        )
        print( glue::glue( "Original: features - {nrow(sce.object)} Cells - {ncol(sce.object)}" ) )

        normalisation_plots$library <- ggplot2::ggplot( as.data.frame( sce.object$scater_qc$all ) ) +
            ggplot2::geom_histogram( ggplot2::aes( total_counts ) ) +
            ggplot2::scale_x_continuous( limits=c( 0, NA ) )

        if( has_spike ) {
            normalisation_plots$spike <- ggplot2::ggplot( as.data.frame( sce.object$scater_qz$feature_control_ERCC ) ) +
                ggplot2::geom_histogram( ggplot2::aes( pct_counts ) ) +
                ggplot2::scale_x_continuous( limits=c( 0, 100 ) )
        }

        if( has_mito ) {
            normalisation_plots$mito <-
                ggplot2::ggplot( as.data.frame( sce.object$scater_qc$feature_control_Mt ) ) +
                ggplot2::geom_histogram( ggplot2::aes( pct_counts ) ) +
                ggplot2::scale_x_continuous( limits = c( 0, 100 ) )
        }


    }

    ## Filter cells ---------------------------------------------------------------------

    if( filter_cells ) {

        # Initialize values
        total_counts <- sce.object$scater_qc$all$log10_total_counts
        total_features <- sce.object$scater_qc$all$log10_total_features_by_counts
        pct_counts_Mt <- sce.object$scater_qc$feature_control_Mt$pct_counts
        pct_counts_ERCC <- sce.object$scater_qc$feature_control_ERCC$pct_counts
        mito_drop <- rep( FALSE, length( total_counts ) )
        spike_drop <- rep( FALSE, length( total_counts ) )

        # Drop cells based if they are outliers in either library size, features, spike_in, mitocondrial genes.
        # We really only care about smaller cells though. Large outliers may be
        # biologically revelent or doublets
        libsize_drop <- scater::isOutlier( total_counts, nmads=nmads, type="both", log=TRUE )
        feature_drop <- scater::isOutlier( total_features, nmads=nmads, type="both", log=TRUE )
        if( has_mito ) mito_drop <- scater::isOutlier( pct_counts_Mt, nmads=nmads, type="higher" )
        if( has_spike ) spike_drop <- scater::isOutlier( pct_counts_ERCC, nmads=nmads, type="higher" )

        if( verbose ) {
            tibble::tibble( sum( mito_drop ), sum( spike_drop ), sum( libsize_drop ), sum( feature_drop ) ) %>% print()
        }

        # Subset cells based on factors determined
        sce.object <- sce.object[,!( libsize_drop | feature_drop | mito_drop | spike_drop )]

        if( verbose ) {
            normalisation_steps <- normalisation_steps %>%
                add_row( type="cell_quality_filtering", nfeatures=nrow( sce.object ), ncells=ncol( sce.object ) )
            print( glue::glue( "Cell filter: features - {nrow(sce.object)} Cells - {ncol(sce.object)}" ) )
        }
    }

    ## Filter features ----------------------------------------------------------------------
    if( filter_features ) {
        ave_counts <- Matrix::rowMeans( SingleCellExperiment::counts( sce.object ) )
        keep <- ave_counts >= min_ave_expression

        sce.object <- sce.object[keep,]

        # Get rid of duplicated features because we're working with symbols. Ideally, we'd
        # rename but I don't know how to do that.
        sce.object <- sce.object[ !duplicated( rownames( sce.object ) ),]

        if( verbose ) {
            normalisation_plots$ave_counts <- tibble::tibble( ave_counts=ave_counts ) %>%
                ggplot2::ggplot() +
                ggplot2::geom_histogram( ggplot2::aes( ave_counts ) ) +
                ggplot2::scale_x_log10() +
                ggplot2::geom_vline( xintercept=min_ave_expression )

            top_features <- ave_counts %>% sort() %>% tail( 20 ) %>% names()

            counts_top_features <- as.matrix( SingleCellExperiment::counts( sce.object[top_features] ) ) %>%
                reshape2::melt( varnames=c( "feature", "cell" ), value.name="count") %>%
                dplyr::mutate( feature=factor( feature, levels=top_features ) )

            avecounts_top_features <- counts_top_features %>%
                group_by( feature ) %>%
                summarise( count=mean( count ) )

            normalisation_plots$top_counts <- counts_top_features %>%
                ggplot2::ggplot( ggplot2::aes( feature, count + 1 ) ) +
                ggplot2::geom_point( shape = "|" ) +
                ggplot2::geom_point( data=avecounts_top_features, color="blue", size=4 ) +
                ggplot2::scale_y_log10() +
                ggplot2::coord_flip()

            normalisation_steps <- normalisation_steps %>%
                add_row( type="feature_expression_filtering", nfeatures=nrow( sce.object ), ncells=ncol( sce.object ) )
            print( glue::glue( "Feature filter: features - {nrow(sce.object)} Cells - {ncol(sce.object)}" ) )
        }
    }

    ## Normalise ----------------------------------------------------------------------------
    if( normalisation == "scran_size_factors" ) {
        clusters <- scran::quickCluster( sce.object, min.size=100, assay.type="counts" )
        sce.object <- scran::computeSumFactors( sce.object, clusters=clusters, min.mean=0.1, positive=TRUE )
        sce.object <- sce.object[,SingleCellExperiment::sizeFactors( sce.object ) > 0]

        if( has_spike ) {
            sce.object <- scran::computeSpikeFactors( sce, type="ERCC", general.use=FALSE )
            if( any( is.na( SingleCellExperiment::sizeFactors( sce.object, type="ERCC" ) ) ) ) {
                warning( "Some cells do not have any spike-ins, this will cause an error later is analysis. Remove spike-ins." )
            }
        }
        #Throw a warning if there is low correlation between sumFactors and library sizes, falling back to library sizes
        if( cor( sce.object$scater_qc$all$total_counts, SingleCellExperiment::sizeFactors( sce.object ) ) < 0.5 ) {
            warning( "Low correlation between sumFactors and library sizes, falling back to library sizes" )
            SingleCellExperiment::sizeFactors( sce.object ) <- scater::librarySizeFactors( sce.object )
        }

        sce.object <- scater::normalize( sce.object )

        if( verbose ) {
            normalisation_plots$size_factors <- tibble::tibble(
                size_factors=SingleCellExperiment::sizeFactors( sce.object ),
                library_size=sce.object$scater_qc$all$total_counts
            ) %>%
                ggplot2::ggplot( ggplot2::aes( size_factors, library_size ) ) +
                ggplot2::geom_point() + ggplot2::geom_smooth( method="lm" )
        }
    } else {
        stop( "Normalisation not supported" )
    }
    if( verbose ) {
        normalisation_steps <- normalisation_steps %>%
            add_row( type="normalisation", nfeatures=nrow( sce.object ), ncells=ncol( sce.object ) )
        print(glue::glue("Normalised: features - {nrow(sce.object)} Cells - {ncol(sce.object)}"))
    }
    ## Select highly variable features ------------------------------------------------------
    if( filter_hvg ) {
        var_fit <- scran::trendVar( sce.object, method="spline", use.spikes=has_spike )
        var_out <- scran::decomposeVar( sce.object, var_fit ) %>% as.data.frame()

        if( verbose ) {
            if( has_spike ) {
                var_out$spike <- SingleCellExperiment::isSpike( sce.object )
            } else {
                var_out$spike <- FALSE
            }

            normalisation_plots$variance_mean <- var_out %>%
                ggplot2::ggplot( ggplot2::aes( mean, total ) ) +
                ggplot2::geom_point( ggplot2::aes( color=spike ) ) +
                ggplot2::geom_smooth()
            normalisation_plots$fdr_bio <- var_out %>%
                ggplot2::ggplot( ggplot2::aes( FDR, bio ) ) +
                ggplot2::geom_point() +
                ggplot2::geom_hline( yintercept=hvg_bio, color="red" ) +
                ggplot2::geom_vline( xintercept=hvg_fdr, color="red" )
        }

        var_out <- var_out[order( var_out$bio, decreasing=TRUE ),]
        hvg_out <- var_out[which( var_out$FDR <= hvg_fdr & var_out$bio >= hvg_bio ),]
        if( nrow( hvg_out ) < min_variable_fraction * num_genes ) {
            n_features <- min( nrow( var_out ), ceiling( min_variable_fraction * num_genes ) )

            hvg_out <- var_out[seq( 1, n_features ),]
        }

        sce.object <- sce.object[rownames( hvg_out ),]

        if( verbose ) {
            normalisation_steps <- normalisation_steps %>%
                add_row( type="feature_variability_filtering", nfeatures=nrow( sce.object ), ncells=ncol( sce.object ) )
            print( glue::glue( "Variable features filtered: features - {nrow(sce.object)} Cells - {ncol(sce.object)}" ) )
        }
    }

    expr_norm_filt <- SingleCellExperiment::logcounts( sce.object ) %>% Matrix::t()
    count_filt <- SingleCellExperiment::counts( sce.object ) %>% Matrix::t()

    ## Iterative filtering on variability ---------------------------------------------------
    repeat {
        feature_sds <- SingleCellExperiment::counts( sce.object ) %>% apply( 1, stats::sd )
        cell_sds <- SingleCellExperiment::counts( sce.object ) %>% apply( 2, stats::sd )

        features_filtered <- which( feature_sds > 0, useNames=TRUE )
        cells_filtered <- which( cell_sds > 0, useNames=TRUE )
        sce.object <- sce.object[features_filtered, cells_filtered]

        if( ( min( c( feature_sds, 1 ), na.rm=TRUE ) > 0 ) && ( min( c( cell_sds, 1 ), na.rm=TRUE ) > 0 ) ) {
            break
        }
    }

    if( verbose ) {
        normalisation_steps <- normalisation_steps %>%
            add_row( type = "final_filtering", nfeatures = nrow( sce.object ), ncells = ncol( sce.object ) )
        print( glue::glue( "Final filtering: features - {nrow(sce.object)} Cells - {ncol(sce.object)}" ) )
    }

    ## Output -------------------------------------------------------------------------------
    if(verbose) {
        type <- NULL # satisfy r cmd check

        normalisation_plots$n_retained <-
            normalisation_steps %>%
            mutate( type = factor( type, levels = rev( type ) ) ) %>%
            tidyr::gather( "dimension", "n", -type ) %>%
            ggplot2::ggplot() +
            ggplot2::geom_bar( ggplot2::aes_string( "type", "n", fill="dimension" ), position="dodge", stat="identity" ) +
            ggplot2::facet_wrap( ~dimension, scales="free_x" ) +
            ggplot2::coord_flip()
    } else {
        normalisation_steps <- NULL
    }

    lst( data=sce.object,
         normalisation_plots,
         info=lst( has_spike,
                   has_mito,
                   normalisation_steps ) )
}
