install.packages( "gam" )

source( "/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/R/ubbanormaliser.R" )
source( "/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/R/ubbadimred.R" )
source( "/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/R/ubbacluster.R" )

pbmc8k <- DropletUtils::read10xCounts( "/Users/natem/codeSandbox/tcr/0_data/pbmc8k/GRCh38", type="sparse", col.names=TRUE )
rownames( pbmc8k ) <- SummarizedExperiment::rowData( pbmc8k )$Symbol
sce.norm <- normalise_filter_counts( pbmc8k, filter_hvg=FALSE )
sce.dimred <- dimensionality_reduction( sce.norm$data, limit_hvg=TRUE, num_hvg=2000, n_dimred="auto", verbose=TRUE )
sce.clust <- basic_cluster( sce.dimred$data, n_dimreds=sce.dimred$info$actual_num_pcs, resolution=0.6, verbose=TRUE )

sce.traj <- slingshot::slingshot( data=sce.clust$data, clusterLabels="cluster", reducedDim="PCA" )
saveRDS( sce.traj, "/Users/natem/codeSandbox/sce.traj.RData" )

summary( sce.traj$slingPseudotime_1 )

colors <- colorRampPalette( RColorBrewer::brewer.pal( 11, "Spectral" )[-6])(100)
plot( SingleCellExperiment::reducedDim( sce.traj, "PCA" ), col=colors[cut( sce.traj$slingPseudotime_1, breaks=100)],
      pch=16, asp=1 )
lines( slingshot::SlingshotDataSet( sce.traj ), lwd=2 )

t <- sce.traj$slingPseudotime_1
Y <- SummarizedExperiment::assay( sce.traj, "logcounts" )

var.fit <- scran::trendVar( sce.traj, method="loess", use.spikes=FALSE )
var.out <- scran::decomposeVar( sce.traj, var.fit )
top.hvg <- head( order( var.out$bio, decreasing=TRUE ), 1000 )
Y <- Y[top.hvg,]

gam.pval <- apply( Y, 1, function( z ) {
    d <- data.frame( z=z, t=t )
    tmp <- gam::gam( z ~ gam::lo( t ), data=d )
    p <- summary( tmp )[4][[1]][1,5]
    p
} )

topgenes <- names( sort( gam.pval, decreasing=FALSE ))[1:100]
heatdata <- SingleCellExperiment::logcounts( sce.traj )[rownames( SingleCellExperiment::logcounts( sce.traj ) ) %in% topgenes, order( t, na.last=NA )]
heatclus <- sce.traj$cluster[order( t, na.last=NA )]

ce <- clusterExperiment::ClusterExperiment( as.matrix( heatdata ), heatclus, transformation=log1p )

png( filename="/Users/natem/Dropbox (Scripps Research)/Personal/Code/R/TCR/ubba/test/heatmap.png", width=2000, height=2000 )
clusterExperiment::plotHeatmap( ce, clusterSamplesData="orderSamplesValue", visualizeData="transformed" )
dev.off()
