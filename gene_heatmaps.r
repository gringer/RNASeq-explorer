#!/usr/bin/Rscript

source("snif_gene_list_heat_map.r");

MIMR.red <- "#DF2726";
MIMR.orange <- "#EC772D";

library(gdata);
genes.df <- read.xls("short lists for heatmaps Nb_DBP-FITC.xlsx", sheet=2);
geneNames <- genes.df$mgi_symbol;
gap.poss <- which(grepl("^:",geneNames));
genes.df$heatmapID <- rep(1:length(gap.poss),diff(c(gap.poss,nrow(genes.df)+1)));

for(id in unique(genes.df$heatmapID)){
    heatmap.in.mat <- as.matrix(genes.df[genes.df$heatmapID==id,grep("log2FoldChange",colnames(genes.df))]);
    heatmap.in.pval <- as.matrix(genes.df[genes.df$heatmapID==id,grep("padj",colnames(genes.df))]);
    rownames(heatmap.in.mat) <- sub("^:","",geneNames)[genes.df$heatmapID==id];
    setName <- rownames(heatmap.in.mat)[1];
    heatmap.in.mat <- heatmap.in.mat[-1,];
    heatmap.in.pval <- heatmap.in.pval[-1,];
    ## replace column names with gene set names
    colnames(heatmap.in.mat) <- sub("C_D.*$","DBP",sub(".C_N.*$",".Nb",sub("_clean","",sub("^.*?_","",colnames(heatmap.in.mat)))));
    ## display heat map
                                        #zcols <- colorRampPalette(c("cyan",rgb(0,0.5,1),"blue","black","red",rgb(1,0,0.5),"magenta"))(64);
    zcols <- colorRampPalette(c("darkblue","chartreuse","yellow",MIMR.orange,MIMR.red))(64);
    ##zlim <- c(median(dds.vstpk.mat)-max(dds.vstpk.mat),max(dds.vstpk.mat));
    zlim <- c(-1,1);
    heat.clipped.mat <- heatmap.in.mat;
    heat.clipped.mat[heat.clipped.mat > zlim[2]] <- zlim[2];
    heat.clipped.mat[heat.clipped.mat < zlim[1]] <- zlim[1];
    heat.clipped.mat <- heat.clipped.mat[nrow(heat.clipped.mat):1,];
    heatmap.in.pval <- heatmap.in.pval[nrow(heatmap.in.pval):1,];
    gringene.heat.map(heat.clipped.mat, pmat = heatmap.in.pval,
                      zcols=zcols, zlim=zlim, filename=sprintf("HeatMap_%s.pdf",gsub(" ","_",setName)), longestLabel="CD11b.DBP");
    gringene.heat.map(heat.clipped.mat,
                      zcols=zcols, zlim=zlim, filename=sprintf("HeatMap_noSig_%s.pdf",gsub(" ","_",setName)), longestLabel="CD11b.DBP");
}

heatmap.in.mat <- as.matrix(genes.df[,grep("log2FoldChange",colnames(genes.df))]);
heatmap.in.pval <- as.matrix(genes.df[,grep("padj",colnames(genes.df))]);
rownames(heatmap.in.mat) <- sub("^:","",geneNames);
rownames(heatmap.in.pval) <- sub("^:","",geneNames);
setName <- rownames(heatmap.in.mat)[1];
heatmap.in.mat <- heatmap.in.mat[-1,];
heatmap.in.pval <- heatmap.in.pval[-1,];
## replace column names with gene set names
colnames(heatmap.in.mat) <- sub("C_D.*$","DBP",sub(".C_N.*$",".Nb",sub("_clean","",sub("^.*?_","",colnames(heatmap.in.mat)))));
## display heat map
##zcols <- colorRampPalette(c("cyan",rgb(0,0.5,1),"blue","black","red",rgb(1,0,0.5),"magenta"))(64);
##zcols <- colorRampPalette(c("darkblue","yellow","red"))(64);
zcols <- colorRampPalette(c("darkblue","chartreuse","yellow",MIMR.orange,MIMR.red))(64);
##zlim <- c(median(dds.vstpk.mat)-max(dds.vstpk.mat),max(dds.vstpk.mat));
zlim <- c(-1,1);
heat.clipped.mat <- heatmap.in.mat;
heat.clipped.mat[heat.clipped.mat > zlim[2]] <- zlim[2];
heat.clipped.mat[heat.clipped.mat < zlim[1]] <- zlim[1];
heat.clipped.mat <- heat.clipped.mat[nrow(heat.clipped.mat):1,];
heatmap.in.pval <- heatmap.in.pval[nrow(heatmap.in.pval):1,];
gringene.heat.map(heat.clipped.mat, pmat = heatmap.in.pval,
                  zcols=zcols, zlim=zlim, filename="HeatMap_All.pdf", longestLabel="CD11b.DBP");
gringene.heat.map(heat.clipped.mat,
                  zcols=zcols, zlim=zlim, filename="HeatMap_noSig_All.pdf", longestLabel="CD11b.DBP");


pdf(file="Key_HeatMap_All.pdf",
    width=1.25, height=1.25, pointsize=8)
snif.heat.map.key(heatmap.in.mat, zlim, zcols, xlab="");
dev.off()

