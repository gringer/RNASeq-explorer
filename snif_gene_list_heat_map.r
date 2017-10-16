 #
#
#

#library(gplots)
#source(file.path(".","_DESeq2_Support","snif_DESeq2_PCA.R")) # currently, to get to.figure.file()



# function to map in.v values (assumed to come from range [in.v.lower;in.v.upper])
# to out.v values in [out.v.lower; out.v.upper]
snif.scale <- function(
  in.v,  
  in.v.lower=min(in.v),  in.v.upper=max(in.v),
  out.v.lower=0,         out.v.upper=1
  ) { 
  return( (out.v.upper-out.v.lower) * (in.v-in.v.lower)/(in.v.upper-in.v.lower) + out.v.lower) 
}


# function to plot (independently) the heatmap key
snif.heat.map.key <- function( z, zlim=range(z, na.rm=TRUE), zcols, xlab="", ... ) {

  # prep
  nbin <- length(zcols)
  
  dlim <- range(z, na.rm=TRUE);
  z.key <- seq(from=zlim[1], to=zlim[2], length.out=nbin)
  plot.range <- c(min(dlim,zlim),
                  max(dlim,zlim));
  d <- density(z[(z>=plot.range[1]) & (z<=plot.range[2])], na.rm=TRUE);
  plot.range <- plot.range + diff(plot.range) * c(-0.1,0.1);
  #cat(sprintf("dlim: %s, zlim: %s\n",paste(dlim,collapse=";"),paste(zlim, collapse=";")));
  
  # background
  par( mar=c(2, 1.25, 0.1, 0.1), mgp=c(1.5,0.5,0));
  
  image( x=z.key, y=0, z=matrix(z.key, ncol=1), col=zcols, xlim=plot.range, zlim=zlim, yaxt="n", bty="n", 
         ylab="", xlab="", cex.axis = 0.6, useRaster=TRUE, ...);
  if(max(plot.range) > max(zlim)){
    rect(xleft=max(zlim),xright=max(plot.range),ybottom=-1.5,ytop=1.5,col=tail(zcols,1), border=NA);
  }
  if(min(plot.range) < min(zlim)){
    rect(xleft=min(plot.range),xright=min(zlim),ybottom=-1.5,ytop=1.5,col=head(zcols,1), border=NA);
  }
  
#  pz <- pretty(seq(from=plot.range[1]*diff(dlim),to=plot.range[2]*diff(dlim),length.out=nbin),10); 
#  axis(1, at=snif.scale(pz, min(dlim), max(dlim), plot.range[1], plot.range[2]), labels=pz )
  title(xlab=xlab);
  
  # density
  par( mgp=c(1.25,0.25,0));
 
  yt <- pretty(d$y) ; yt <- yt[yt>0]; 
  title(ylab="Density", line=0.25);
  
  dx <- d$x
  dy <- snif.scale(d$y, 0, 1.01*max(d$y), -1, 1)
  lines(dx, dy, lwd=3, col="#80808080")
  lines(dx, dy, lwd=1, col="black")
  
  invisible(d)
  
}

signif.cols <- function(fcVals, cols, zlim, ...){
  extArgs <- list(...);
  fcBreaks <- 
    if("breaks" %in% names(extArgs)){
      extArgs$breaks;
    } else {
      seq(zlim[1],zlim[2],length.out=length(zcols));
    }
  fcCols <- cols[cut(fcVals, breaks = fcBreaks, include.lowest=TRUE)];
  ## Colour weighting taken from Darel Rex's colour function
  ## http://alienryderflex.com/saturation.html
  fcGreys <- colSums(col2rgb(fcCols)/255 * c(0.299,0.587,0.114));
  return(ifelse(fcGreys < 0.5, "#D0D0D0", "#202020"));
}

signif.stars <- function(pvals){
    res <- matrix("",nrow=nrow(pvals), ncol=ncol(pvals));
    res[pvals < 0.05] <- "*";
    res[pvals < 0.01] <- "**";
    res[pvals < 0.001] <- "***";
    return(res);
}

gringene.heat.map <- function( 
  data.mat,
  pmat=NULL,
  label.pad.nchar=1,
  filename=NULL, pdf.pointsize=8,
  zlim=NULL,
  zcols=NULL,
  order.by="hc",
  longestLabel="***", post.plot=NULL, ...) {
  
  ### Derived from Alex's snif.gene.list.heat.map, but made more 'stupid'
  
  ### data.mat: matrix of data to display
  ###   it is assumed that this has been previously cleaned for diplay
  ###   [column names as display names]
  ### pmat: matrix of [adjusted] probability scores
  
  if(class(data.mat) != "matrix"){
    stop("Report must be a matrix");
  }
  
  ## reverse matrix Y order to fit with expected display
  data.mat <- data.mat[nrow(data.mat):1,,drop=FALSE];
  if(!is.null(pmat)){
    pmat <- pmat[nrow(pmat):1,,drop=FALSE];
  }
  
  # pad dimnames of mat
  label.pad.str <- paste0(rep(" ", label.pad.nchar), collapse="")
  rownames(data.mat) <- paste0(label.pad.str, rownames(data.mat), label.pad.str, " ") # + extra space to accommodate italics
  colnames(data.mat) <- paste0(label.pad.str, colnames(data.mat), label.pad.str)
  
  
  # make a dummy pdf to get font sizes
  if (!is.null( filename)) 
    pdf( file      = file.path(tempdir(),".dummy.pdf"), pointsize = pdf.pointsize )  # dummy PDF to get font size
  
  max.font.height  <- strheight( "M", units="inches")
  max.font.width   <- strwidth(  "M", units="inches")
  col.label.height <- max( strwidth(colnames(data.mat), units="inches") )
  col.width        <- strwidth(sprintf("XX%sXX",longestLabel), units="inches") # want to fit in anything from "" to "***" with 2 chars on either side
  row.label.width  <- max( strwidth(rownames(data.mat) , units="inches"))

  if (!is.null( filename))
    dev.off()
  # graphics setup
  op <- par(no.readonly=TRUE)

  if (!is.null( filename)) {
    if (length(grep( paste0(".","pdf","$"), filename))<=0) {
      filename <- paste0(filename,".","pdf")
    }
    pdf( file      = filename,
         width     = col.width       * ncol(data.mat)     + row.label.width,
         height    = max.font.height * (nrow(data.mat) * 2 + 4) + col.label.height,  # 2: relative size of row vs font height
         pointsize = pdf.pointsize )
    #cat("Creating:",filename,"\n")
    on.exit(dev.off())
  }
  par( mai=c( col.label.height , 0, max.font.height*4, row.label.width ) , mgp=c(0,0, 0));
  ## plot!
  if (is.null(zlim)){
    warning("zlim not specified, using data range");
    zlim <- range(c(data.mat), na.rm=TRUE);
  }
  data.mat[data.mat < zlim[1]] <- zlim[1];
  data.mat[data.mat > zlim[2]] <- zlim[2];
  if (is.null(zcols)){
    zcols <- heat.colors(64);
  }
  # Plot the heatmap, if breaks are declared - use them! (Sam)
  image( t(data.mat), xaxt="n", yaxt="n", bty='n',useRaster=FALSE, zlim=zlim, col=zcols, ...);
  #cat("Plotting...");
  
  #cat(" done!\n");
  xPoss <- (c(col(data.mat))-1)/(ncol(data.mat)-1);
  if(ncol(data.mat) == 1){
    xPoss <- 0;
  }
  yPoss <- (c(row(data.mat))-1)/(nrow(data.mat)-1);
  if(nrow(data.mat) == 1){
    yPoss <- 0;
  }
  if(!is.null(pmat) && all(dim(pmat) == dim(data.mat))){
    text(x = xPoss, y = yPoss,
         labels = signif.stars(pmat),
         col = signif.cols(data.mat, zlim=zlim, cols=zcols, ...));
  } else if(!is.null(pmat)){
    stop("Probability matrix has the wrong size");
  }
  lasType <- ifelse(!is.null(filename),1,3);
  axAts <- if(ncol(data.mat) == 1){
    0;
  } else {
      (0:(ncol(data.mat)-1))/(ncol(data.mat)-1);
  }
  ayAts <- if(nrow(data.mat) == 1){
    0;
  } else {
    (0:(nrow(data.mat)-1))/(nrow(data.mat)-1);
  }
  axis(1, at=axAts, labels = colnames(data.mat), tick=FALSE, las=lasType);
  axis(4, at=ayAts, labels = rownames(data.mat), tick=FALSE, las=1);
  if(!is.null(post.plot)){
    do.call(post.plot, args=list());
  }
  if (!is.null( filename))
    dev.off();
  par(op);
}




# wrapper for the previous
snif.make.gene.list.heat.maps <- function( all.gene.sets, report, comp.labels, dest.dir=tempdir(), ...) {
  
  z <- c(as.matrix(report[,paste0(names(comp.labels),".", "log2FoldChange")]))
  z <- z[!is.na(z)]
  zlim <- max(abs(z))
  zlim <- c(-1*zlim, zlim)
  zcols <- colorRampPalette(c("cyan",rgb(0,0.5,1),"blue","black","red",rgb(1,0,0.5),"magenta"))(64)
  
  # - key
  pdf(file=file.path(dest.dir, paste0("GeneSet_L2FC_Heatmap_KEY", ".pdf" )),
      width=1.25, height=1.25, pointsize=8)
  snif.heat.map.key( z, zlim, zcols, xlab="all L2FCs", ...) 
  dev.off()
  
  # - 1 heatmap / gene set
  
  for (gene.set.name in unique(all.gene.sets$fullName)) {
    
    gene.set <-  unique(subset(all.gene.sets, (fullName == gene.set.name) & 
                                 (ensembl_gene_id %in% rownames(report)))) # drops repeated rows, in case someone made a mistake somewhere...
    if (nrow(gene.set) > 2) {
      cat("Creating gene set heatmap for",gene.set.name,"\n")
      snif.gene.list.heat.map( 
        gene.set, myreport=report, 
        comps=names(comp.labels),
        comp.names=comp.labels,
        filename=file.path(dest.dir, paste0("GeneSet_L2FC_Heatmap-", make.names(gene.set.name) )),
        zlim=zlim, zcols=zcols,
        order="row.names"     # default: hc
      )
    } else {
      warning("! SKIPPING gene set heatmap for",gene.set.name,": not enough unique genes!\n")
    }
  }
  
}



# function to make Excel file of gene set reports
snif.make.gene.list.Excel <- function( xl.gene.sets, report,
                                       annotation, project, dest.dir=tempdir()) {

  if (!require(XLConnect)) stop("Requires XLConnect")
  
  # report column selection
  report.df <- data.frame(signif(report,4));
  merged.df <- merge(annotation, report.df, by.x="ensembl_gene_id", by.y="row.names");
  rownames(merged.df) <- merged.df$ensembl_gene_id;
  merged.df <- merged.df[,!(colnames(merged.df) %in% c("Row.names","mgi_id","entrezgene"))]; # exclude ID cruft
  
  setNames.duplicated <- names(which(table(sub("^.*:","",unique(all.gene.sets$fullName))) > 1));
  
  # file handling: must remove previous if exists otherwise won't overwrite
  xl.file.name <- sprintf("%s/GeneSet_Reports_%s.xls", dest.dir, project);
  if (file.exists(xl.file.name))
    file.remove(xl.file.name)
  
  wb <- loadWorkbook(filename=xl.file.name , create = TRUE);
  done.gene.set.names <- c()
  
  # 1 sheet per gene set
  for (gene.set.name in unique(all.gene.sets$fullName)) {
    
    gene.set <-  unique(subset(all.gene.sets, (fullName == gene.set.name) & 
                                 (ensembl_gene_id %in% rownames(report)))) # drops repeated rows, in case someone made a mistake somewhere...
    
    xl.name <- sub("^.*:","",gene.set.name);
    if(xl.name %in% setNames.duplicated){
      xl.name <- sub(":","_",gene.set.name);
    }
    
    if (xl.name %in% done.gene.set.names) {
      warning("! WARNING: gene set named '",gene.set.name,"' has already been written, ignoring !")
      next
    }
    
    if (is.null(gene.set) || (nrow(gene.set) <= 2)) {
      warning("! SKIPPING gene set report for",gene.set.name,": not enough unique genes!\n")
      next
    }
    
    cat("Creating report for",gene.set.name,"\n")
    createSheet(wb , xl.name);
    sub.report <- merged.df[as.character(gene.set$ensembl_gene_id),];
    sub.report <- sub.report[order(as.character(sub.report$mgi_symbol) ),];
#    print(head(sub.report));
#  return();
    writeWorksheet( wb , 
                    data=sub.report, 
                    sheet=xl.name , startRow = 1, startCol = 1,
                    header = TRUE );
    setColumnWidth(wb, sheet = xl.name, column = 1:(dim(sub.report)[2]));
    done.gene.set.names <- c(done.gene.set.names, xl.name);

  }
  
  saveWorkbook(wb)
  cat("Created",xl.file.name,"\n")

}


# EOF
