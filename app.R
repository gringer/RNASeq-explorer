## RNASeq Data Expolorer
library(shiny);
values <- reactiveValues();
## Shiny App for MM_Th2_DC_RNAseq

## Notes on editing this Shiny App
# * Make the output as dynamic as possible
# * Consider different input situations, will the function "do what I mean" for each situation?
# ** single transcript, single condition
# ** single transcript, multiple conditions
# ** multiple transcript, single condition
# ** multiple transcript, two conditions
# ** multiple transcript, multiple conditions

sortConds <- function(condNames){
  ## Produces a properly-sorted version
  condPops <- sub("\\..*$","",condNames);
  condConds <- sub("\\..*$","",sub("^.*?\\.","",condNames));
  return(sort(condNames));
  #return(condNames[order(-xtfrm(condConds),condPops)]);
}

#    values$allSubNames <- isolate(values$allSubNames)[grep("(CD)?(11b|103|TN)(\\+?)\\.(Nb|PBS|DBP.FITC|UT)",isolate(values$allSubNames))];


dataDirs <- list.dirs("RNASeq_Files/DESeq2/", recursive=FALSE, full.names=FALSE);

DE.comparisons <- NULL;

for(dDir in dataDirs){
  dataFile <- list.files(sprintf("RNASeq_Files/DESeq2/%s", dDir),
                         pattern=sprintf("^(David|Alex)_%s_All.*.csv.gz$", dDir), full.names=TRUE);
  data.topLine <- scan(dataFile, what=character(), nlines=1, quiet=TRUE, sep=",")[-1];
  DE.comparisons <- rbind(DE.comparisons,data.frame(experiment=dDir, test=sortConds(data.topLine)));
}
DE.comparisons$key <- paste0(DE.comparisons$experiment,":",DE.comparisons$test);

outliers <- c("CD11b+.UT.1","CD11b+.UT.2","CD11b+.UT.r1");

geneLists.df <- read.csv("RNASeq_Files/GeneLists_2016-Jul-29.csv", stringsAsFactors=FALSE);
ensemblData.df <- read.csv("RNASeq_Files/ensembl_gene_data.csv.gz", stringsAsFactors=FALSE);

clickedPoints <- list();
nextPointID <- 1;

## rescale data to fit within a given range, setting missing values to just outside the range
rescale <- function(values, newMin, newMax=NULL){
  if(is.null(newMax)){
    newMax <- newMin[2];
    newMin <- newMin[1];
  }
  values <- (values - min(values, na.rm=TRUE)) / (max(values, na.rm=TRUE)-min(values, na.rm=TRUE));
  values[is.na(values)] <- -0.05;
  return(unlist(values * (newMax-newMin) + newMin));
}

## generate points with multidimensional visual properties
rectPoints <- function(x, y, cex=rep(1,length(x)), len=rep(0,length(x)), ang=rep(0,length(x)), 
                       vis=rep(TRUE, length(x)), col=rep("#000000A0", length(x))){
  xAsp <- (diff(par("usr")[1:2])/par("pin")[1]) / 25.4;
  yAsp <- (diff(par("usr")[3:4])/par("pin")[2]) / 25.4;
  r90 <- (pi/2);
  for(ip in 1:length(x)){
    x1 <- x[ip]; y1 <- y[ip]; cip <- cex[ip];
    points(x=x1, y=y1, cex=cip, pch=ifelse(vis[ip],16,1), col=col[ip]);
    if(len[ip] > 0){
      aR <- ang[ip] * (pi/180);
      x2=x1+(len[ip])*xAsp*cos(aR);
      y2=y1+(len[ip])*yAsp*sin(aR);
      aP <- (ang[ip]+90) * (pi/180);
      aPo <- (ang[ip]+270) * (pi/180);
      px <- c(x1+(xAsp*cip)*cos(aP),x2+(xAsp*cip)*cos(aP),x2+(xAsp*cip)*cos(aPo),x1+(xAsp*cip)*cos(aPo));
      py <- c(y1+(yAsp*cip)*sin(aP),y2+(yAsp*cip)*sin(aP),y2+(yAsp*cip)*sin(aPo),y1+(yAsp*cip)*sin(aPo));
      polygon(x=px, y=py, density=ifelse(vis[ip],NA,0), col=col[ip], border=col[ip]);
    }
  }
}

nextCustom <- if("Custom1" %in% unique(geneLists.df$listName)){
  maxCustom <- max(as.numeric(sub("^Custom","",grep("^Custom[0-9]+$",geneLists.df$listName, value = TRUE))));
  nextCustom <- paste0("Custom",maxCustom+1);
} else {
  "Custom1";
  };

## Define UI for dataset viewer application
ui <- pageWithSidebar(
  
  # Application title
  headerPanel("FR RNASeq data explorer"),
  
  sidebarPanel(
    fluidRow(
      column(6,radioButtons("datasetName", "Data Set", choices=c("VSTPk","ImmGen"), selected = "VSTPk")),
      column(6,
             conditionalPanel(condition = "(input.outputPanel == 'heatmap')", tags$div(title="Merge heatmap results into a single column",checkboxInput("mergeHM", "Merge")))
      )
    ),
    verticalLayout(
      fluidRow(column(6,conditionalPanel(condition = "(input.outputPanel == 'heatmap') || (input.outputPanel == 'venn')",radioButtons("colourScheme", "Colour Scheme",
                                     c("white/red","blue/white/red","blue/cyan/red/brown","blue/yellow/red")))),
               column(6,conditionalPanel(condition = "(input.datasetName != 'VSTPk') && (input.datasetName != 'ImmGen') && (input.outputPanel == 'heatmap')",checkboxInput("showPadj","Show Padj Values")))),
      conditionalPanel(condition = "(input.outputPanel == 'heatmap')",
                       fluidRow(column(6,radioButtons("sampleOrder", "Sample Order",
                                                      c("Default","Treatment","Population","Cluster"))),
                                column(6,radioButtons("geneOrder", "Gene Order",
                                                      c("Default","Fold Change","Cluster")))))),
    radioButtons("geneSetClass", "Gene Selection", choices=c("All", "List", "Custom"), selected="All", inline=TRUE),
    conditionalPanel(condition = "input.geneSetClass == 'All'",
                     numericInput("expressionThreshold", "Expression Threshold", value = -3, min = -3, max = 21, width="6em")),
    conditionalPanel(condition = "input.geneSetClass == 'Custom'",
                     fluidRow(
                       column(6,textInput("customGeneSet", "Gene Names (comma separated)", value = "")),
                       column(6,textInput("customSetName", "Set Name", value = "")))),
    conditionalPanel(condition = "input.geneSetClass == 'List'",
      verticalLayout(
        selectInput("geneCat", "Gene List Category:", 
                    choices = unique(geneLists.df$listCategory)),
        selectInput("geneList", "Gene List:",
                    choices = "Loading..."))),
    conditionalPanel(condition = "input.datasetName == 'VSTPk'",
        radioButtons("replicateJoin", "Replicates", c("Combined","Split"), inline = TRUE)),
    fluidRow(
      column(8,selectizeInput("dataSubsets", "Data Subsets", choices=NULL, selected=NULL, multiple=TRUE)),
      column(4,actionButton("selectAll", label="Select All"))),
    fluidRow(
      column(2,downloadButton('downloadCSV', 'Download CSV')))
    ),
  
  # Main panel display for the multiple plot types.
  mainPanel(
    tabsetPanel(id = "outputPanel",
      tabPanel("Heat Maps", value = "heatmap",
               verticalLayout(
                 fluidRow(
                   column(12,downloadButton("pdfHeatMap", label = "View PDF"),
                          downloadButton("heatMapKey", label = "key PDF"))),
                 fluidRow(column(12,plotOutput("heatmapView", height="800px"))))
      ),
      tabPanel("Volcano Plots", value = "volcano",
               verticalLayout(
                 fluidRow(
                   column(3,downloadButton("pdfVolcano", label = "View PDF")),
                   column(3,downloadButton("pdfSmallVolcano", label = "Small PDF"))),
                 fluidRow(
                   column(3, numericInput("volcanoYLimit", label="Y axis limit", value=60, min=0)),
                   column(3, numericInput("volcanoXLimit", label="X axis limit", value=4, min=0)),
                   column(3, numericInput("ablineLog10Pvalue", label="-Log10 p-value",value=2, min=0)),
                   column(3, numericInput("ablineFoldChange", label="FC Threshold", value=1.5, min=0, step=0.1))),
                 fluidRow(column(12,plotOutput("volcanoView", brush=brushOpts(id="plot_brush"), height="600px"))),
                 h2("Selected point(s)"),
                 tableOutput("volcanoInfo"))
      ),
      tabPanel("MA Plots", value = "MA",
               verticalLayout(
                 fluidRow(
                   column(3,downloadButton("pdfMA", label = "View PDF")),
                   column(4, numericInput("MALog10Pvalue", label="-Log10 p-value Threshold",value=2, min=0)),
                   column(5, numericInput("MAFoldChange", label="Fold change Threshold", value=1.5, min=0, step=0.1))),
                 fluidRow(column(12,plotOutput("MAView", brush=brushOpts(id="plot_brush"), height="600px"))),
                 h2("Selected point(s)"),
                 tableOutput("MAInfo"))
      ),
      tabPanel("Venn Diagrams", value = "venn",
               verticalLayout(
                 fluidRow(
                   column(3,downloadButton("pdfVenn", label = "View PDF")),
                   column(4, numericInput("vennLog10Pvalue", label="-Log10 p-value Threshold",value=2, min=0)),
                   column(5, numericInput("vennFoldChange", label="Fold change Threshold", value=1.5, min=0, step=0.1))),
                 fluidRow(column(12,plotOutput("vennView", hover = "venn_hover", brush=brushOpts(id="venn_brush"), height="600px"))),
                 h2("Diagram Summary"),
                 tableOutput("vennSummary"),
                 h2("Selected point(s)"),
                 tableOutput("vennInfo"))
      ),
      tabPanel("Scatter Plots", value = "scatter",
               verticalLayout(
                 fluidRow(
                   column(3,downloadButton("pdfScatter", label = "View PDF")),
                   column(9,selectizeInput("scatterGeneNames", "Gene Name(s)", choices="[choose data subsets]", multiple = TRUE))),
                 fluidRow(
                   column(9,plotOutput("scatterView", hover = "scatter_hover", brush=brushOpts(id="scatter_brush"), height="600px")),
                   column(3,plotOutput("scatterKey", height="600px"))
                 ),
                 h1("Selected point(s)"),
                 tableOutput("scatterInfo")
               )
      ),
      tabPanel("PCA", value = "pca",
               verticalLayout(
                 fluidRow(
                   column(2, downloadButton("pdfPca", label = "View PCA")),
                   column(4, selectInput("pcaDESelector", "DE Set Selector", 
                                         choices=c("All", DE.comparisons$key), selected="All", selectize=FALSE, multiple=FALSE)),
                   column(3, numericInput("pcaLog10Pvalue", label="-Log10 p-value Threshold",value=2, min=0)),
                   column(3, numericInput("pcaFoldChange", label="Fold change Threshold", value=1.5, min=0, step=0.1))),
                 fluidRow(
                   column(3, downloadButton("pdfSmallPca", label = "Small PCA")),
                   column(4, selectInput("pcaInput1real", "x Principal Component", 
                              choices=c(paste0("Component ",1:5)), selected="Component 1")),
                   column(5,selectInput("pcaInput2real", "y Principal Component", 
                                        choices=c(paste0("Component ",1:5)), selected="Component 2"))),
                 plotOutput("pcaView", height="800px"))
      ),
      tabPanel("Log", value = "log",
               verticalLayout(
                 h1("Debug Output"),
                 verbatimTextOutput("logStdout")
               )
      )
    ),
    h2("Data Sources"),
    uiOutput("sourceFileOutput")
  )
);

logOutput <- function(input, outputName){
  timeStr <- as.character(Sys.time());
  cat(file = "accesslog.txt", append=TRUE, sprintf("Output requested from '%s' on %s:\n", outputName, timeStr));
  for(n in names(input)){
    if(is.character(input[[n]]) || (is.numeric(input[[n]]) && (length(input[[n]]) == 1))){
      cat(file = "accesslog.txt", append=TRUE, sprintf("  [%s] - %s\n", n, substring(paste(input[[n]], collapse=";"),1,100)));
    }
  }
}

logDebug <- function(logObject){
  isolate({
    res <- if((typeof(logObject) == "character") && (length(logObject) == 1)){
      logObject;
    } else {
      paste(capture.output(logObject), collapse="\n     ");
    }
    cat("     ", res, "\n", sep="");
    values$logText[values$logPos] <- paste0(res,"\n");
    values$logPos <- values$logPos + 1;
    if(values$logPos > 100){
      values$logPos <- 1;
    }
  })
}

showSourceFiles <- function(grid=FALSE){
  if(!grid){
    mtext(paste(unique(sub("^.*/","",values$fileName)), collapse=";"), side=3, cex=1/length(values$fileName));
  } else {
    grid.text(paste(unique(sub("^.*/","",values$fileName)), collapse=";"), y=unit(0.9,"npc"), gp=gpar(cex=1/length(values$fileName)));
  }
}

getInputMatrix <- function(input, dsName = input$datasetName, subNames = input$dataSubsets, exprThresh = input$expressionThreshold){
  if(is.null(dsName) || is.null(subNames) || 
       ((dsName == "VSTPk" || dsName == "ImmGen") && any(grepl("\\.vs\\.",subNames))) ||
       ((dsName != "VSTPk" && dsName != "ImmGen") && !any(grepl("\\.vs\\.",subNames)))){
    logDebug("Returning from getInputMatrix...");
    values$fileName <- NULL;
    return(NULL);
  }
  if(!is.null(values[[dsName]])){
    logDebug(paste0("Retrieving cached input matrix: ",dsName));
  } else {
    if (dsName == "VSTPk"){
      dataFile <- values$vstpkFileName;
      dataColNames <- tail(scan(dataFile, what=character(), sep=",", nlines=1, quiet=TRUE),-1);
      subData <- as.matrix(read.csv(dataFile, row.names=1, stringsAsFactors=FALSE));
      colnames(subData) <- dataColNames;
      values[[dsName]] <- subData;
      values[[paste0("fileName.",dsName)]] <- dataFile;
    } else if(dsName == "ImmGen"){
      dataFile <- values$immgenFileName;
      dataColNames <- tail(scan(dataFile, what=character(), sep=",", nlines=1, quiet=TRUE),-1);
      subData <- as.matrix(read.csv(dataFile, row.names=1, stringsAsFactors=FALSE));
      colnames(subData) <- dataColNames;
      values[[dsName]] <- subData;
      values[[paste0("fileName.",dsName)]] <- dataFile;
    } else {
      dataFile <- list.files(sprintf("RNASeq_Files/DESeq2/%s", dsName),
                             pattern=sprintf("^(David|Alex)_%s_All.*.csv.gz$", dsName),
                             full.names=TRUE);
      dataColNames <- tail(scan(dataFile, what=character(), sep=",", nlines=1, quiet=TRUE),-1);
      subData <- as.matrix(read.csv(dataFile, row.names=1, stringsAsFactors=FALSE));
      ## only keep rows with unique names
      colnames(subData) <- dataColNames;
      values[[dsName]] <- subData;
      values[[paste0("fileName.",dsName)]] <- dataFile;
    }
  } 
  values$fileName <- values[[paste0("fileName.",dsName)]];
  if(!is.null(exprThresh) && (exprThresh > -3)){
    maxExpression <- values$maxExpression;
    if(is.null(maxExpression)){
      values$maxFile <- sub("DESeq2_blind_local_all","maxExpression",values$vstpkFileName);
      maxExpression.mat <- read.csv(values$maxFile, row.names=1);
      maxExpression <- maxExpression.mat[,1]; names(maxExpression) <- rownames(maxExpression.mat);
      values$maxExpression <- maxExpression;
    }
    values$fileName <- isolate(c(values$fileName,values$maxFile));
    threshNames <- names(maxExpression)[maxExpression > exprThresh];
    retVal <- (values[[dsName]])[rownames(values[[dsName]]) %in% threshNames,,drop=TRUE];
  } else {
    retVal <- values[[dsName]];
  }
  colsToInvert = values$invert[values$invert %in% colnames(retVal)];
  if(length(colsToInvert) > 0){
    retVal[,colsToInvert] <- -retVal[,colsToInvert];
  }
  return(retVal);
}

getPadjMatrix <- function(input, refMatrix, dsName = input$datasetName, subNames = input$dataSubsets){
  if((dsName %in% c("VSTPk", "ImmGen")) || (is.null(refMatrix))){
    return(NULL);
  }
  pData <- NULL;
  padjName <- paste0("padj.",dsName);
  if(!is.null(values[[padjName]])){
    logDebug(paste0("Retrieving cached p-value matrix:",padjName));
    pData <- (values[[padjName]]);
  } else {
    dataFile <- list.files(sprintf("RNASeq_Files/DESeq2/%s", dsName),
                           pattern=sprintf("^padj_(David|Alex)_%s_All.*.csv.gz$", dsName),
                           full.names=TRUE);
    values$fileName <- c(values$fileName,dataFile);
    pColNames <- tail(scan(dataFile, what=character(), sep=",", nlines=1, quiet=TRUE),-1);
    pData <- as.matrix(read.csv(dataFile, row.names=1, stringsAsFactors=FALSE));
    pData[is.na(pData)] <- 1.1; # set NA p-values to 1.1 [i.e. a easily detected value with low probability]
    pData[pData == 0] <- 1.1; # hack... set 0 p-value also to 1.1 [fix garbage input data generation]
    colnames(pData) <- pColNames; ## replace with actual names, rather than R's transformed names
    values[[padjName]] <- pData;
  }
  return(pData[match(toupper(rownames(refMatrix)),toupper(rownames(pData))), colnames(refMatrix), drop=FALSE]);
}

getMatrixRange <- function(input){
  rangeDSName <- paste0("range",input$datasetName);
  if(!is.null(values[[rangeDSName]])){
    return(values[[rangeDSName]]);
  }
  if(input$datasetName %in% c("VSTPk","ImmGen")){
    resRange <- c(0,16);
    values[[rangeDSName]] <- resRange;
    return(resRange);
  }
  subData <- getInputMatrix(input);
  if(is.null(subData)){
    return(c(-1,1));
  }
  dataRange <- range(subData[abs(subData)>0], na.rm=TRUE);
  quantRange <- quantile(subData[abs(subData)>0],c(0.01,0.99), na.rm=TRUE);
  resRange <- colMeans(rbind(dataRange,quantRange));
  if(resRange[1] < (-0.5*resRange[2])){
    resRange <- c(-max(abs(resRange)),max(abs(resRange)));
  } else {
    resRange <- c(0,resRange[2]);
  }
  values[[rangeDSName]] <- resRange;
  return(resRange);
}

getRampCols <- function(colourScheme){
  MIMR.red <- "#DF2726";
  MIMR.orange <- "#EC772D";
  rampCols <- if(colourScheme == "blue/white/red"){
    c("darkblue","#D0D0FF","white","#FFD0D0",MIMR.red);
  } else if(colourScheme == "blue/cyan/red/brown"){
    c("darkblue","cyan3","aliceblue","white","lightpink",MIMR.red,"brown4");
  } else if(colourScheme == "blue/yellow/red"){
    c("darkblue","chartreuse","yellow",MIMR.orange,MIMR.red);
  } else if(colourScheme == "white/red"){
    c("white","#FFD0D0",MIMR.red);
  }
  return(rampCols);
}

getColBreaks <- function(limits=NULL, nBreaks=65, dataName){
  # 'getBreaks function' - Returns scaled colour breaks for different heatmaps (line 200);

  if(is.null(limits)){
    return(NULL);
  }
  if(!is.numeric(nBreaks)){
    return(NULL);
  } else {
    if( nBreaks %% 2 == 0 ) { 
      sideBreaks = (nBreaks/2) 
      } else {
        sideBreaks = ((nBreaks-1)/2);
  }}
  if(dataName %in% c("VSTPk","ImmGen")){
    breaks=c(seq(min(limits), max(limits), length.out=nBreaks));
    return(breaks);
  } else {
    breaks = c(seq(min(limits), 0-10^-4, length.out=sideBreaks), 0, seq(0+10^-4, max(limits), length.out=sideBreaks))^2;
    breaks = (breaks/max(breaks))*max(limits);
    breaks[1:sideBreaks] = breaks[1:sideBreaks]*-1
    return(breaks);
  }
}
  
getGeneSetMatrix <- function(input, geneLimit=100, forceSplit=FALSE){
  dsName <- input$datasetName;
  subNames <- input$dataSubsets;
  if(is.null(dsName) || is.null(subNames) || (length(subNames) < 1)){
    return(NULL);
  }
  plotNames <- NULL;
  geneList <- if(input$geneSetClass == "Custom"){
    if(is.null(input$geneSetClass)){
      return(NULL);
    };
    geneNames <- unlist(strsplit(gsub(", +",",",input$customGeneSet),","));
    if(any(grepl("^ENSMUSG",geneNames))){
      plotNames <- geneNames;
      geneNames;
    } else {
      plotNames <- geneNames;
      ensemblData.df$ensembl_gene_id[match(toupper(geneNames),
                                           toupper(ensemblData.df$mgi_symbol))];
    }
  } else if(input$geneSetClass == "List") {
    plotNames <- geneLists.df[(geneLists.df$listCategory %in% input$geneCat) 
                              & (geneLists.df$listName %in% input$geneList),"plotName"];
    geneLists.df[(geneLists.df$listCategory %in% input$geneCat) 
                 & (geneLists.df$listName %in% input$geneList),"ensembl_gene_id"];
  } else {
    "All";
  }
  if((input$geneSetClass != "All") && ((length(subNames) < 1) || (length(geneList) < 1))){
    return(NULL);
  }
  subData <- getInputMatrix(input);
  if(is.null(subData)){
    return(NULL);
  }
  colLookups <- NULL;
  if (dsName %in% c("VSTPk", "ImmGen")){
    colLookups <- sub("^X([0-9])","\\1",colnames(subData));
  } else {
    colLookups <- colnames(subData);
  }
  ## remove replicates
  if(input$replicateJoin == "Combined"){
    colLookups[colLookups %in% outliers] <- "**OUTLIER**";
    colLookups <- sub("(\\.r)?\\.r?[0-9]+$","",colLookups);
  }
  if(any(!subNames %in% colLookups)){
    logDebug("-- At least one SubName not found in colLookups --\n");
    logDebug(head(subNames[!subNames %in% colLookups]));
    logDebug(head(subNames));
    logDebug(head(colLookups));
    logDebug(head(subData));
    return(NULL);
  }
  if(input$geneSetClass != "All"){
    ens.ids <- sub("^.*?_","",rownames(subData));
    subData <- as.matrix(subData[match(geneList, ens.ids), colLookups %in% subNames, drop=FALSE]);
    ## only keep rows with unique names
    subData <- subData[match(unique(rownames(subData)),rownames(subData)),,drop=FALSE];
  } else {
    tmpRowOrder <- 1:nrow(subData);
    if(geneLimit != -1) {
      # only show first few lines
      tmpRowOrder <- head(order(-apply(subData[,colLookups %in% subNames,drop=FALSE],1,sd)),geneLimit);
    }
    ## only keep rows with unique names
    subData <- subData[match(unique(rownames(subData)[tmpRowOrder]),rownames(subData)),colLookups %in% subNames,drop=FALSE];
  }
  colnames(subData) <- colLookups[colLookups %in% subNames];
  if((dsName %in% c("VSTPk")) && (input$replicateJoin == "Combined") && (!forceSplit)){
    ## combine replicate columns
    colIDs <- colnames(subData);
    subData <- sapply(subNames, function(x){apply(subData[,colIDs == x, drop=FALSE],1,mean)});
  } else {
    # ERROR
    if(any(!subNames %in% colnames(subData))){
      logDebug(head(subNames));
      logDebug(head(colnames(subData)));
      logDebug(str(subData));
    }
    subData <- subData[, subNames, drop=FALSE];
  }
  if(is.null(subData)){
    return(NULL);
  }
  rowOrder <- 1:nrow(subData);
  ## Apply additional reordering as necessary
  if(input$geneOrder == "Cluster" && (nrow(subData) > 1)){
    ## remove rows with NA values (e.g. label rows)
    rowOrder <- rowOrder[complete.cases(subData[rowOrder,])];
    ## cluster by hclust method
    rowOrder <- rowOrder[hclust(dist(subData[rowOrder,]))$order];
  } else if((input$geneOrder == "Expression") || (input$geneOrder == "Fold Change")){
    rowOrder <- rowOrder[order(rowSums(subData[rowOrder,]))];
  }
  colOrder <- 1:ncol(subData);
  if(input$sampleOrder == "Cluster" && (ncol(subData) > 1)){
    colOrder <- colOrder[hclust(dist(t(subData[,colOrder, drop=FALSE])))$order];
  } else if(input$sampleOrder == "Population"){
    colPops <- sub("\\..*$","",colnames(subData));
    colConds <- sub("\\..*$","",sub("^.*?\\.","",colnames(subData)));
    colOrder <- colOrder[order(colPops,-xtfrm(colConds))];
  } else if(input$sampleOrder == "Treatment"){
    colPops <- sub("\\..*$","",colnames(subData));
    colConds <- sub("\\..*$","",sub("^.*?\\.","",colnames(subData)));
    colOrder <- colOrder[order(-xtfrm(colConds),colPops)];
  }
  if((input$geneSetClass == "All") && (geneLimit != -1)) {
    subData <- rbind(subData[rowOrder,colOrder,drop=FALSE],"Abridged..." = rep(NA,length(colOrder)));
  } else {
    subData <- subData[rowOrder,colOrder,drop=FALSE];
  }
  return(subData);
}


# Define server logic required to summarize and view the selected dataset
server <- function(input, output, session) {

  source("snif_gene_list_heat_map.r");
  values$selectButtonText = "Select All";

  output$heatmapView <- renderPlot({
    ## If this is the first time for plotting, change the panel ID and return quickly
    values$outputTable <- NULL;
    subNames <- input$dataSubsets;
    if(is.null(subNames) || (length(subNames) < 1) || is.null(input$geneSetClass)){
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), ann=FALSE, frame.plot=FALSE, axes=FALSE);
      text(0,0,"[Please Select Some Data Sets]", adj=c(0.5,0.5), offset=0);
    } else {
      subData <- getGeneSetMatrix(input);
      if(is.null(subData)){
        NULL;
      } else {
        subNames <- colnames(subData);
        zlim <- getMatrixRange(input);
        zcols <- colorRampPalette(getRampCols(input$colourScheme))(64);
        breaks <- getColBreaks(zlim, 65, input$datasetName);
        pMat <- NULL;
        if(!is.null(input$showPadj) && input$showPadj){
          pMat <- getPadjMatrix(input, subData);
          ## remove Ensembl codes
          rownames(pMat) <- sub("_.*$","",rownames(pMat));
        }
        ## remove Ensembl codes
        rownames(subData) <- sub("_.*$","",rownames(subData));
        if(input$mergeHM){
          gsmGeneNames <- rownames(subData);
          gsmSampNames <- colnames(subData);
          subData <- matrix(apply(subData,1,mean), ncol=1);
          rownames(subData) <- gsmGeneNames;
          colnames(subData) <- sprintf("Merged_%d", length(gsmSampNames));
          if(!is.null(pMat)){
            pMat <- matrix(apply(pMat, 1, min), ncol=1);
            rownames(pMat) <- gsmGeneNames;
            colnames(pMat) <- sprintf("Merged_%d", length(gsmSampNames));
          }
        }
        gringene.heat.map(subData, pmat = pMat, 
                          longestLabel=subNames[which(nchar(subNames) == max(nchar(subNames)))[1]], 
                          zlim=zlim, zcols=zcols, breaks=breaks);
      }
    }
  });

  ## multi-dimensional scatter plot
  plotScatter <- function(plotData, plotType=1){
    if(is.null(plotData) || (nrow(plotData) < 1)){
      plot(NA,xlim=c(0,1), ylim=c(0,1), axes=FALSE, ann=FALSE);
      text(0.5,0.5,"Please choose some genes to display");
    } else {
      if(plotType == 3){
        layout(matrix(2:1, ncol=2), widths=c(0.8,0.2));
      }
      if(bitwAnd(plotType,2) != 0){
        par(mar=rep(0.5,4));
        plot(NA, xlim=c(-0.1,1.1), ylim=c(0,17), axes=FALSE, ann=FALSE);
        for(ig in 1:min(7,nrow(plotData))){
          text(0.5,18-(ig*2),rownames(plotData)[ig], xpd=TRUE);
          text(x=c(0,1),17-(ig*2),range(plotData[ig,], na.rm=TRUE), xpd=TRUE);
          if(ig==1){
            arrows(x0=c(0.45,0.55), x1=c(0.25,0.75), y0=17-(ig*2), length=0.1, lwd=2);
          }
          if(ig==2){
            arrows(x0=c(1/3,2/3), y0=17-(ig*2)+c(0.5,-0.5), y1=17-(ig*2)-c(0.5,-0.5), length=0.1, lwd=2);
          }
          if(ig==3){
            rectPoints(x=c(1/3,1/2,2/3), y=rep(17-(ig*2),3), cex=c(0.75,2.375,4));
          }
          if(ig==4){
            rectPoints(x=c(1/3,1/2,2/3), y=rep(17-(ig*2),3), cex=rep(2,3), col=hsv(h=0.5,v=c(0.1, 0.55, 1), alpha=0.75));
          }
          if(ig==5){
            rectPoints(x=c(1/3,1/2,2/3), y=rep(17-(ig*2),3), cex=rep(2,3), col=hsv(h=c(1/3,1/2,2/3), v=1, alpha=0.75));
          }
          if(ig==6){
            rectPoints(x=c(1/3,1/2,2/3), y=rep(17-(ig*2),3), cex=rep(2,3), len=c(0.75,2.375,4));
          }
          if(ig==7){
            rectPoints(x=c(1/3,1/2,2/3), y=rep(17-(ig*2),3), cex=rep(2,3), len=rep(2.375,3), ang=c(0,135,270));
          }
        }
      }
      if(bitwAnd(plotType,1) != 0){
        par(mar=c(5,4,1,0)+0.1);
        if(nrow(plotData) == 1){
          plotOrder <- order(plotData[1,]);
          plot.default(x=1:ncol(plotData),y=rescale(plotData[1,plotOrder], range(plotData[1,], na.rm=TRUE)), 
                       xlab="Sample", ylab=rownames(plotData), cex=2,
                       pch=ifelse(is.na(apply(plotData,2,max)),1,16), col="#000000A0");
        } else if(nrow(plotData) == 2){
          plot.default(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)), 
                       y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                       xlab=rownames(plotData)[1], ylab=rownames(plotData)[2], cex=2,
                       pch=ifelse(is.na(apply(plotData,2,max)),1,16), col="#000000A0");
        } else if(nrow(plotData) == 3){
          plot.default(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)), 
                       y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                       xlab=rownames(plotData)[1], ylab=rownames(plotData)[2], type="n");
          rectPoints(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)),
                     y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                     cex=rescale(plotData[3,], 0.75, 4), 
                     vis=!is.na(apply(plotData,2,max)));
        } else if(nrow(plotData) == 4){
          plot.default(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)), 
                       y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                       xlab=rownames(plotData)[1], ylab=rownames(plotData)[2], type="n");
          rectPoints(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)),
                     y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                     cex=rescale(plotData[3,], 0.75, 4),
                     col=hsv(h=0.5,v=rescale(plotData[4,], 0.1, 1), alpha=0.75),
                     vis=!is.na(apply(plotData,2,max)));
        } else if(nrow(plotData) == 5){
          plot.default(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)), 
                       y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                       xlab=rownames(plotData)[1], ylab=rownames(plotData)[2], type="n");
          rectPoints(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)),
                     y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                     cex=rescale(plotData[3,], 0.75, 4),
                     col=hsv(h=rescale(plotData[5,], 1/3, 2/3), v=rescale(plotData[4,], 0.1, 1), alpha=0.75),
                     vis=!is.na(apply(plotData,2,max)));
        } else if(nrow(plotData) == 6){
          plot.default(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)), 
                       y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                       xlab=rownames(plotData)[1], ylab=rownames(plotData)[2], type="n");
          rectPoints(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)),
                     y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                     cex=rescale(plotData[3,], 0.75, 4),
                     col=hsv(h=rescale(plotData[5,], 1/3, 2/3), v=rescale(plotData[4,], 0.1, 1), alpha=0.75),
                     len=rescale(plotData[6,], 0.75, 4),
                     vis=!is.na(apply(plotData,2,max)));
        } else if(nrow(plotData) >= 7){
          plot.default(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)), 
                       y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                       xlab=rownames(plotData)[1], ylab=rownames(plotData)[2], type="n");
          rectPoints(x=rescale(plotData[1,], range(plotData[1,], na.rm=TRUE)),
                     y=rescale(plotData[2,], range(plotData[2,], na.rm=TRUE)),
                     cex=rescale(plotData[3,], 0.75, 4),
                     col=hsv(h=rescale(plotData[5,], 1/3, 2/3), v=rescale(plotData[4,], 0.1, 1), alpha=0.75),
                     len=rescale(plotData[6,], 0.75, 4),
                     ang=rescale(plotData[7,], 0, 270),
                     vis=!is.na(apply(plotData,2,max)));
        }
        if(!is.null(input$scatter_brush)){
          sub.df <- data.frame(t(values$scatterData));
          if(prod(dim(sub.df)) > 0){
            if(ncol(sub.df) == 1){
              sub.df <- sub.df[order(sub.df[,1]),,drop=FALSE];
              sub.df$ID = 1:nrow(sub.df);
              sub.df <- sub.df[,2:1];
            }
            ## use nearPoints() for clicked points, brushedPoints() for brushed points
            res <- brushedPoints(sub.df, input$scatter_brush, colnames(sub.df)[1], colnames(sub.df)[2], allRows=TRUE);
            if(sum(res$selected_, na.rm=TRUE) > 0){
              res <- subset(res, selected_);
              text(x=res[,1], y=res[,2], labels=rownames(res), pos=1, xpd=NA, cex=0.75, offset=1);
            }
          }
        }
        showSourceFiles();
      }
    }
  }
  
  output$sourceFileOutput <- renderUI({
      HTML(paste(unique(sub("^.*/","",values$fileName)), collapse="<br />"));
  });
  
  output$scatterView <- renderPlot({
    dsSelected <- input$datasetName;
    subNames <- input$dataSubsets;
    if(length(subNames) == 0){
      values$scatterData <- NULL;
      values$scatterWholeData <- NULL;
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), ann=FALSE, frame.plot=FALSE, axes=FALSE);
      text(0,0,"[Please Select Some Data Sets]", adj=c(0.5,0.5), offset=0);
    } else {
      if(is.null(values$scatterDataSet) || (values$scatterDataSet != dsSelected) || is.null(values$scatterWholeData)){
        subData <- getGeneSetMatrix(input, geneLimit=-1);
        values$scatterDataSet <- dsSelected;
        values$scatterWholeData <- subData;
        ## remove Ensembl IDs from graphical display
        rownames(values$scatterWholeData) <- sub("_.*$","",rownames(values$scatterWholeData));
        values$scatterData <- NULL;
        return(NULL);
      }
      if(length(input$scatterGeneNames) < 1){
        plot(NA,xlim=c(0,1), ylim=c(0,1), axes=FALSE, ann=FALSE);
        text(0.5,0.5,"Please choose some genes to display");
        values$scatterData <- NULL;
        return(NULL);
      }
      values$scatterData <- round(values$scatterWholeData[input$scatterGeneNames,,drop=FALSE],2);
      plotScatter(values$scatterData);
    }
  });
  output$scatterKey <- renderPlot({
    if(!is.null(values$scatterData) && (prod(dim(values$scatterData)) > 0)){
      plotScatter(values$scatterData, plotType=2);
    }
  });
  output$scatterInfo <- renderTable({
    if(is.null(values$scatterData)){
      return(NULL);
    }
    sub.df <- data.frame(t(values$scatterData));
    if(prod(dim(sub.df)) == 0){
      return(NULL);
    }
    if(ncol(sub.df) == 1){
      sub.df <- sub.df[order(sub.df[,1]),,drop=FALSE];
      sub.df$ID = 1:nrow(sub.df);
      sub.df <- sub.df[,2:1];
    }
    ## use nearPoints() for clicked points, brushedPoints() for brushed points
    res <- brushedPoints(sub.df, input$scatter_brush, colnames(sub.df)[1], colnames(sub.df)[2], allRows=TRUE);
    if(sum(res$selected_, na.rm=TRUE) == 0){
      return(NULL);
    }
    res <- subset(res, selected_);
    res$selected_ <- NULL; # remove 'selected' column
    return(values$outputTable <- res);
  });
  output$pdfScatter <- downloadHandler(
    filename = function() {
      dateStr <- format(Sys.Date(),"%Y-%b-%d");
      paste0('scatter_', dateStr, '.pdf');
    },
    content = function(file) {
      logOutput(input, "scatter");
      pdf(file=file, width=10, height=8);
      plotScatter(values$scatterData, plotType=3);
      dev.off();
    },
    contentType = "application/pdf"
  );
  ## Volcano plot
  output$volcanoView <- renderPlot({
    ## If this is the first time for plotting, change the panel ID and return quickly
    dsSelected <- input$datasetName;
    subNames <- input$dataSubsets;
    values$outputTable <- NULL;
    if(is.null(subNames) || (length(subNames) < 1)){
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), ann=FALSE, frame.plot=FALSE, axes=FALSE);
      text(0,0,"[Please Select A Data Set]", adj=c(0.5,0.5), offset=0);
      return(); 
    }
  
    subData <- getGeneSetMatrix(input, geneLimit=-1);
    if(is.null(subData)){
      return(NULL);
    }
    zlim <- getMatrixRange(input);
    zcols <- colorRampPalette(getRampCols(input$colourScheme))(64);
    pMat <- getPadjMatrix(input, subData);
    if(is.null(pMat)){
      return(NULL);
    }
    ## aggregate results to generate one value per gene
    if(ncol(subData) > 1){
      subData <- apply(subData,1,mean, na.rm=TRUE);
      pMat <- apply(pMat,1,function(x){10^mean(log10(x), na.rm=TRUE)});
    }
    xVals <- subData;
    yVals <- -log10(pMat);
    full.cases <- !is.na(xVals) & !is.na(yVals);
    xVals <- xVals[full.cases];
    yVals <- yVals[full.cases];
    over.limit.y <- (yVals > input$volcanoYLimit);
    over.limit.x <- (abs(xVals) > input$volcanoXLimit);
    yVals[over.limit.y] <- input$volcanoYLimit;
    xVals[over.limit.x] <- input$volcanoXLimit * sign(xVals[over.limit.x]);
    over.limit <- over.limit.y | over.limit.x;
    if(!is.null(pMat)){
      plot(x=xVals, y=yVals, pch=ifelse(over.limit,24,16), cex=2,
           xlim=c(-input$volcanoXLimit,input$volcanoXLimit), ylim=c(0,input$volcanoYLimit),
           col=ifelse(abs(xVals)>=log2(input$ablineFoldChange),
                      ifelse(yVals>=input$ablineLog10Pvalue, "#DF272680", "#000000A0"),
                      ifelse(yVals>=input$ablineLog10Pvalue, "#F0909080", "#000000A0")), 
           xlab="log2FC", ylab=expression(-'log'[10](pAdj)));
      count.down <- sum((xVals <= -log2(input$ablineFoldChange)) & (yVals >= input$ablineLog10Pvalue));
      count.up <- sum((xVals >= log2(input$ablineFoldChange)) & (yVals >= input$ablineLog10Pvalue));
      text(x=min(xVals)*0.75, y=input$volcanoYLimit, labels=count.down, col="#DF2726", cex=3, pos=1);
      text(x=max(xVals)*0.75, y=input$volcanoYLimit, labels=count.up, col="#DF2726", cex=3, pos=1);
      abline(h=input$ablineLog10Pvalue, v=c(log2(input$ablineFoldChange),-log2(input$ablineFoldChange)), col="gray60", lty=2);
    } else {
      warning("pMat is null");
    }
  });
  
  output$volcanoInfo <- renderTable({
    dsSelected <- input$datasetName;
    if(dsSelected == "VSTPk" || dsSelected == "ImmGen"){
      return(NULL);
    }
    subNames <- input$dataSubsets;
    subData <- getGeneSetMatrix(input, geneLimit=-1);
    pMat <- getPadjMatrix(input, subData);
    if(is.null(pMat)){
      return(NULL);
    }
    colnames(subData) <- paste0(colnames(subData),".L2FC");
    colnames(pMat) <- paste0(colnames(pMat),".pAdj");
    ## aggregate results to generate one value per gene
    subAggr <- if(ncol(subData) > 1){
      apply(subData,1,mean);
    } else {
      subData;
    }
    pAggr <- if(ncol(pMat) > 1){
      apply(pMat,1,function(x){-mean(log10(x), na.rm=TRUE)});
    } else {
      -log10(pMat);
    }
    sub.df <- data.frame(cbind(round(subData,2), signif(pMat,4)));
    colnames(sub.df) <- c(colnames(subData),colnames(pMat));
    sub.df$meanL2FC <- subAggr;
    sub.df$meanAdjP <- pAggr;
    ## use nearPoints() for clicked points, brushedPoints() for brushed points
    res <- brushedPoints(sub.df, input$plot_brush, "meanL2FC", "meanAdjP", allRows=TRUE);
    if(sum(res$selected_, na.rm=TRUE) == 0){
      return(NULL);
    }
    res <- subset(res, selected_);
    res$selected_ <- NULL; # remove 'selected' column
    return(values$outputTable <- res);
  });

  output$MAView <- renderPlot({
    ## If this is the first time for plotting, change the panel ID and return quickly
    dsSelected <- input$datasetName;
    subNames <- input$dataSubsets;
    values$outputTable <- NULL;
    if(is.null(subNames) || (length(subNames) < 1)){
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), ann=FALSE, frame.plot=FALSE, axes=FALSE);
      text(0,0,"[Please Select A Data Set]", adj=c(0.5,0.5), offset=0);
      return(); 
    }
    subData <- getGeneSetMatrix(input, geneLimit=-1);
    if(is.null(subData)){
      return();
    }
    zlim <- getMatrixRange(input);
    # get median gene Values
    geneMedian <- apply(getInputMatrix(input, dsName="VSTPk", subNames=""),1,median);
    subData <- subData[rownames(subData) %in% names(geneMedian),,drop=FALSE];
    # retrieve p values
    pMat <- getPadjMatrix(input, subData);
    if(is.null(pMat)){
      return();
    }
    ## aggregate results to generate one value per gene
    if(ncol(subData) > 1){
      subData <- matrix(apply(subData,1,mean, na.rm=TRUE), ncol = 1, dimnames=list(rownames(subData),NULL));
      pMat <- matrix(apply(pMat,1,function(x){10^mean(log10(x), na.rm=TRUE)}), ncol=1, dimnames=list(rownames(subData),NULL));
    }
    #     dataCols <- ifelse(subData < zlim[1],1,
    #                        ifelse(subData > zlim[2],length(zcols),
    #                               (length(zcols)-1)*(subData-zlim[1])/(diff(zlim)) + 1));
    yVals <- subData;
    xVals <- matrix(geneMedian[rownames(subData)],ncol = 1);
    pVals <- -log10(pMat);
    full.cases <- !is.na(xVals) & !is.na(yVals) & !is.infinite(xVals) & !is.infinite(yVals);
    xVals <- xVals[full.cases];
    yVals <- yVals[full.cases];
    pVals <- pVals[full.cases];
    if(!is.null(pMat)){
      plot(x=xVals, y=yVals, pch=16, cex=2, 
           col=ifelse(abs(yVals)>=log2(input$ablineFoldChange),
                      ifelse(pVals>=input$ablineLog10Pvalue, "#DF272680", "#000000A0"),
                      ifelse(pVals>=input$ablineLog10Pvalue, "#F0909080", "#000000A0")), 
           xlab="Median VSTPk Expression (all samples)", ylab="log2FC");
      abline(h=c(log2(input$ablineFoldChange),-log2(input$ablineFoldChange)), col="gray60", lty=2)
    } else {
      warning("pMat is null");
    }
  });

  output$MAInfo <- renderTable({
    dsSelected <- input$datasetName;
    if(dsSelected == "VSTPk" || dsSelected == "ImmGen"){
      return(NULL);
    }
    subNames <- input$dataSubsets;
    subData <- getGeneSetMatrix(input, geneLimit=-1);
    if(is.null(subData)){
      return();
    }
    # get median gene Values
    geneMedian <- apply(getInputMatrix(input, dsName="VSTPk", subNames=""),1,median);
    ## filter subData to only include genes with expression values
    subData <- subData[rownames(subData) %in% names(geneMedian),,drop=FALSE];
    pMat <- getPadjMatrix(input, subData);
    if(is.null(pMat)){
      return(NULL);
    }
    colnames(subData) <- paste0(colnames(subData),".L2FC");
    colnames(pMat) <- paste0(colnames(pMat),".pAdj");
    geneMedian <- geneMedian[rownames(subData)];
    ## aggregate results to generate one value per gene
    subAggr <- if(ncol(subData) > 1){
      apply(subData,1,mean);
    } else {
      subData;
    }
    pAggr <- if(ncol(pMat) > 1){
      apply(pMat,1,function(x){-mean(log10(x), na.rm=TRUE)});
    } else {
      -log10(pMat);
    }
    sub.df <- data.frame(cbind(geneMedian,round(subData,2), signif(pMat,4)));
    colnames(sub.df) <- c("VSTPkMed",colnames(subData),colnames(pMat));
    sub.df$meanL2FC <- subAggr;
    sub.df$meanAdjP <- pAggr;
    ## use nearPoints() for clicked points, brushedPoints() for brushed points
    res <- brushedPoints(sub.df, input$plot_brush, "VSTPkMed", "meanL2FC", allRows=TRUE);
    if(sum(res$selected_, na.rm=TRUE) == 0){
      return(NULL);
    }
    res <- subset(res, selected_);
    res$selected_ <- NULL; # remove 'selected' column
    return(values$outputTable <- res);
  });

  output$vennView <- renderPlot({
    ## If this is the first time for plotting, change the panel ID and return quickly
    dsSelected <- input$datasetName;
    if(dsSelected == "VSTPk" || dsSelected == "ImmGen"){
      return(NULL);
    }
    values$outputTable <- NULL;
    subNames <- input$dataSubsets;
    if(is.null(subNames) || (length(subNames) < 2)){
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), ann=FALSE, frame.plot=FALSE, axes=FALSE);
      text(0,0,"[Please Select At Least Two Data Sets]", adj=c(0.5,0.5), offset=0);
    } else {
      source("venn_diagrams.R");
      subData <- getGeneSetMatrix(input, geneLimit=-1);
      pMat <- getPadjMatrix(input, subData);
      cat(file="out.txt", dim(subData), dim(pMat),"\n");
      if(is.null(pMat) || identical(all.equal(length(subData),0),TRUE) || identical(all.equal(ncol(subData),1),TRUE) ){
        return(NULL);
      }
      #identical(all.equal(length(subData),0),TRUE) || identical(all.equal(ncol(subData)
      valid.mat <- matrix((abs(subData) >= log2(input$vennFoldChange)) & (-log10(pMat) >= input$vennLog10Pvalue), nrow = nrow(subData), ncol=ncol(subData));
      colnames(valid.mat) <- colnames(subData);
      rownames(valid.mat) <- rownames(subData);
      valid.mat[is.na(valid.mat)] <- FALSE;
      zcols <- colorRampPalette(getRampCols(input$colourScheme))(ncol(subData));
      values$vennResult <- draw.venn(lapply(colnames(valid.mat),function(x){which(valid.mat[,x])}), euler.d = FALSE, scaled=FALSE,
                                     fill=zcols, add=TRUE, category=colnames(subData),cat.default.pos='text', percent=FALSE);
      cat(file="out.txt", dim(valid.mat), "\n", append=TRUE);
    }
  });

  output$vennInfo <- renderTable({
    subNames <- input$dataSubsets;
    if(is.null(subNames) || (length(subNames) < 2)){
      return(NULL);
    }
    subData <- getGeneSetMatrix(input, geneLimit=-1);
    pMat <- getPadjMatrix(input, subData);
    if(is.null(pMat) || identical(all.equal(length(subData),0),TRUE) || identical(all.equal(ncol(subData),1),TRUE) ){
      return(NULL);
    }
    valid.mat <- (abs(subData) >= log2(input$vennFoldChange)) & (-log10(pMat) >= input$vennLog10Pvalue);
    valid.mat[is.na(valid.mat)] <- FALSE;
    zcols <- colorRampPalette(getRampCols(input$colourScheme))(ncol(subData));
    sub.df <- values$vennResult;
    if(is.null(sub.df)){
      return(NULL);
    }
    res <- brushedPoints(sub.df, input$venn_brush, "x", "y", allRows=TRUE);
    genes <- as.data.frame(valid.mat);
    final <- NULL;
    for( colStr in res[res$selected_,"groups"] ) {
      colNumeric <- as.numeric(unlist(strsplit(colStr, ",")));
      if(length(colNumeric) < ncol(genes)){
        runningNames <- unique(rownames(genes)[which(apply(genes[,(colNumeric),drop=FALSE], 1, all)
                                                     & apply(!genes[,(-colNumeric),drop=FALSE], 1, all))]);
      } else{
        runningNames <- unique(rownames(genes)[which(apply(genes[,(colNumeric),drop=FALSE], 1, all))]);
      }
      final <- c(runningNames,final);
    }
    if(is.null(final)){
      return();
    }
    lFilt <- round(subData[final,,drop=FALSE],2);
    colnames(lFilt) <- paste0(colnames(lFilt),".L2FC");
    pFilt <- signif(pMat[match(toupper(final),toupper(rownames(pMat))),,drop=FALSE],4);
    colnames(pFilt) <- paste0(colnames(pFilt),".padj");
    lDir <- apply(lFilt,1,function(x){paste(ifelse(sign(x)==-1,"D","U"), collapse="")});
    pSig <- apply(pFilt,1,function(x){-log10(x) >= input$vennLog10Pvalue});
    lSig <- apply(lFilt,1,function(x){abs(x) >= log2(input$vennFoldChange)});
    bothSig <- apply(pSig & lSig,2,function(x){paste(ifelse(x,"S","-"), collapse="")});
    lFilt <- cbind(lFilt, pFilt, dir.sig=paste(lDir,bothSig,sep="."));
    lFilt <- lFilt[order(lFilt[,"dir.sig"], rownames(lFilt)),];
    values$outputTable <- lFilt;
    rownames(lFilt) <- sub("_.*$","",rownames(lFilt)); # remove Ensembl names from display
    return(lFilt);
  });

  output$vennSummary <- renderTable({
    subNames <- input$dataSubsets;
    if(is.null(subNames) || (length(subNames) < 2)){
      return(NULL);
    }
    subData <- getGeneSetMatrix(input, geneLimit=-1);
    pMat <- getPadjMatrix(input, subData);
    if(is.null(pMat) || identical(all.equal(length(subData),0),TRUE) || identical(all.equal(ncol(subData),1),TRUE) ){
      return(NULL);
    }
    valid.mat <- (abs(subData) >= log2(input$vennFoldChange)) & (-log10(pMat) >= input$vennLog10Pvalue);
    valid.mat[is.na(valid.mat)] <- FALSE;
    zcols <- colorRampPalette(getRampCols(input$colourScheme))(ncol(subData));
    sub.df <- values$vennResult;
    if(is.null(sub.df)){
      return(NULL);
    }
    anySig <- rownames(subData)[(apply(valid.mat,1,any))];
    lFilt <- round(subData[anySig,,drop=FALSE],2);
    colnames(lFilt) <- paste0(colnames(lFilt),".L2FC");
    logDebug(head(anySig));
    logDebug(head(pMat));
    pFilt <- signif(pMat[match(toupper(anySig),toupper(rownames(pMat))),,drop=FALSE],4);
    colnames(pFilt) <- paste0(colnames(pFilt),".padj");
    lDirSign <- apply(lFilt,1,sign);
    pSig <- apply(pFilt,1,function(x){-log10(x) >= input$vennLog10Pvalue});
    lSig <- apply(lFilt,1,function(x){abs(x) >= log2(input$vennFoldChange)});
    lDirSign[!(pSig & lSig)] <- 0;
    lDir <- apply(lDirSign,2,function(x){paste(ifelse(x==-1,"D",ifelse(x==1,"U","-")), collapse="")});
    bothSig <- apply(pSig & lSig,2,function(x){paste(ifelse(x,"S","-"), collapse="")});
    lFilt <- cbind(lFilt, pFilt, sig.dir=paste(bothSig,lDir,sep="."));
    lFilt <- lFilt[order(lFilt[,"sig.dir"], rownames(lFilt)),];
    tSig <- tapply(lFilt[,"sig.dir"],lFilt[,"sig.dir"],length);
    return(data.frame(row.names = names(tSig), count = tSig));
  });

  output$pcaView <- renderPlot({
    dsSelected <- input$datasetName;
    if(dsSelected != "VSTPk" && dsSelected != "ImmGen"){
      return(NULL);
    }
    values$outputTable <- NULL;
    subNames <- input$dataSubsets;
    if(is.null(subNames) || (length(subNames) < 3) ||
         is.null(input$pcaInput1real) || is.null(input$pcaInput2real)){
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), ann=FALSE, frame.plot=FALSE, axes=FALSE);
      text(0,0,"[Please select at least 3 Data Sets]", adj=c(0.5,0.5), offset=0);
      return(); 
    }
    subData <- as.data.frame(getGeneSetMatrix(input, geneLimit=-1));
    if(is.null(subData) || (length(subData) == 0) || (ncol(subData) < 3)){
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), ann=FALSE, frame.plot=FALSE, axes=FALSE);
      text(0,0,"[Please select at least 3 Data Sets]", adj=c(0.5,0.5), offset=0);
      return(); 
    }
    library(FactoMineR);
    subData <- subData[complete.cases(subData),];
    if(input$pcaDESelector != "All"){
      selectDSExp <- sub(":.*$","",input$pcaDESelector);
      selectDSTest <- sub("^.*?:","",input$pcaDESelector);
      ## get L2FC and P values associated with selected dataset
      DEData <- getInputMatrix(input, dsName=selectDSExp, subNames=selectDSTest)[,selectDSTest, drop=FALSE];
      DEpMat <- getPadjMatrix(input, DEData, dsName=selectDSExp, subNames=selectDSTest);
      filter.names <- rownames(DEData)[(abs(DEData[,1]) >= log2(input$pcaFoldChange)) &
                                        (-log10(DEpMat[,1]) >= input$pcaLog10Pvalue)];
      ## hack to work around not using full names
      filter.names <- unique(sub("_.*$","",filter.names));
      ## subset data on filtered list
      subData <- subData[toupper(rownames(subData)) %in% toupper(filter.names),];
    }
    ## only choose 500 most different genes
    #     if(input$geneSetClass == "All") {
    #       subData <- subData[head(order(-apply(subData,1,sd)),500),];
    #     }
    subData <- data.frame(t(subData));
    subData$Condition <- sub("(\\.r)?\\.r?[0-9]+$", "", row.names(subData));
    zcols <- vector();
    library(PKI);
    for( i in 1:length(subData$Condition)){
      zcols[i] <- hcl(h=round(as.numeric(PKI.digest(charToRaw(subData$Condition[i]))[1])/256*360), c=85, l=35);
    }
    res.pca <- PCA(subData, quali.sup=which(colnames(subData) %in% c("Condition")), graph=FALSE);
    # arrange table to be printed as a CSV
    values$outputTable <- res.pca$var$coord[order(apply(abs(res.pca$var$coord),1,max),decreasing=TRUE),]
    title <- ifelse(input$geneSetClass == "All", 
                    c("Individuals factor map (PCA) - All genes"), 
                    paste("Individuals factor map(PCA) -", input$geneCat, "-", input$geneList))
    if(input$replicateJoin == 'Combined'){
      plot(res.pca, choix="ind", axes=c(as.numeric(gsub("[^0-9]","",input$pcaInput1real)),
                                        as.numeric(gsub("[^0-9]","",input$pcaInput2real))), 
           habillage="ind", col.hab=zcols, invisible="quali", title = title);
    } else {
      plot(res.pca, choix="ind", axes=c(as.numeric(gsub("[^0-9]","",input$pcaInput1real)),
                                        as.numeric(gsub("[^0-9]","",input$pcaInput2real))), 
           habillage="ind", col.hab=zcols, title = title);
    }
  });

  output$heatMapKey <- downloadHandler(
    filename = function() {
      paste0('heatmap_',ifelse(input$customSetChosen,input$customSetName,input$geneList), '_', Sys.Date(), '.pdf');
    },
    content = function(file) {
      subData <- getGeneSetMatrix(input);
      if(is.null(subData)){
        return(NULL);
      }
      zlim <- getMatrixRange(input);
      zcols <- colorRampPalette(getRampCols(input$colourScheme))(64);
      breaks <- getColBreaks(zlim, 65, input$datasetName);
      pdf(file=file, width=1.25, height=1.25, pointsize=8);
      snif.heat.map.key(subData, zlim=zlim, zcols=zcols, xlab="", breaks=breaks);
      dev.off();
    },
    contentType = "application/pdf"
  );
  
  output$pdfVolcano <- downloadHandler(
    filename = function() {
      paste0('volcano_',paste(input$dataSubsets,collapse="_"), '_', Sys.Date(), '.pdf');
    },
    content = function(file) {
      logOutput(input, "volcano");
      subData <- getGeneSetMatrix(input, geneLimit=-1);
      if(is.null(subData)){
        return(NULL);
      }
      zlim <- getMatrixRange(input);
      zcols <- colorRampPalette(getRampCols(input$colourScheme))(64);
      pMat <- getPadjMatrix(input, subData);
      if(is.null(pMat)){
        return(NULL);
      }
      ## aggregate results to generate one value per gene
      if(ncol(subData) > 1){
        subData <- apply(subData,1,mean, na.rm=TRUE);
        pMat <- apply(pMat,1,function(x){10^mean(log10(x), na.rm=TRUE)});
      }
      xVals <- subData;
      yVals <- -log10(pMat);
      full.cases <- !is.na(xVals) & !is.na(yVals);
      xVals <- xVals[full.cases];
      yVals <- yVals[full.cases];
      over.limit.y <- (yVals > input$volcanoYLimit);
      over.limit.x <- (abs(xVals) > input$volcanoXLimit);
      yVals[over.limit.y] <- input$volcanoYLimit;
      xVals[over.limit.x] <- input$volcanoXLimit * sign(xVals[over.limit.x]);
      over.limit <- over.limit.y | over.limit.x;
      if(!is.null(pMat)){
        pdf(file=file, width=6, height=6, pointsize=9);
        plot(x=xVals, y=yVals, pch=ifelse(over.limit,24,16), cex=2,
             xlim=c(-input$volcanoXLimit,input$volcanoXLimit), ylim=c(0,input$volcanoYLimit),
             col=ifelse(abs(xVals)>=log2(input$ablineFoldChange),
                        ifelse(yVals>=input$ablineLog10Pvalue, "#DF272680", "#000000A0"),
                        ifelse(yVals>=input$ablineLog10Pvalue, "#F0909080", "#000000A0")), 
             xlab="log2FC", ylab=expression(-'log'[10](pAdj)));
        count.down <- sum((xVals <= -log2(input$ablineFoldChange)) & (yVals >= input$ablineLog10Pvalue));
        count.up <- sum((xVals >= log2(input$ablineFoldChange)) & (yVals >= input$ablineLog10Pvalue));
        text(x=min(xVals)*0.75, y=input$volcanoYLimit, labels=count.down, col="#DF2726", cex=3, pos=1);
        text(x=max(xVals)*0.75, y=input$volcanoYLimit, labels=count.up, col="#DF2726", cex=3, pos=1);
        abline(h=input$ablineLog10Pvalue, v=c(log2(input$ablineFoldChange),-log2(input$ablineFoldChange)), col="gray60", lty=2);
        showSourceFiles();
        dev.off();
      } else {
        warning("pMat is null");
      }      
    },
    contentType = "application/pdf"
  );

  output$pdfSmallVolcano <- downloadHandler(
    filename = function() {
      paste0('volcano_',paste(input$dataSubsets,collapse="_"), '_', Sys.Date(), '.pdf');
    },
    content = function(file) {
      logOutput(input, "small_volcano");
      subData <- getGeneSetMatrix(input, geneLimit=-1);
      if(is.null(subData)){
        return(NULL);
      }
      zlim <- getMatrixRange(input);
      zcols <- colorRampPalette(getRampCols(input$colourScheme))(64);
      pMat <- getPadjMatrix(input, subData);
      if(is.null(pMat)){
        return(NULL);
      }
      graphLab <- colnames(subData)[1];
      ## aggregate results to generate one value per gene
      if(ncol(subData) > 1){
        subData <- apply(subData,1,mean, na.rm=TRUE);
        pMat <- apply(pMat,1,function(x){10^mean(log10(x), na.rm=TRUE)});
        graphLab <- paste(colnames(subData), collapse=";");
      }
      xVals <- subData;
      yVals <- -log10(pMat);
      full.cases <- !is.na(xVals) & !is.na(yVals);
      xVals <- xVals[full.cases];
      yVals <- yVals[full.cases];
      over.limit.y <- (yVals > input$volcanoYLimit);
      over.limit.x <- (abs(xVals) > input$volcanoXLimit);
      yVals[over.limit.y] <- input$volcanoYLimit;
      xVals[over.limit.x] <- input$volcanoXLimit * sign(xVals[over.limit.x]);
      over.limit <- over.limit.y | over.limit.x;
      if(!is.null(pMat)){
        pdf(file=file, width=6, height=6, pointsize=9);
        par(mar=c(7,9,6,2), mgp = c(5,1.5,0), lwd=2, las=1);
        plot(x=xVals, y=yVals, pch=ifelse(over.limit,24,16), cex=2, main=graphLab, cex.main = 3,
             xlim=c(-input$volcanoXLimit,input$volcanoXLimit), ylim=c(0,input$volcanoYLimit),
             col=ifelse(abs(xVals)>=log2(input$ablineFoldChange),
                        ifelse(yVals>=input$ablineLog10Pvalue, "#DF272680", "#000000A0"),
                        ifelse(yVals>=input$ablineLog10Pvalue, "#F0909080", "#000000A0")), 
             xlab=expression('log'[2](FC)), ylab=expression(-'log'[10](pAdj)), cex.lab=3, cex.axis=2);
        #axis(side=1, at=seq(-20,20,by=2));
        count.down <- sum((xVals <= -log2(input$ablineFoldChange)) & (yVals >= input$ablineLog10Pvalue));
        count.up <- sum((xVals >= log2(input$ablineFoldChange)) & (yVals >= input$ablineLog10Pvalue));
        text(x=-input$volcanoXLimit*0.75, y=input$volcanoYLimit, labels=count.down, col="#DF2726", cex=3, pos=1, offset=1);
        text(x=input$volcanoXLimit*0.75, y=input$volcanoYLimit, labels=count.up, col="#DF2726", cex=3, pos=1, offset=1);
        abline(h=input$ablineLog10Pvalue, v=c(log2(input$ablineFoldChange),-log2(input$ablineFoldChange)),
               col="#80808080", lty=2);
        showSourceFiles();
        dev.off();
      } else {
        warning("pMat is null");
      }      
    },
    contentType = "application/pdf"
  );

  output$pdfMA <- downloadHandler(
    filename = function() {
      paste0('MA_',ifelse(input$customSetChosen,input$customSetName,input$geneList), '_', Sys.Date(), '.pdf');
    },
    content = function(file) {
      logOutput(input, "MA");
      subData <- getGeneSetMatrix(input, geneLimit=-1);
      if(is.null(subData)){
        return(NULL);
      }
      zlim <- getMatrixRange(input);
      # get median gene Values
      geneMedian <- apply(getInputMatrix(input, dsName="VSTPk", subNames=""),1,median);
      subData <- subData[rownames(subData) %in% names(geneMedian),,drop=FALSE];
      # retrieve p values
      pMat <- getPadjMatrix(input, subData);
      if(is.null(pMat)){
        return();
      }
      ## aggregate results to generate one value per gene
      if(ncol(subData) > 1){
        subData <- matrix(apply(subData,1,mean, na.rm=TRUE), ncol = 1, dimnames=list(rownames(subData),NULL));
        pMat <- matrix(apply(pMat,1,function(x){10^mean(log10(x), na.rm=TRUE)}), ncol=1, dimnames=list(rownames(subData),NULL));
      }
      #     dataCols <- ifelse(subData < zlim[1],1,
      #                        ifelse(subData > zlim[2],length(zcols),
      #                               (length(zcols)-1)*(subData-zlim[1])/(diff(zlim)) + 1));
      if(!is.null(pMat)){
        pdf(file=file, width=6, height=6, pointsize=9);
        yVals <- subData;
        xVals <- matrix(geneMedian[rownames(subData)],ncol = 1);
        pVals <- -log10(pMat);
        full.cases <- !is.na(xVals) & !is.na(yVals);
        xVals <- xVals[full.cases];
        yVals <- yVals[full.cases];
        pVals <- pVals[full.cases];
        plot(x=xVals, y=yVals, pch=16, cex=2, 
             col=ifelse(abs(yVals)>=log2(input$ablineFoldChange),
                        ifelse(pVals>=input$ablineLog10Pvalue, "#DF272680", "#000000A0"),
                        ifelse(pVals>=input$ablineLog10Pvalue, "#F0909080", "#000000A0")), 
             xlab="Median VSTPk Expression (all samples)", ylab="log2FC");
        abline(h=c(log2(input$ablineFoldChange),-log2(input$ablineFoldChange)), col="gray60", lty=2)
        showSourceFiles();
        dev.off();
      } else {
        warning("pMat is null");
      }
      
    },
    contentType = "application/pdf"
  );

  output$pdfVenn <- downloadHandler(
    filename = function() {
      paste0('venn_',ifelse(input$customSetChosen,input$customSetName,input$geneList), '_', Sys.Date(), '.pdf');
    },
    content = function(file) {
      logOutput(input, "venn");
      subData <- getGeneSetMatrix(input, geneLimit=-1);
      if(is.null(subData)){
        return(NULL);
      }
      pMat <- getPadjMatrix(input, subData);
      if(is.null(pMat)){
        return(NULL);
      }
      valid.mat <- (abs(subData) >= log2(input$vennFoldChange)) & (-log10(pMat) >= input$vennLog10Pvalue);
      zcols <- colorRampPalette(getRampCols(input$colourScheme))(ncol(subData));
      
      pdf(file=file,family='Helvetica', pointsize=8, useDingbats=FALSE);
      draw.venn(apply(valid.mat,2,which), fill=zcols, euler.d = FALSE,
                scaled=FALSE, add=TRUE, category=colnames(subData), percent=FALSE);
      showSourceFiles(grid=TRUE);
      dev.off();
    },
    contentType = "application/pdf"
  );  

  output$pdfHeatMap <- downloadHandler(
    filename = function() {
      paste0('heatmap_out_', Sys.Date(), '.pdf');
    },
    content = function(file) {
      logOutput(input, "heatmap");
      subData <- getGeneSetMatrix(input);
      if(is.null(subData)){
        return(NULL);
      }
      subNames <- colnames(subData);
      zlim <- getMatrixRange(input);
      zcols <- colorRampPalette(getRampCols(input$colourScheme))(64);
      breaks <- getColBreaks(zlim, 65, input$datasetName);
      pMat <- NULL;
      if(!is.null(input$showPadj) && input$showPadj){
        pMat <- getPadjMatrix(input, subData);
        ## remove Ensembl IDs
        rownames(pMat) <- sub("_.*$","", rownames(pMat));
      }
      ## remove Ensembl IDs from subData
      rownames(subData) <- sub("_.*$","", rownames(subData));
      if(input$mergeHM){
        gsmGeneNames <- rownames(subData);
        gsmSampNames <- colnames(subData);
        subData <- matrix(apply(subData,1,mean), ncol=1);
        rownames(subData) <- gsmGeneNames;
        colnames(subData) <- sprintf("Merged_%d", length(gsmSampNames));
        if(!is.null(pMat)){
          pMat <- matrix(apply(pMat, 1, min), ncol=1);
          rownames(pMat) <- gsmGeneNames;
          colnames(pMat) <- sprintf("Merged_%d", length(gsmSampNames));
        }
      }
      gringene.heat.map(subData, pmat = pMat,
                        longestLabel=subNames[which(nchar(subNames) == max(nchar(subNames)))[1]], 
                        filename=file, zlim=zlim, zcols=zcols, breaks=breaks, post.plot=showSourceFiles);
    },
    contentType = "application/pdf"
  );
  
  output$pdfPca <- downloadHandler(
    filename=function() {
      paste0('PCA_',ifelse(input$customSetChosen,input$customSetName,input$geneList),'_',Sys.Date(),'.pdf')
    },
    content=function(file) {
      logOutput(input, "PCA");
      require(FactoMineR);
      subData <- as.data.frame(getGeneSetMatrix(input, geneLimit=-1));
      if(is.null(subData) || (length(subData) == 0) || (ncol(subData) < 3)){
        return(); 
      }
      subData <- subData[complete.cases(subData),]
      if(input$pcaDESelector != "All"){
        selectDSExp <- sub(":.*$","",input$pcaDESelector);
        selectDSTest <- sub("^.*?:","",input$pcaDESelector);
        ## get L2FC and P values associated with selected dataset
        DEData <- getInputMatrix(input, dsName=selectDSExp, subNames=selectDSTest)[,selectDSTest, drop=FALSE];
        DEpMat <- getPadjMatrix(input, DEData, dsName=selectDSExp, subNames=selectDSTest);
        filter.names <- rownames(DEData)[(abs(DEData[,1]) >= log2(input$pcaFoldChange)) &
                                           (-log10(DEpMat[,1]) >= input$pcaLog10Pvalue)];
        ## hack to work around not using full names
        filter.names <- unique(sub("_.*$","",filter.names));
        ## subset data on filtered list
        subData <- subData[toupper(rownames(subData)) %in% toupper(filter.names),];
      }
      subData <- data.frame(t(subData));
      subData$Condition <- sub("(\\.r)?\\.r?[0-9]+$", "", row.names(subData));
      zcols <- vector();
      require(PKI);
      for( i in 1:length(subData$Condition)){
        zcols[i] <- hcl(h=round(as.numeric(PKI.digest(charToRaw(subData$Condition[i]))[1])/256*360), c=85, l=35);
      }
      res.pca <- PCA(subData, quali.sup=which(colnames(subData) %in% c("Condition")), graph=FALSE);
      title <- ifelse(input$geneSetClass == "All", 
                      c("Individuals factor map (PCA) - All genes"), 
                      paste("Individuals factor map(PCA) -", input$geneCat, "-", input$geneList))
      pdf(file=file, width=6, height=6, pointsize=9);
      if(input$replicateJoin == 'Combined'){
        plot(res.pca, choix="ind", axes=c(as.numeric(gsub("[^0-9]","",input$pcaInput1real)),
                                          as.numeric(gsub("[^0-9]","",input$pcaInput2real))), 
             habillage="ind", col.hab=zcols, invisible="quali", title=title);
      } else {
        plot(res.pca, choix="ind", axes=c(as.numeric(gsub("[^0-9]","",input$pcaInput1real)),
                                          as.numeric(gsub("[^0-9]","",input$pcaInput2real))), 
             habillage="ind", col.hab=zcols, title=title);
      }
      showSourceFiles();
      dev.off();      
    },
    contentType = "application/pdf"
  );

  output$pdfSmallPca <- downloadHandler(
    filename=function() {
      paste0('PCA_',ifelse(input$customSetChosen,input$customSetName,input$geneList),'_',Sys.Date(),'.pdf')
    },
    content=function(file) {
      logOutput(input, "smallPCA");
      require(FactoMineR);
      subData <- as.data.frame(getGeneSetMatrix(input, geneLimit=-1));
      if(is.null(subData) || (length(subData) == 0) || (ncol(subData) < 3)){
        return(); 
      }
      subData <- subData[complete.cases(subData),]
      if(input$pcaDESelector != "All"){
        selectDSExp <- sub(":.*$","",input$pcaDESelector);
        selectDSTest <- sub("^.*?:","",input$pcaDESelector);
        ## get L2FC and P values associated with selected dataset
        DEData <- getInputMatrix(input, dsName=selectDSExp, subNames=selectDSTest)[,selectDSTest, drop=FALSE];
        DEpMat <- getPadjMatrix(input, DEData, dsName=selectDSExp, subNames=selectDSTest);
        filter.names <- rownames(DEData)[(abs(DEData[,1]) >= log2(input$pcaFoldChange)) &
                                           (-log10(DEpMat[,1]) >= input$pcaLog10Pvalue)];
        ## hack to work around not using full names
        filter.names <- unique(sub("_.*$","",filter.names));
        ## subset data on filtered list
        subData <- subData[toupper(rownames(subData)) %in% toupper(filter.names),];
      }
      subData <- data.frame(t(subData));
      axes = axes=c(as.numeric(gsub("[^0-9]","",input$pcaInput1real)),
                    as.numeric(gsub("[^0-9]","",input$pcaInput2real)));
      ## Remove Ensembl ID
      colnames(subData) <- sub("_.*$","",colnames(subData));
      ## Clean data names
      rownames(subData) <- sub("^CD","",rownames(subData));
      rownames(subData) <- sub("^DBP\\+FITC\\.","DBP.",rownames(subData));
      rownames(subData) <- sub("(pool(_no)?_dep)\\.","\\1",rownames(subData));
      ## Determine data categories
      subData$pop <- sub("\\..*$","",rownames(subData));
      subData$rep <- sub("^.*\\.","",rownames(subData));
      subData$sample <- sub("^(.*)(\\.r)?\\..*$","\\1",rownames(subData));
      subData$allTreat <- sub("^.*?\\.(.*?)(\\.r)?\\.[^\\.]*$","\\1",rownames(subData));
      subData$marker <- sub("^.*\\.","", subData$allTreat);
      subData$treat <- sub("\\..*$","", subData$allTreat);
      subData$allTreat <- factor(subData$allTreat);
      ## carry out PCA
      out.pca <- PCA(subData,
                     quali.sup=which(colnames(subData)
                                     %in% c("pop","rep","allTreat",
                                            "marker","treat","sample")),
                     graph=FALSE);
      colorLookup <- if(length(unique(subData$pop)) > 1){
        unlist(list("DBP+FITC" = "#00B000B0",
                    "Nb" = "#00B900B0",
                    "UT" = "#008CE5B0",
                    "PBS" = "#008CE5B0"));
      } else {
        unlist(list("DBP+FITC" = "#00B000B0",
                    "UT" = "#00B900B0",
                    "Nb" = "#008CE5B0",
                    "PBS" = "#008CE5B0"));
      }
      colorLookup[unique(subData$pop)[!(unique(subData$pop) %in% colorLookup)]] <- "#808080";
      styleLookup <- unlist(list("Nb" = 21, "DBP+FITC" = 21,
                                 "Nb_AFPos" = 21, "NB_AFNeg" = 21,
                                 "PBS" = 22, "UT" = 22,
                                 "TN" = 24,
                                 "103" = 22, "CD103" = 22,
                                 "103Neg" = 22, "CD103Neg" = 22,
                                 "326" = 22, "CD326" = 22,
                                 "11b" = 21, "CD11b" = 21));
      styleLookup[unique(subData$pop)[!(unique(subData$pop) %in% styleLookup)]] <- 21;
      styleLookup[unique(subData$treat)[!(unique(subData$treat) %in% styleLookup)]] <- 21;
      ## produce file output
      pdf(file=file, width=6, height=6, pointsize=9);
      par(mfrow=c(2,2), mgp = c(1.5,1.5,0), mar=c(5,5,0.5,0.5));
      for(yf in c(-1, 1)){
        for(xf in c(-1, 1)){
          plot(out.pca$ind$coord[,axes[1]] * xf,
               out.pca$ind$coord[,axes[2]] * yf,
               xlim=sort(range(out.pca$ind$coord[,axes[1]]) * 1.2 * xf),
               ylim=sort(range(out.pca$ind$coord[,axes[2]]) * 1.2 * yf),
               xlab = sprintf("PC#%d (%0.2f%%)", axes[1],
                              out.pca$eig[["percentage of variance"]][axes[1]]),
               ylab = sprintf("PC#%d (%0.2f%%)", axes[2],
                              out.pca$eig[["percentage of variance"]][axes[2]]),
               pch = if(length(unique(subData$pop)) > 1){
                 ## distinguish populations (symbol) and treatments (colour)
                 styleLookup[subData[rownames(out.pca$ind$coord),"pop"]]
               } else {
                 ## distinguish preparation (symbol) and treatments (colour)
                 styleLookup[subData[rownames(out.pca$ind$coord),"treat"]]},
               bg = colorLookup[subData[rownames(out.pca$ind$coord),"treat"]],
               cex = 5, cex.lab = 3, axes=FALSE, frame.plot=TRUE,
               panel.first = abline(h=0,v=0,lty="dotted")
          );
          showSourceFiles();
          roundXLim <- round(signif(sort(range(out.pca$ind$coord[,axes[1]]) * 1.2), 2));
          roundYLim <- round(signif(sort(range(out.pca$ind$coord[,axes[2]]) * 1.2), 2));
          #axis(1,at=roundXLim * xf, labels=roundXLim);
          #axis(2,at=roundYLim * yf, labels=roundYLim);
          arrows(x0=sort(range(out.pca$ind$coord[,axes[1]]) * 0.9 * xf)[((xf+1)/2)+1],
                 x1=sort(range(out.pca$ind$coord[,axes[1]]) * 1.2 * xf)[((xf+1)/2)+1],
                 y0=0, lwd=2, length = 0.15);
          arrows(y0=sort(range(out.pca$ind$coord[,axes[2]]) * 0.9 * yf)[((yf+1)/2)+1],
                 y1=sort(range(out.pca$ind$coord[,axes[2]]) * 1.2 * yf)[((yf+1)/2)+1],
                 x0=0, lwd=2, length = 0.15);
          if(length(unique(subData$pop)) > 1){
            for(tv in unique(subData$pop)){
              text(out.pca$quali.sup$coord[tv,axes[1]] * xf,
                   out.pca$quali.sup$coord[tv,axes[2]] * yf, tv, cex=2);
            }
          } else {
            tv <- unique(subData$pop);
            text(out.pca$quali.sup$coord[tv,axes[1]] * xf,
                 out.pca$quali.sup$coord[tv,axes[2]] * yf, paste0("[",tv,"]"), cex=2);
          }
          if(length(unique(subData$treat)) > 1){
            for(tv in unique(subData$treat)){
              text(out.pca$quali.sup$coord[tv,axes[1]] * xf,
                   out.pca$quali.sup$coord[tv,axes[2]] * yf, tv, cex=2);
            }
          } else {
            tv <- unique(subData$treat);
            text(out.pca$quali.sup$coord[tv,axes[1]] * xf,
                 out.pca$quali.sup$coord[tv,axes[2]] * yf, paste0("[",tv,"]"), cex=2);
          }
        }
      }
      dev.off();      
    },
    contentType = "application/pdf"
  );

  output$downloadCSV <- downloadHandler(
    filename = function() { 
      paste0('selected_genes_', input$datasetName, '.csv', sep='') 
    },
    content = function(file) {
      logOutput(input, "CSVdownload");
      if(input$outputPanel == "heatmap"){
        subData <- getGeneSetMatrix(input, geneLimit = -1);
        if(input$datasetName %in% c("VSTPk", "ImmGen")){
          write.csv(subData, file);
        } else {
          pMat <- getPadjMatrix(input, subData);
          colnames(subData) <- paste0(colnames(subData),".L2FC");
          colnames(pMat) <- paste0(colnames(pMat),".padj");
          combined <- cbind(round(subData,2), signif(pMat,4));
          write.csv(combined[,rep(1:ncol(subData),each=2)+c(0,ncol(subData))], file);
        }
      } else {
        outData <- values$outputTable;
        if(any(grepl("_",rownames(outData)))){
          outData <- cbind(geneName=sub("_.*$","",rownames(outData)),ensemblID=sub("^.*?_","",rownames(outData)),outData);
          write.csv(outData, file, row.names=FALSE);
        } else {
          write.csv(outData, file);
        }
      }
    },
    contentType = "text/csv"
  )

#   output$selectButton <- renderUI({
#     actionButton("selectAll", label=values$selectButtonText);
#   });

  output$pcaInput1 <- renderUI({
    if(!is.null(values$outputTable)) {
      dsName <- input$datasetName;
      if(is.null(dsName)){
        return(NULL);
      }     
      subNames <- input$dataSubsets;
      if( length(subNames) > 2) {
        selectInput("pcaInput1real", "x Principal Component", 
                    choices=c(paste0('Component ',1:max(5,length(subNames)-1))), selected=input$pcaInput1real);
      }
    }
  })
  
  output$pcaInput2 <- renderUI({
    if(!is.null(values$outputTable)) {
      dsName <- input$datasetName;
      if(is.null(dsName)){
        return(NULL);
      }     
      subNames <- input$dataSubsets;
      if( length(subNames) > 2) {
        choice <- c(paste0('Component ',1:max(5,length(subNames)-1)));
        selectInput("pcaInput2real", "y Principal Component", choices=choice, 
                    selected=ifelse(is.null(input$pcaInput2real), choice[2], input$pcaInput2real));
      }
    }
  })

  output$logStdout <- renderText({
    vOut <- 
    if(values$logPos > 1){
      c(values$logText[values$logPos:100], values$logText[1:(values$logPos-1)]);
    } else {
      values$logText;
    }
    vOut <- vOut[vOut != ""];
    return(paste(sprintf("%3d: %s",seq_along(vOut), vOut), collapse=""));
  })

  ################ UI Updating functions

  updateDataSets <- function(input){
    dsSelected <- input$datasetName;
    ## data set name
    if(input$outputPanel %in% c("volcano", "MA", "venn")){
      if(is.null(dsSelected) || (dsSelected %in% c("VSTPk", "ImmGen"))){
        dsSelected <- dataDirs[1];
        values$selectButtonText <- "Select All";
      }
      updateRadioButtons(session, "datasetName", "Data Set", choices=dataDirs, selected=dsSelected);
    }
    if (input$outputPanel == "pca"){
      if(is.null(dsSelected) || (dsSelected != "VSTPk" && dsSelected != "ImmGen")){
        if(is.null(dsSelected)){
          dsSelected <- "VSTPk";
        }
        values$selectButtonText <- "Select All";
      }
      updateRadioButtons(session, "datasetName", "Data Set", choices=c("VSTPk","ImmGen"), selected=dsSelected);
    }
    if (input$outputPanel == "heatmap"){
      if(is.null(dsSelected)){
        dsSelected <- "VSTPk";
      }
      updateRadioButtons(session, "datasetName", "Data Set", choices=c("VSTPk", "ImmGen", dataDirs), selected=dsSelected);
      if(!(dsSelected %in% c("VSTPk","ImmGen"))){
        if(input$colourScheme == "white/red"){
          updateRadioButtons(session, "colourScheme", selected="blue/white/red");
        }
        updateRadioButtons(session, "geneOrder", choices=c("Default","Fold Change","Cluster"));
      }
      if(dsSelected %in% c("VSTPk","ImmGen")){
        if(input$colourScheme != "white/red"){
          updateRadioButtons(session, "colourScheme", selected="white/red");
        }
        updateRadioButtons(session, "geneOrder", choices=c("Default","Expression","Cluster"));
      }
    }
    ## also update data subsets (if necessary)
    subNames <- input$dataSubsets;
    if (dsSelected == "VSTPk" || dsSelected == "ImmGen"){
      subChoices <- values$allSubNames;
      if(dsSelected == "ImmGen"){
        subChoices <- values$allImmgenSubNames;
      }
      joinReplicates <- input$replicateJoin;
      if(!is.null(joinReplicates) && (joinReplicates == "Combined") && (dsSelected == "VSTPk")){
        subChoices <- subChoices[!(subChoices %in% outliers)];
        subChoices <- unique(sub("(\\.r)?\\.r?[0-9]+$","",subChoices));
      }
      if(!all(subNames %in% subChoices)){
        subNames <- NULL;
        values$selectButtonText <- "Select All";
      }
      updateSelectizeInput(session, "dataSubsets", choices=sortConds(subChoices), selected = subNames);
    } else {
      dataFile <- list.files(sprintf("RNASeq_Files/DESeq2/%s", dsSelected),
                             pattern=sprintf("^(David|Alex)_%s_All.*.csv.gz$", dsSelected), full.names=TRUE);
      data.topLine <- scan(dataFile, what=character(), nlines=1, quiet=TRUE, sep=",")[-1];
      if(!all(subNames %in% data.topLine)){
        subNames <- NULL;
        values$selectButtonText <- "Select All";
      }
      updateSelectizeInput(session, "dataSubsets", choices=sortConds(data.topLine), selected=subNames);
    }
  }

  ############# Observer events

  ## Change data sets whenever tab is changed
  ## Change data subsets whenever data set name is changed
  ## Change data subsets whenever replicate join is changed
  observeEvent(paste(input$outputPanel,input$datasetName,input$replicateJoin), updateDataSets(input));

# select/deselect all using action button
  observeEvent((input$selectAll), {
    buttonName <- values$selectButtonText;
    if(!buttonName %in% c("Select All", "Deselect")){
      return;
    }
    # Update selection after click
    if (buttonName == "Select All"){
      values$selectButtonText <- "Deselect";
    } else {
      values$selectButtonText <- "Select All";
    }
    dsName <- input$datasetName;
    if (dsName == "VSTPk" || dsName == "ImmGen"){
      subNames <- values$allSubNames;
      if(dsName == "ImmGen"){
        subNames <- values$allImmgenSubNames;
      }
      if(!is.null(input$replicateJoin) && (input$replicateJoin == "Combined") && (dsName == "VSTPk")){
        subNames <- subNames[!(subNames %in% outliers)];
        subNames <- unique(sub("(\\.r)?\\.r?[0-9]+$","",subNames));
      }
      if (buttonName == "Select All"){
        if(input$replicateJoin == "Split"){
          subNames <- subNames[!(subNames %in% outliers)];
        }
        updateSelectizeInput(session=session, inputId="dataSubsets", selected = subNames);
      } else {
        updateSelectizeInput(session=session, inputId="dataSubsets", selected = character(0));
      }
    } else {
      dataFile <- list.files(sprintf("RNASeq_Files/DESeq2/%s", input$datasetName),
                             pattern=sprintf("^(David|Alex)_%s_All.*.csv.gz$",input$datasetName), full.names=TRUE);
      data.topLine <- scan(dataFile, what=character(), nlines=1, quiet=TRUE, sep=",")[-1];
      if (buttonName == "Select All"){
        updateSelectizeInput(session=session, inputId="dataSubsets", selected = data.topLine);
      } else {
        updateSelectizeInput(session=session, inputId="dataSubsets", selected = character(0));
      }
    }
  });

  observeEvent(values$scatterWholeData, {
    if(!is.null(values$scatterWholeData)){
      updateSelectizeInput(session, inputId="scatterGeneNames", choices = rownames(values$scatterWholeData), server=TRUE, selected=NULL);
    } else {
      updateSelectizeInput(session, inputId="scatterGeneNames", choices = "[choose data subsets]", selected=NULL);
    }
  });

  ## observer to update gene list choices when list category is changed
  observeEvent(input$geneCat, {
    updateSelectInput(session=session, inputId="geneList", 
                      choices=unique(geneLists.df[geneLists.df$listCategory == input$geneCat,"listName"]));    
  });

  ## Set up variables for application
  values$vstpkFileName <- 
    "RNASeq_Files/DESeq2/published_lowFiltered_vstpk_DESeq2_blind_local_all_2016-May-05.csv";
  values$immgenFileName <- "RNASeq_Files/immgen/immgen_DC_geneLevel_transformed_2016-Jul-21.csv";
  values$allSubNames <- tail(scan(isolate(values$vstpkFileName), what=character(), nlines=1, quiet=TRUE, sep=","),-1);
  values$allImmgenSubNames <- tail(scan(isolate(values$immgenFileName), what=character(), nlines=1, quiet=TRUE, sep=","),-1);

  values$logText <- rep("",100); ## set up log buffer of last 100 messages
  values$logPos <- 1;
}

shinyApp(ui = ui, server = server);
