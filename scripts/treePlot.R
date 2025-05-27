library(ape)
args                        <- commandArgs()
name_summaryFile            <- args[6]      # 100_1stAnalysisSummary.txt
Gene_tree_newick            <- args[7]      # 115_1st
Rearranged_gene_tree_newick <- args[8]
groupName_for_highlight     <- args[9]
outfileName                 <- args[10]      # 115_1st

#print("name_summaryFile")
#print(name_summaryFile)
#print("outfileName")
#print(outfileName)
#q()

#print("groupName_for_highlight")
#print(groupName_for_highlight)

################# Node name
#nodeNameLabel_change_swich <- "on"
nodeNameLabel_change_swich <- "off"
#########################

greenPrefixes <- c(); purplePrefixes <- c(); orangePrefixes <- c(); magentaPrefixes <-c(); bluePrefixes <- c()
redPrefixes <- c();
OrthogroupBasalNode <- "";

#######
line.picker <- function(keyWord)
{
  container <- c()
  frag <- 0
  lineStock <- c()
  for(line in infile$V1) {
    #print(line)

    ### Collect lines 
    if (regexpr('^>', line) < 0) {
      if(frag == 1) {
        lineStock <- c(lineStock, line)
      }
    }
  
   if (regexpr('^>', line) > 0)
   {
      if(frag == 1)
      {
        container <- c(container, lineStock)
        lineStock <- c()
        break
      }

      #keyWord1 <- paste('>',　keyWord,　sep='')
      if (regexpr(keyWord, line) > 0)
      {
        container <- c(container, line)
        frag <- 1
      }

    }
  }
  container <- c(container, lineStock)    
  container <- container[-1]
  return(container)
}


preab.sub <- function (keyWordA)
{
  #print("keyWordA")
  #print(keyWordA)
  keyWordA <- paste(keyWordA, "( |$)", sep = "")
  keyWordA <- paste(">", keyWordA, sep = "")
  #print(keyWordA)
  #print("### infile$V1 START ###")
  #print(infile$V1)
  #print("### infile$V1 END ###")

  if(any(i <- grep(keyWordA, infile$V1)))
  {
    #print("Found keyWordA")
    #print(keyWordA)
    #print("")
    containerA <- line.picker(keyWordA)
  } else {
    #print("Not Found keyWordA")
    #print(keyWordA)
    #print("")
    # print (paste('  ',keyWordA,' does not exist.', sep=''))
    containerA  <- NULL
  }
  return(containerA)
}


fontNumChange <- function (tr)
{
  tipFontNums              <- rep(1, length(tr$tip.label))
  #fontNum[queryIDNum]     <- 4
  return (tipFontNums)
}


queryNameInversion <- function (tr, queryNames)
{
  queryTipNums <- c()
  for(i in 1:length(tr$tip.label)){
    for (queryName in queryNames){
      if(regexpr(queryName, tr$tip.label[i]) > 0){
        #print(queryName)
        #print(tr$tip.label[i])
        #print("")
        queryTipNums <- c(queryTipNums, i)
      }
    }
  }
  return(queryTipNums)
}


pickUp_leafNum <- function (tr, queryNames)
{
  queryTipNums <- c()
  for(i in 1:length(tr$tip.label)){
    for (queryName in queryNames){
      if(regexpr(queryName, tr$tip.label[i]) > 0){
        #print(queryName)
        #print(tr$tip.label[i])
        #print("")
        queryTipNums <- c(queryTipNums, i)
      }
    }
  }
  return(queryTipNums)
}


tipColorChange <- function(tr)
{
  orangeNum <- c(); blueNum <- c(); redNum <- NULL; greenNum <- c(); magentaNum <- c(); purpleNum <- c(); humanNum <- c()

  for(i in 1:length(tr$tip.label)){

    for (redPrefix in redPrefixes){
      if(regexpr(redPrefix, tr$tip.label[i]) > 0){
        redNum <- c(redNum,i)
      }
    }

    for (greenPrefix in greenPrefixes){
      if(regexpr(greenPrefix, tr$tip.label[i]) > 0){
        greenNum <- c(greenNum,i)
      }
    }

    for (purplePrefix in purplePrefixes){
      if(regexpr(purplePrefix, tr$tip.label[i]) > 0){
        purpleNum <- c(purpleNum,i)
      }
    }
  
    for (orangePrefix in orangePrefixes){
      if(regexpr(orangePrefix, tr$tip.label[i]) > 0){
        orangeNum <- c(orangeNum,i)
      }
    }
  
    for (magentaPrefix in magentaPrefixes){
      if(regexpr(magentaPrefix, tr$tip.label[i]) > 0){
        magentaNum <- c(magentaNum,i)
      }
    }  
  
    for (bluePrefix in bluePrefixes){
      #print(bluePrefix)
      if(regexpr(bluePrefix, tr$tip.label[i]) > 0){
        blueNum <- c(blueNum,i)
      }
    }
  
  }

  tipColorNums             <- rep("black",length(tr$tip.label))
  tipColorNums[redNum]     <- "red"
  tipColorNums[greenNum]   <- "darkgreen"
  tipColorNums[purpleNum]  <- "purple"
  tipColorNums[orangeNum]  <- "darkorange1"
  tipColorNums[magentaNum] <- "hotpink2"
  tipColorNums[blueNum]    <- "blue"

  return(tipColorNums)
}


make_colorPrefixes <- function (taxonSampling_colorTMP, colorFN)
{
  colorPrefixes <- c()
  for (line in taxonSampling_colorTMP)
  {
    #print("line")
    #print(line)
    #print("")
    #spPrefixFN <- sub(' +.*$', "", line)   
    spPrefixFN <- sub('_.*$', "", line)
    if(regexpr(colorFN, line) > 0){
        colorPrefixes <- c(colorPrefixes, spPrefixFN)
    }
  }
  return(colorPrefixes)
}


nodeNameLabel_change <- function(tr)
{
  for(p in 1:length(tr$node.label)){
    #print(tr$node.label[p])
    if(regexpr('D=N', tr$node.label[p])> 0){
      tr$node.label[p] <- sub('_.*$', "", tr$node.label[p])
    } else {
      tr$node.label[p] <- sub('_.*$', "D", tr$node.label[p])
    }
    #print(tr$node.label[p])
  }
  return(tr)
}


BScolorChange <- function(tr)
{
  #print("Species_tree")
  #print(Species_tree)

  startChr <- "r"

  BSvalueColors <- NULL
  if (Gene_tree_newick == "SpeciesTree"){
      #print("OrthogroupBasalNode")
      #print(OrthogroupBasalNode)
      for(p in 1:length(tr$node.label)){
        if (regexpr(OrthogroupBasalNode, tr$node.label[p]) > 0){
          BSvalueColors <- c(BSvalueColors, "blue")
        } else {
          BSvalueColors <- c(BSvalueColors, 1)    
        }
      }
  } else {
      for(p in 1:length(tr$node.label)){
        if (regexpr("r", tr$node.label[p]) > 0){
          BSvalueColors <- c(BSvalueColors, 2)
        #} else if (as.numeric(tr$node.label[p]) < as.numeric(Rearrangement_BS_value_threshold)){
        #  BSvalueColors <- c(BSvalueColors, 2)
        } else {
          BSvalueColors <- c(BSvalueColors, 1)    
        }
      }
  }
  

#  if (is.null(Rearrangement_BS_value_threshold)){
#  } else {
#    for(p in 1:length(tr$node.label)){
#      if (tr$node.label[p] == "r"){
#      } else {   
#        #print(tr$node.label[p])
#        #print(Rearrangement_BS_value_threshold)
#        #print("\n")
#        #print("tr$node.label[p]")
#        #print(tr$node.label[p])
#        if (as.numeric(tr$node.label[p]) < as.integer(Rearrangement_BS_value_threshold)){
#          BSvalueColors <- c(BSvalueColors, 2) 
#        }
#      }
#    }
#  }


  return(BSvalueColors)
}


numbering_edgeWidth <- function (tr)
{

  nums_thickBranchLeaves <- c()
  for(i in 1:length(tr$tip.label)){
    for(leafName_for_thickBranch in leafNames_for_thickBranch){
      if (tr$tip.label[i] == leafName_for_thickBranch){
        nums_thickBranchLeaves <- c(nums_thickBranchLeaves, i)
      }
    }
  }

  edgeWidths_for_eachBranch <- NULL
  if(is.null(nums_thickBranchLeaves)){
    edgeWidths_for_eachBranch <- rep(1.5, dim(tr$edge)[1])
  } else {
    edgeWidths_for_eachBranch <- rep(1.5, dim(tr$edge)[1])
    nums_thickBranchLeaves <- c(nums_thickBranchLeaves, Num_1stQuery)
    wh_members <- which.edge(tr, nums_thickBranchLeaves)
    edgeWidths_for_eachBranch[wh_members] <- 4
  }
  return(edgeWidths_for_eachBranch)
}


PNG_treeDrawing <- function (tr, prefix)
{  
  png.file <- paste(outfileName, prefix, sep = "")
  pngWidth <- NULL

  pngHeight <- NULL
  if (length(tr$tip.label) > 200) {
    pngWidth <- 1500
    pngHeight = 2700
  } else if (length(tr$tip.label) > 100) {
    pngWidth <- 1200
    pngHeight = 1800
  } else if (length(tr$tip.label) > 50) {
    pngWidth <- 1000
    pngHeight = 1200
  } else if (length(tr$tip.label) > 10) {
    pngWidth <- 1000
    pngHeight = 900
  } else {
    pngWidth <- 800
    pngHeight = 600
  }

  #print("edgeWidths_for_eachBranch")
  #print(edgeWidths_for_eachBranch)
  png(png.file, width = pngWidth, height = pngHeight)
  plot(tr,                   no.margin=TRUE, underscore = TRUE, use.edge.length=TRUE, cex = 0.9, font = tipFontNums, tip.col = tipColorNums, edge.width = edgeWidths_for_eachBranch)
  if(regexpr("Rearranged_geneTree", prefix) > 0){
      #plot(tr, label.offset = 2, no.margin=TRUE, underscore = TRUE, use.edge.length=TRUE, cex = 0.9, font = tipFontNums, tip.col = tipColorNums, edge.width = edgeWidths_for_eachBranch)
      #tiplabels(pch = tipLabelPCH, col = tipLabelColor, adj = -0.01, cex = 1.5)
      tiplabels(pch = tipLabelPCH, col = tipLabelColor, cex = 1.5)
  } else {
      #plot(tr,                   no.margin=TRUE, underscore = TRUE, use.edge.length=TRUE, cex = 0.9, font = tipFontNums, tip.col = tipColorNums, edge.width = edgeWidths_for_eachBranch)
      add.scale.bar()
  }
  
  #####nodelabels(tr$node.label, adj = c(1.2,-0.5), frame = "n", font = nodeLabelFontNums, cex=nodeLabelFontSizeNums, col = nodeLabelFontColorNums)
  #adj_nodelabel = 0
  #if(regexpr("Rearranged_geneTree", prefix) > 0){
  #    adj_nodelabel = -0.18
  #}

  nodelabels(tr$node.label, adj = c(-0.05,0.5), frame = "n", col = nodeLabelFontColorNums)

  #if(!is.null(Num_allQueries)){
  #  #tiplabels(tr$tip.label[Num_allQueries], Num_allQueries, cex=1.0, adj = adj_nodelabel, bg = "gray40", col="white")
  #  tiplabels(tr$tip.label[Num_allQueries], Num_allQueries, cex=1.0, adj = 0, bg = "gray40", col="white")
  #}
  if(!is.null(Num_1stQuery)){
    #tiplabels(tr$tip.label[Num_1stQuery],   Num_1stQuery,   cex=1.0, adj = adj_nodelabel, bg = "navyblue", col="white")
    tiplabels(tr$tip.label[Num_1stQuery],   Num_1stQuery,   cex=1.0, adj = 0, bg = "navyblue", col="white")
  }
  dev.off()
}

PDF_treeDrawing <- function (tr, prefix)
{
  pdf.file <- paste(outfileName, prefix, sep = "")
  pdfWidth  <- NULL
  pdfHeight <- NULL
  if (length(tr$tip.label) > 200) {
    pdfWidth  = 38
    pdfHeight = 28
  } else if (length(tr$tip.label) > 100) {
    pdfWidth  = 31
    pdfHeight = 21
  } else if (length(tr$tip.label) > 50) {
    pdfWidth  = 24
    pdfHeight = 14
  } else if (length(tr$tip.label) > 10) {
    if (Gene_tree_newick == "SpeciesTree"){
        pdfWidth  = 7
        pdfHeight = 7
    } else {
        pdfWidth  = 25
        pdfHeight = 10
    }
  } else {
    if (Gene_tree_newick == "SpeciesTree"){
        pdfWidth  = 7
        pdfHeight = 15
    } else {
        pdfWidth  = 15
        pdfHeight = 7
    }
  }

  pdf(pdf.file, width = pdfWidth, height = pdfHeight)
  plot      (tr, no.margin=TRUE, underscore = TRUE, use.edge.length=TRUE, cex = 0.9, font = tipFontNums, tip.col = tipColorNums, edge.width = edgeWidths_for_eachBranch)
  nodelabels(tr$node.label, adj = c(-0.05,0.5), frame = "n", col = nodeLabelFontColorNums)
  #if(!is.null(Num_allQueries)){
  #  tiplabels (tr$tip.label[Num_allQueries], Num_allQueries, cex=1.0, adj = 0, bg = "gray40", col="white")
  #}
  if(!is.null(Num_1stQuery)){
    tiplabels (tr$tip.label[Num_1stQuery],   Num_1stQuery, cex=1.0,   adj = 0, bg = "navyblue", col="white")
  }
  if(regexpr("gene_tree", prefix) > 0){
      add.scale.bar()
  }
  dev.off()
}

##################################################################


infile                      <- read.table(name_summaryFile, na.strings = FALSE, sep = '\t')
Querys_used_in_the_analysis <- preab.sub("QuerySequence")
Rooting                     <- preab.sub("Rooting")

taxonSampling_color <- c()
if (Gene_tree_newick == "SpeciesTree"){
    taxonSampling_color <- preab.sub("TaxonSampling")
}else{
    taxonSampling_color <- preab.sub("TaxonSampling_color")
}
greenPrefixes   = make_colorPrefixes(taxonSampling_color, "Green")
purplePrefixes  = make_colorPrefixes(taxonSampling_color, "Purple")
orangePrefixes  = make_colorPrefixes(taxonSampling_color, "Orange")
magentaPrefixes = make_colorPrefixes(taxonSampling_color, "Magenta")
bluePrefixes    = make_colorPrefixes(taxonSampling_color, "Blue")
redPrefixes     = make_colorPrefixes(taxonSampling_color, "Red")

queryNames <- c()
for (line in Querys_used_in_the_analysis)
{
  line <- sub(' +.*$', "", line)   
  queryNames <- c(queryNames, line)
}

leafNames_for_thickBranch <- preab.sub(groupName_for_highlight)

Rearrangement_BS_value_threshold <- c()
Rearrangement_BS_value_threshold <- preab.sub("Rearrangement_BS_value_threshold")


####
if (Gene_tree_newick == "SpeciesTree"){
    #print("sssss")
    #print("Gene_tree_newick")
    #print(Gene_tree_newick)
    #q()
    infile <- read.table(name_summaryFile, na.strings = FALSE, sep = '\t')
    Species_tree <- preab.sub(Gene_tree_newick)
    Species_tree <- read.tree(text = Species_tree)
    #print("Species_tree")
    #print(Species_tree)
    if (is.null(Species_tree)){
        print("No Species_tree")
        q()
    }

    Querys_used_in_the_analysis <- preab.sub("QuerySequence")
    Rooting <- preab.sub("Rooting")
    Species_tree <- ladderize(Species_tree, TRUE)
    leafNames_for_thickBranch <- preab.sub(groupName_for_highlight)
    edgeWidths_for_eachBranch <- numbering_edgeWidth(Species_tree)
    tipFontNums  <- fontNumChange(Species_tree)
    tipColorNums <- tipColorChange(Species_tree)
    #print("tipColorNums")
    #print(tipColorNums)
    #q()

    nodeLabelFontColorNums <- c()
    OrthogroupBasalNode <- preab.sub("OrthogroupBasalNode")
    nodeLabelFontColorNums   <- BScolorChange(Species_tree)

    querySpecies <- preab.sub("QuerySpecies")
    Num_allQueries <- c()
    Num_1stQuery <- queryNameInversion(Species_tree, querySpecies)

    PDF_treeDrawing(Species_tree, prefix=".pdf")

    q()
}
######################################################################################################

Gene_tree <- preab.sub(Gene_tree_newick)
if (is.null(Gene_tree)){
    print("No Gene_tree")
    q()
}

Gene_tree <- read.tree(text = Gene_tree)
Gene_tree <- ladderize(Gene_tree, TRUE)
Gene_tree$edge.length[Gene_tree$edge.length<0]<-0   ### nagative branch length, replace with 0

tipFontNums             <- fontNumChange(Gene_tree)
tipColorNums            <- tipColorChange(Gene_tree)
#nodeLabelFontNums      <- rep(1,length(Gene_tree$tip.label))
#nodeLabelFontSizeNums  <- rep(0.9, length(Gene_tree$tip.label))
nodeLabelFontColorNums  <- BScolorChange(Gene_tree)

Num_allQueries         <- queryNameInversion(Gene_tree, queryNames)
Num_1stQuery            <- queryNameInversion(Gene_tree, queryNames[1])

edgeWidths_for_eachBranch <- numbering_edgeWidth(Gene_tree)

PNG_treeDrawing(Gene_tree, prefix="GeneTree.png")
PDF_treeDrawing(Gene_tree, prefix="GeneTree.pdf")


######################################################################################################
Rearranged_gene_tree <- preab.sub(Rearranged_gene_tree_newick)
if(is.null(Rearranged_gene_tree)) {
  q()
}

Rearranged_gene_tree <- read.tree(text = Rearranged_gene_tree)
Rearranged_gene_tree <- ladderize(Rearranged_gene_tree, TRUE)

tipFontNums  <- fontNumChange(Rearranged_gene_tree)
tipColorNums <- tipColorChange(Rearranged_gene_tree)

nodeLabelFontColorNums <- c()
if (nodeNameLabel_change_swich == "on")
{
  Rearranged_gene_tree     <- nodeNameLabel_change(Rearranged_gene_tree)
  nodeLabelFontColorNums   <- BScolorChange(Rearranged_gene_tree)
}

Num_rooting                <- pickUp_leafNum(Rearranged_gene_tree, Rooting)
tipLabelPCH                <- rep(1, length(Rearranged_gene_tree$tip.label))
tipLabelPCH[Num_rooting]   <- 16
tipLabelColor              <- rep("white", length(Rearranged_gene_tree$tip.label))
tipLabelColor[Num_rooting] <- "Black"

#Num_allQueries <- queryNameInversion(Rearranged_gene_tree, queryNames)
Num_1stQuery   <- queryNameInversion(Rearranged_gene_tree, queryNames[1])

edgeWidths_for_eachBranch <- numbering_edgeWidth(Rearranged_gene_tree)

#txt.file <- paste(outfileName, "Rearranged_geneTree.txt", sep = "")
#write.tree(Rearranged_gene_tree, file=txt.file)

PNG_treeDrawing(Rearranged_gene_tree, prefix="Rearranged_geneTree.png")
PDF_treeDrawing(Rearranged_gene_tree, prefix="Rearranged_geneTree.pdf")
