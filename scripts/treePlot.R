library(ape)
args                             <- commandArgs()
name_summaryFile                 <- args[6]      # 100_1stAnalysisSummary.txt
Gene_tree_newick                 <- args[7]      # SpeciesTree 1st_gene_tree_newick 2nd_rearranged_gene_tree_newick
Rearranged_gene_tree_newick      <- args[8]      # dummy_rearranged_species_tree_newick 1st_rearranged_gene_tree_newick 2nd_rearranged_gene_tree_newick
thick_branch_species_name_line   <- args[9]      # Rooting_4_2ndTree Rooting_4_1stTree Orthogroup
root_species_name_line           <- args[10]     # Rooting_4_2ndTree Rooting_4_1stTree
outfileName                      <- args[11]     # speciesTree 115_1st 240_2nd

#print("# thick_branch_species_name_line")
#print(thick_branch_species_name_line)
#print("# root_species_name_line")
#print(root_species_name_line)
#print("outfileName")
#print(outfileName)
#q()

#print("thick_branch_species_name_line")
#print(thick_branch_species_name_line)

################# Node name
#nodeNameLabel_change_swich <- "on"
nodeNameLabel_change_swich <- "off"
#########################

greenPrefixes <- c(); purplePrefixes <- c(); orangePrefixes <- c(); magentaPrefixes <-c(); bluePrefixes <- c()
redPrefixes <- c();
OrthogroupBasalNode <- "";

#######
# df_infile$V1 に全行が入っている前提
parse_records <- function(lines) {
  records <- list()
  current_header <- NULL
  current_lines  <- character()

  for (ln in lines) {
    if (grepl("^>", ln)) {
      # 直前のレコードを保存
      if (!is.null(current_header)) {
        records[[current_header]] <- current_lines
      }
      # 新しいレコード開始
      current_header <- sub("^>", "", ln)  # ">" を外す
      current_lines  <- character()
    } else {
      # 本文行を追加（空行も許容したければ条件調整）
      current_lines <- c(current_lines, ln)
    }
  }
  # 最後のレコードを保存
  if (!is.null(current_header)) {
    records[[current_header]] <- current_lines
  }
  records
}

# 使い方

extract_section_lines <- function(keyword_header_line)
{
    #print("### extract_section_lines () ###")
    vector_record <- c()
    frag <- 0
    vector_section_lines <- c()
    for(line in df_infile$V1) {
        #print(line)
    
        ### Collect lines 
        if (regexpr('^>', line) < 0) {
            if(frag == 1) {
                vector_section_lines <- c(vector_section_lines, line)
            }
        }
      
       if (regexpr('^>', line) > 0)
       {
          if(frag == 1)
          {
              vector_record <- c(vector_record, vector_section_lines)
              vector_section_lines <- c()
              break
            }
    
          #keyWord1 <- paste('>',　keyword_header_line,　sep='')
          if (regexpr(keyword_header_line, line) > 0)
          {
              vector_record <- c(vector_record, line)
              frag <- 1
            }
    
        }
      }

    vector_record <- c(vector_record, vector_section_lines)
    #print("# vector_record")
    #print(vector_record)
    vector_section_lines <- vector_record[-1]
    #print("# vector_section_lines")
    #print(vector_section_lines)
    #print("# q70")
    #q()
    return(vector_section_lines)
}


get_section_lines_by_header <- function (keyword_name_line)
{
    #print("### get_section_lines_by_header() ###")
    #print("# keyword_name_line")
    #print(keyword_name_line)
    keyword_name_line <- paste(keyword_name_line, "( |$)", sep = "")
    keyword_name_line <- paste(">", keyword_name_line, sep = "")
    #print(keyword_name_line)
    #print("### df_infile$V1 START ###")
    #print(df_infile$V1)
    #print("### df_infile$V1 END ###")
  
    if(any(i <- grep(keyword_name_line, df_infile$V1)))
    {
        #print("# Found keyword_name_line")
        #print(keyword_name_line)
        #print("")
        value_matched <- extract_section_lines(keyword_name_line)
        #print("# value_matched")
        #print(value_matched)
        #print("### q line 90")
        #q()
    } else {
        #print("# Not Found keyword_name_line")
        #print(keyword_name_line)
        #print("")
        # print (paste('  ',keyword_name_line,' does not exist.', sep=''))
        value_matched  <- NULL
    }
    
    #print("value_matched")
    #print(value_matched)

    return(value_matched)
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
      #for(p in 1:length(tr$node.label)){
      #  if (regexpr(OrthogroupBasalNode, tr$node.label[p]) > 0){
      #    BSvalueColors <- c(BSvalueColors, "blue")
      #  } else {
      #    BSvalueColors <- c(BSvalueColors, 1)    
      #  }
      #}
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
  #print("### numbering_edgeWidth() ###")
  #print("# leafNames_for_thickBranch")
  #print(leafNames_for_thickBranch)
  nums_thickBranchLeaves <- c()
  for(i in 1:length(tr$tip.label)){
    #print("# tr$tip.label[i]")
    #print(tr$tip.label[i])
    for(leafName in leafNames_for_thickBranch){
      #print("# leafName")
      #print(leafName)
      if (tr$tip.label[i] == leafName){
        #print("hit")
        nums_thickBranchLeaves <- c(nums_thickBranchLeaves, i)
        break  # ← 内側 for を抜ける（外側は続く）
      } else {
        NULL
        #print("not")
      }
    }
    #print("")
  }

  #print("nums_thickBranchLeaves")
  #print(nums_thickBranchLeaves)
  #q()

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
  #print("### PNG_treeDrawing() ###")
  #print("prefix")
  #print(prefix)
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
      #print("1111")
      #plot(tr, label.offset = 2, no.margin=TRUE, underscore = TRUE, use.edge.length=TRUE, cex = 0.9, font = tipFontNums, tip.col = tipColorNums, edge.width = edgeWidths_for_eachBranch)
      #tiplabels(pch = tipLabelPCH, col = tipLabelColor, adj = -0.01, cex = 1.5)
      tiplabels(pch = tipLabelPCH, col = tipLabelColor, cex = 1.5)
  } else {
      #print("22222")
      #plot(tr,                   no.margin=TRUE, underscore = TRUE, use.edge.length=TRUE, cex = 0.9, font = tipFontNums, tip.col = tipColorNums, edge.width = edgeWidths_for_eachBranch)
      #tiplabels(pch = tipLabelPCH, col = tipLabelColor, cex = 1.5)
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

PDF_treeDrawing <- function(tr, prefix)
{
  pdf.file <- paste(outfileName, prefix, sep = "")
  #print("pdf.file")
  #print(pdf.file)
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

df_infile <- read.table(name_summaryFile, na.strings = FALSE, sep = '\t')
#print("is.data.frame(df_infile)")
#is.data.frame(df_infile)

#vector_queries <- get_section_lines_by_header("QuerySequence")
#print("# vector_queries")
#print(vector_queries)
#print("# q 494")

records <- parse_records(df_infile$V1)
# いま生成済みの records からヘッダ一覧を確認
headers <- names(records)
print("# headers (names(records))")
print(headers)

vector_queries <- records[["QuerySequence"]]
#print("# vector_queries")
#print(vector_queries)

Rooting_species = ""
#if (thick_branch_species_name_line == "Rooting_4_2ndTree"){
#    Rooting_species <- get_section_lines_by_header("Rooting_4_2ndTree")
#} else if (thick_branch_species_name_line == "Rooting_4_1stTree"){
#    Rooting_species <- get_section_lines_by_header("Rooting_4_1stTree")
#}

#print("# root_species_name_line")
#print(root_species_name_line)
#Rooting_species <- get_section_lines_by_header(root_species_name_line)
#print("# Rooting_species")
#print(Rooting_species)
Rooting_species <- records[[root_species_name_line]]
#print("# Rooting_species")
#print(Rooting_species)
#print("q507")
#q()


taxonSampling_color <- c()
if (Gene_tree_newick == "SpeciesTree"){
    #taxonSampling_color <- get_section_lines_by_header("TaxonSampling")
    taxonSampling_color <- records[["TaxonSampling"]]
}else{
    #taxonSampling_color <- get_section_lines_by_header("TaxonSampling_color")
    taxonSampling_color <- records[["TaxonSampling_color"]]
}
#print("# taxonSampling_color")
#print(taxonSampling_color)
greenPrefixes   = make_colorPrefixes(taxonSampling_color, "Green")
purplePrefixes  = make_colorPrefixes(taxonSampling_color, "Purple")
orangePrefixes  = make_colorPrefixes(taxonSampling_color, "Orange")
magentaPrefixes = make_colorPrefixes(taxonSampling_color, "Magenta")
bluePrefixes    = make_colorPrefixes(taxonSampling_color, "Blue")
redPrefixes     = make_colorPrefixes(taxonSampling_color, "Red")
#print("magentaPrefixes")
#print(magentaPrefixes)
#q()

queryNames <- c()
for (line in vector_queries)
{
  line <- sub(' +.*$', "", line)   
  queryNames <- c(queryNames, line)
}

#leafNames_for_thickBranch <- get_section_lines_by_header(thick_branch_species_name_line)
leafNames_for_thickBranch <- records[[thick_branch_species_name_line]]
#print("# leafNames_for_thickBranch")
#print(leafNames_for_thickBranch)
#print("q 555")
#q()

Rearrangement_BS_value_threshold <- c()
#Rearrangement_BS_value_threshold <- get_section_lines_by_header("Rearrangement_BS_value_threshold")
Rearrangement_BS_value_threshold <- records[["Rearrangement_BS_value_threshold"]]

####
#print("481")
#q()

if (Gene_tree_newick == "SpeciesTree"){
    #print("sssss")
    #print("Gene_tree_newick")
    #print(Gene_tree_newick)
    #q()
    df_infile <- read.table(name_summaryFile, na.strings = FALSE, sep = '\t')
    #Species_tree <- get_section_lines_by_header(Gene_tree_newick)
    Species_tree <- records[[Gene_tree_newick]]
    Species_tree <- read.tree(text = Species_tree)
    #print("Species_tree")
    #print(Species_tree)
    if (is.null(Species_tree)){
        print("No Species_tree")
        q()
    }

    #vector_queries <- get_section_lines_by_header("QuerySequence")
    vector_queries <- records[["QuerySequence"]]

    #Rooting <- get_section_lines_by_header("Rooting")
    Species_tree <- ladderize(Species_tree, TRUE)
    #leafNames_for_thickBranch <- get_section_lines_by_header(thick_branch_species_name_line)
    edgeWidths_for_eachBranch <- numbering_edgeWidth(Species_tree)
    tipFontNums  <- fontNumChange(Species_tree)
    tipColorNums <- tipColorChange(Species_tree)
    #print("tipColorNums")
    #print(tipColorNums)

    nodeLabelFontColorNums <- c()
    #OrthogroupBasalNode <- get_section_lines_by_header("OrthogroupBasalNode")
    OrthogroupBasalNode <- records[["OrthogroupBasalNode"]]

    nodeLabelFontColorNums   <- BScolorChange(Species_tree)

    querySpecies <- get_section_lines_by_header("QuerySpecies")
    querySpecies <- records[["QuerySpecies"]]

    Num_allQueries <- c()
    Num_1stQuery <- queryNameInversion(Species_tree, querySpecies)

    PDF_treeDrawing(Species_tree, prefix=".pdf")

    q()
}
#print("1111")
#q()
######################################################################################################
#print("529")

Gene_tree <- get_section_lines_by_header(Gene_tree_newick)
Gene_tree <- records[[Gene_tree_newick]]

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
#print("555")
Rearranged_gene_tree <- get_section_lines_by_header(Rearranged_gene_tree_newick)
Rearranged_gene_tree <- records[[Rearranged_gene_tree_newick]]
if(is.null(Rearranged_gene_tree)) {
  print("# No Rearranged_gene_tree. Stopped.")
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

Num_rooting_species                <- pickUp_leafNum(Rearranged_gene_tree, Rooting_species)
#print("# Num_rooting_species")
#print(Num_rooting_species)
tipLabelPCH                        <- rep(1, length(Rearranged_gene_tree$tip.label))
tipLabelPCH[Num_rooting_species]   <- 16
tipLabelColor                      <- rep("white", length(Rearranged_gene_tree$tip.label))
tipLabelColor[Num_rooting_species] <- "Black"

#Num_allQueries <- queryNameInversion(Rearranged_gene_tree, queryNames)
Num_1stQuery <- queryNameInversion(Rearranged_gene_tree, queryNames[1])

edgeWidths_for_eachBranch <- numbering_edgeWidth(Rearranged_gene_tree)

#txt.file <- paste(outfileName, "Rearranged_geneTree.txt", sep = "")
#write.tree(Rearranged_gene_tree, file=txt.file)

PNG_treeDrawing(Rearranged_gene_tree, prefix="Rearranged_geneTree.png")
PDF_treeDrawing(Rearranged_gene_tree, prefix="Rearranged_geneTree.pdf")
