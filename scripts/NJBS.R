library(ape)

args            <- commandArgs()
infile          <- args[6]
outGroup        <- args[7]
outfile.nwk     <- args[8]

set.seed(123)

#print ('  ## 085_NJBSa.R starts.')

#print('infile')
#print(infile)
#print("")
#print('outGroup')
#print(outGroup)
#print("")
#print('outfile.nwk')
#print(outfile.nwk)
#print("")

#  print (paste('  infile', infile))

#outfile.nwk <- paste(outfile, '.nwk', sep = '')

## Open the 'infile' file
if(file.access(infile) != 0){
  print (paste(infile,' does not exist.'))
} else {
  infile.phy  <- read.dna(infile)
}

## Estimate distance
dist.TN93 <- dist.dna(infile.phy, model = "TN93", pairwise.deletion=TRUE, gamma=5)

## #Estimate NJ tree
nj.TN93 <- njs(dist.TN93)

## Reroot by outGroup
nj.TN93 <- root(nj.TN93,outGroup,r=T)

## ladderize from bottom
nj.TN93 <- ladderize(nj.TN93,TRUE)

## BS analysis with partitioning data (1st+2nd)
#nj.boot.nj.TN93 <- boot.phylo(nj.TN93, infile.phy,     function(xx) root(nj(dist.dna(xx, model = "TN93", pairwise.deletion = TRUE, gamma = 5)),outGroup,r=T), 100, 2)
#nj.boot.nj.TN93 <- boot.phylo(nj.TN93, infile.phy,     function(xx) root(njs(dist.dna(xx, model = "TN93", pairwise.deletion = TRUE)), outGroup,r=T), 100, 2)
nj.boot.nj.TN93 <- boot.phylo(nj.TN93, infile.phy,     function(xx) root(njs(dist.dna(xx, model = "TN93", pairwise.deletion = TRUE)), outGroup,r=T), 100)
nj.TN93$node.label <- nj.boot.nj.TN93

## Change bs calues: NA -> 0
#for(p in 1:length(nj.TN93$node.label)){
#  print(p)
#  print(nj.TN93$node.label[p])
#}
#cat('\n\n')

nj.TN93$node.label[is.na(nj.TN93$node.label)]<-0

#for(p in 1:length(nj.TN93$node.label)){
#  print(p)
#  print(nj.TN93$node.label[p])
#}
#cat('\n')

## Write tree
write.tree(nj.TN93, file = outfile.nwk)

