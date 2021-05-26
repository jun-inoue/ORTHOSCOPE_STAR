library(ape)

args   <- commandArgs()
MyTree <- args[6]
outfile <- args[7]

MyTree   <- read.tree(MyTree)

rootBS100 <- MyTree
rootBS100$node.label[2]="100"
write.tree(rootBS100, file=outfile)
