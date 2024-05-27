library(ape)

args <- commandArgs()
MyTree <- args[6]
outfile <- args[7]
#print("MyTree")
#print(MyTree)
#q()

## ladderize from bottom
tr <- read.tree(MyTree)
tr.ladderrize <- ladderize(tr, FALSE)
#tr.ladderrize <- ladderize(tr)
write.tree(tr.ladderrize, file=outfile)

