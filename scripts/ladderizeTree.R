library(ape)

args <- commandArgs()
MyTree <- args[6]
outfile <- args[7]
#print("MyTree")
#print(MyTree)
#q()

## ladderize from bottom
tr <- read.tree(MyTree)
### Rooting species branch off from bottom right
tr.ladderrize <- ladderize(tr, FALSE)
### Rooting species branch off from top left
#tr.ladderrize <- ladderize(tr)
write.tree(tr.ladderrize, file=outfile)

