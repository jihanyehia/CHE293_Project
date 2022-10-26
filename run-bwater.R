source("water.R")
#------------------------------
# Read and process a .hb2 file
#------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	#struct <- '1n78.hb2'
  struct <- '78-r5-90k-140k-merged-1032'
} else if (length(args)>=1) {
	struct <- args[1]
}
print(struct)
hbonds <- process_hb2(struct)
print(hbonds[1:5,])
#---------------------------------------------------------------
# Retrieve all hydrogen bonds that involve single-water bridges
#---------------------------------------------------------------
swb <- find.swb(hbonds)
print(swb)

dwb <- find.dwb(hbonds)
print(dwb)

# print(hb.network(hbonds))
sum(nucle.aacid(hbonds))
