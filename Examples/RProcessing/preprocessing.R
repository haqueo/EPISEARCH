source("./random.genetic.dataset.R")
require(cluster)
require(mcclust)
require(infotheo)

# Create the dataset and write it 
D <- random.genetic.dataset(300,100,3)
write.table(D,"./genetic.txt",sep="\t",row.names=FALSE,col.names = FALSE,quote=FALSE)


# Create the VI Filtering file
viDists <- c()
for (i in 1:99){
  viDists <- c(viDists,vi.dist(D[,i],D[,"Class"]))
}
sorted.dat <- dat[,order(viDists)]
write(x = sort(viDists), file = "../vifilter/filterfile.txt",sep="\n")
write.table(sorted.dat,"../vifilter/sortedviGenetic.txt",sep="\t",row.names=FALSE,col.names = FALSE,quote=FALSE)

# Perform PAM clustering, then create the clutering files
diss <- matrix(nrow=99,ncol=99)
for (i in 1:99){
  for (j in 1:i){
    diss[i,j] = vi.dist(D[,i],D[,j])
    diss[j,i] = diss[i,j]
  }
}
results <- pam(x=diss,diss=TRUE,k=10)
clusterIDs <- sort(results$clustering) - 1
clusterSizes <- c()
for (i in 0:9){
  clusterSizes <- c(clusterSizes,sum(clusterIDs == i))
}
write(x = clusterIDs, file = "../clusterfilter/clusterIDs.txt",sep="\n")
write(x = clusterSizes, file = "../clusterfilter/clusterSizes.txt",sep="\n")
dat.clust.sort <- D[,order(results$clustering)]
write.table(sorted.dat,"../clusterfilter/sortedclusterGenetic.txt",sep="\t",row.names=FALSE,col.names = FALSE,quote=FALSE)



