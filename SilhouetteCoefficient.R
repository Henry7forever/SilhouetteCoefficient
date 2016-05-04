#find which cluster dose the element belong
findCluster <- function(clusters , n)
{
  for (i in 1:length(clusters)) {
    if ((-n) %in% clusters[[i]])
    {
      return(i)
    }
  }
}


caculateCoupling <- function(clusterTable , data)
{
  clusters = clusterTable[, 1]
  #"couplingVector" for collecting the coupling values
  couplingVector = NULL
  for (theNode in 1:nrow(data)) {
    belongCluster = findCluster(clusters, theNode)
    otherCluster = clusters[-belongCluster]
    #distance for collecting corresponding distance between each clusters another 
    #and the node
    distance = NULL
    
    if (length(otherCluster) > 0)
    {
      #iterate all clusters
      for (i_otherC in 1:length(otherCluster)) {
        deCluster = otherCluster[[i_otherC]]
        #iterate all elements in deCluster
        #and record the distance between theNode and each node in deCluster
        distanceACluster = NULL
        for (index in 1:length(deCluster)) {
          thatNode = -deCluster[index]
          distanceACluster = c(distanceACluster , data[theNode, thatNode])
        }
        #do average in each cluster distance record
        avgDistanceACluster = mean(distanceACluster , na.rm = T)
        distance = c(distance , avgDistanceACluster)
      }
      theCoupling = min(distance , na.rm = T)
    }
    else   #only one cluster 
    {
      theCoupling = 0
    }
    couplingVector = c(couplingVector , theCoupling)
  }
  return(couplingVector)
}


caculateCohesion <- function(clusterTable , data)
{
  clusters = clusterTable[, 1]
  cohesionVector = NULL
  for (theNode in 1:nrow(data)) {
    belongClusterI = findCluster(clusters , theNode) #belonging cluster index
    belongCluster = clusters[[belongClusterI]]
    
    distance = NULL
    if (length(belongCluster) > 1)
    {
      for (i in 1:length(belongCluster)) {
        thatNode = -belongCluster[i]
        distance = c(distance , data[theNode, thatNode])
      }
      theCohesion = mean(distance , na.rm = T)
    }
    else    #only one element in the cluster
    {
      theCohesion = 0
    }
    cohesionVector = c(cohesionVector , theCohesion)
  }
  return(cohesionVector)
}

#read data from file
data_original = read.table("input.csv" , sep = ",", header = T)
data_original
#transfer to dist format for hcluster algorithm
data = as.dist(data_original)
data

hc = hcluster(data)
hc
plot(hc)


threshold = 15
#clustering under the threshold
sliceLv =  which((hc$height) < threshold)
sliceLv

hc$merge
obj = hc$labels
obj = as.factor(obj)

#'cluster' for collecting each cluster by hcluster's result
#last merge for recording the each clusters' last merge step
clusters = list()
lastMerge = NULL

#initialization : each object are in the cluster of single themselve
for (i in 1:length(hc$labels)) {
  clusters = c(clusters , list(-i))
  lastMerge[i] = 0
}

#step by step merging
for (i in 1:length(sliceLv)) {
  merging = (hc$merge)[i,]
  if (merging[1] > 0)
  {
    whichClusterA = which(lastMerge == merging[1])
  }
  else if (merging[1] < 0)
  {
    #should be modified by using findCluster()
    whichClusterA = which(clusters %in% merging[1])
  }
  if (merging[2] > 0)
  {
    whichClusterB = which(lastMerge == merging[2])
  }
  else if (merging[2] < 0)
  {
    #should be modified by using findCluster()
    whichClusterB = which(clusters %in% merging[2])
  }
  #refactoring the record
  clusters[[whichClusterA]] = c(clusters[[whichClusterA]] , clusters[[whichClusterB]])
  clusters = clusters[-whichClusterB]
  lastMerge[whichClusterA] = i
  lastMerge = lastMerge[-whichClusterB]
}
#combine the clusters and lastMerge into one table
clusterTable = cbind(clusters , lastMerge)

clusterTable
hc$merge
clusterTable[, 1]

b = caculateCoupling(clusterTable , data_original)
a = caculateCohesion(clusterTable , data_original)

larger = c()
for (i in 1:length(a)) {
  if (a[i] > b[i]) {
    larger[i] = a[i]
  }
  else {
    larger[i] = b[i]
  }
}
larger

#doing the fomula
sc = (b - a) / larger
sc[sc %in% c(1, -1)] = 0
sc

#prepare the data and write into the file
output = NULL
for (i in 1:length(sc)) {
  output = cbind(output , sc[i])
}
colnames(output) = colnames(data_original)
output = data.frame(output)
output

write.table(
  output ,
  file = "output.csv" ,
  append = F ,
  sep = "," ,
  row.names = F
)



