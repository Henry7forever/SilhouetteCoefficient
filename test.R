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
  clusters = clusterTable[,1]
  couplingVector = NULL
  for (theNode in 1:nrow(data)) {
    belongCluster = findCluster(clusters,theNode)
    otherCluster = clusters[-belongCluster]
    distance = NULL
    
    if(length(otherCluster) > 0)
    {
      for (i_otherC in 1:length(otherCluster)) {
        deCluster = otherCluster[[i_otherC]]
        distanceACluster = NULL
        for (index in 1:length(deCluster)) {
          thatNode = -deCluster[index] 
          distanceACluster = c(distanceACluster ,data[theNode,thatNode])
          
        }
        avgDistanceACluster = mean(distanceACluster ,na.rm = T)
        distance = c(distance , avgDistanceACluster)
      }
      theCoupling = min(distance ,na.rm = T)
    }
    else
    {
      theCoupling = 0
    }
    couplingVector = c(couplingVector ,theCoupling)
  }
  return(couplingVector)
}
caculateCohesion <- function(clusterTable , data)
{
  clusters = clusterTable[,1]
  cohesionVector = NULL
  for (theNode in 1:nrow(data)) {
    belongClusterI = findCluster(clusters , theNode)
    belongCluster = clusters[[belongClusterI]]
    
    distance = NULL
    
    if(length(belongCluster) > 1)
    {
      for (i in 1:length(belongCluster)) {
        thatNode = -belongCluster[i]
        distance = c(distance , data[theNode,thatNode])
      }
      theCohesion = mean(distance , na.rm = T)
    }
    else
    {
      theCohesion = 0
    }
    cohesionVector = c(cohesionVector , theCohesion)
  }
  return(cohesionVector)
}


data_original = read.table("input.csv" , sep = ",", header = T)
data_original
data = as.dist(data_original)
data
hc = hcluster(data)
hc
plot(hc)
threshold = 15

sliceLv =  which((hc$height) < threshold)
sliceLv

hc$merge
obj = hc$labels
obj = as.factor(obj)

clusters = list()
lastMerge = NULL

#initialization
for (i in 1:length(hc$labels)) {
  clusters = c(clusters , list(-i))
  lastMerge[i]= 0
}

for (i in 1:length(sliceLv)) {
  merging = (hc$merge)[i, ]
  
  
  if (merging[1] > 0)
  {
    whichClusterA = which(lastMerge == merging[1])
  }
  else if (merging[1] < 0)
  {
    whichClusterA = which(clusters %in% merging[1] )
  }
  if (merging[2] > 0)
  {
    whichClusterB = which(lastMerge == merging[2])
  }
  else if (merging[2] < 0)
  {
    whichClusterB = which(clusters %in% merging[2] )
  }
  #whichClusterA = which(lastMerge == merging[1])
  #whichClusterB = which(lastMerge == merging[2])
  clusters[[whichClusterA]] = c(clusters[[whichClusterA]] , clusters[[whichClusterB]])
  clusters = clusters[-whichClusterB]
  lastMerge[whichClusterA] = i
  lastMerge = lastMerge[-whichClusterB]
}


clusterTable = cbind(clusters , lastMerge)





clusterTable

hc$merge

clusterTable[,1]

b=caculateCoupling(clusterTable ,data_original)
a=caculateCohesion(clusterTable ,data_original)

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
sc = (b-a) / larger
sc[sc %in% c(1,-1)] = 0
sc

output=NULL
for (i in 1:length(sc)) {
  output = cbind( output , sc[i] )
}
colnames(output) = colnames(data_original) 
output = data.frame(output)
output

write.table(output , file = "output.csv" , append = F , sep = "," , row.names = F)