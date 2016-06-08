#######################################################################################
## We create a third function EstClust which takes a similarity matrix Sim and runs
## MaxIt consecutive regulated random walks (1000 by default) with the node removal
## parameter m (5 by default).  EstClust reads the profile of recurrence times stored
## in Recur for each walk and selects the top 5% of recurrence times as Spikes which
## signify probably transitions into previously unexplored clusters of nodes.  The
## record of the random walk stored in I is then segmented based on the occurrence of
## the spikes, with the locations of the segment boundaries stored in the vector Clust.
## Each node is assigned to a segment if the walk visited that node at least m/2
## times while the walk was within that segment.  The assignments are stored in the
## vector Assign.  The symmetric matrix EmpSim stores in element (i,j) the proportion
## of times nodes i and j were assigned to the same segment over MaxIt random walks.
## The function returns EmpSim.
#######################################################################################


EstClust <- function(Sim, MaxIt=1000, m=5) {
  N <- nrow(Sim)  ## N is the number of nodes in the data set

  ## Initialize the N x N matrix Ensemble with 0s for every element.  Initialize
  ## Iter at 1.

  Ensemble <- matrix(0, nrow=N, ncol = N)

  Iter <- 1

  ## Continue until MaxIt iterations are completed.

  while(Iter <= MaxIt) {

    ## Conduct a regulated random walk based on Sim with M=m.  Store the resulting
    ## recurrence time profile as Rec.

    Series <- MakeSeries(Sim, M=m)
    Rec <- Series$Recur

    ## Select the locations of the outliers of recurrence times and store them as the
    ## vector Spikes.  Make sure the first spike is placed at the beginning, and make
    ## sure it is not redundant.  Then pass through the elements of Spikes and remove
    ## any spike which occurs immediately after the preceding spike, so that there is
    ## always some space between spikes in the recurrence times.
    ## Concatenate Spikes with the location of the end of the recurrence time
    ## profile, i.e., put a final spike at the end of the profile.


    Spikes <- unique(c(1, which(Rec >= quantile(Rec, 0.95))))
    i <- 2
    while(i <= length(Spikes)) {
      if((Spikes[i] - Spikes[i-1]) == 1)  Spikes <- Spikes[-i]
      i <- i + 1
    }
    Spikes <- c(Spikes, (length(Series$Recur)+1))

    ## Create the vector Clust the same length as Spikes, and set its first element = 1.
    ## For each successive element of Clust, compute the number of steps in the random
    ## walk that occurred prior to the corresponding spike.  The number of steps is the
    ## sum of the recurrence times from one spike to the next, plus the number of steps
    ## prior to that (stored in Clust[i-1]).  Clust now stores the number of steps in the
    ## random walk between each successive spike in the recurrence times.

    Clust <- rep(0, length(Spikes))
    Clust[1] <- 1
    for(i in 2:length(Spikes))
      Clust[i] <- sum(Series$Recur[Spikes[(i-1)]:(Spikes[i]-1)]) + Clust[i-1]

    ## Create the vector Assign of length N.  For each data node i, check each segment j of
    ## the random walk, as delineated by Clust, and compute the number of times node i
    ## was visited during the random walk within segment j.  Store that number in Temp.
    ## If this number is at least m/2, assign node i to segment j

    SpikesRmv <- cumsum(Rec)[Spikes]
    SpikesRmv[length(Spikes)] <- length(Series$I)

    Assign <- rep(0, N)

    for(k in 1:N){
      FirstVst <- which(Series$I == k)[1]
      Assign[k] <- which(SpikesRmv >= FirstVst)[1]
    }

    ## For each pair of nodes i and j, if both nodes are assigned to the same segment,
    ## increment elements (i,j) and (j,i) of the Ensemble matrix by 1.  Then increment
    ## Iter by 1.

    for(i in 2:N) {
      for(j in 1:(i-1)) {
        if(Assign[i] == Assign[j]) {
          Ensemble[i,j] <- Ensemble[i,j] + 1
          Ensemble[j,i] <- Ensemble[j,i] + 1
        }
      }
    }

    Iter <- Iter + 1

  } ## end of while loop

  ## Divide all elements of Ensemble by MaxIt, so that element (i,j) is the proportion
  ## of times nodes i and j were assigned to the same segment of the random walk over
  ## a total of MaxIt random walks.  Return this matrix.

  Ensemble <- Ensemble/MaxIt
  rownames(Ensemble) <- rownames(Sim)
  colnames(Ensemble) <- colnames(Sim)
  return(Ensemble)

} ## end of function EstClust
