#############################################################################
## We create another function MakeSeries which carries out a regulated random
## walk on the nodes of the data set based on the similarity matrix W and a
## parameter M (whose default is 5), which is the maximum number of times
## any node can be visited during the random walk.
## After the random walk, the function returns a vector of recurrence
## times called Recur, which in its rth position gives the number of time
## steps between the (r-1)st node removal and the rth node removal (the
## 0th node removal is just the beginning of the random walk).
## The function also returns a vector I which is a record of the nodes
## visited during the random walk in sequential order.
#############################################################################

# called in EstClust

MakeSeries <- function(W, M=5) {

  N <- nrow(W)  ## N is the number of data nodes

  ## Compute the degree of each node as the sum of its row in W, and
  ## store the results in the vector Degree.
  ## Create the diagonal matrix D^{-1} whose diagonal elements are the
  ## reciprocals of the degrees of the nodes, and store it as Dinv.
  ## Create the transition matrix Q = Dinv * W.  All entries in Q are
  ## between 0 and 1, with Q_{ij} being the probability of moving
  ## from node i to node j during the random walk.

  Degree <- rowSums(W)
  Dinv <- diag(1/Degree)
  Q <- Dinv %*% W

  ## Vstd is a vector whose vth element is the number of times node v
  ## has been visited during the random walk.  Initially all elements
  ## equal zero.
  ## Vec is a vector which shows the nodes that are available for
  ## visitation at any point in the walk.  Initially nodes 1:N are
  ## available.

  Vstd <- rep(0,N)
  Vec <- 1:N

  ## We select the starting node from among those nodes which have the
  ## greatest degree, as follows:
  ## First, sort the degrees in descending order and store the sorted
  ## vector as Deg.
  ## Next find the largest i such that the sum of the first i elements
  ## of Deg is less than half of the sum of all the degrees.
  ## Then choose only those nodes whose degrees are no less than the
  ## ith element of Deg, and store their indices in the vector Top.
  ## Next, create the corresponding vector Probs as the ratio of the degree
  ## of each node in Top to the sum of all the degrees of the nodes in Top.
  ## Then, from among those nodes in Top randomly select one as the
  ## starting node, with the probability of choosing the jth node of Top
  ## given by Probs_j, and store this as the first element of I.
  ## Finally, the element of Vstd corresponding to the starting node is
  ## updated to equal 1, since that node has now been visited once.

  Deg <- sort(Degree, decreasing=TRUE)
  Sum <- 0
  i <- 1
  while(Sum < 0.5*sum(Degree)) {  # the goal is to find elements of top degrees.
    Sum <- sum(Deg[1:i])
    i <- i + 1
  }
  Top <- which(Degree >= Deg[i])
  Probs <- Degree[Top]/sum(Degree[Top])
  I <- sample(x=Vec[Top], size=1)
  Vstd[I[1]] <- 1

  ## Initialize L (the length of I), Recur (the vector of recurrence times),
  ## and Count (the number of steps since the last node removal).

  L <- 1
  Recur <- vector("numeric")
  Count <- 1

  ## Begin regulated random walk.  Continue walk until every node has been
  ## visited at least once, except one.

  while(length(which(Vstd == 0)) > 1) {

    ## Select those nodes which have not yet been eliminated from the walk.  Store
    ## their indices in Wh.  These are the nodes eligible for visitation.

    Wh <- which(Vec > 0)

    ## If the row of the transition matrix Q corresponding to the current node in
    ## the walk has at least one nonzero entry corresponding to the nodes eligible
    ## for visitation, select that row of Q, and multiply it entry-wise by the
    ## factor exp(Vstd[Wh]/M), so that the probability of moving to a node which
    ## has been previously visited is enhanced based on the number of previous
    ## visits.  Then normalize the resulting vector Row by dividing its elements
    ## by their sum, creating a probability vector Prob.  Then randomly select an
    ## eligible node as the next node, with the probability of choosing the jth node
    ## given by Prob_j, and store its index as the next element of I.
    ## Otherwise, select the next node uniformly from those eligible, and store its
    ## index as the next element of I.

    if(sum(Q[I[L],Wh], na.rm=TRUE) > 0 ) {
      Row <- Q[I[L],Wh]
      Prob <- Row / sum(Row)
      I <- c(I, sample(x=Vec[Wh], size=1, prob=Prob))
    } else {
      I <- c(I, sample(x=Vec[Wh], size=1))
    }

    ## Update L to the new length of vector I.  Update the number of visits made to
    ## the node selected as the next node, increasing the number by one.
    ## If that node has now been visited M times, its index is replaced with 0,
    ## thus removing it from the random walk.  Then the vector Recur is updated by
    ## concatenating it with the number of steps since the previous node removal,
    ## which is stored in Count.  Count is reset to 0.
    ## No matter what, Count is increased by one.

    L <- length(I)
    Vstd[I[L]] <- Vstd[I[L]] + 1
    if(Vstd[I[L]] == M) {
      Vec[I[L]] <- 0
      Recur <- c(Recur, Count)
      Count <- 0
    }
    Count <- Count + 1

  }  ## end of while loop

  ## Concatenate I with the index of the last node which has not been visited.
  ## Return the vectors Recur and I.

  I <- c(I, which(Vstd == 0))
  return(list(I=I, Recur=Recur))

} ## end of function MakeSeries
