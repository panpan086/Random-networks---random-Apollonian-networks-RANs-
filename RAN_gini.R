## Compute the Gini index of a simulated RAN of index k at time n >= 2
rm(list=ls(all=TRUE))
## install.packages("igraph")
library(igraph)
giniApoll <- function(n, k){
  
  ## Generate the initial adjacent matrix of a RAN of index k
  adj <- matrix(1, k, k) - diag(k)
  
  ## Label the newcomer of the network
  d <- k + 1
  
  ## Collect all active k-cliques at the point
  activecliques <- matrix(seq(1, k, by = 1), k, 1)
  repeat{
    
    ## Count the number of active k-cliques
    actcliqnum <- ncol(activecliques)
    
    ## Choose an active k-clique at random (all active cliques being equally likely)
    luckytri <- sample(actcliqnum, 1)
    
    ## Add new active cliques upon the newcomer joining into the network
    activecliques <- cbind(activecliques,rbind(combn(activecliques[, luckytri], k - 1), rep(d, k)))
    
    ## Adjust the addjacency matrix after the newcomer joins into the network
    adj_add <- matrix(0, d - 1, 1)
    adj_add[activecliques[,luckytri]] <- 1
    adj <- cbind(adj, adj_add)
    adj <- rbind(adj, c(adj_add, 0))
    
    ## Deactivate the clique that recruits the newcomer
    activecliques <- activecliques[,-luckytri]
    d <- d + 1
   if(d > n + k - 1){
      break
    }
  }
  
  ## Sort the degrees of the vertices
  degseq <- sort(colSums(adj))
  
  ## Calculate the proportions of simulated degrees with respect to the maximum degree
  propseq <- degseq/((n + k - 1)*(n + k - 2))
  
  ## Compute cumulative degrees
  cumpropseq <- cumsum(propseq)
  
  ## Compute the revised Gini index 
  ## Reference: Goswami, S., Murthy, C.A. and Das, A.K. (2018). Sparsity measure of a network graph: Gini index. Information Sciences, 462, 16--39.
  gini <- 1 - 2*((2*sum(cumpropseq) - 1)/(2*(n + k - 1)))
  
  ## Return the revised Gini index
  return(gini)
}
