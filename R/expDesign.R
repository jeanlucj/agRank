#' Allocation of varieties to farmers
#'
#' Generate experimental designs to allocate vartieties to farmers.
#' 17 experimental designs are included, namely they are AP, RD, KM1~KM6, GD1~GD9 (see the vignette for details)
#'
#' @param nvar the number of varieties
#' @param nobs the number of farmers (observers), must be greater than or equal to nvar / 3
#' @param method must be one of "alpha", "random", "KM1", "KM2", "KM3", "KM4", "KM5", "KM6",
#' "GD1", "GD2", "GD3", "GD4", "GD5", "GD6", "GD7", "GD8", "GD9"
#' @param marker for KM1~KM6 and GD1~GD9, a nvar * nloci matrix should be specified, which is coded in (-1, 0, 1)
#' and is allowed to have NA entry
#'
#'
#'
#'
#' @return an nobs * 3 matrix, with each row being the label of varieties assigned to that farmer
#'
#'
#'
#' @examples
#' #use alpha design
#' expDesign(120, 245, method = 'alpha')
#'
#' #use GD2
#' data(marker)
#' M = marker[sample(368, 120), ] #the marker matrix
#' expDesign(120, 245, method = 'GD2', marker = M)
#'
#'
#'
#'
#' @export


expDesign = function(nvar, nobs, method, marker = NA){
  #marker is a nvar * nloci matrix coded in {-1, 0, 1}
  library(doParallel)
  library(stats)
  library(rrBLUP)

  if(is.matrix(marker)){
    #impute the marker matrix and calculate additive relationship matrix
    temp.res = A.mat(marker, shrink = T, return.imputed = T)
    K = temp.res$A
    marker = temp.res$imputed
  }

  ###FUNCTIONS TO GENERATE ONE REPLICA

  #alpha designs
  library(agricolae)
  alpha = function(nvar){
    design = design.alpha(1:nvar, 3, 3)$sketch
    #concatenate three replicates
    design = matrix(as.numeric(rbind(design[[1]], design[[2]], design[[3]])), nvar, 3)
    return(design)
  }

  #random allocations
  rand = function(nvar){
    nobs = nvar / 3

    #random allocation of varieties
    design = matrix(sample(nvar), nobs, 3)

    return(design)
  }

  #KM1
  KM1 = function(marker){
    #marker is a nvar * nloci matrix coded in {-1, 0, 1}.
    nvar = nrow(marker)
    nobs = nvar / 3
    rownames(marker) = 1:nvar

    #store the results
    design = matrix(0, nobs, 3) #the allocation of varieties is row-wise


    #initialize
    varCluster = kmeans(marker, 3)
    size = varCluster$size
    cluster = varCluster$cluster
    niter = 0

    while(1){
      niter = niter + 1

      #varieties are labeled from 1 to nvar
      cluster1 = as.numeric(names(cluster[cluster == 1]))
      cluster2 = as.numeric(names(cluster[cluster == 2]))
      cluster3 = as.numeric(names(cluster[cluster == 3]))

      #first assign min(size) varieties to farmers
      chosen_size = min(size)
      #cluster with minimum size
      min_cluster = which(size == chosen_size)

      chosen1 = sample(size[1], chosen_size)
      chosen2 = sample(size[2], chosen_size)
      chosen3 = sample(size[3], chosen_size)


      startrow = match(0, design[, 1])
      design[startrow:(chosen_size + startrow - 1), 1] = cluster1[chosen1]
      design[startrow:(chosen_size + startrow - 1), 2] = cluster2[chosen2]
      design[startrow:(chosen_size + startrow - 1), 3] = cluster3[chosen3]

      #get the remaining varieties
      remain = c(cluster1[-chosen1], cluster2[-chosen2], cluster3[-chosen3])

      #break
      if(length(remain) == 3){
        #browser()
        design[nobs, ] = remain
        break

      }

      if(length(remain) == 0){
        #browser()
        break
      }



      #re-cluster the remaining varieties
      varCluster = kmeans(marker[remain, ], 3)
      size = varCluster$size
      cluster = varCluster$cluster

    }

    return(design)

  }

  #KM2
  KM2 = function(marker, K){
    #marker is a nvar * nloci matrix coded in {-1, 0, 1}.
    #K is the nvar * nvar relationship matrix

    #browser()

    nvar = nrow(marker)
    nobs = nvar / 3
    rownames(marker) = 1:nvar

    #store the results
    design = matrix(0, nobs, 3) #the allocation of varieties is row-wise


    #initialize
    varCluster = kmeans(marker, 3)
    size = varCluster$size
    cluster = varCluster$cluster
    niter = 0

    while(1){
      niter = niter + 1

      #varieties are labeled from 1 to nvar
      cluster1 = as.numeric(names(cluster[cluster == 1]))
      cluster2 = as.numeric(names(cluster[cluster == 2]))
      cluster3 = as.numeric(names(cluster[cluster == 3]))

      cluster_list = list(cluster1, cluster2, cluster3)

      #first assign min(size) varieties to farmers
      chosen_size = min(size)
      #cluster with minimum size
      min_cluster = which(size == chosen_size)


      #for each cluster, choose the most related chosen_size varieties
      #chosen is a list with i-th component being varieties chosen from cluster i
      #remain is a list with i-th component being varieties remaining from cluster i
      chosen = list(c(), c(), c())
      remain = c()
      for(i in 1:3){

        if(i %in% min_cluster){ #if this is the cluster with the smallest size

          chosen[[i]] = cluster_list[[i]][sample(chosen_size)]



        } else{
          cluster_each = cluster_list[[i]]
          #compute row-wise sum (diagonal elements deducted)
          similarity = apply(K[cluster_each, cluster_each], 1, FUN = sum) - diag(K[cluster_each, cluster_each])
          names(similarity) = cluster_each #assign labels to the similarity vector

          #order the similarity in decreasing order
          ordered_var = as.numeric(names(sort(similarity, decreasing = T)))

          #choose top chosen_size varieties and randomize them
          chosen[[i]] = ordered_var[1:chosen_size][sample(chosen_size)]

          #remaining varieties
          remain = c(remain, ordered_var[-(1:chosen_size)])

        }

      }

      startrow = match(0, design[, 1])
      design[startrow:(chosen_size + startrow - 1), 1] = chosen[[1]]
      design[startrow:(chosen_size + startrow - 1), 2] = chosen[[2]]
      design[startrow:(chosen_size + startrow - 1), 3] = chosen[[3]]



      #break
      if(length(remain) == 3){
        #browser()
        design[nobs, ] = remain
        break

      }

      if(length(remain) == 0){
        #browser()
        break
      }



      #re-cluster the remaining varieties
      varCluster = kmeans(marker[remain, ], 3)
      size = varCluster$size
      cluster = varCluster$cluster

    }

    return(design)

  }

  #KM3
  KM3 = function(marker, K){

    #marker is a nvar * nloci matrix coded in {-1, 0, 1}.
    #K is the additive relationship matrix computed by A.mat function in rrBLUP
    nvar = nrow(marker)
    nobs = nvar / 3
    rownames(marker) = 1:nvar

    #store the results
    design = matrix(0, nobs, 3) #the allocation of varieties is row-wise

    #cluster varieties into 3 clusters
    init_cluster = kmeans(marker, 3)

    #a list with i-th component being the varieties in i-th cluster
    cluster = list()
    for(i in 1:3){
      cluster[[i]] = as.numeric(names(init_cluster$cluster[init_cluster$cluster == i]))
    }

    #choose the most related nobs items for cluster with size > nobs
    remain = c() #varieties remained (and should be assigned to other clusters)
    for(i in which(init_cluster$size > nobs)){

      cluster_each = cluster[[i]]

      #compute row-wise sum (diagonal elements deducted)
      similarity = apply(K[cluster_each, cluster_each], 1, FUN = sum) - diag(K[cluster_each, cluster_each])
      names(similarity) = cluster_each #assign labels to the similarity vector

      #order the similarity in decreasing order
      ordered_var = as.numeric(names(sort(similarity, decreasing = T)))

      #choose top nobs varieties
      cluster[[i]] = ordered_var[1:nobs]

      #remaining varieties
      remain = c(remain, ordered_var[-(1:nobs)])

    }

    #for clusters with size < nobs, they each choose varieties (from remaining varieties)
    #that are most similar to them, until they reach the size of nobs

    size_less_than_nobs = which(init_cluster$size < nobs) #the label of clusters whose size < nobs
    new_size_less_than_nobs = size_less_than_nobs

    while(1){

      for(i in size_less_than_nobs){

        cluster_each = cluster[[i]]

        #remaining varieties are column-wise
        K_each = K[cluster_each, remain]

        #column-wise mean correlation
        #now the each element in similarity is the similarity between that variety and cluster_each
        if(length(cluster_each) == 1){ #if there's only one variety in that cluster

          similarity = K_each

        } else{

          K_each = as.matrix(K_each) #in case there's only one variety in remain, then K_each is treated as a column matrix
          similarity = apply(K_each, 2, mean)

        }
        names(similarity) = remain #assign labels to similarity

        #add the most similar variety to that cluster
        ordered_var = as.numeric(names(sort(similarity, decreasing = T)))
        cluster[[i]] = c(cluster_each, ordered_var[1])



        #delete the variety just chosen from remain
        remain = remain[remain != ordered_var[1]]

        #if this cluster reaches the size of nobs
        if(length(cluster[[i]]) == nobs){

          #delete this cluster from clusters with size < 3
          new_size_less_than_nobs = new_size_less_than_nobs[new_size_less_than_nobs != i]

        }

      }

      size_less_than_nobs = new_size_less_than_nobs

      #if all the clusters reach the size of nobs
      if(length(size_less_than_nobs) == 0){

        break

      }


    }


    #each farmer choose one variety from each cluster
    design[, 1] = cluster[[1]][sample(nobs)]
    design[, 2] = cluster[[2]][sample(nobs)]
    design[, 3] = cluster[[3]][sample(nobs)]


    return(design)
  }

  #KM4
  KM4 = function(marker, K){

    #marker is a nvar * nloci matrix coded in {-1, 0, 1}.
    #K is the additive relationship matrix computed by A.mat function in rrBLUP

    nvar = nrow(marker)
    nobs = nvar / 3
    rownames(marker) = 1:nvar

    #store the results
    design = matrix(0, nobs, 3) #the allocation of varieties is row-wise

    #cluster varieties into 3 clusters
    init_cluster = kmeans(marker, 3)

    #a list with i-th component being the varieties in i-th cluster
    cluster = list()
    for(i in 1:3){
      cluster[[i]] = as.numeric(names(init_cluster$cluster[init_cluster$cluster == i]))
    }

    #choose the most related nobs items for cluster with size > nobs
    remain = c() #varieties remained (and should be assigned to other clusters)
    for(i in which(init_cluster$size > nobs)){

      cluster_each = cluster[[i]]

      #compute row-wise sum (diagonal elements deducted)
      similarity = apply(K[cluster_each, cluster_each], 1, FUN = sum) - diag(K[cluster_each, cluster_each])
      names(similarity) = cluster_each #assign labels to the similarity vector

      #order the similarity in decreasing order
      ordered_var = as.numeric(names(sort(similarity, decreasing = T)))

      #choose top nobs varieties
      cluster[[i]] = ordered_var[1:nobs]

      #remaining varieties
      remain = c(remain, ordered_var[-(1:nobs)])

    }

    #for the remaining varieties, each variety choose its favorate cluster to go (according to correlation)
    size_less_than_nobs = which(init_cluster$size < nobs) #the label of clusters whose size < nobs


    for(i in 1:length(remain)){

      remain_var = remain[i]

      #similarity measure, each element is the similarity between remain_var and that cluster
      similarity = rep(0, length(size_less_than_nobs))
      names(similarity) = size_less_than_nobs #assign labels

      for(j in 1:length(size_less_than_nobs)){

        cluster_label = size_less_than_nobs[j]
        var_in_cluster = cluster[[cluster_label]]
        similarity[j] = mean(K[remain_var, var_in_cluster])

      }

      #order the cluster according to similarity
      ordered_cluster = as.numeric(names(sort(similarity, decreasing = T)))

      #remain_var enters the most similar cluster with size < nobs
      cluster[[ordered_cluster[1]]] = c(cluster[[ordered_cluster[1]]], remain_var)

      #if this cluster reaches the size of nobs
      if(length(cluster[[ordered_cluster[1]]]) == nobs){

        #delete this cluster from size_less_than_nobs
        size_less_than_nobs = size_less_than_nobs[size_less_than_nobs != ordered_cluster[1]]

      }



    }


    #each farmer choose one variety from each cluster
    design[, 1] = cluster[[1]][sample(nobs)]
    design[, 2] = cluster[[2]][sample(nobs)]
    design[, 3] = cluster[[3]][sample(nobs)]


    return(design)
  }

  #KM5
  KM5 = function(marker, K){
    #marker is a nvar * nloci matrix coded in {-1, 0, 1}.
    #K is the additive relationship matrix computed by A.mat function in rrBLUP
    nvar = nrow(marker)
    nobs = nvar / 3
    rownames(marker) = 1:nvar

    #store the results
    design = matrix(0, nobs, 3) #the allocation of varieties is row-wise

    #cluster varieties into nobs clusters
    init_cluster = kmeans(marker, nobs)

    #a list with i-th component being the varieties in i-th cluster
    cluster = list()
    for(i in 1:nobs){
      cluster[[i]] = as.numeric(names(init_cluster$cluster[init_cluster$cluster == i]))
    }

    #choose the most related 3 items for cluster with size > 3
    remain = c() #varieties remained (and should be assigned to other clusters)
    for(i in which(init_cluster$size > 3)){

      cluster_each = cluster[[i]]

      #compute row-wise sum (diagonal elements deducted)
      similarity = apply(K[cluster_each, cluster_each], 1, FUN = sum) - diag(K[cluster_each, cluster_each])
      names(similarity) = cluster_each #assign labels to the similarity vector

      #order the similarity in decreasing order
      ordered_var = as.numeric(names(sort(similarity, decreasing = T)))

      #choose top 3 varieties
      cluster[[i]] = ordered_var[1:3]

      #remaining varieties
      remain = c(remain, ordered_var[-(1:3)])

    }

    #for clusters with size < 3, they each choose varieties (from remaining varieties)
    #that are most similar to them, until they reach the size of 3

    size_less_than3 = which(init_cluster$size < 3) #the label of clusters whose size < 3
    new_size_less_than3 = size_less_than3

    while(1){

      for(i in size_less_than3){

        cluster_each = cluster[[i]]

        #remaining varieties are column-wise
        K_each = K[cluster_each, remain]

        #column-wise mean correlation
        #now the each element in similarity is the similarity between that variety and cluster_each
        if(length(cluster_each) == 1){ #if there's only one variety in that cluster

          similarity = K_each

        } else{

          K_each = as.matrix(K_each) #in case there's only one variety in remain, then K_each is treated as a column matrix
          similarity = apply(K_each, 2, mean)

        }
        names(similarity) = remain #assign labels to similarity

        #add the most similar variety to that cluster
        ordered_var = as.numeric(names(sort(similarity, decreasing = T)))
        cluster[[i]] = c(cluster_each, ordered_var[1])



        #delete the variety just chosen from remain
        remain = remain[remain != ordered_var[1]]

        #if this cluster reaches the size of 3
        if(length(cluster[[i]]) == 3){

          #delete this cluster from clusters with size < 3
          new_size_less_than3 = new_size_less_than3[new_size_less_than3 != i]

        }

      }

      size_less_than3 = new_size_less_than3

      #if all the clusters reach the size of 3
      if(length(size_less_than3) == 0){

        break

      }


    }


    #each cluster represents one farmer
    for(i in 1:nobs){
      design[i, ] = cluster[[i]]
    }



    return(design)
  }

  #KM6
  KM6 = function(marker, K){
    #marker is a nvar * nloci matrix coded in {-1, 0, 1}.
    #K is the additive relationship matrix computed by A.mat function in rrBLUP

    nvar = nrow(marker)
    nobs = nvar / 3
    rownames(marker) = 1:nvar

    #store the results
    design = matrix(0, nobs, 3) #the allocation of varieties is row-wise

    #cluster varieties into nobs clusters
    init_cluster = kmeans(marker, nobs)

    #a list with i-th component being the varieties in i-th cluster
    cluster = list()
    for(i in 1:nobs){
      cluster[[i]] = as.numeric(names(init_cluster$cluster[init_cluster$cluster == i]))
    }

    #choose the most related 3 items for cluster with size > 3
    remain = c() #varieties remained (and should be assigned to other clusters)
    for(i in which(init_cluster$size > 3)){

      cluster_each = cluster[[i]]

      #compute row-wise sum (diagonal elements deducted)
      similarity = apply(K[cluster_each, cluster_each], 1, FUN = sum) - diag(K[cluster_each, cluster_each])
      names(similarity) = cluster_each #assign labels to the similarity vector

      #order the similarity in decreasing order
      ordered_var = as.numeric(names(sort(similarity, decreasing = T)))

      #choose top 3 varieties
      cluster[[i]] = ordered_var[1:3]

      #remaining varieties
      remain = c(remain, ordered_var[-(1:3)])

    }

    #for the remaining varieties, each variety choose its favorate cluster to go (according to correlation)
    size_less_than3 = which(init_cluster$size < 3) #the label of clusters whose size < 3


    for(i in 1:length(remain)){

      remain_var = remain[[i]]

      #similarity measure, each element is the similarity between remain_var and that cluster
      similarity = rep(0, length(size_less_than3))
      names(similarity) = size_less_than3 #assign labels

      for(j in 1:length(size_less_than3)){

        cluster_label = size_less_than3[j]
        var_in_cluster = cluster[[cluster_label]]
        similarity[j] = mean(K[remain_var, var_in_cluster])

      }

      #order the cluster according to similarity
      ordered_cluster = as.numeric(names(sort(similarity, decreasing = T)))

      #remain_var enters the most similar cluster with size < 3
      cluster[[ordered_cluster[1]]] = c(cluster[[ordered_cluster[1]]], remain_var)

      #if this cluster reaches the size of 3
      if(length(cluster[[ordered_cluster[1]]]) == 3){

        #delete this cluster from size_less_than3
        size_less_than3 = size_less_than3[size_less_than3 != ordered_cluster[1]]

      }



    }


    #each cluster represents one farmer
    for(i in 1:nobs){
      design[i, ] = cluster[[i]]
    }


    return(design)
  }

  #GD1
  GD1 = function(K){
    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #loop over all 3 clusters
    for(i in 1:3){

      remain = remain[sample(length(remain))] #randomize

      var_each = remain[1] #choose one variety to enter cluster i

      #compute similarity between var_each and all the other remaining varieties
      similarity = K[var_each, remain[-1]]
      names(similarity) = remain[-1] #assign labels


      #choose the most related (nobs - 1) varieties to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = TRUE)))

      cluster[[i]] = c(var_each, ordered_var[1:(nobs - 1)])

      #delete those varieties from remain
      remain = setdiff(remain, cluster[[i]])

    }

    #assign varieties: each farmer chooses one from each cluster
    design = matrix(0, nobs, 3)
    design[, 1] = sample(cluster[[1]])
    design[, 2] = sample(cluster[[2]])
    design[, 3] = sample(cluster[[3]])


    return(design)

  }

  #GD2
  GD2 = function(K){
    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #initialize three clusters
    remain = remain[sample(length(remain))] #randomize
    cluster[[1]] = remain[1] #this variety goes to the i-th cluster
    #delete this variety from remain
    remain = setdiff(remain, cluster[[1]])
    #varieties alrealdy in clusters
    var_in_cluster = cluster[[1]]
    #loop over the other 2 clusters
    for(i in 2:3){

      #compute mean similarity between varieties in cluster[[i-1]] and all the other remaining varieties
      K_each = K[var_in_cluster, remain]

      if(length(var_in_cluster) == 1){ #if there's only 1 variety in the previous cluster

        similarity = K_each

      } else{

        similarity = apply(K_each, 2, mean)

      }

      names(similarity) = remain #assign labels


      #choose the most unrelated variety to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = FALSE)))

      #the most unrelated variety goes to cluster i
      cluster[[i]] = ordered_var[1]
      var_in_cluster = c(var_in_cluster, cluster[[i]])

      #delete this variety from remain

      remain = setdiff(remain, cluster[[i]])

    }

    #for clusters with size < nobs, they each choose varieties (from remaining varieties)
    #that are most similar to them, until they reach the size of nobs

    size_less_than_nobs = c(1, 2, 3) #the label of clusters whose size < nobs
    new_size_less_than_nobs = size_less_than_nobs

    while(1){

      for(i in size_less_than_nobs){

        cluster_each = cluster[[i]]

        #remaining varieties are column-wise
        K_each = K[cluster_each, remain]

        #column-wise mean correlation
        #now the each element in similarity is the similarity between that variety and cluster_each
        if(length(cluster_each) == 1){ #if there's only one variety in that cluster

          similarity = K_each

        } else{

          K_each = as.matrix(K_each) #in case there's only one variety in remain, then K_each is treated as a column matrix
          similarity = apply(K_each, 2, mean)

        }

        names(similarity) = remain #assign labels to similarity

        #add the most similar variety to that cluster
        ordered_var = as.numeric(names(sort(similarity, decreasing = T)))
        cluster[[i]] = c(cluster_each, ordered_var[1])



        #delete the variety just chosen from remain
        remain = remain[remain != ordered_var[1]]

        #if this cluster reaches the size of nobs
        if(length(cluster[[i]]) == nobs){

          #delete this cluster from clusters with size < 3
          new_size_less_than_nobs = new_size_less_than_nobs[new_size_less_than_nobs != i]

        }

      }

      size_less_than_nobs = new_size_less_than_nobs

      #if all the clusters reach the size of nobs
      if(length(size_less_than_nobs) == 0){

        break

      }


    }


    #each farmer choose one variety from each cluster
    design = matrix(0, nobs, 3)
    design[, 1] = cluster[[1]][sample(nobs)]
    design[, 2] = cluster[[2]][sample(nobs)]
    design[, 3] = cluster[[3]][sample(nobs)]



    return(design)


  }

  #GD3
  GD3 = function(K){
    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #initialize three clusters
    remain = remain[sample(length(remain))] #randomize
    cluster[[1]] = remain[1] #this variety goes to the i-th cluster
    #delete this variety from remain
    remain = setdiff(remain, cluster[[1]])
    #varieties alrealdy in clusters
    var_in_cluster = cluster[[1]]
    #loop over the other 2 clusters
    for(i in 2:3){

      #compute mean similarity between varieties in cluster[[i-1]] and all the other remaining varieties
      K_each = K[var_in_cluster, remain]

      if(length(var_in_cluster) == 1){ #if there's only 1 variety in the previous cluster

        similarity = K_each

      } else{

        similarity = apply(K_each, 2, mean)

      }

      names(similarity) = remain #assign labels


      #choose the most unrelated variety to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = FALSE)))

      #the most unrelated variety goes to cluster i
      cluster[[i]] = ordered_var[1]
      var_in_cluster = c(var_in_cluster, cluster[[i]])

      #delete this variety from remain

      remain = setdiff(remain, cluster[[i]])

    }

    #for the remaining varieties, each variety choose its favorate cluster to go (according to correlation)
    size_less_than_nobs = c(1, 2, 3) #the label of clusters whose size < nobs


    for(i in 1:length(remain)){

      remain_var = remain[[i]]

      #similarity measure, each element is the similarity between remain_var and that cluster
      similarity = rep(0, length(size_less_than_nobs))
      names(similarity) = size_less_than_nobs #assign labels

      for(j in 1:length(size_less_than_nobs)){

        cluster_label = size_less_than_nobs[j]
        var_in_cluster = cluster[[cluster_label]]
        similarity[j] = mean(K[remain_var, var_in_cluster])

      }

      #order the cluster according to similarity
      ordered_cluster = as.numeric(names(sort(similarity, decreasing = T)))

      #remain_var enters the most similar cluster with size < nobs
      cluster[[ordered_cluster[1]]] = c(cluster[[ordered_cluster[1]]], remain_var)

      #if this cluster reaches the size of nobs
      if(length(cluster[[ordered_cluster[1]]]) == nobs){

        #delete this cluster from size_less_than_nobs
        size_less_than_nobs = size_less_than_nobs[size_less_than_nobs != ordered_cluster[1]]

      }



    }




    #each farmer choose one variety from each cluster
    design = matrix(0, nobs, 3)
    design[, 1] = cluster[[1]][sample(nobs)]
    design[, 2] = cluster[[2]][sample(nobs)]
    design[, 3] = cluster[[3]][sample(nobs)]


    return(design)

  }

  #GD4
  GD4 = function(K){

    #K: additive relationship matrix


    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #loop over all nobs clusters
    for(i in 1:nobs){

      remain = remain[sample(length(remain))] #randomize

      var_each = remain[1] #choose one variety to enter cluster i

      #compute similarity between var_each and all the other remaining varieties
      similarity = K[var_each, remain[-1]]
      names(similarity) = remain[-1] #assign labels


      #choose the most related 2 varieties to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = TRUE)))

      cluster[[i]] = c(var_each, ordered_var[1:2])

      #delete those varieties from remain
      remain = setdiff(remain, cluster[[i]])

    }

    #assign varieties: each farmer chooses one from each cluster
    design = matrix(0, nobs, 3)
    for(i in 1:nobs){
      design[i, ] = cluster[[i]]
    }

    return(design)

  }

  #GD5
  GD5 = function(K){

    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #initialize nobs clusters
    remain = remain[sample(length(remain))] #randomize
    cluster[[1]] = remain[1] #this variety goes to the first cluster
    #delete this variety from remain
    remain = setdiff(remain, cluster[[1]])
    #varieties alrealdy in clusters
    var_in_cluster = cluster[[1]]
    #loop over the other (nobs - 1) clusters
    for(i in 2:nobs){

      #compute mean similarity between varieties in cluster[[i-1]] and all the other remaining varieties
      K_each = K[var_in_cluster, remain]

      if(length(var_in_cluster) == 1){ #if there's only 1 variety in the previous cluster

        similarity = K_each

      } else{

        similarity = apply(K_each, 2, mean)

      }

      names(similarity) = remain #assign labels


      #choose the most unrelated variety to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = FALSE)))

      #the most unrelated variety goes to cluster i
      cluster[[i]] = ordered_var[1]
      var_in_cluster = c(var_in_cluster, cluster[[i]])

      #delete this variety from remain

      remain = setdiff(remain, cluster[[i]])

    }
    #for clusters with size < 3, they each choose varieties (from remaining varieties)
    #that are most similar to them, until they reach the size of 3

    size_less_than3 = 1:nobs #the label of clusters whose size < 3
    new_size_less_than3 = size_less_than3

    while(1){

      for(i in size_less_than3){

        cluster_each = cluster[[i]]

        #remaining varieties are column-wise
        K_each = K[cluster_each, remain]

        #column-wise mean correlation
        #now the each element in similarity is the similarity between that variety and cluster_each
        if(length(cluster_each) == 1){ #if there's only one variety in that cluster

          similarity = K_each

        } else{

          K_each = as.matrix(K_each) #in case there's only one variety in remain, then K_each is treated as a column matrix
          similarity = apply(K_each, 2, mean)

        }
        names(similarity) = remain #assign labels to similarity

        #add the most similar variety to that cluster
        ordered_var = as.numeric(names(sort(similarity, decreasing = T)))
        cluster[[i]] = c(cluster_each, ordered_var[1])



        #delete the variety just chosen from remain
        remain = remain[remain != ordered_var[1]]

        #if this cluster reaches the size of 3
        if(length(cluster[[i]]) == 3){

          #delete this cluster from clusters with size < 3
          new_size_less_than3 = new_size_less_than3[new_size_less_than3 != i]

        }

      }

      size_less_than3 = new_size_less_than3

      #if all the clusters reach the size of 3
      if(length(size_less_than3) == 0){

        break

      }


    }


    #each farmer choose one variety from each cluster
    design = matrix(0, nobs, 3)
    for(i in 1:nobs){
      design[i, ] = cluster[[i]]
    }

    return(design)


  }

  #GD6
  GD6 = function(K){

    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #initialize three clusters
    remain = remain[sample(length(remain))] #randomize
    cluster[[1]] = remain[1] #this variety goes to the i-th cluster
    #delete this variety from remain
    remain = setdiff(remain, cluster[[1]])
    #varieties alrealdy in clusters
    var_in_cluster = cluster[[1]]
    #loop over the other (nobs - 1) clusters
    for(i in 2:nobs){

      #compute mean similarity between varieties in cluster[[i-1]] and all the other remaining varieties
      K_each = K[var_in_cluster, remain]

      if(length(var_in_cluster) == 1){ #if there's only 1 variety in the previous cluster

        similarity = K_each

      } else{

        similarity = apply(K_each, 2, mean)

      }

      names(similarity) = remain #assign labels


      #choose the most unrelated variety to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = FALSE)))

      #the most unrelated variety goes to cluster i
      cluster[[i]] = ordered_var[1]
      var_in_cluster = c(var_in_cluster, cluster[[i]])

      #delete this variety from remain

      remain = setdiff(remain, cluster[[i]])

    }


    #for the remaining varieties, each variety choose its favorate cluster to go (according to correlation)
    size_less_than3 = 1:nobs #the label of clusters whose size < 3


    for(i in 1:length(remain)){

      remain_var = remain[[i]]

      #similarity measure, each element is the similarity between remain_var and that cluster
      similarity = rep(0, length(size_less_than3))
      names(similarity) = size_less_than3 #assign labels

      for(j in 1:length(size_less_than3)){

        cluster_label = size_less_than3[j]
        var_in_cluster = cluster[[cluster_label]]
        similarity[j] = mean(K[remain_var, var_in_cluster])

      }

      #order the cluster according to similarity
      ordered_cluster = as.numeric(names(sort(similarity, decreasing = T)))

      #remain_var enters the most similar cluster with size < 3
      cluster[[ordered_cluster[1]]] = c(cluster[[ordered_cluster[1]]], remain_var)

      #if this cluster reaches the size of 3
      if(length(cluster[[ordered_cluster[1]]]) == 3){

        #delete this cluster from size_less_than3
        size_less_than3 = size_less_than3[size_less_than3 != ordered_cluster[1]]

      }



    }


    #each cluster represents one farmer
    design = matrix(0, nobs, 3)
    for(i in 1:nobs){
      design[i, ] = cluster[[i]]
    }

    return(design)

  }

  #GD7
  GD7 = function(K){

    #K: additive relationship matrix


    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #loop over all nobs clusters
    for(i in 1:nobs){

      remain = remain[sample(length(remain))] #randomize

      var_each = remain[1] #choose one variety to enter cluster i

      #compute similarity between var_each and all the other remaining varieties
      similarity = K[var_each, remain[-1]]
      names(similarity) = remain[-1] #assign labels


      #choose the most unrelated 2 varieties to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = FALSE)))

      cluster[[i]] = c(var_each, ordered_var[1:2])

      #delete those varieties from remain
      remain = setdiff(remain, cluster[[i]])

    }

    #assign varieties: each farmer chooses one from each cluster
    design = matrix(0, nobs, 3)
    for(i in 1:nobs){
      design[i, ] = cluster[[i]]
    }


    return(design)


  }

  #GD8
  GD8 = function(K){

    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #initialize nobs clusters
    remain = remain[sample(length(remain))] #randomize
    cluster[[1]] = remain[1] #this variety goes to the first cluster
    #delete this variety from remain
    remain = setdiff(remain, cluster[[1]])
    #varieties alrealdy in clusters
    var_in_cluster = cluster[[1]]
    #loop over the other (nobs - 1) clusters
    for(i in 2:nobs){

      #compute mean similarity between varieties in cluster[[i-1]] and all the other remaining varieties
      K_each = K[var_in_cluster, remain]

      if(length(var_in_cluster) == 1){ #if there's only 1 variety in the previous cluster

        similarity = K_each

      } else{

        similarity = apply(K_each, 2, mean)

      }

      names(similarity) = remain #assign labels


      #choose the most related variety to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = TRUE)))

      #the most unrelated variety goes to cluster i
      cluster[[i]] = ordered_var[1]
      var_in_cluster = c(var_in_cluster, cluster[[i]])

      #delete this variety from remain

      remain = setdiff(remain, cluster[[i]])

    }
    #for clusters with size < 3, they each choose varieties (from remaining varieties)
    #that are most similar to them, until they reach the size of 3

    size_less_than3 = 1:nobs #the label of clusters whose size < 3
    new_size_less_than3 = size_less_than3

    while(1){

      for(i in size_less_than3){

        cluster_each = cluster[[i]]

        #remaining varieties are column-wise
        K_each = K[cluster_each, remain]

        #column-wise mean correlation
        #now the each element in similarity is the similarity between that variety and cluster_each
        if(length(cluster_each) == 1){ #if there's only one variety in that cluster

          similarity = K_each

        } else{

          K_each = as.matrix(K_each) #in case there's only one variety in remain, then K_each is treated as a column matrix
          similarity = apply(K_each, 2, mean)

        }
        names(similarity) = remain #assign labels to similarity

        #add the most dissimilar variety to that cluster
        ordered_var = as.numeric(names(sort(similarity, decreasing = FALSE)))
        cluster[[i]] = c(cluster_each, ordered_var[1])



        #delete the variety just chosen from remain
        remain = remain[remain != ordered_var[1]]

        #if this cluster reaches the size of 3
        if(length(cluster[[i]]) == 3){

          #delete this cluster from clusters with size < 3
          new_size_less_than3 = new_size_less_than3[new_size_less_than3 != i]

        }

      }

      size_less_than3 = new_size_less_than3

      #if all the clusters reach the size of 3
      if(length(size_less_than3) == 0){

        break

      }


    }


    #each farmer choose one variety from each cluster
    design = matrix(0, nobs, 3)
    for(i in 1:nobs){
      design[i, ] = cluster[[i]]
    }


    return(design)

  }

  #GD9
  GD9 = function(K){

    #K: additive relationship matrix
    nvar = dim(K)[1]
    nobs = nvar / 3

    #varieties that haven't been assigned to any clusters
    remain = 1:nvar

    #a list where the i-th component is the varieties in i-th cluster
    cluster = list()

    #initialize three clusters
    remain = remain[sample(length(remain))] #randomize
    cluster[[1]] = remain[1] #this variety goes to the i-th cluster
    #delete this variety from remain
    remain = setdiff(remain, cluster[[1]])
    #varieties alrealdy in clusters
    var_in_cluster = cluster[[1]]
    #loop over the other (nobs - 1) clusters
    for(i in 2:nobs){

      #compute mean similarity between varieties in cluster[[i-1]] and all the other remaining varieties
      K_each = K[var_in_cluster, remain]

      if(length(var_in_cluster) == 1){ #if there's only 1 variety in the previous cluster

        similarity = K_each

      } else{

        similarity = apply(K_each, 2, mean)

      }

      names(similarity) = remain #assign labels


      #choose the most unrelated variety to enter cluster i
      ordered_var = as.numeric(names(sort(similarity, decreasing = TRUE)))

      #the most unrelated variety goes to cluster i
      cluster[[i]] = ordered_var[1]
      var_in_cluster = c(var_in_cluster, cluster[[i]])

      #delete this variety from remain

      remain = setdiff(remain, cluster[[i]])

    }


    #for the remaining varieties, each variety choose its favorate cluster to go (according to correlation)
    size_less_than3 = 1:nobs #the label of clusters whose size < 3


    for(i in 1:length(remain)){

      remain_var = remain[[i]]

      #similarity measure, each element is the similarity between remain_var and that cluster
      similarity = rep(0, length(size_less_than3))
      names(similarity) = size_less_than3 #assign labels

      for(j in 1:length(size_less_than3)){

        cluster_label = size_less_than3[j]
        var_in_cluster = cluster[[cluster_label]]
        similarity[j] = mean(K[remain_var, var_in_cluster])

      }

      #order the cluster according to similarity
      ordered_cluster = as.numeric(names(sort(similarity, decreasing = FALSE)))

      #remain_var enters the most similar cluster with size < 3
      cluster[[ordered_cluster[1]]] = c(cluster[[ordered_cluster[1]]], remain_var)

      #if this cluster reaches the size of 3
      if(length(cluster[[ordered_cluster[1]]]) == 3){

        #delete this cluster from size_less_than3
        size_less_than3 = size_less_than3[size_less_than3 != ordered_cluster[1]]

      }



    }


    #each cluster represents one farmer
    design = matrix(0, nobs, 3)
    for(i in 1:nobs){
      design[i, ] = cluster[[i]]
    }

    return(design)

  }


  ratio = ceiling(nobs / (nvar / 3)) #how many replicas needed

  #require nvar can be divided by 3
  if(nvar %% 3 != 0){
    stop('nvar must be a multiple of 3')
  }

  #require nvar / 3 <= nobs so that each variety is compared at least once

  if(nvar / 3 > nobs){
    stop('Require nvar / 3 <= nobs so that each variety is compared at least once')
  }

  #if(method %in% c('alpha', 'random', 'KM1', 'KM2', 'KM3', 'KM4', 'KM5', 'KM6',
   #                'GD1', 'GD2', 'GD3', 'GD4', 'GD5', 'GD6', 'GD7', 'GD8', 'GD9') == FALSE){

  #}

  if(method == 'alpha'){
    
    ratio = ceiling(nobs / nvar)

    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      alpha(nvar)

    }
  } else if(method == 'random'){

    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      rand(nvar)

    }

  } else if(method == 'KM1'){
    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }

    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      KM1(marker)

    }
  } else if(method == 'KM2'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }

    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      KM2(marker, K)

    }

  } else if(method == 'KM3'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }
    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      KM3(marker, K)

    }

  } else if(method == 'KM4'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }

    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      KM4(marker, K)

    }

  } else if(method == 'KM5'){
    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }


    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      KM5(marker, K)

    }

  } else if(method == 'KM6'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }


    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      KM6(marker, K)

    }

  } else if(method == 'GD1'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }

    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD1(K)

    }

  } else if(method == 'GD2'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }

    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD2(K)

    }

  } else if(method == 'GD3'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }
    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD3(K)

    }
  } else if(method == 'GD4'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }


    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD4(K)

    }

  } else if(method == 'GD5'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }

    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD5(K)

    }
  } else if(method == 'GD6'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }
    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD6(K)

    }
  } else if(method == 'GD7'){
    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }
    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD7(K)

    }

  } else if(method == 'GD8'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }
    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD8(K)

    }
  } else if(method == 'GD9'){

    if(!is.matrix(marker)){
      stop('marker matrix not specified')
    }
    design = foreach(k = 1:ratio, .combine = 'rbind') %do% {

      GD9(K)

    }
  } else{

    stop('Method must be one of alpha, random, KM1, KM2, KM3, KM4, KM5, KM6, GD1, GD2, GD3, GD4, GD5, GD6, GD7,
         GD8, GD9')

  }



  #choose the first nobs rows to ensure the number of occurrence of each variety is "equalized"
  design = design[1:nobs, ]

  #randomize the design again
  design = design[sample(nobs), ]

  return(design)

}
