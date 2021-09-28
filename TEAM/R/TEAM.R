#' @title Testing on an Aggregation Tree Method
#' @description This function performs multiple testing embedded in an aggregation tree structure in order to identify local differences between two probability density functions
#' @importFrom stats ks data.table ggplot2 plyr dplyr
#' @param partition_info Partition for first layer of aggregation tree
#' @param alpha Target false discovery rate (FDR) level
#' @param L Number of layers in the aggregation tree
#' @return A \code{\link{list}} containing the number of pooled observations in each bin (n), the number of bins/leaves at each layer (m.l), the discoveries (S.list) in each layer and the estimated layer-specific thresholds (c.hats)
#' @references Pura J, Li X, Chan C, Xie J. TEAM: A Multiple Testing Algorithm on the Aggregation Tree for Flow Cytometry Analysis \url{https://arxiv.org/abs/1906.07757}
#' @examples
#' ## Example with 1D pdfs: find where case density is higher than control density
#' set.seed(1)
#' # Simulate local shift difference for each sample from mixture of normals
#' N1 <- N2 <- 1e6
#' require(ks) #loads rnorm.mixt function
#' #Controls
#' x1 <- rnorm.mixt(N1,mus=c(0.2,0.89),sigmas=c(0.04,0.01),props=c(0.97,0.03))
#' #Cases
#' x2 <- rnorm.mixt(N2,mus=c(0.2,0.88),sigmas=c(0.04,0.01),props=c(0.97,0.03))
#' #Create partition
#' #Create data frames for data
#' x1.df <- data.frame(X=x1)
#' x2.df <- data.frame(X=x2)
#' #Choose number of bins
#' m = 2^14 #Close to ((N1+N2))^(2/3)
#' partition1D <- create_partition_info(df1=x1.df,df2=x2.df,m=m)
#' res <- TEAM(partition_info=partition1D,alpha=0.05,L=3)
#' #Indices of layer 1 bins rejected at each aggregation layer (not run)
#' #res$S.list
#' #Map rejected bins to corresponding regions (not run)
#' #partition1D$bin.df$x.bd[res$S.list[[3]]]
#' @export TEAM
TEAM = function(partition_info,alpha,L){

  #Local copy
  bin_df_ <- setDT(partition_info$bin.df)

  #Create layer 1 indices and sort by index
  bin_df_ <- bin_df_[, L1 := indx]
  #bin_df_ <- setorder(bin_df_,indx)

  #number of pooled observations per bin
  N1 = partition_info$N1
  N2 = partition_info$N2

  #Global/Initial variables
  c.hats = vector("integer",L) #set of critical values
  S = NULL #set of rejections regions
  #m.excl = NULL #set of excluded regions
  S.list=vector("list",L)
  m.l.vec=vector("list",L)
  #Label bins across row (1,m+1,2*m+1,...)

  theta0 = N2/(N1+N2)
  n = round((N1+N2)/(partition_info$m^partition_info$nd))

  #Loop through layers
  for (l in seq(L)){

    if(l > 1){
      #Create bin info for next layers
      Ll <- paste0("L",l)
      bin_df_ <-  bin_df_[, (Ll) := {
        idx <- ceiling(rleid(L1) / 2^(l-1))
        numL1 <- uniqueN(L1)
        if(numL1 > 1 && between(numL1 %% 2^(l-1),1,2^(l-1)-1)) {
          idx[idx == idx[.N]] <- NA
        }
        idx
      }]

      c_hat_prev = c.hats[l-1]
    } else{
      c_hat_prev = NULL
    }

    #Run procedure
    proc_res = run_proc(l=l,
                        n=n,
                        partition_info_l=bin_df_,
                        theta0=theta0,
                        c_hat_prev=c_hat_prev,
                        alpha=alpha)

    S = proc_res$rej

    S.list[l] = list(unique(c(unlist(S.list),S)))

    m.l.vec[l] = proc_res$m.l

    #append c hat vector
    c.hats[l] = proc_res$c.hat

    #Update bin_df_ based on previous layer discoveries
    bin_df_ <- bin_df_[!(bin_df_$L1 %in% S), ]

  }


  return(list("n"=n,
              "m.l"=m.l.vec,
              "S.list"=S.list,
              "c.hats"=c.hats))
}

##### Helper Functions #####

#' create_partition_info
#' Creates a partition of the 1D or 2D sample space of the pooled observations into mutually disjoint bins.
#' The bins form the leaves or finest-level regions for layer 1 of the aggregation tree.
#' @param df1 A \code{\link{data.frame}} with 1 or 2 columns ("X","Y") corresponding to the reference sample
#' @param df2 A \code{\link{data.frame}} with 1 or 2 columns ("X","Y") corresponding to the non-reference sample (want to find regions enriched for these observations)
#' @param m A positive integer specifying the number of bins in layer 1
#' @import dplyr ggplot2
# @return A \code{\link{list}} containing the the pooled observation \code{\link{data.frame}} (dat),
# a \code{\link{data.frame}} containing the segments/rectangles that define each bin and their layer 1 indices (bin.df),
# the breaks along each dimension for the bins
#' @export create_partition_info
create_partition_info <- function(df1,df2,m){

  N1 = nrow(df1)
  N2 = nrow(df2)

  df <- rbind(df1,df2)
  df$lab <- rep(c(0,1),times=c(N1,N2))

  nd <- ncol(df1)#length(a)

  if(nd==1){ #1D setting

    df$quant = with(df, cut_number(df$X, n = m))

    x_quant = aggregate(lab~quant,df,sum)
    n = as.vector(table(df$quant))
    breaks = quantile(df$X,probs = seq(0,1,length.out = m+1))

    out <- data.frame(indx = seq(m^nd),
                      x.bd = x_quant$quant,
                      #y.bd = NULL,
                      n = n,
                      x = x_quant$lab)
    first.axis = NA

  } else{ #2D setting

    #Determine dim with highest variance and cut first along that dim
    var.x <- var(df$X)
    var.y <- var(df$Y)
    first.axis <- ifelse(var.x>var.y,1,2)

    #Create breaks
    breaks = create.xy.breaks(dat=df,m=m,first.axis=first.axis)


    a <- as.list(rep(m,nd))
    b <- seq(m^nd)

    mat.indx <- matrix(NA,nrow=length(b),ncol=length(a))
    tmp.b <- b
    for(d in seq_along(a)){
      mat.indx[,length(a)-d+1] <- bitwShiftR(abs(bitwOr(1L,bitwShiftL(a[[d]]+tmp.b-1L,1L))%%bitwShiftL(a[[d]],2L)-
                                                   bitwShiftL(a[[d]],1L)),1L)+1L
      tmp.b <- ((tmp.b-1L) %/% a[[d]])+1L
    }

    indx_snake <- array(NA,unlist(a))
    indx_snake[mat.indx[b,]] <- b

    if(first.axis==1){
      x.bds <- NULL
      for(i in seq(length(breaks$x)-1)){
        if(i == 1){
          bd <- paste0("(",breaks$x[i],",",breaks$x[i+1],"]")
        } else{
          bd <- paste0("[",breaks$x[i],",",breaks$x[i+1],"]")
        }
        x.bds <- c(x.bds,bd)
      }

      x.bd <- rep(x.bds,each=length(breaks$y))
      y.bd <- unlist(lapply(breaks$y,function(x) names(x[["n"]])))
      n = unlist(lapply(breaks$y, '[[', "n"))[c(t(indx_snake))]
      x = unlist(lapply(breaks$y, '[[', "x"))[c(t(indx_snake))]

      out <- data.frame(indx = seq(m^nd),
                        x.bd = x.bd,
                        y.bd = y.bd[c(t(indx_snake))],
                        n = n,
                        x = x
      ) %>%
        dplyr::arrange(indx)
    } else{
      y.bds <- NULL
      for(i in seq(length(breaks$y)-1)){
        if(i == 1){
          bd <- paste0("(",breaks$y[i],",",breaks$y[i+1],"]")
        } else{
          bd <- paste0("[",breaks$y[i],",",breaks$y[i+1],"]")
        }
        y.bds <- c(y.bds,bd)
      }

      x.bd <- unlist(lapply(breaks$x,function(x) names(x[["n"]])))
      y.bd <- rep(y.bds,each=length(breaks$x))
      n = unlist(lapply(breaks$x, '[[', "n"))[c(t(indx_snake))]
      x = unlist(lapply(breaks$x, '[[', "x"))[c(t(indx_snake))]

      out <- data.frame(indx = seq(m^nd),
                        x.bd = x.bd[c(t(indx_snake))],
                        y.bd = y.bd,
                        n = n,
                        x = x
      ) %>%
        dplyr::arrange(indx)
    }

  }

  list(dat = df,
       N1=N1,
       N2=N2,
       nd=nd,
       m=m,
       bin.df = out,
       breaks=breaks,
       first.axis=first.axis)

}

#' create.xy.breaks
#' @import stats dplyr
#' @param dat input \code{\link{data.frame}} of observations
#' @param m number of bins
#' @param first.axis first axis to partition
create.xy.breaks <- function(dat,m,first.axis){

  if(first.axis==1){
    brk.x <- quantile(dat$X,seq(0,1,length.out = m+1))
    brk.y <- vector("list",m)

    for(i in seq(m)){
      x.subset = which(between(dat$X,brk.x[i],brk.x[i+1]))
      brk.y[[i]]$breaks = quantile(dat$Y[x.subset],seq(0,1,length.out = m+1))
      brk.y[[i]]$binnedpts = cut(dat$Y[x.subset],brk.y[[i]]$breaks,include.lowest = TRUE,dig.lab=13)
      brk.y[[i]]$x = aggregate(lab[x.subset]~brk.y[[i]]$binnedpts,data=dat,FUN=sum)$lab
      brk.y[[i]]$n = table(brk.y[[i]]$binnedpts)
    }
  }
  else{
    brk.y <- quantile(dat$Y,seq(0,1,length.out = m+1))
    brk.x <- vector("list",m)

    for(j in seq(m)){
      y.subset = which(between(dat$Y,brk.y[j],brk.y[j+1]))
      brk.x[[j]]$breaks = quantile(dat$X[y.subset],seq(0,1,length.out = m+1))
      brk.x[[j]]$binnedpts = cut(dat$X[y.subset],brk.x[[j]]$breaks,include.lowest = TRUE)
      brk.x[[j]]$x = aggregate(lab[y.subset]~brk.x[[j]]$binnedpts,data=dat,FUN=sum)$lab
      brk.x[[j]]$n = table(brk.x[[j]]$binnedpts)
    }
  }

  return(list(x=brk.x,y=brk.y,first.axis=first.axis))
}

#' matsplitter
#' @param M matrix
#' @param r number of rows to split
#' @param c number of columns to split
matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1 # %/% = divide and round up
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/(r*c)
  lapply(1:N, function(x) M[rci==x])
}

#' run_proc
#' @importFrom stats
#' @param l current layer
#' @param n bin size
#' @param partition_info_l updated partition information at current layer
#' @param theta0 proportion of cohort 2 to total
#' @param c_hat_prev count that controls FDR at previous layer l-1
#' @param alpha nominal FDR level
run_proc <- function(l,n,partition_info_l,theta0,c_hat_prev,alpha){

  Ll <- paste0("L",l)

  #Only keep Ll-1 and Ll

  #Update bin_df each time
  #Aggregate layer 1 x's to get layer (l) x's
  xlform = as.formula(paste("x ~", Ll))
  x_l = aggregate(xlform, data=partition_info_l, sum)$x
  n_l = 2^(l-1)*n
  m_l = length(x_l)

  a_l =  1/(m_l*log(m_l))

  #min_c = qbinom(alpha,n.l,theta0,lower.tail=FALSE) #Indices that are immediately ignored - add this to calculate
  x_l_valid = sort(unique(x_l[x_l>n_l*theta0]),decreasing = TRUE)

  #Compute p-values only for unique values of valid x's
  x_l_sort <- sort(x_l,decreasing = FALSE)

  #Get the rank based on sum(I(Xi>x)) = sum(G0i(X_i) <= c)
  rank_x_l_valid = vapply(x_l_valid,function(z) sum(x_l > z),numeric(1))

  pval_valid = compute_pval(l=l,
                            xs=x_l_valid,
                            theta0 = theta0,
                            n_l = n_l,
                            a_l = a_l,
                            c_prev = c_hat_prev)

  fdr.est = m_l*pval_valid/pmax(rank_x_l_valid,1L)
  c.hat = ifelse(length(which(fdr.est <= alpha))>0,
                 pval_valid[max(which(fdr.est <= alpha))],
                 a_l)

  if(length(which(fdr.est <= alpha))>0){
    #Map discoveries to L1
    rej_l = which(x_l > x_l_valid[max(which(fdr.est <= alpha))])
    rej_L1 = unique(partition_info_l$L1[partition_info_l[,get(Ll)]%in%rej_l])
  } else{
    rej_L1 = NULL
  }

  return(list(c.hat=c.hat,rej=rej_L1,m.l=m_l))
}

#' compute_pval
#' @importFrom stats
#' @param l current layer
#' @param xs vector of counts
#' @param theta0 proportion of cohort 2 to total
#' @param n_l bin size at current layer
#' @param a_l max rejection count at current layer
#' @param c_prev count that controls FDR at previous layer l-1
compute_pval <- function(l,xs,theta0,n_l,a_l,c_prev){

  f1 <- function(z){
    #sums over the dbinoms for vector of a_{1,r}'s...
    #need row products
    return(ifelse(!is.matrix(z),
                  sum(prod(dbinom(x=z,size=n_l/2,prob=theta0)),na.rm = TRUE),
                  sum(apply(dbinom(x=z,size=n_l/2,prob=theta0),1,prod),na.rm = TRUE)))
  }
  f2 <- function(w,cprev){
    #sum over the values of c_k's
    ##print(head(valid.counts(w,c.prev=c.prev)))
    val = ifelse(is.list(valid.counts(w,c.prev=cprev)),
                 sum(sapply(valid.counts(w,c.prev=cprev),f1),na.rm=TRUE),
                 ifelse(is.null(valid.counts(w,c.prev=cprev)),0,
                        f1(valid.counts(w,c.prev=cprev))))
    return(val)
  }

  if(l==1){
    out = pbinom(xs,n_l,theta0,lower.tail = FALSE)
    return(out)
  } else{

    z_prev_l = qbinom(c_prev,n_l/2,theta0,lower.tail = FALSE)-1L
    za_prev_l = qbinom(a_l,n_l,theta0,lower.tail = FALSE)-1L

    max.x = min(2*z_prev_l,za_prev_l)

    out = NULL

    for(j in seq_along(xs)){
      min.x = xs[j]
      #print(min.x)

      if(min.x + 1 > max.x){
        #print("prob = 0")
        prob = 0
        #print(prob)
      } else{
        #print("prob =/= 0")
        tmp_z = seq(min.x+1,max.x)
        prob = f2(tmp_z,z_prev_l)/pbinom(q=z_prev_l,size=n_l/2,prob=theta0)^2
        #print(prob)
      }
      out = c(out,prob)
      #print(out)
    }

    return(out)
  }

}

#' expand.mat
#' @param mat matrix of counts
#' @param vec vector of counts
expand.mat = function(mat, vec) {
  #### print(paste(nrow(mat),length(vec)))
  #### print(paste("nrow=",nrow(mat)))
  out = matrix(0, nrow = as.numeric(nrow(mat)) * as.numeric(length(vec)),
               ncol = as.numeric(ncol(mat) + 1)) #deal with integer overflow
  for (i in 1:ncol(mat)) out[, i] = mat[, i]
  out[, ncol(mat) + 1] = rep(vec, each = nrow(mat))
  return(out)
}

#' valid.counts
#' @param x count
#' @param c.prev count that controls FDR at previous layer l-1
valid.counts = function(x,c.prev){
  #Outputs a 2-column matrix of valid counts
  vec = seq(from = 0, to = c.prev)
  mat = matrix(vec, ncol = 1)

  mat  = expand.mat(mat, vec)
  mat = mat[rowSums(mat) %in% x, ]

  if(!is.matrix(mat)){
    return(mat)
  }
  else if(is.matrix(mat) & nrow(mat)<1){
    return(NULL)
  }
  else{
    #mat = mat[rowSums(mat) %in% x, ]
    idx = split(seq(nrow(mat)), rowSums(mat))
    return(lapply(idx, function(i, x) x[i, ], mat))
  }
}

