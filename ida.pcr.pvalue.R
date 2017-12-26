ida.pcr.pvalue <- function(x.pos,y.pos,data,graphEst,method="local",
                     y.notparent = FALSE, verbose=FALSE, all.dags=NA)
{
  ## Purpose: Estimate the causal effect of x on y; the graphEst and correlation
  ## matrix have to be precomputed; all DAGs can be precomputed;
  ## Orient undirected edges at x in a way so that no new collider
  ## is introduced
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - mcov: Covariance matrix that was used to estimate graphEst
  ## - graphEst: Fit of PC Algorithm (semidirected)
  ## - method: "local" - local (all combinations of parents in regr.)
  ##           "global" - all DAGs
  ## - y.notparent: if TRUE, the effect of x <- y is ignored;
  ##                (remove y from all parents set pa1 or pa2)
  ##                if FALSE, the effect of x <- y is set to zero
  ## - verbose: if TRUE, details on regressions that were used
  ## - all.dags: All DAGs in the format of function allDags; if this is
  ##   available, no new function call allDags is done
  ## ----------------------------------------------------------------------
  ## Value: causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 7 Jan 2010, 11:18
  ## Replaced lm.cov with pcr (Partial Component Regression)  
  ## Added code to select regression coefficient for the max variance component and compute its p-value

  tmpColl <- FALSE
  idx.name <- names(data)[x.pos]
  
  ## prepare adjMatrix and skeleton
  amat <- wgtMatrix(graphEst)
  amat[which(amat!=0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel!=0] <- 1
  
  if (method=="local") {
    ##############################
    ## local method
    ## Main Input: mcov, graphEst
    ##############################
    ## find unique parents of x
    wgt.est <- (wgtMatrix(graphEst)!=0)
    if (y.notparent) {
      ## Direct edge btw. x.pos and y.pos towards y.pos
      wgt.est[x.pos, y.pos] <- FALSE
    }
    tmp <- wgt.est-t(wgt.est)
    tmp[which(tmp<0)] <- 0
    wgt.unique <- tmp
    pa1 <- which(wgt.unique[x.pos,]!=0)
    if (y.pos %in% pa1) {
      ## y is parent of x -> zero effect
      beta.hat <- 0
    } else { ## y.pos not in pa1
      ## find ambiguous parents of x
      wgt.ambig <- wgt.est-wgt.unique
      pa2 <- which(wgt.ambig[x.pos,]!=0)
      if (verbose) {
        cat("\n\nx=",x.pos,"y=",y.pos,"\n")
        cat("pa1=",pa1,"\n")
        cat("pa2=",pa2,"\n")
      }
      
      ## estimate beta
      if (length(pa2)==0) {
        # beta.hat <- lm.cov(mcov,y.pos,c(x.pos,pa1))
        pcr.model <- pcr(Rainfall~., data=data[,c(x.pos,pa1,y.pos)])
        max.var.comp <-  which.max(explvar(pcr.model))
        var.hat <- explvar(pcr.model)[max.var.comp]
        beta.hat <- data.frame(pcr.model$coefficients)[idx.name, max.var.comp]
        # Compute p-value of the causal effect
        beta.hat <- computePValue(idx.name, beta.hat, data, c(x.pos, pa1), y.pos, 1)
        if(is.na(beta.hat)) var.hat <- NA

        if (verbose) cat("Fit - y:",y.pos,"x:",c(x.pos,pa1),
                         "|b.hat=",beta.hat,"\n")
      } else {
        ## at least one undirected parent
        beta.hat <- NA
        var.hat <- NA
        ii <- 1
        
        ## no member of pa2
        pa2.f <- pa2
        pa2.t <- NA
        tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
        if (!tmpColl) {
          # beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
          pcr.model <- pcr(Rainfall~., data=data[,c(x.pos,pa1,y.pos)])
          max.var.comp <-  which.max(explvar(pcr.model))
          beta.hat[ii] <- data.frame(pcr.model$coefficients)[idx.name, max.var.comp]
          var.hat[ii] <- explvar(pcr.model)[max.var.comp]
          # Compute p-value of the causal effect
          beta.hat[ii] <- computePValue(idx.name, beta.hat[ii], data, c(x.pos, pa1), y.pos, 1)
          if(is.na(beta.hat[ii])) var.hat[ii] <- NA
          if (verbose) cat("Fit1 - y:",y.pos,"x:",c(x.pos,pa1),
                           "|b.hat=",beta.hat[ii],"\n")
        }
        ## exactly one member of pa2
        for (i2 in 1:length(pa2)) {
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
          if (!tmpColl) {
            ii <-  ii+1
            if (y.pos %in% pa2.t) {
              beta.hat[ii] <- 0
            } else {
              # beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa2[i2]))
              pcr.model <- pcr(Rainfall~., data=data[,c(x.pos,pa1,pa2[i2],y.pos)])
              max.var.comp <-  which.max(explvar(pcr.model))
              beta.hat[ii] <- data.frame(pcr.model$coefficients)[idx.name, max.var.comp]
              var.hat[ii] <- explvar(pcr.model)[max.var.comp]
              # Compute p-value of the causal effect
              beta.hat[ii] <- computePValue(idx.name, beta.hat[ii], data, c(x.pos,pa1,pa2[i2]), y.pos, 1)
              if(is.na(beta.hat[ii])) var.hat[ii] <- NA
              if (verbose) cat("Fit2 - y:",y.pos,"x:",c(x.pos,pa1,pa2[i2]),
                               "|b.hat=",beta.hat[ii],"\n")
            } ## if (y.pos %in% pa2.t)
          } ## if (!tmpColl)
        } ## for (i2 in 1:length(pa2))
        
        ## higher order subsets of pa2
        if (length(pa2)>1) {
          for (i in 2:length(pa2)) {
            pa.tmp <- combn(pa2,i,simplify=TRUE)
            n.comb <- ncol(pa.tmp)
            for (j in 1:n.comb) {
              pa2.f <- setdiff(pa2,pa.tmp[,j])
              pa2.t <- pa.tmp[,j]
              tmpColl <- check.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)
              if (!tmpColl) {
                ii <- ii+1
                if (y.pos %in% pa2.t) {
                  beta.hat[ii] <- 0
                } else {
                  # beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa.tmp[,j]))
                  pcr.model <- pcr(Rainfall~., data=data[,c(x.pos,pa1,pa.tmp[,j],y.pos)])
                  max.var.comp <-  which.max(explvar(pcr.model))
                  beta.hat[ii] <- data.frame(pcr.model$coefficients)[idx.name, max.var.comp]
                  var.hat[ii] <- explvar(pcr.model)[max.var.comp]
                  # Compute p-value of the causal effect
                  beta.hat[ii] <- computePValue(idx.name, beta.hat[ii], data, c(x.pos,pa1,pa.tmp[,j]), y.pos, 1)
                  if(is.na(beta.hat[ii])) var.hat[ii] <- NA
                  if (verbose) {
                    cat("Fit3 - y:",y.pos,"x:",c(x.pos,pa1,pa.tmp[,j]),
                        "|b.hat=",beta.hat[ii],"\n")
                  } ## if (verbose)
                } ## if (y.pos %in% pa2.t)
              } ## if (!tmpColl)
            } ## for (j in 1:n.comb)
          } ## for (i in 2:length(pa2))
        } ## if (length(pa2)>1)
      } ## if (length(pa2) == 0)
    } ## if (y.pos %in% pa1)
    
  } else {
    ##############################
    ## global method
    ## Main Input: mcov, graphEst
    ##############################
    p <- numNodes(graphEst)
    am.pdag <- wgtMatrix(graphEst)
    am.pdag[am.pdag!=0] <- 1
    if (y.notparent) {
      ## Direct edge btw. x.pos and y.pos towards y.pos
      am.pdag[x.pos, y.pos] <- 0
    }
    
    ## find all DAGs if not provided externally
    if (is.na(all.dags)) {
      ad <- allDags(am.pdag,am.pdag,NULL)
    } else {
      ad <- all.dags
    }
    n.dags <- nrow(ad)
    beta.hat <- rep(NA,n.dags)
    for (i in 1:n.dags) {
      ## compute effect for every DAG
      gDag <- as(matrix(ad[i,],p,p),"graphNEL")
      ## path from y to x
      ## rev.pth <- sp.between(gDag,as.character(y.pos),
      ##                    as.character(x.pos))[[1]]$path
      ## if (length(rev.pth)>1) {
      ## if reverse path exists, beta=0
      ##  beta.hat[i] <- 0
      ## } else {
      ## path from x to y
      ##       pth <- sp.between(gDag,as.character(x.pos),
      ##                       as.character(y.pos))[[1]]$path
      ##   if (length(pth)<2) {
      ## sic! There is NO path from x to y
      ##   beta.hat[i] <- 0
      ## } else {
      ## There is a path from x to y
      wgt.unique <- t(matrix(ad[i,],p,p)) ## wgt.est is wgtMatrix of DAG
      pa1 <- which(wgt.unique[x.pos,]!=0)
      if (y.pos %in% pa1) {
        beta.hat[i] <- 0
      } else {
        # beta.hat[i] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
        pcr.model <- pcr(Rainfall~., data=data[,c(x.pos,pa1,y.pos)])
        max.var.comp <-  which.max(explvar(pcr.model))
        beta.hat[ii] <- data.frame(pcr.model$coefficients)[idx.name, max.var.comp]
        # Compute p-value of the causal effect
        beta.hat[ii] <- computePValue(idx.name, beta.hat[ii], data, c(x.pos,pa1), y.pos, 1)
        if (verbose) cat("Fit - y:",y.pos,"x:",c(x.pos,pa1),
                         "|b.hat=",beta.hat[i],"\n")
      } ## if (y.pos %in% pa1)
      ##  } ## if length(pth)
      ## } ## if rev.pth
    } ## for n.dags
  } ## if (method = "local")
  # print(beta.hat)
  if(length(beta.hat)==1 && is.na(beta.hat))
    beta.hat
  else{
    # idx <- which.max(var.hat)
    beta.hat
  }
}