# Functions for HMM domain determination

runBioHMM.2 <- function (mval, datainfo, useCloneDists = TRUE, criteria = "AIC", 
                       delta = NA, var.fixed = FALSE, epsilon = 1e-06, numiter = 30000) 
{
  crit = TRUE
  if (criteria == "AIC") {
    aic = TRUE
    bic = FALSE
  } else if (criteria == "BIC") {
    bic = TRUE
    aic = FALSE
  } else crit = FALSE
  if ((crit == 1) || (crit == 2)) {
    if (criteria == "BIC") {
      if (is.na(delta)) {
        delta <- c(1)
      }
    }
    res <- try(fit.model.2(
      obs = mval, datainfo = datainfo, useCloneDists = useCloneDists, 
      aic = aic, bic = bic, 
      delta = delta, var.fixed = var.fixed, epsilon = epsilon, 
      numiter = numiter
    ))$out.list$state
  }
  else {
    cat("You must enter AIC or BIC for the criteria argument\n")
  }
}


fit.model.2 <- function (obs, datainfo = NULL, useCloneDists = TRUE, 
                       aic = TRUE, bic = FALSE, delta = 1, var.fixed = FALSE, 
                       epsilon = 1e-06, numiter = 30000) 
{
  library(cluster)
  kb <- datainfo$start
  if (useCloneDists) {
    dists.pre = kb[2:length(kb)] - kb[1:(length(kb) - 1)]
    dists = dists.pre/(max(dists.pre))
  } else {
    dists <- rep(1, length(kb))
  }
  covars <- as.matrix(dists)
  obs.ord <- obs[order(kb)]
  kb.ord <- kb[order(kb)]
  ind.nonna <- which(!is.na(obs.ord))
  data <- obs.ord[ind.nonna]
  kb <- kb.ord[ind.nonna]
  numobs <- length(data)
  if (numobs > 5) {
    temp2 <- clara(data, 2)
    init.mean.two <- temp2$medoids
    init.var.two <- vector()
    if (var.fixed == FALSE) {
      for (i in 1:2) {
        if (length(temp2$data[temp2$clustering == i]) > 
            1) 
          init.var.two[i] <- log(sqrt(var(temp2$data[temp2$clustering == 
                                                       i])))
        else init.var.two[i] <- log(0.5)
      }
    } else {
      init.var.two[1:2] <- log(sqrt(var(data)))
    }
    z2.init <- c(init.mean.two[, 1], init.var.two, -1, -3.6, 
                 -3.6, 0)
    z.pre <- run.nelder(numobs, z2.init, data, 
                        covars, var.fixed, epsilon, numiter, i)
    if (!is.nan(z.pre$x[1])) {
      z2 <- find.param.two(z.pre, var.fixed)
    } else {
      z2 <- NULL
    }
    if (aic) {
      factor <- 2
    } else if (bic) {
      factor <- log(numobs) * delta
    } else {
      stop("No criteria selected")
    }
    z <- z2
    nstates <- 2
    trans.mat <- list()
    for (j in 1:(length(data) - 1)) {
      trans.mat[[j]] <- z$LH.trans + exp(-(covars[j, 
                                                  1]^(z$rate1)) * prod(covars[j, -1])) * z$RH.trans
    }
    Vit.seg <- Viterbi.two(data, 
                           z, trans.mat)
    maxstate.unique <- unique(Vit.seg)
    mean <- rep(0, length(data))
    var <- rep(0, length(data))
    for (m in 1:length(maxstate.unique)) {
      mean[Vit.seg == maxstate.unique[m]] <- mean(data[Vit.seg == 
                                                         maxstate.unique[m]])
      var[Vit.seg == maxstate.unique[m]] <- var(data[Vit.seg == 
                                                       maxstate.unique[m]])
    }
    out <- cbind(matrix(Vit.seg, ncol = 1), matrix(mean, 
                                                   ncol = 1), matrix(var, ncol = 1))
    out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
    out.all[ind.nonna, 1:3] <- out
    out.all[, 4] <- obs.ord
    out.all <- as.data.frame(out.all)
    dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
    numstates <- length(unique(Vit.seg))
  } else {
    out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
    out.all[ind.nonna, 1] <- c(rep(1, numobs))
    out.all[ind.nonna, 2] <- c(rep(mean(obs.ord), numobs))
    out.all[ind.nonna, 3] <- c(rep(var(obs.ord), numobs))
    out.all[ind.nonna, 4] <- obs.ord
    out.all <- as.data.frame(out.all)
    dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
    numstates = 1
  }
  list(out.list = out.all, nstates.list = numstates)
}


runBioHMM.3 <- function (mval, datainfo, useCloneDists = TRUE,
                         criteria = "AIC", delta = NA, var.fixed = FALSE,
                         epsilon = 1e-06, numiter = 30000)
{
  crit = TRUE
  if (criteria == "AIC") {
    aic = TRUE
    bic = FALSE
  } else if (criteria == "BIC") {
    bic = TRUE
    aic = FALSE
  } else crit = FALSE
  if ((crit == 1) || (crit == 2)) {
    if (criteria == "BIC") {
      if (is.na(delta)) {
        delta <- c(1)
      }
    }
    res <- try(fit.model.3(obs = mval, datainfo = datainfo, useCloneDists = useCloneDists, aic = aic, bic = bic, delta = delta, var.fixed = var.fixed, epsilon = epsilon, numiter = numiter))$out.list$state
  }
  else {
    cat("You must enter AIC or BIC for the criteria argument\n")
  }
}

fit.model.3 <- function (obs, datainfo = NULL, useCloneDists = TRUE,
                         aic = TRUE, bic = FALSE, delta = 1, var.fixed = FALSE,
                         epsilon = 1e-06, numiter = 30000)
{
  kb <- datainfo$start 
  if (useCloneDists) {
    dists.pre = kb[2:length(kb)] - kb[1:(length(kb) - 1)] 
    dists = dists.pre/(max(dists.pre)) 
  } else {
    dists <- rep(1, length(kb))
  }
  covars <- as.matrix(dists) 
  obs.ord <- obs[order(kb)] 
  kb.ord <- kb[order(kb)] 
  ind.nonna <- which(!is.na(obs.ord)) 
  data <- obs.ord[ind.nonna] 
  kb <- kb.ord[ind.nonna] 
  numobs <- length(data) 
  if (numobs > 5) { 
    temp3 <- clara(data, 3)
    init.mean.three <- vector(length = 3)
    init.mean.three <- temp3$medoids
    init.var.three <- vector(length = 3)
    if (var.fixed == FALSE) {
      for (i in 1:3) {
        if (length(temp3$data[temp3$clustering == i]) > 1) {
          init.var.three[i] <- log(sqrt(var(temp3$data[temp3$clustering == i])))
        } else {
          init.var.three[i] <- log(0.5)
        }
      }
    } else {
      init.var.three[1:3] <- log(sqrt(var(data)))
    }
    z3.init <- c(init.mean.three[, 1], init.var.three, -0.7, 
                 -0.7, -3.6, -3.6, -3.6, -3.6, -3.6, -3.6, 0)
    z.pre <- run.nelder(numobs, z3.init, data, covars, var.fixed, epsilon, numiter, i)
    if (!is.nan(z.pre$x[1])) {
      z3 <- find.param.three(z.pre, var.fixed)
    } else {
      z3 <- NULL
    }
    z <- z3
    trans.mat <- list()
    for (j in 1:(length(data) - 1)) {
      trans.mat[[j]] <- z$LH.trans + exp(-(covars[j, 1]^(z$rate1)) * prod(covars[j, -1])) * z$RH.trans
    }
    Vit.seg <- Viterbi.three(data, z, trans.mat)
    maxstate.unique <- unique(Vit.seg)
    mean <- rep(0, length(data))
    var <- rep(0, length(data))
    for (m in 1:length(maxstate.unique)) {
      mean[Vit.seg == maxstate.unique[m]] <- mean(data[Vit.seg == 
                                                         maxstate.unique[m]])
      var[Vit.seg == maxstate.unique[m]] <- var(data[Vit.seg == maxstate.unique[m]])
    }
    out <- cbind(matrix(Vit.seg, ncol = 1), matrix(mean, ncol = 1), matrix(var, ncol = 1))
    out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
    out.all[ind.nonna, 1:3] <- out
    out.all[, 4] <- obs.ord
    out.all <- as.data.frame(out.all)
    dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
    numstates <- length(unique(Vit.seg))
  } else {
    out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
    out.all[ind.nonna, 1] <- c(rep(1, numobs))
    out.all[ind.nonna, 2] <- c(rep(mean(obs.ord), numobs))
    out.all[ind.nonna, 3] <- c(rep(var(obs.ord), numobs))
    out.all[ind.nonna, 4] <- obs.ord
    out.all <- as.data.frame(out.all)
    dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
    numstates = 1
  }
  list(out.list = out.all, nstates.list = numstates)
}


delete.het <- function(data, r6 = dm.r6){
  if (!r6){
    euc.coords <- GRanges(
      seqnames = Rle(c("2L", "2R", "3L", "3R", "X")),
      ranges = IRanges(
        start = c(1, 1600000, 1, 1, 1),
        end = c(22000000, 21146700, 22900000, 27900000, 22422800)
      )
    )
  } else {
    euc.coords <- GRanges(
      seqnames = Rle(c("2L", "2R", "3L", "3R", "X")),
      ranges = IRanges(
        start = c(1, 5712495, 1, 4174279, 103614),
        end = c(22000000, 25259177, 22906900, 32074278, 23020964)
      )
    )
  }
  
  data.coords <- GRanges(
    seqnames = Rle(data$chr),
    ranges = IRanges(
      start = data$start,
      end = data$end
    )
  )
  
  data.x.euc <- subsetByOverlaps(data.coords, euc.coords)
  return(merge(data, 
               data.frame(chr = seqnames(data.x.euc),
                          start = start(data.x.euc),
                          end = end(data.x.euc)),
               by = c("chr", "start", "end")) %>% arrange(chr, start))
}

hmm.2.wrapper <- function(DATAs, by.chr = T){
  # A function to perform HMM and assign domain score (0 or 1)
  hmm2 <- function(exp_list) {
    data.domains <- lapply(exp_list, function(chromdf){
    BioHMM.output <- runBioHMM.2(mval = chromdf$DamID.value, datainfo = chromdf, useCloneDists = T)
    classifier <- aggregate(x=chromdf$DamID.value, by=list(BioHMM.output), FUN=quantile, probs = 0.75)
    print(classifier)
    if (nrow(classifier) == 2){
      signals <- BioHMM.output == classifier[which.max(classifier$x),1]
      nonsignals <- BioHMM.output == classifier[which.min(classifier$x),1]
      BioHMM.output[signals] <- 1
      BioHMM.output[nonsignals] <- 0
      chromdf$domain <- BioHMM.output
    } else {
      print("There's a problem with classification")
      return()
    }
    if (by.chr){
      return(chromdf)
    } else{
      chromdf <- merge(DATA[, 1:lastcol], chromdf[, -c(2:4)], by = "ID")
      return(chromdf)
    }
  })
    if (!by.chr){
      return(data.domains)
    } else {
      return(do.call(rbind, data.domains))
    }
  }
  
  if (by.chr){
    DATA.HMM <- lapply(DATAs, hmm2)
  } else {
    DATA.HMM <- hmm2(DATAs)
  }
  # remove null elements from lists
  return(lapply(DATA.HMM, function(lll) Filter(Negate(function(x) is.null(unlist(x))), lll)))
  
}

hmm.3.wrapper <- function(DATAs, by.chr = T){
  # A function to perform HMM and assign domain score (0 or 1)
  hmm3 <- function(exp_list){
    data.domains <- lapply(exp_list, function(chromdf){
      BioHMM.output <- runBioHMM.3(mval = chromdf$DamID.value, datainfo = chromdf, useCloneDists = T)
      classifier <- tapply(chromdf$DamID.value, BioHMM.output, FUN = median)
      amounts <- tapply(chromdf$DamID.value, BioHMM.output, FUN = length)
      print(classifier)
      print(amounts)
      print("");
      print("-------------------------")
      
      #sort classifier from min to max and chooses max to represent true signals
      mvp <- order(classifier)[3]
      signals <- BioHMM.output == mvp
      nonsignals <- BioHMM.output != mvp
      BioHMM.output[signals] <- 1
      BioHMM.output[nonsignals] <- 0
      
      chromdf$domain <- BioHMM.output
      
      if (by.chr){
        return(chromdf)
      } else{
        chromdf <- merge(DATA[, 1:lastcol], chromdf[, -c(2:4)], by = "ID")
        return(chromdf)
      }
      print("=========================")
    })
    
    if (!by.chr){
      return(data.domains)
    } else {
      return(do.call(rbind, data.domains))
    }
  }
  
  if (by.chr){
    DATA.HMM <- lapply(DATAs, hmm3)
  } else {
    DATA.HMM <- hmm3(DATAs)
  }
}



