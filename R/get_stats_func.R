# -mean age at 50% maturity
# -RDI, RDD 
# -mean size of catch 
require(mizer)
require(tidyverse)

# Jon, I saw your comment on temperature about TBM. 
# Since temperature is not varying during sim, the effect of temperature are included in the rates and there is no need to input temperature anywhere

# RDI/RDD/Rmax
getSpawnStats <- function(object)
{

rddDat <- getRDI(object@params,object@n[dim(object@n)[1],,],object@n_pp[dim(object@n_pp)[1],])
rdiDat <- getRDD(object@params,object@n[dim(object@n)[1],,],object@n_pp[dim(object@n_pp)[1],])
rmaxDat <- object@params@species_params$R_max

return(rbind(rddDat,rdiDat,rmaxDat))
}

# Predation mortality
predMortDat <- getPredMort(object@params,object@n[dim(object@n)[1],,],object@n_pp[dim(object@n_pp)[1],])


#Function that returns the age of individuals at 50% maturity 

getAgeAtMatM2<-function (object, 
                         max_age = 200
                        ) {
  
  sim <- object
  species <- dimnames(sim@n)$sp
  
  idx_last<- dim(sim@n)[1]
  
  
  idx <- which(dimnames(sim@n)$sp %in% species)
  species <- dimnames(sim@n)$sp[idx]
  age <- seq(0, max_age, length.out = max_age*10)
  ws <- array(dim = c(length(species), length(age)), 
              dimnames = list(Species = species, Age = age))
  

    g <- getEGrowth(sim@params, sim@n[dim(sim@n)[1], , ], 
                           sim@n_pp[dim(sim@n)[1], ])
    
  # this loop fills ws with the size at age (could change age cells as they are linearyly spaced so not much resolution on fast growing species)
  for (j in 1:length(species)) {
    i <- idx[j]
    g_fn <- stats::approxfun(sim@params@w, g[i, ])
    myodefun <- function(t, state, parameters) {
      return(list(g_fn(state)))
    }
    ws[j, ] <- deSolve::ode(y = sim@params@species_params$w_min[i], 
                            times = age, func = myodefun)[, 2]
    
  }
  
  # Calculate size at maturation 
  # Get weight at maturation 
  # plug weight into age function to get age at maturation  
    
  #' RF Not 100% sure how this works
  #' psi_fn output the size when 50% of the species has matured, we should find the age at this size in ws with age_fn
  #' getting back warnings though
  #' when increasing the cell res of ws, getting different results for some species
  
  
  agemat<-rep(NA, length(species))
  names(agemat)<-species
  
  for (SpIdx in 1:(length(species))) {
    # i <- idx[j]
    psi_fn <- stats::approxfun(sim@params@psi[SpIdx,], dimnames(sim@params@psi)$w)
    
    age_fn<- stats::approxfun( ws[SpIdx, ], as.numeric(colnames(ws)))
    agemat[SpIdx]<-age_fn(psi_fn(.5))
    
  }
  
  
  return(agemat)
}


#' mean size of catch
getMeanSizeCatch <- function(object)
{
catchData <- getFMort(sim) # this function's unit is 1/year so ind caught per year
# but it doesn't check is there is fish biomass available to actually fish
# adding 0 to every cell with size > w_inf
w_inf <- sim@params@species_params$w_inf
for(SpIdx in 1:(length(w_inf)))
    catchData[,SpIdx, which(sim@params@w > w_inf[SpIdx])] <- 0

indCatchYear <- plyr::aaply(catchData,c(1,2),sum) # num of ind caught per year each year per species
sumCatchYear <- plyr::aaply(catchData,c(1,2), function(x){x%*%sim@params@w}) # ind * weigth caught summed over species per species per year

meanCatchData <- sumCatchYear/indCatchYear # mean size per species per year caught

return(meanCatchData[dim(meanCatchData)[1],]) # only sending last time step to stick with other function but since all time steps are saved you can tweak as you like
}


