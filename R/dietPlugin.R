# Used to produce 4D diets arrays and then plot
library(ggplot2)

# give it a mizer sim and it creates an array used in plotDietComp to get the diet data of the simulation's last time step
getDietComp<- function(sim, n = sim@n[dim(sim@n)[1],,],  n_pp = sim@n_pp[dim(sim@n_pp)[1],], 
                       feedinglevel=getFeedingLevel(object, n = n,n_pp = n_pp)){
  
  # initialisation to get it working with latest Mizer version
  object <- sim@params
  no_sp <- dim(object@species_params)[1]
  no_w <- length(object@w)
  no_w_full <- length(object@w_full)
  
  
  pred_kernel <- array(0, dim = c(no_sp, no_w, no_w_full),
                       dimnames = list(sp = object@species_params$species,
                                       w_pred = signif(object@w, 3),
                                       w_prey = signif(object@w_full, 3)))
  diet_comp<-array(0, c(no_sp, no_w, no_sp + 1, no_w_full),
                   dimnames=list( predator=as.character(object@species_params$species), pred_size = object@w,
                                  prey = c(as.character(object@species_params$species), "background"),
                                  prey_size = object@w_full))
  
  # Filling pred_kernel  
  
  Beta <- log(object@species_params$beta)
  sigma <- object@species_params$sigma
  # w_full has the weights from the smallest relevant plankton, to the largest fish
  x_full <- log(object@w_full)
  # We choose the origin of the x axis to be at the smallest plankton size
  x_full <- x_full - x_full[1]
  dx <- x_full[2] - x_full[1]
  # rr is the maximal log predator/prey mass ratio
  rr <- Beta + 3 * sigma
  ri <- floor(rr / dx)
  # might need to add some initialisation here
  # res@ft_pred_kernel_e <- matrix(0, nrow = no_sp, ncol = length(x_full))
  for (i in 1:no_sp) {
    # print(i)
    # We compute the feeding kernel terms and their fft.
    phi <- exp(-(x_full - Beta[i])^2 / (2 * sigma[i]^2))
    phi[x_full > rr[i]] <- 0
    phi[1] <- 0
    # Fourier transform of feeding kernel for evaluating available energy
    # res@ft_pred_kernel_e[i, ] <- fft(phi)
    # Fourier transform of feeding kernel for evaluating predation rate
    phi_p <- rep(0, no_w_full)
    phi_p[(no_w_full - ri[i] + 1):no_w_full] <- phi[(ri[i] + 1):2]
    # res@ft_pred_kernel_p[i, ] <- fft(phi_p)
    # Full feeding kernel array
    
    min_w_idx <- no_w_full - no_w + 1
    for (k in seq_len(no_w)) 
      pred_kernel[i, k, (min_w_idx - 1 + k):1] <- phi[1:(min_w_idx - 1 + k)]
    
  }
  
  # Starting the function
  #Biomass by species;
  n_total_in_size_bins<- sweep(n, 2, object@dw , "*")
  b_tot <- sweep(n_total_in_size_bins, 2, object@w , "*")
  
  #Biomass of resource as prey; scaled to reflect pred size kernel; might have to change if we start using interaction with resource spectrum like Hartvig et al. 2011
  
  #Note that we multiply the amount available by the availability parameter in the species parameter file 
  b_background <- (sweep( pred_kernel[,,], c(3), object@dw_full*object@w_full*n_pp, "*")) 
  # b_benthos <- sweep( object@pred_kernel[,,], c(3), object@dw_full*object@w_full*n_bb, "*") * object@species_params$avail_BB
  # b_algae <- sweep( object@pred_kernel[,,], c(3), object@dw_full*object@w_full*n_aa, "*") * object@species_params$avail_AA
  
  
  #Search rate *  feeding level * predator biomass
  b_background<- sweep(b_background, c(1,2), object@search_vol,"*") #Scale up by search volume
  b_background<- sweep(b_background, c(1,2), feedinglevel,"*") # Scale according to feeding level. Prey eaten: g prey / year / g predator
  b_background_tot<-sweep(b_background,c(1,2), b_tot, "*") # Prey eaten: total g prey/ year  (given predator biomass density)
  
  # b_benthos <- sweep(b_benthos, c(1,2), object@search_vol,"*") #Scale up by search volume
  # b_benthos <- sweep(b_benthos, c(1,2), feedinglevel,"*") # Scale according to feeding level. Prey eaten: g prey / year / g predator
  # b_benthos_tot <- sweep(b_benthos,c(1,2), b_tot, "*") # Prey eaten: total g prey/ year  (given predator biomass density)
  # 
  # b_algae <- sweep(b_algae, c(1,2), object@search_vol,"*") #Scale up by search volume
  # b_algae <- sweep(b_algae, c(1,2), feedinglevel,"*") # Scale according to feeding level. Prey eaten: g prey / year / g predator
  # b_algae_tot<-sweep(b_algae,c(1,2), b_tot, "*") # Prey eaten: total g prey/ year  (given predator biomass density)
  
  # Handy indices 
  # no_w<- dim(n)[2]
  # no_sp<- dim(n)[1]
  
  # Index of predator size classes 
  idx_sp<- object@w_full %in% object@w
  
  for(i in 1:no_w){
    for(j in 1:no_sp){    
      diet_comp[j,i,1:no_sp,idx_sp]<- sweep(sweep( b_tot, c(1), object@interaction[j, 1:no_sp], "*"), c(2), pred_kernel[j,i,idx_sp], "*")
    }
  }
  
  # Search rate *  feeding level * predator biomass
  diet_comp[,,1:no_sp,]<- sweep(sweep(sweep(diet_comp[,,1:no_sp,], c(1,2), object@search_vol,"*"), c(1,2),feedinglevel,"*"), c(1,2),b_tot,"*")  # Prey eaten: total g prey/ year  (given predator biomass density)
  
  # Store background eaten 
  diet_comp[,,no_sp+1,]<- b_background_tot
  # diet_comp_all[,,no_sp+2,]<- b_benthos_tot
  # diet_comp_all[,,no_sp+3,]<- b_algae_tot
  
  return(diet_comp)
  #Save in sim object; divide by the number of time steps, and add time up to get average across time 
  
} 

# give it a mizer sim and associated diet_comp
plotDietComp<-function(object,diet_comp = diet_comp, prey=dimnames(diet_comp)$prey, min_w=.001,
                       predator=dimnames(diet_comp)$predator, timeaverage=FALSE, print_it = T){
  
  prey_nam<-prey
  pred_nam<-predator
  
  out<-diet_comp 
  
  prey<-apply(out, c(1,2,3), FUN=sum) #Sum across size classess with in prey 
  tot<-apply(prey, c(1,2), FUN=sum) #Sum across prey species 
  
  prey_prop<-sweep(prey, c(1,2), tot, "/") # Get proportion of diet for each species
  
  no_pred<- length(dimnames(prey_prop)[[1]])
  no_pred_w<- length(dimnames(prey_prop)[[2]])
  no_prey<- length(dimnames(prey_prop)[[3]])
  
  #Stacked  bar chart 
  plot_dat<-expand.grid(dimnames(prey_prop)[[1]], dimnames(prey_prop)[[2]], dimnames(prey_prop)[[3]])
  colnames(plot_dat)<-c("predator","predsize","prey")
  plot_dat$predsize<-as.numeric(as.character(log10(as.numeric(as.character(plot_dat$predsize)))))
  
  plot_dat$value<- as.vector(prey_prop)
  
  
  species<-object@params@species_params$species
  wmin<-object@params@w[object@params@species_params$w_min_idx]
  wmax<-object@params@w[object@params@species_params$w_max_idx]
  
  for ( i in 1:length(species)){
    plot_dat$value[plot_dat$predator==species[i] & plot_dat$predsize < log10(wmin[i])]<- 0
    plot_dat$value[plot_dat$predator==species[i] & plot_dat$predsize > log10(wmax[i])]<- 0
  }
  
  
  plot_dat<-subset(plot_dat,predsize> log10(min_w))
  
  dsub<-plot_dat[ plot_dat$prey %in% prey_nam, ]
  dsub<-dsub[ dsub$predator %in% pred_nam, ]
  dsub[is.na(dsub)]<-0
  
  #nice colors
  no_sp <- length(unique(dsub$predator))
  colfunc <- colorRampPalette(c("firebrick3","darkturquoise", "orange"))
  colGrad <- colfunc(no_sp)
  backgroundColor <- colorRampPalette("salmon2")
  colGrad <- c(colGrad,backgroundColor(1))
  names(colGrad) <- c(1:no_sp, "background")
  
  p<-  ggplot(dsub) +
    geom_area(aes(x = predsize, y = value, fill = prey), position = 'stack') + 
    facet_wrap(~predator, ncol = 5) + 
    scale_fill_manual(name = "Species", values = colGrad)+ # color gradient
    scale_x_continuous(name = "log10 predator mass (g)") + 
    scale_y_continuous(name = "Proportion of diet by mass (g)") + 
    theme(legend.title=element_text(),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), 
          legend.position="bottom",
          strip.background = element_blank(),
          legend.key = element_rect(fill = "white"))+
    guides(fill=guide_legend(nrow = 1))+
    ggtitle(NULL)

  if(print_it) return(p)
  
}  