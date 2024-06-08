#!/usr/bin/env Rscript

#Virtual Colon
#
#The front-end function that simulates a mammalian virtual colon

#libraries
require(cplexAPI)
require(sybil)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
require(BacArena)

#set working directory
setwd("/***/Virtual_Colon_Publication")
nameslist <- readRDS("InputData/vmh_translation_curated.RDS")

#A data frame that contains the name, description,
#input and diffusion values of each compound to be added in the envinroment.
substances  <-  readRDS("InputData/diet_diffusion.RDS")

#A matrix that corresponds to the blueprints of arena. The values correspond
#to each layer of the colon
#(1: lumen, 2: outer mucus layer, 3: inner mucus layer,
#4: host cell layer 5: host space)
m_setup <-  readRDS("InputData/matrix_colon_layers.RDS")

#A matrix that corresponds to the position of host cells (value 1, otherwise 0).
m_p <-  readRDS("InputData/matrix_colonocytes_position.RDS")

#A matrix that corresponds to the inner mucus layer
#where a diffusion barrier will be introduced (value 1, otherwise 0).
m_occ  <-  readRDS("InputData/matrix_diffusion_barrier.RDS")

#A list of host models in sybil format.
colonocyte  <-  readRDS("InputData/colonocytes.RDS")

#A list of bacterial models in sybil format.
bacteria  <-  readRDS("InputData/SIHUMIx_for_virtual_colon.RDS")

#A data frame similar to the substances one, 
#where the diffusion parameters of a great selection of compounds
#were set to the default diffusion value of BacArena (glucose)
generic_diffusion  <-  readRDS("InputData/generic_diffusion.RDS")

# A data frame that contains the combination of host and bacterial species
simulation_sequence <- readRDS("InputData/simulations_sequence.RDS")

# reading array number to be used to select models 
j <-  as.numeric(commandArgs(trailingOnly = TRUE))

#A number giving the amount of iterations.
timeloops  <-  6

#A number giving the time (in h) per iteration (as in BacArena)
tstep  <-  1

#A number giving the number of replications.
replications  <-  10

# A binary argument. Set to TRUE if no host models are needed.
no_host <- F

# A binary argument. Set to TRUE if no bacterial models are needed.
no_bacterium <- F



#Diffusion functions - START
  
get_Dgrid_list <- function(sub_dat, M_setup){
  Dgrid_list <- list()
  for (i in 1:nrow(sub_dat)) { # create matrix with spatial diffusion constant for every substance
    M_diff <- M_setup
    D_fec <- sub_dat$D_fec[i]*1e-5; D_omu <- sub_dat$D_omu[i]*1e-5 # numbers are in 10^5 cm^2/s
    D_imu <- sub_dat$D_imu[i]*1e-5; D_hum <- sub_dat$D_hum[i]*1e-5
    D_inter <- sub_dat$D_inter[i]*1e-5
    M_diff[which(M_diff == 1)] <- D_fec; M_diff[which(M_diff == 2)] <- D_omu
    M_diff[which(M_diff == 3)] <- D_imu; M_diff[which(M_diff == 5)] <- D_hum
    M_diff[which(M_diff == 4)] <- D_inter
    prop2d <- list(); class(prop2d) <- "prop.2D"
    prop2d$x.mid <- M_diff; prop2d$y.mid <- M_diff
    prop2d$x.int <- rbind(M_diff, rep(0, ncol(M_diff)))
    prop2d$y.int <- cbind(M_diff, rep(0, nrow(M_diff)))
    Dgrid_list[i] <- list(prop2d)
  }
  return(Dgrid_list)
}

get_Vgrid_list <- function(sub_dat, M_setup){
  Vgrid_list <- list()
  for (i in 1:nrow(sub_dat)) { # create matrix with spatial diffusion constant for every substance
    M_adv <- M_setup
    M_zero <- matrix(0, nrow = nrow(M_adv) , ncol = ncol(M_adv))
    A_fec <- sub_dat$A_fec[i]; A_omu <- sub_dat$A_omu[i] # numbers are in cm/s
    A_imu <- sub_dat$A_imu[i]; A_hum <- sub_dat$A_hum[i]
    A_inter <- sub_dat$A_inter[i]*1e-5
    M_adv[which(M_adv == 1)] <- A_fec; M_adv[which(M_adv == 2)] <- A_omu
    M_adv[which(M_adv == 3)] <- A_imu; M_adv[which(M_adv == 5)] <- A_hum
    M_adv[which(M_adv == 4)] <- A_inter
    prop2d <- list(); class(prop2d) <- "prop.2D"
    prop2d$x.mid <- M_zero; prop2d$y.mid <- M_adv
    prop2d$x.int <- rbind(M_zero, rep(0, ncol(M_adv)))
    prop2d$y.int <- cbind(M_adv, M_adv[,ncol(M_adv)])
    Vgrid_list[i] <- list(prop2d)
  }
return(Vgrid_list)
}

get_Diffmat_list <- function(sub_dat, M_setup){
  Diffmat_list <- list()
  for (i in 1:nrow(sub_dat)) {
    if (!is.na(sub_dat$init_conc[i])) {
      init_diffmat <- sub_dat$init_area[i]
      area <- as.numeric(unlist(strsplit(",", x = init_diffmat, fixed = T)))
      init_diffmat <- M_setup; init_diffmat[which(!init_diffmat %in% area)] <- 0
    }else{
      init_diffmat <- matrix(0, nrow = nrow(M_setup), ncol = ncol(M_setup))
    }
    Diffmat_list[[i]] <- init_diffmat
  }
  return(Diffmat_list)
}

nutrientPDE <- function(t, y, parms)  {
  with(as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- ReacTran::tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid)$dC
    return(list(dCONC))
  })
}

mucusPDE <- function(t, y, parms)  {
  with(as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- ReacTran::tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid)$dC
    return(list(dCONC))
  })
}

#Diffusion functions - END

# Diffusion of compounds
# For specific substances
dgrid.list <- get_Dgrid_list(substances, m_setup)  # diffusion constants
vgrid.list <- get_Vgrid_list(substances, m_setup)  # advection constants
diffmat.list <- get_Diffmat_list(substances, m_setup)

# Generic diffusion as glucose. U only in 4,5, E only in 1,2,3,4
dgrid.diff.list <- get_Dgrid_list(generic_diffusion, m_setup)  # diffusion constants
vgrid.diff.list <- get_Vgrid_list(generic_diffusion, m_setup)  # advection constants
diffmat.diff.list <- get_Diffmat_list(generic_diffusion, m_setup)

# Dimensions of arena
nrow <- dim(m_setup)[1] # x axis of matrix
ncol <- dim(m_setup)[2] # y axis of matrix

# Computing bacterial abundances
# Total corrected Caecum-Colon-Rectum Volume of host
# in cm(3) PMID: 20007641
volume <- 0.37

# Lumen area : 60 * 24 over 179 * 24 (total  nrow * ncol)
lumen <- (60 * 24)/(179 * 24)
out <- (37 * 24)/(179 * 24)
inn <- (19 * 24)/(179 * 24)
# Microbial input 200μl with 10^8 cfu/ml 
# (personal communication with PI of Geva-Zatorsky2017)  	 
input <- 1e8*0.2

simlist <- vector("list",replications)
savemain <- ls()
savemain <- c(savemain,"p","savemain")

for (p in 1:replications) {
  
  # Initialization of arena
  arena <- BacArena::Arena(m = nrow, n = ncol)
  arena@tstep <- tstep
  # OCCUPYM
  arena@occupyM <- m_occ
  # Volume of arena assuming 3D
  myvolume <- sqrt((arena@Lx * arena@Ly)) * arena@Lx * arena@Ly

  # Adding generic diffusion
  arena <- BacArena::addSubs(arena,
                             template = T,
                             mediac = generic_diffusion$Exchange[1:1500],
                             diffmat = diffmat.diff.list[1:1500],
                             Dgrid = dgrid.diff.list[1:1500],
                             Vgrid = vgrid.diff.list[1:1500],
                             pde = as.list(generic_diffusion$pde[1:1500]),
                             unit = 'mM',
                             addAnyway = T,
                             smax = as.numeric(generic_diffusion$Input_mM[c(1:1500)]),
                             add = T)

  # Adding dietary compounds to arena
  arena <- BacArena::addSubs(arena,
                             template = T,
                             mediac = substances$Exchange[1:92],
                             diffmat = diffmat.list[1:92],
                             Dgrid = dgrid.list[1:92],
                             Vgrid = vgrid.list[1:92],
                             pde = as.list(substances$pde[1:92]),
                             unit = 'mM',
                             addAnyway = T,
                             smax = as.numeric((substances$Input_mol[1:92]*0.1)/volume),
                             add = F)

  # Adding blood, mucus and extra compounds to arena
  arena <- BacArena::addSubs(arena,
                             template = T,
                             mediac = substances$Exchange[c(93:135)],
                             diffmat = diffmat.list[c(93:135)],
                             Dgrid = dgrid.list[c(93:135)],
                             Vgrid = vgrid.list[c(93:135)],
                             pde = as.list(as.character(substances$pde[c(93:135)])),
                             unit = 'mM',
                             addAnyway = T,
                             smax = as.numeric(substances$Input_mM[c(93:135)]*0.1),
                             add = T)

  # Simulations - only diffusion simulated
  sim <- BacArena::simEnv(arena, time = 1, continue = T, sec_obj = "mtf", with_shadow = F, verbose = F)


  # Introduction of host model
  # 70 cells should do the work of 1000 Bionumbers ID 110648
  # One cell starts with standard weight but its max weight can be up to
  # (1000/70) more than the standard max weight.
  # Scaling of human size based on bacterial one
  # Calculations based on Bionumbers ID 108890
  # 13 micrometer diameter => 132.73229 micrometer^2 (circle)
  # bacterial: 4.42 m^2 (default)
  # Simulations - adding host models
  if (no_host == F) {
    recon <- colonocyte[[which(names(colonocyte) == simulation_sequence$Sample[j])]] 	 
    t.col.ratio = (132.73229/4.42)
    HUM <- BacArena::Human(model = recon, growtype = "linear", speed = 0,
                           setAllExInf = T,
                           cellweight_sd =  0.132 * t.col.ratio,
                           cellweight_mean = 0.4890 * t.col.ratio,
                           maxweight =  1.172 * t.col.ratio * (1000/70))
    sim <- BacArena::addOrg(object = sim, specI = HUM, posmat = m_p)
  }

  # Introduction of bacterial models
  for (mybac in 1:length(bacteria)) {
    # selecting bacterium
    bac <- bacteria[[mybac]][[1]]
    bac@mod_name <- names(bacteria)[mybac]
    bac@mod_id <- names(bacteria)[mybac]
    abundance_factor <- bacteria[[mybac]][[2]]
    # Computing bacterial abundances
    final <- input * myvolume * abundance_factor
    finallumen <- final * lumen
    finalout <- final * out
    finalinn <- final * inn
    for (m in 1:length(nameslist$met_id)) {
      if (length(grep(pattern = nameslist$reac_id_gapseq[[m]], x = bac@react_id, value = F, fixed = T)) > 0) { 
        bac@react_id[grep(pattern = nameslist$reac_id_gapseq[[m]], x = bac@react_id, value = F, fixed = T)] <- nameslist$reac_id_vmh_bigg[[m]]}}
                
    ex <- findExchReact(bac)
    ex <- ex[grep("^EX_", ex@react_id),]
      for (i in seq_along(ex@react_id)) {
      	   rea <- gsub("\\(e\\)","\\(d\\)",ex@react_id[i])
      	   met <- ex@met_id[i]
      	   ub  <- 0
      	   lb  <- ex@lowbnd[i]
      	   name <- paste(react_name(bac)[ex@react_pos[i]],"(diet side)")
      	   bac <- addReact(bac, id = rea, reactName = name,met = met, Scoef = -1, reversible = T, ub = ub, lb = lb)
      	   bac@react_attr[bac@react_num,] <- react_attr(bac)[ex@react_pos[i],]
      	}
      	
    BAC_l <- BacArena::Bac(model = bac, type = paste0("lumen_",bac@mod_id),
                           setAllExInf = T )
    
    BAC_i <- BacArena::Bac(model = bac, type = paste0("inner_",bac@mod_id),
                           setAllExInf = T )
    
    BAC_o <- BacArena::Bac(model = bac, type = paste0("outer_",bac@mod_id),
                           setAllExInf = T )
    
    sim <- BacArena::addOrg(sim,BAC_l, m = 60, m0 = 1, amount = round(finallumen,0))
    sim <- BacArena::addOrg(sim,BAC_o, m = 97, m0 = 61, amount = round(finalout,0))
    sim <- BacArena::addOrg(sim,BAC_i, m = 116, m0 = 98, amount = round(finalinn,0))
  }  
  # Set by chance a point in outer mucus as reference
  reference_outer <- sample(grep(pattern = "2",x = m_setup, fixed = T),size = 1)

  # Vector, where the details of reference outer mucus are saved
  rf <- vector()

  # Call of position points of inner mucus
  inner_points <- grep(pattern = "3",x = m_setup,fixed = T)

  # Start of sim loop
  loops <- timeloops

  for (t in 1:loops)  {

      # if t=2 1. the reference outer mucus is measured
      # if t>2 Check if the composition of at least one inner compound
      # is equal or less to the respective outer, if yes no restriction movement

      if (t == 2) {
        for (i in substances$Exchange[97:104]) {
          rf[[i]] <- (unlist(BacArena::extractMed(object = sim,
                                                  time = 1,
                                                  mediac = i)))[reference_outer]}}
      if (t > 2) {

        for (z in inner_points) {
          y = NULL
          x = NULL
          for (i in substances$Exchange[97:104]) {

            if ((unlist(BacArena::extractMed(object = sim,time = t - 1,mediac = i)))[z] <= rf[[i]]) {y[i] = TRUE}
            else {y[i] = FALSE}}

          if (TRUE  %in% y) {

            if (z > nrow(m_setup) & (z <= ((nrow(m_setup)*ncol(m_setup)) - (nrow(m_setup))))) {
              x[5] = (-nrow(m_setup) + z);x[3] = (-nrow(m_setup) + z - 1);x[4] = (-nrow(m_setup) + z + 1) ;
              x[6] = (nrow(m_setup) + z);x[7] = (nrow(m_setup) + z + 1);x[8] = (nrow(m_setup) + z - 1)}

            if (z > ((nrow(m_setup)*ncol(m_setup)) - (nrow(m_setup)))) {
              x[5] = (-nrow(m_setup) + z);x[3] = (-nrow(m_setup) + z - 1);x[4] = (-nrow(m_setup) + z + 1)}

            if (z < nrow(m_setup)) {x[5] = (nrow(m_setup) + z);x[3] = (nrow(m_setup) + z - 1);x[4] = (nrow(m_setup) + z + 1)}

            x[1] = (1 + z);x[2] = (-1 + z);
            for (u in x) {if (m_setup[u] == 3) {sim@occupyM[u] <- 0;
            #print("Grids whoose occupyM is zero: ");#print(u); 
		}}
            #print("Inner mucus destruction in the matrix point(s): "); #print(z);
            #print("Inner mucus substance(s) at outer mucus levels: "); #print(y)
            #print("---------------------------------------------------------------")
		}}
      }
      #print("---------------------------------------------------------------")
      #print("---------------------------------------------------------------")
      #print("Iteration: ")
      #print(t)

      sim <- BacArena::simEnv(sim, time = 1, continue = T, sec_obj = "mtf", with_shadow = F, verbose = F)
      }
  # Return results
  simlist[[p]] <- sim
  rm(list = ls()[which((ls() %in% savemain) == F)])    
}
saveRDS(simlist, file = paste("OutputData/SIHUMIx/",simulation_sequence$Sample[j],
                                ".RDS",sep = ""), compress = T)
warnings()
