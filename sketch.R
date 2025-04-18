# Alberto Rovellini
# 04/16/2025
# This is a collection of tools that help calibrate Atlantis GOA and evaluate its output
# This is important for inclusion of this model as part of the ensemble work, to explore HCRs, etc.
# Some criteria explored here are taken or adapted from Kaplan et al. (2016) and include:

# General red face tests from biom.txt and catch.txt files
# Persistence - no species go extinct
# Relative biomass over time
# Relative WAA and NAA over time
# Trends in WAA and NAA

# Hindcast validation 
# In cases where biomass exceeds / is below the initial values, is it still within the fistorical bounds for a stock?
# The difficulty here is that this will require data, e.g. catch and biomass streams or other model output
# We could use Stock SMART for some stocks (though units are a headache, e.g. total vs Sb vs female SSB vs other)
# Raw biomass compared to historical values
# Raw catch compared to historcical values

# you have a lot of this code so do not reinvent it especially for the NC files
library(here)
library(tidyverse)
library(tidync)
library(ncdf4)

run_no <- 1961

outpath <- "../Parametrization/output_files/data/"
outdir <- paste0(outpath, "out_", run_no)

# read in files
biomfile <- paste0(outdir, "/outputGOA0", run_no, "_testBiomIndx.txt")

# read in functional groups file
grp <- read.csv("data/GOA_Groups.csv", header = T)

# age structured vs biomass pools
verts <- grp %>% filter(NumCohorts > 1) %>% pull(Code)
pools <- grp %>% filter(NumCohorts == 1) %>% pull(Code)

ncfile <- paste0(outdir, "/outputGOA0", run_no, "_test.nc")

fl <- 'data/GOA_WGS84_V4_final.bgm'
bgm <- rbgm::read_bgm(fl)
goa_sf <- rbgm::box_sf(bgm)
boundary_boxes <- goa_sf %>% sf::st_set_geometry(NULL) %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes


# Relative biomass -------------------------------------------------------

plot_relbiom <- function(biomfile){
  
  biom <- read.delim(biomfile, header = T, sep = " ")
  
  # relative biomass and catch calculations
  relbiom <- biom %>% 
    select(Time, RelKWT:RelDR) %>% 
    pivot_longer(-Time, names_to = "Code") %>%
    mutate(Code = gsub("Rel","",Code))
  
  
  # persistence at the end of the run
  thres <- 0.01 # what do we mean by extinct? <1%?
  relbiom %>%
    filter(Time == max(Time),
           value < 0.01)
  
  # decide on a burn-in, important for the relative plots I think
  burnin <- 30
  
  # relative biomass over time plots
  plotls <- list(verts = verts, pools = pools)
  for (i in 1:2){
    
    grp_name <- names(plotls)[i]  # Get the name as a string
    grp_plot <- plotls[[i]]
    
    rel_p <- relbiom %>%
      #filter(Time >= burnin * 365) %>%
      filter(Code %in% grp_plot) %>%
      ggplot(aes(x = Time / 365, y = value))+
      geom_line()+
      geom_hline(yintercept = 1, color = "blue", linetype = "dashed")+
      geom_hline(yintercept = 0.5, color = "green", linetype = "dashed")+
      geom_hline(yintercept = 2, color = "green", linetype = "dashed")+
      geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed")+
      geom_hline(yintercept = 4, color = "orange", linetype = "dashed")+
      geom_hline(yintercept = 0.125, color = "red", linetype = "dashed")+
      geom_hline(yintercept = 8, color = "red", linetype = "dashed")+
      scale_y_continuous(limits = c(0,NA))+
      theme_bw()+
      facet_wrap(~Code, scales = "free", ncol = 4)
    
    ggsave(paste0("plots/relbiom_", grp_name,".png"), rel_p, width = 7, height = 14)
  }
  
}

plot_relbiom(biomfile)


# Relative numbers at age -------------------------------------------------

# utility functions
# function sum over depth layers in each array slice
collapse_array <- function(mat){
  mat2 <- apply(mat, 3, colSums)
  mat3 <- data.frame(t(mat2))
  colnames(mat3) <- 0:108
  mat3
}


# function to set values in the boundary boxes to NA
setNA <- function(mat) {
  mat2 <- mat
  if(length(dim(mat2))==3) mat2[,(boundary_boxes+1),]<-NA
  if(length(dim(mat2))==2) mat2[(boundary_boxes+1),] <- NA
  mat2
}

# get names, drop migratory species as the relative numbers are not going to work here
# TODO: figure out a way to handle migratory spp, e.g. averages
vert_nomig <- grp %>% 
  filter(Code %in% verts,
         NumMigrations == 0) %>% 
  pull(Name) 

extract_naa <- function(ncfile){
  
  this_ncfile <- tidync(ncfile)
  this_ncdata <- nc_open(ncfile)
  
  ts <- ncdf4::ncvar_get(this_ncdata,varid = "t") %>% as.numeric
  tyrs <- ts/(60*60*24*365)
  
  # do one fg at a time, then bring them back together
  naa_frame <- data.frame()
  for (i in 1:length(vert_nomig)){
    
    fg <- vert_nomig[i]
    
    # Get numbers by box
    abun_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums",name)) %>% # filter for abundance variables
      filter(grepl(fg,name)) # filter for specific functional group
    
    abun1 <- purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this_ncdata) %>% 
      lapply(setNA) %>%
      purrr::map(apply,MARGIN=3,FUN=sum,na.rm=T) %>% 
      bind_cols() %>% 
      suppressMessages() %>% 
      set_names(abun_vars$name) %>% 
      mutate(t=tyrs)
    
    abun2 <- abun1 %>%
      pivot_longer(cols = -t,names_to = 'age_group',values_to = 'abun') %>%
      mutate(age=parse_number(age_group)) %>%
      mutate(year = ceiling(t)) %>%
      group_by(year, age_group, age) %>%
      summarise(abun = mean(abun)) %>%
      ungroup() %>%
      mutate(Name = vert_nomig[i]) %>%
      dplyr::select(year, Name, age, abun)
    
    # get relative values from t0
    abun3 <- abun2 %>%
      group_by(Name, age) %>%
      # Get the base value (year 0) for each group
      mutate(base_abun = first(abun[year == 0])) %>%
      # Calculate the normalized abundance
      mutate(rel_abun = abun / base_abun) %>%
      # Optionally remove the temporary base_abun column
      select(-base_abun) %>%
      ungroup()
    
    naa_frame <- rbind(naa_frame, abun3) %>%
      arrange(year,Name,age)
    
  }
  
  return(naa_frame)
  
}

# apply function to the nc files
naa <- extract_naa(ncfile)
# add column for age as factor
naa$`Age class` <- factor(naa$age)

# plot
naa_plot <- naa %>%
  ggplot(aes(x = year, y = rel_abun, color = `Age class`))+
  geom_line()+
  geom_hline(yintercept = 1, color = "blue", linetype = "dashed")+
  geom_hline(yintercept = 0.5, color = "green", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "green", linetype = "dashed")+
  geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed")+
  geom_hline(yintercept = 4, color = "orange", linetype = "dashed")+
  geom_hline(yintercept = 0.125, color = "red", linetype = "dashed")+
  geom_hline(yintercept = 8, color = "red", linetype = "dashed")+
  scale_color_grey()+
  theme_bw()+
  labs(x = "Year", y = 'Relative abundance')+
  facet_wrap(~Name, scales = 'free', ncol = 4)
naa_plot

# save
ggsave(paste0("plots/relNAA.png"), naa_plot, width = 7, height = 14)

# linear trends over end of run, tiles: how is it changing?
# This will plot the annual percentage change of the relative abundance
trend_yrs <- 10
trends <- naa %>%
  filter(year > (max(year)-trend_yrs)) %>%
  group_by(Name, age, `Age class`) %>%
  nest() %>%
  mutate(model = map(data, ~lm(rel_abun ~ year, data = .)),
         slope = map_dbl(model, ~coef(.)[2])) %>%
  select(Name, age, `Age class`, slope) %>%
  mutate(slope = slope * 100) %>% # turn to annual percentage change
  ungroup()

# Bubble plot
trend_bubbles <- trends %>%
  ggplot()+
  geom_point(aes(x = `Age class`, y = Name, size = abs(slope), 
                 color = slope > 0),
             shape = 1)+
  scale_color_manual(values = c("blue", "red"), labels = c("Negative", "Positive"))+
  scale_size_continuous(range = c(1, 20))+
  theme_bw()+
  labs(title = "Annual % change in relative abundance over last 10 years")
trend_bubbles
# large values in either direction make everything else difficult to assess but I think that is a good feature - encourages you to address the biggest issues first
ggsave(paste0("plots/slope_naa.png"), trend_bubbles, width = 7, height = 7)


# Relative weight at age --------------------------------------------------
# WAA - do the same set of plots:
# relative WAA over time
# trends in last 10 years
# function that extracts weight at age from NepCDF files
extract_waa <- function(ncfile){
  
  this_ncfile <- tidync(ncfile)
  this_ncdata <- nc_open(ncfile)
  
  ts <- ncdf4::ncvar_get(this_ncdata,varid = "t") %>% as.numeric
  tyrs <- ts/(60*60*24*365)
  
  # do one fg at a time, then bring them back together
  waa_frame <- data.frame()
  for (i in 1:length(vert_nomig)){
    
    fg <- vert_nomig[i]
    print(paste("doing", fg))
    
    # get the attributes associated with each functional group
    fg_atts <- grp %>% filter(Name==fg)
    
    #Extract from the output .nc file the appropriate reserve N time series variables
    resN_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_ResN",name)) %>% # filter for reserve N
      filter(grepl(fg,name)) # filter for specific functional group
    
    #Extract from the output .nc file the appropriate structural N time series variables
    strucN_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_StructN",name)) %>% # filter for structural N
      filter(grepl(fg,name)) # filter for specific functional group
    
    # Get numbers by box
    abun_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums",name)) %>% # filter for abundance variables
      filter(grepl(fg,name)) # filter for specific functional group
    
    # A note on handling boundary boxes and sediments
    # Previous uses of this code have weighted WAA by NAA in the boundary boxes, which is 0
    # This worked for the purpose of excluding BB from mean WAA calcs, but this is technically inaccurate
    # same applies to sediments
    # We should use NAs instead
    # The good news is that the outcome does not change (I tried both ways), but I think the NA method makes more sense
    
    if(nrow(resN_vars)==0) {return("no data.")}
    else {
      # # Actually pull the data from the .nc
      resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this_ncdata) 
      strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this_ncdata)
      nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this_ncdata) %>% lapply(setNA) #numbers by age group,box,layer,time
      totnums <-nums %>% purrr::map(apply,MARGIN=3,FUN=function (x) sum(x, na.rm = T)) # total numbers by age group, time
      relnums <- purrr::map2(nums,totnums,sweep,MARGIN=3,FUN=`/`)
     
       # add the two matrices to get total nitrogen weight
      rnsn <- purrr::map2(resN,strucN,`+`)
      
      # WAA in the sediments remains at initial values throughout the run
      # normally this is not an issue because we use NAA to weight the average and the NAA value in the sed is E-16
      # It does become a problem when the group is extinct though, as then this non-0 value brings WAA up!
      rnsn <- purrr::map(rnsn, function(x) {
        x[7, , ] <- NA
        return(x)
      })
      
      rnsn_summ <- purrr::map2(rnsn,relnums,`*`) %>% 
        purrr::map(apply,MARGIN=3,FUN=function(x) sum(x, na.rm = T)) %>% # mean total N by time
        bind_cols() %>% # bind age groups elements together
        suppressMessages() %>% 
        set_names(resN_vars$name) %>% 
        mutate(t=tyrs) %>%
        # pivot to long form
        pivot_longer(cols = contains(fg_atts$Name),names_to = 'age_group',values_to = 'totN') %>%
        mutate(age=parse_number(age_group)) %>% 
        mutate(weight=totN*20*5.7/1000000) %>%   # convert totN to weight/individual in kg
        #dplyr::filter(t>0) %>%
        mutate(year = ceiling(t)) %>%
        group_by(year, age_group, age) %>%
        summarise(weight = mean(weight)) %>%
        ungroup() %>%
        mutate(Name = fg_atts$Name) %>%
        select(-age_group)
      
      relwgt <- rnsn_summ %>%
        group_by(Name, age) %>%
        # Get the base value (year 0) for each group
        mutate(base_weight = first(weight[year == 0])) %>%
        # Calculate the normalized abundance
        mutate(rel_weight = weight / base_weight) %>%
        # Do this as relative percentage change, positive or negative
        # rowwise() %>%
        # prop = ifelse(abun / init >= 1, 
        #               (abun - init) / init * 100, 
        #               -1*((init - abun) / abun * 100)) %>%
        select(-base_weight) %>%
        ungroup()
      
      waa_frame <- rbind(waa_frame, relwgt) %>%
        arrange(year,Name,age)
      
    }
    
  }
  
  return(waa_frame)
}

# apply function to the nc files
waa <- extract_waa(ncfile)
# add column for age as factor
waa$`Age class` <- factor(waa$age)

# plot
waa_plot <- waa %>%
  ggplot(aes(x = year, y = rel_weight, color = `Age class`))+
  geom_line()+
  geom_hline(yintercept = 1, color = "blue", linetype = "dashed")+
  geom_hline(yintercept = 0.5, color = "green", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "green", linetype = "dashed")+
  geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed")+
  geom_hline(yintercept = 4, color = "orange", linetype = "dashed")+
  # geom_hline(yintercept = 0.125, color = "red", linetype = "dashed")+
  # geom_hline(yintercept = 8, color = "red", linetype = "dashed")+
  scale_color_grey()+
  #scale_y_continuous(limits = c(-8,8))+ # if it's worse than that it does not matter how bad it is
  theme_bw()+
  labs(title = "Weight at age", x = "Year", y = 'kg')+
  facet_wrap(~Name, scales = 'free', ncol = 4)
waa_plot

# save
ggsave(paste0("plots/relWAA.png"), waa_plot, width = 7, height = 14)

# linear trends over end of run, tiles: how is it changing?
# This will plot the annual percentage change of the relative abundance
trend_yrs <- 10
trends <- waa %>%
  filter(year > (max(year)-trend_yrs)) %>%
  group_by(Name, age, `Age class`) %>%
  nest() %>%
  mutate(model = map(data, ~lm(rel_weight ~ year, data = .)),
         slope = map_dbl(model, ~coef(.)[2])) %>%
  select(Name, age, `Age class`, slope) %>%
  mutate(slope = slope * 100) %>% # turn to annual percentage change
  ungroup()

# Option 2: Bubble plot
trend_bubbles <- trends %>%
  ggplot()+
  geom_point(aes(x = `Age class`, y = Name, size = abs(slope), 
                 color = slope > 0),
             shape = 1)+
  scale_color_manual(values = c("blue", "red"), labels = c("Negative", "Positive"))+
  scale_size_continuous(range = c(1, 20))+
  theme_bw()+
  labs(title = "Annual % change in relative weigth at age over last 10 years")
trend_bubbles
# large values in either direction make everything else difficult to assess but I think that is a good feature - encourages you to address the biggest issues first
ggsave(paste0("plots/slope_waa.png"), trend_bubbles, width = 7, height = 7)
