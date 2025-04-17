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
biom <- read.delim(paste0(outdir, "/outputGOA0", run_no, "_testBiomIndx.txt"), header = T, sep = " ")
catch <- read.delim(paste0(outdir, "/outputGOA0", run_no, "_testCatch.txt"), header = T, sep = " ")

# read in functional groups file
grp <- read.csv("data/GOA_Groups.csv", header = T)

# age structured vs biomass pools
verts <- grp %>% filter(NumCohorts > 1) %>% pull(Code)
pools <- grp %>% filter(NumCohorts == 1) %>% pull(Code)

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

# WAA and NAA
ncfile <- paste0(outdir, "/outputGOA0", run_no, "_test.nc")

# function sum over depth layers in each array slice
collapse_array <- function(mat){
  mat2 <- apply(mat, 3, colSums)
  mat3 <- data.frame(t(mat2))
  colnames(mat3) <- 0:108
  mat3
}

fl <- 'data/GOA_WGS84_V4_final.bgm'
bgm <- rbgm::read_bgm(fl)
goa_sf <- rbgm::box_sf(bgm)
boundary_boxes <- goa_sf %>% sf::st_set_geometry(NULL) %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes
# function to set values in the boundary boxes to NA
setNA <- function(mat) {
  mat2 <- mat
  if(length(dim(mat2))==3) mat2[,(boundary_boxes+1),]<-NA
  if(length(dim(mat2))==2) mat2[(boundary_boxes+1),] <- NA
  mat2
}

# get names, drop migratory species as the relative numbers are not going to work here
# TODO: figure out a way to handle migratory spp, e.g. averages
vert_names <- grp %>% 
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
  for (i in 1:length(vert_names)){
    
    fg <- vert_names[i]
    
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
      mutate(Name = vert_names[i]) %>%
      dplyr::select(year, Name, age, abun)
    
    # get relative values from t0
    abun3 <- abun2 %>%
      group_by(Name, age) %>%
      # Get the base value (year 0) for each group
      mutate(base_abun = first(abun[year == 0])) %>%
      # Calculate the normalized abundance
      mutate(norm_abun = abun / base_abun) %>%
      # Optionally remove the temporary base_abun column
      select(-base_abun) %>%
      ungroup()
    
    naa_frame <- rbind(naa_frame, abun3) %>%
      arrange(year,age)
    
  }
  
  return(naa_frame)
  
}

# apply function to the nc files
naa <- extract_naa(ncfile)
# add column for age as factor
naa$`Age class` <- factor(naa$age)

# plot
naa_plot <- naa %>%
  ggplot(aes(x = year, y = norm_abun, color = `Age class`))+
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

# make a figure
ggsave(paste0("plots/relNAA.png"), naa_plot, width = 7, height = 14)


# some ideas for alternative figures, let's see what's most useful
# terminal year rel abundance, tiles: how different from the beginning is it?
# linear trends over end of run, tiles: how is it changing?


