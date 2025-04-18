#' Atlantis Model Evaluation Functions
#' 
#' This script contains functions for evaluating Atlantis GOA model outputs
#' and producing standardized validation plots.

#' Plot relative biomass over time
#'
#' Reads in biomass index file and creates plots of relative biomass over time
#' for both age-structured and biomass pool species.
#'
#' @param biomfile Path to the biomass index file
#' @param grp_data Dataframe containing functional group information
#' @param verts Vector of vertebrate species codes
#' @param pools Vector of biomass pool species codes
#' @param run_no Run number for including in output file names
#' @param outdir Directory to save plots to
#'
#' @return Invisibly returns the processed biomass data
plot_relbiom <- function(biomfile, grp_data, verts, pools, run_no, outdir) {
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  biom <- read.delim(biomfile, header = T, sep = " ")
  
  # relative biomass calculations
  relbiom <- biom %>% 
    select(Time, RelKWT:RelDR) %>% 
    pivot_longer(-Time, names_to = "Code") %>%
    mutate(Code = gsub("Rel", "", Code)) %>%
    left_join(grp_data %>% select(Code,Name), by = "Code")
  
  # check for extinction at the end of the run
  thres <- 0.01 # what do we mean by extinct? <1%?
  extinct_groups <- relbiom %>%
    filter(Time == max(Time),
           value < thres) %>%
    pull(Name)
  
  if(length(extinct_groups) > 0) {
    warning(paste("The following groups went extinct (biomass < 1% of initial):", 
                  paste(extinct_groups, collapse = ", ")))
  }
  
  # relative biomass over time plots
  plotls <- list(verts = verts, pools = pools)
  for (i in 1:2) {
    
    grp_name <- names(plotls)[i]  
    grp_plot <- plotls[[i]]
    
    rel_p <- relbiom %>%
      filter(Code %in% grp_plot) %>%
      ggplot(aes(x = Time / 365, y = value)) +
      geom_line() +
      geom_hline(yintercept = 1, color = "blue", linetype = "dashed") +
      geom_hline(yintercept = 0.5, color = "green", linetype = "dashed") +
      geom_hline(yintercept = 2, color = "green", linetype = "dashed") +
      geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed") +
      geom_hline(yintercept = 4, color = "orange", linetype = "dashed") +
      geom_hline(yintercept = 0.125, color = "red", linetype = "dashed") +
      geom_hline(yintercept = 8, color = "red", linetype = "dashed") +
      scale_y_continuous(limits = c(0, NA)) +
      theme_bw() +
      labs(x = "Year", y = "Relative biomass") +
      facet_wrap(~Name, scales = "free", ncol = 4)
    
    ggsave(paste0(outdir, "/relbiom_", grp_name, "_", run_no, ".png"), rel_p, width = 7, height = 14)
  }
  
  # Return the relbiom data invisibly
  invisible(relbiom)
}

#' Analyze biomass trends
#'
#' Analyzes trends in biomass over the last n years of the simulation
#' and creates bubble plots to visualize annual percentage changes.
#'
#' @param biomfile Path to the biomass index file
#' @param grp_data Dataframe containing functional group information
#' @param trend_yrs Number of years for trend analysis
#' @param run_no Run number for including in output file names
#' @param outdir Directory to save plots to
#'
#' @return Invisibly returns the calculated trends
analyze_biomass_trends <- function(biomfile, grp_data, trend_yrs = 10, run_no, outdir) {
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  biom <- read.delim(biomfile, header = T, sep = " ")
  
  # Get relative biomass
  relbiom <- biom %>% 
    select(Time, RelKWT:RelDR) %>% 
    pivot_longer(-Time, names_to = "Code") %>%
    mutate(Code = gsub("Rel", "", Code)) %>%
    # Convert time to years
    mutate(Year = Time / 365)
  
  # Get vertebrates and pools
  verts <- grp_data %>% filter(NumCohorts > 1) %>% pull(Code)
  pools <- grp_data %>% filter(NumCohorts == 1) %>% pull(Code)
  
  # Add group type
  relbiom <- relbiom %>%
    mutate(GroupType = case_when(
      Code %in% verts ~ "Vertebrates",
      Code %in% pools ~ "Biomass Pools",
      TRUE ~ "Other"
    ))
  
  # Calculate trends over the end of the simulation
  max_year <- max(relbiom$Year)
  cutoff_year <- max_year - trend_yrs
  
  trends <- relbiom %>%
    filter(Year > cutoff_year) %>%
    group_by(Code, GroupType) %>%
    nest() %>%
    mutate(
      model = map(data, ~lm(value ~ Year, data = .)),
      slope = map_dbl(model, ~coef(.)[2]),
      # Convert to annual percentage change
      pct_change = slope * 100 / mean(map_dbl(data, ~mean(.$value))),
      # Get final value for determining size of effects
      final_value = map_dbl(data, ~.$value[which.max(.$Year)])
    ) %>%
    select(Code, GroupType, pct_change, final_value) %>%
    ungroup()
  
  # Add group names for better plot labels
  trends <- trends %>%
    left_join(grp_data %>% select(Code, Name), by = "Code")
  
  # Create a sorted bar plot of all groups
  all_groups_plot <- ggplot(trends) +
    geom_col(aes(x = reorder(Code, pct_change), y = pct_change, fill = pct_change > 0)) +
    scale_fill_manual(values = c("blue", "red"), 
                      labels = c("Decreasing", "Increasing"),
                      name = "Direction") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    facet_wrap(~GroupType, scales = "free_y") +
    labs(
      title = paste0("Annual biomass change (%) over last ", trend_yrs, " years"),
      x = "Functional Group",
      y = "Annual Change (%)"
    ) +
    theme_bw()
  
  ggsave(
    paste0(outdir, "/biomass_trends_all_", run_no, ".png"), 
    all_groups_plot, 
    width = 12, 
    height = 10
  )
  
  # Return the trends data invisibly
  invisible(trends)
}

#' Set values in boundary boxes to NA
#'
#' Helper function to set values in boundary boxes to NA
#'
#' @param mat Matrix or array to modify
#' @param boundary_boxes Vector of boundary box IDs
#'
#' @return Modified matrix or array with NA values in boundary boxes
setNA <- function(mat, boundary_boxes) {
  mat2 <- mat
  if(length(dim(mat2)) == 3) mat2[, (boundary_boxes + 1), ] <- NA
  if(length(dim(mat2)) == 2) mat2[(boundary_boxes + 1), ] <- NA
  mat2
}

#' Collapse array over depth layers
#'
#' Helper function to sum over depth layers in each array slice
#'
#' @param mat 3D array to collapse
#'
#' @return 2D dataframe with summed values
collapse_array <- function(mat) {
  mat2 <- apply(mat, 3, colSums)
  mat3 <- data.frame(t(mat2))
  colnames(mat3) <- 0:108
  mat3
}

#' Extract and plot numbers at age
#'
#' Extracts numbers at age from NetCDF files and creates plots of:
#' 1. Relative numbers at age over time
#' 2. Trends in numbers at age over the last n years
#'
#' @param ncfile Path to the NetCDF file
#' @param grp_data Dataframe containing functional group information
#' @param boundary_boxes Vector of boundary box IDs
#' @param run_no Run number for including in output file names
#' @param outdir Directory to save plots to
#' @param trend_yrs Number of years for calculating trends (default 10)
#'
#' @return Invisibly returns a list with the processed data
extract_plot_naa <- function(ncfile, grp_data, boundary_boxes, run_no, outdir, trend_yrs = 10) {
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  this_ncfile <- tidync(ncfile)
  this_ncdata <- nc_open(ncfile)
  
  ts <- ncdf4::ncvar_get(this_ncdata, varid = "t") %>% as.numeric
  tyrs <- ts / (60 * 60 * 24 * 365)
  
  # get names, drop migratory species as the relative numbers are not going to work here
  vert_nomig <- grp_data %>% 
    filter(NumCohorts > 1,
           NumMigrations == 0) %>% 
    pull(Name)
  
  # Extract numbers at age
  naa_frame <- data.frame()
  for (i in 1:length(vert_nomig)) {
    
    fg <- vert_nomig[i]
    
    # Get numbers by box
    abun_vars <- hyper_vars(this_ncfile) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums", name)) %>% # filter for abundance variables
      filter(grepl(fg, name)) # filter for specific functional group
    
    if (length(abun_vars$name) == 0) {
      warning(paste("No abundance data found for", fg))
      next
    }
    
    abun1 <- purrr::map(abun_vars$name, ncdf4::ncvar_get, nc = this_ncdata) %>% 
      lapply(function(x) setNA(x, boundary_boxes)) %>%
      purrr::map(apply, MARGIN = 3, FUN = sum, na.rm = T) %>% 
      bind_cols() %>% 
      suppressMessages() %>% 
      set_names(abun_vars$name) %>% 
      mutate(t = tyrs)
    
    abun2 <- abun1 %>%
      pivot_longer(cols = -t, names_to = 'age_group', values_to = 'abun') %>%
      mutate(age = parse_number(age_group)) %>%
      mutate(year = ceiling(t)) %>%
      group_by(year, age_group, age) %>%
      summarise(abun = mean(abun), .groups = "drop") %>%
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
      arrange(year, Name, age)
  }
  
  if (nrow(naa_frame) == 0) {
    warning("No data found for numbers at age")
    return(NULL)
  }
  
  # add column for age as factor
  naa_frame$`Age class` <- factor(naa_frame$age)
  
  # plot relative numbers at age over time
  naa_plot <- naa_frame %>%
    ggplot(aes(x = year, y = rel_abun, color = `Age class`)) +
    geom_line() +
    geom_hline(yintercept = 1, color = "blue", linetype = "dashed") +
    geom_hline(yintercept = 0.5, color = "green", linetype = "dashed") +
    geom_hline(yintercept = 2, color = "green", linetype = "dashed") +
    geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = 4, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = 0.125, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 8, color = "red", linetype = "dashed") +
    scale_color_grey() +
    theme_bw() +
    labs(x = "Year", y = 'Relative abundance') +
    facet_wrap(~Name, scales = 'free', ncol = 4)
  
  # save plot
  ggsave(paste0(outdir, "/relNAA_", run_no, ".png"), naa_plot, width = 7, height = 14)
  
  # Calculate trends over the last n years
  trends <- naa_frame %>%
    filter(year > (max(year) - trend_yrs)) %>%
    group_by(Name, age, `Age class`) %>%
    nest() %>%
    mutate(model = map(data, ~lm(rel_abun ~ year, data = .)),
           slope = map_dbl(model, ~coef(.)[2])) %>%
    select(Name, age, `Age class`, slope) %>%
    mutate(slope = slope * 100) %>% # turn to annual percentage change
    ungroup()
  
  # Bubble plot of trends
  trend_bubbles <- trends %>%
    ggplot() +
    geom_point(aes(x = `Age class`, y = Name, size = abs(slope), 
                   color = slope > 0),
               shape = 1) +
    scale_color_manual(values = c("blue", "red"), 
                       labels = c("Negative", "Positive"),
                       name = "Direction") +
    scale_size_continuous(range = c(1, 20), name = "Annual % change") +
    theme_bw() +
    labs(title = paste("Annual % change in relative abundance over last", trend_yrs, "years"))
  
  # save trend plot
  ggsave(paste0(outdir, "/slope_naa_", run_no, ".png"), trend_bubbles, width = 7, height = 7)
  
  # Return the data invisibly
  invisible(list(naa = naa_frame, trends = trends))
}

#' Extract and plot weight at age
#'
#' Extracts weight at age from NetCDF files and creates plots of:
#' 1. Relative weight at age over time
#' 2. Trends in weight at age over the last n years
#'
#' @param ncfile Path to the NetCDF file
#' @param grp_data Dataframe containing functional group information
#' @param boundary_boxes Vector of boundary box IDs
#' @param run_no Run number for including in output file names
#' @param outdir Directory to save plots to
#' @param trend_yrs Number of years for calculating trends (default 10)
#'
#' @return Invisibly returns a list with the processed data
extract_plot_waa <- function(ncfile, grp_data, boundary_boxes, run_no, outdir, trend_yrs = 10) {
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  this_ncfile <- tidync(ncfile)
  this_ncdata <- nc_open(ncfile)
  
  ts <- ncdf4::ncvar_get(this_ncdata, varid = "t") %>% as.numeric
  tyrs <- ts / (60 * 60 * 24 * 365)
  
  # get names, drop migratory species
  vert_nomig <- grp_data %>% 
    filter(NumCohorts > 1,
           NumMigrations == 0) %>% 
    pull(Name)
  
  # Extract weight at age
  waa_frame <- data.frame()
  for (i in 1:length(vert_nomig)) {
    
    fg <- vert_nomig[i]
    message(paste("Processing weight at age for", fg))
    
    # get the attributes associated with each functional group
    fg_atts <- grp_data %>% filter(Name == fg)
    
    # Extract reserve N time series variables
    resN_vars <- hyper_vars(this_ncfile) %>%
      filter(grepl("_ResN", name)) %>%
      filter(grepl(fg, name))
    
    # Extract structural N time series variables
    strucN_vars <- hyper_vars(this_ncfile) %>%
      filter(grepl("_StructN", name)) %>%
      filter(grepl(fg, name))
    
    # Get numbers by box
    abun_vars <- hyper_vars(this_ncfile) %>%
      filter(grepl("_Nums", name)) %>%
      filter(grepl(fg, name))
    
    if (nrow(resN_vars) == 0 || nrow(strucN_vars) == 0 || nrow(abun_vars) == 0) {
      warning(paste("Missing data for", fg))
      next
    }
    
    # Pull the data from the .nc
    resN <- purrr::map(resN_vars$name, ncdf4::ncvar_get, nc = this_ncdata) 
    strucN <- purrr::map(strucN_vars$name, ncdf4::ncvar_get, nc = this_ncdata)
    nums <- purrr::map(abun_vars$name, ncdf4::ncvar_get, nc = this_ncdata) %>% 
      lapply(function(x) setNA(x, boundary_boxes))
    
    totnums <- nums %>% purrr::map(apply, MARGIN = 3, FUN = function(x) sum(x, na.rm = T))
    relnums <- purrr::map2(nums, totnums, sweep, MARGIN = 3, FUN = `/`)
    
    # Add the two matrices to get total nitrogen weight
    rnsn <- purrr::map2(resN, strucN, `+`)
    
    # Set sediments to NA
    rnsn <- purrr::map(rnsn, function(x) {
      x[7, , ] <- NA
      return(x)
    })
    
    rnsn_summ <- purrr::map2(rnsn, relnums, `*`) %>% 
      purrr::map(apply, MARGIN = 3, FUN = function(x) sum(x, na.rm = T)) %>%
      bind_cols() %>%
      suppressMessages() %>% 
      set_names(resN_vars$name) %>% 
      mutate(t = tyrs) %>%
      # pivot to long form
      pivot_longer(cols = contains(fg_atts$Name), names_to = 'age_group', values_to = 'totN') %>%
      mutate(age = parse_number(age_group)) %>% 
      mutate(weight = totN * 20 * 5.7 / 1000000) %>%   # convert totN to weight/individual in kg
      mutate(year = ceiling(t)) %>%
      group_by(year, age_group, age) %>%
      summarise(weight = mean(weight), .groups = "drop") %>%
      mutate(Name = fg_atts$Name) %>%
      select(-age_group)
    
    relwgt <- rnsn_summ %>%
      group_by(Name, age) %>%
      # Get the base value (year 0) for each group
      mutate(base_weight = first(weight[year == 0])) %>%
      # Calculate the normalized weight
      mutate(rel_weight = weight / base_weight) %>%
      select(-base_weight) %>%
      ungroup()
    
    waa_frame <- rbind(waa_frame, relwgt) %>%
      arrange(year, Name, age)
  }
  
  if (nrow(waa_frame) == 0) {
    warning("No data found for weight at age")
    return(NULL)
  }
  
  # Add column for age as factor
  waa_frame$`Age class` <- factor(waa_frame$age)
  
  # Plot relative weight at age over time
  waa_plot <- waa_frame %>%
    ggplot(aes(x = year, y = rel_weight, color = `Age class`)) +
    geom_line() +
    geom_hline(yintercept = 1, color = "blue", linetype = "dashed") +
    geom_hline(yintercept = 0.5, color = "green", linetype = "dashed") +
    geom_hline(yintercept = 2, color = "green", linetype = "dashed") +
    geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed") +
    geom_hline(yintercept = 4, color = "orange", linetype = "dashed") +
    scale_color_grey() +
    theme_bw() +
    labs(title = "Relative weight at age", x = "Year", y = 'Relative weight') +
    facet_wrap(~Name, scales = 'free', ncol = 4)
  
  # Save plot
  ggsave(paste0(outdir, "/relWAA_", run_no, ".png"), waa_plot, width = 7, height = 14)
  
  # Calculate trends over the last n years
  trends <- waa_frame %>%
    filter(year > (max(year) - trend_yrs)) %>%
    group_by(Name, age, `Age class`) %>%
    nest() %>%
    mutate(model = map(data, ~lm(rel_weight ~ year, data = .)),
           slope = map_dbl(model, ~coef(.)[2])) %>%
    select(Name, age, `Age class`, slope) %>%
    mutate(slope = slope * 100) %>% # turn to annual percentage change
    ungroup()
  
  # Bubble plot of trends
  trend_bubbles <- trends %>%
    ggplot() +
    geom_point(aes(x = `Age class`, y = Name, size = abs(slope), 
                   color = slope > 0),
               shape = 1) +
    scale_color_manual(values = c("blue", "red"), 
                       labels = c("Negative", "Positive"),
                       name = "Direction") +
    scale_size_continuous(range = c(1, 20), name = "Annual % change") +
    theme_bw() +
    labs(title = paste("Annual % change in relative weight at age over last", trend_yrs, "years"))
  
  # Save trend plot
  ggsave(paste0(outdir, "/slope_waa_", run_no, ".png"), trend_bubbles, width = 7, height = 7)
  
  # Return the data invisibly
  invisible(list(waa = waa_frame, trends = trends))
}

#' Generate a summary report of model evaluation
#'
#' Creates a text summary of model evaluation results
#'
#' @param biomfile Path to the biomass index file
#' @param grp_data Dataframe containing functional group information
#' @param trend_yrs Number of years for trend analysis
#' @param run_no Run number for including in output file name
#' @param outdir Directory to save report to
#'
#' @return Invisibly returns the report content
generate_summary_report <- function(biomfile, grp_data, trend_yrs = 10, run_no, outdir) {
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # Read biomass data
  biom <- read.delim(biomfile, header = T, sep = " ")
  
  # Get relative biomass
  relbiom <- biom %>% 
    select(Time, RelKWT:RelDR) %>% 
    pivot_longer(-Time, names_to = "Code") %>%
    mutate(Code = gsub("Rel", "", Code)) %>%
    # Convert time to years
    mutate(Year = Time / 365)
  
  # Check for extinction at the end of the run
  thres <- 0.01 # what do we mean by extinct? <1%?
  extinct_groups <- relbiom %>%
    filter(Time == max(Time),
           value < thres) %>%
    pull(Code)
  
  # Get groups with extreme changes
  extreme_high <- relbiom %>%
    filter(Time == max(Time),
           value > 8) %>%
    pull(Code)
  
  extreme_low <- relbiom %>%
    filter(Time == max(Time),
           value < 0.125, 
           value >= thres) %>%
    pull(Code)
  
  # Calculate trends over the end of the simulation
  max_year <- max(relbiom$Year)
  cutoff_year <- max_year - trend_yrs
  
  trends <- relbiom %>%
    filter(Year > cutoff_year) %>%
    group_by(Code) %>%
    nest() %>%
    mutate(
      model = map(data, ~lm(value ~ Year, data = .)),
      slope = map_dbl(model, ~coef(.)[2]),
      # Convert to annual percentage change
      pct_change = slope * 100 / map_dbl(data, ~mean(.$value)),
      # Get final value
      final_value = map_dbl(data, ~.$value[which.max(.$Year)])
    ) %>%
    select(Code, pct_change, final_value) %>%
    ungroup()
  
  # Identify groups with concerning trends
  strong_trends <- trends %>%
    filter(abs(pct_change) > 5) %>%
    arrange(desc(abs(pct_change)))
  
  # Build the report
  report <- c(
    paste("# Atlantis Model Evaluation Summary - Run", run_no),
    paste("Date:", Sys.Date()),
    paste("Simulation length:", max_year, "years"),
    "",
    "## Persistence Check",
    if(length(extinct_groups) > 0) {
      c(
        "WARNING: The following groups went extinct (biomass < 1% of initial):",
        paste("-", extinct_groups, collapse = "\n")
      )
    } else {
      "All groups persisted throughout the simulation."
    },
    "",
    "## Extreme Biomass Changes",
    if(length(extreme_high) > 0) {
      c(
        "Groups with extreme increases (> 8x initial biomass):",
        paste("-", extreme_high, collapse = "\n")
      )
    } else {
      "No groups showed extreme increases in biomass."
    },
    "",
    if(length(extreme_low) > 0) {
      c(
        "Groups with extreme decreases (< 12.5% of initial biomass, but not extinct):",
        paste("-", extreme_low, collapse = "\n")
      )
    } else {
      "No groups showed extreme decreases in biomass."
    },
    "",
    paste("## Strong Trends Over Last", trend_yrs, "Years"),
    if(nrow(strong_trends) > 0) {
      c(
        "The following groups show strong trends (>5% annual change):",
        paste("-", strong_trends$Code, ":", 
              sprintf("%+.1f", strong_trends$pct_change), "% per year,", 
              sprintf("%.2f", strong_trends$final_value), "x initial biomass",
              collapse = "\n")
      )
    } else {
      "No groups showed strong trends in the final years of the simulation."
    }
  )
  
  # Write the report to a file
  report_file <- paste0(outdir, "/model_summary_", run_no, ".md")
  writeLines(report, report_file)
  message(paste("Summary report written to", report_file))
  
  # Return the report content invisibly
  invisible(report)
}

#' Master function to evaluate Atlantis model output
#'
#' This function runs selected evaluation functions on the specified model run
#'
#' @param run_no Run number to evaluate
#' @param trend_yrs Number of years to use for trend analysis (default: 10)
#' @param plot_biomass Whether to plot relative biomass (default: TRUE)
#' @param plot_naa Whether to plot numbers at age (default: TRUE)
#' @param plot_waa Whether to plot weight at age (default: TRUE)
#' @param plot_biomass_trends Whether to plot biomass trends (default: FALSE)
#' @param generate_report Whether to generate a summary report (default: TRUE)
#' @param output_base_path Base path to the output files (default: "../Parametrization/output_files/data/")
#' @param data_path Path to the model data files (default: "data/")
#' @param bgm_file Path to the BGM file (default: "data/GOA_WGS84_V4_final.bgm")
#' @param plot_output_dir Base directory for saving plots (default: "plots/")
#'
#' @return NULL (saves plots and report as side effects)
#' 
#' @examples
#' # Run all analyses
#' evaluate_atlantis_run(1961)
#' 
#' # Only plot biomass with 15-year trend analysis
#' evaluate_atlantis_run(1961, trend_yrs = 15, plot_naa = FALSE, plot_waa = FALSE)
#' 
#' # Only plot numbers at age
#' evaluate_atlantis_run(1961, plot_biomass = FALSE, plot_waa = FALSE)
evaluate_atlantis_run <- function(run_no, 
                                  trend_yrs = 10,
                                  plot_biomass = TRUE,
                                  plot_naa = TRUE,
                                  plot_waa = TRUE,
                                  plot_biomass_trends = FALSE,
                                  generate_report = TRUE,
                                  output_base_path = "../Parametrization/output_files/data/",
                                  data_path = "data/",
                                  bgm_file = "data/GOA_WGS84_V4_final.bgm",
                                  plot_output_dir = "plots/") {
  
  # Create output directory for this run
  run_output_dir <- paste0(plot_output_dir, "run_", run_no, "/")
  if (!dir.exists(run_output_dir)) {
    dir.create(run_output_dir, recursive = TRUE)
  }
  
  message(paste("Evaluating Atlantis run", run_no))
  
  # Paths to input files
  outdir <- paste0(output_base_path, "out_", run_no)
  biomfile <- paste0(outdir, "/outputGOA0", run_no, "_testBiomIndx.txt")
  ncfile <- paste0(outdir, "/outputGOA0", run_no, "_test.nc")
  
  # Check if input files exist
  if (!file.exists(biomfile)) {
    stop("Biomass index file not found: ", biomfile)
  }
  
  if ((plot_naa || plot_waa) && !file.exists(ncfile)) {
    warning("NetCDF file not found: ", ncfile)
    plot_naa <- FALSE
    plot_waa <- FALSE
  }
  
  # Read in functional groups file
  grp_file <- paste0(data_path, "GOA_Groups.csv")
  if (!file.exists(grp_file)) {
    stop("Functional groups file not found: ", grp_file)
  }
  grp <- read.csv(grp_file, header = TRUE)
  
  # Age structured vs biomass pools
  verts <- grp %>% filter(NumCohorts > 1) %>% pull(Code)
  pools <- grp %>% filter(NumCohorts == 1) %>% pull(Code)
  
  # Check if we need boundary boxes info (for NAA or WAA)
  if (plot_naa || plot_waa) {
    # Read BGM file to get boundary boxes
    fl <- bgm_file
    if (!file.exists(fl)) {
      stop("BGM file not found: ", fl)
    }
    
    bgm <- rbgm::read_bgm(fl)
    goa_sf <- rbgm::box_sf(bgm)
    boundary_boxes <- goa_sf %>% 
      sf::st_set_geometry(NULL) %>% 
      filter(boundary == TRUE) %>% 
      pull(box_id)
  } else {
    boundary_boxes <- NULL
  }
  
  # Run the selected analyses
  if (plot_biomass) {
    message("Analyzing relative biomass")
    plot_relbiom(biomfile, grp, verts, pools, run_no, run_output_dir)
  }
  
  if (plot_biomass_trends) {
    message("Analyzing biomass trends")
    analyze_biomass_trends(biomfile, grp, trend_yrs, run_no, run_output_dir)
  }
  
  if (plot_naa) {
    message("Analyzing numbers at age")
    extract_plot_naa(ncfile, grp, boundary_boxes, run_no, run_output_dir, trend_yrs)
  }
  
  if (plot_waa) {
    message("Analyzing weight at age")
    extract_plot_waa(ncfile, grp, boundary_boxes, run_no, run_output_dir, trend_yrs)
  }
  
  if (generate_report) {
    message("Generating summary report")
    generate_summary_report(biomfile, grp, trend_yrs, run_no, run_output_dir)
  }
  
  message(paste("Evaluation complete. Results saved to", run_output_dir))
}