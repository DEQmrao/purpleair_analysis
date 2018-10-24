library(tidyverse)  ## reverse factors, amongst other things
library(lubridate)  ## date manipulation
library(stringr)    ## string manipulation (str_pad)
library(zeallot)    ## %<-% operator for multiple assignments
library(corrplot)   ## correlation matrix plot
library(geosphere)  ## distance calculations based on lat-lon

#library(ggmap)
#library(mapdata)

# lat-lon of reference nephelometers used in the vignette
ref_lat_lon <- data.frame( site = c("SEL", "TBC","BHP", "HHF"),
                           lat = c(45.4966, 45.3997, 45.4702, 45.5285),
                           lon = c(-122.6029, -122.7451, -122.8162, -122.9724))


##--------------------------------------------------------------------------------
## read_sensor_data: function to read in Purple Air data
## inputs: sensor_dir, default = "data/"
##         sensor_names, default = NULL
## 
## The first parameter to the function is the directory with the Purple Air datafiles.
## The function creates 'in_files' - a list of all files in the directory that contain the word "Primary".
## The code assumes that the first 3 letters of the file name are unique for the sensor.
## Alternately, supply a list of unique names to use to identify the sensors in the optional paramter sensor_names.
## Each of the files in in_files is read in, and the time converted from UTC to PST,
##     and only the pm2.5, temp, and rh fields extracted.
## A list of dataframes that can be referenced using the sensor unique names is returned.
## For example sensor[["alo1"]] or sensor[["rcp2"]]
##---------------------------------------------------------------------------------
read_sensor_data <- function(sensor_dir = "data/", sensor_names = NULL) {
  
  in_files <- as.list(list.files(sensor_dir, pattern = "*Primary*"))
  if (is.null(sensor_names)) {
    col_name <- sapply(in_files, function(x) substring(x, 1, 4))
  }
  else colname == sensor_names
  
  last <- length(in_files)
  sensor <- lapply(1:last, function(i) {
    print(paste0(i, " started......", col_name[[i]]))
    idf <- read.csv(paste0("data/", in_files[[i]]),stringsAsFactors = FALSE)
    idf$pst <- ymd_hms(idf$created_at) - hms("8:00:00")
    idf$day <- as.Date(idf$pst)
    idf$hour <- hour(idf$pst)
    idf$minute <- minute(idf$pst)
    idf1 <- idf %>% group_by(day, hour, minute) %>% 
      summarize(pm25 = mean(PM2.5_CF_1_ug.m3, na.rm = TRUE),
                temp = mean(Temperature_F, na.rm = TRUE),
                rh = mean(Humidity_., na.rm = TRUE))
    idf1$pm25[idf1$pm25 > 800] <- NA
    
    idf1$pst <- paste0(idf1$day, " ", 
                       formatC(idf1$hour, width =2, flag = "0"), ":", 
                       formatC(idf1$minute, width = 2, flag = "0"), ":00")
    startt <- min(ymd_hms(idf1$pst), na.rm = TRUE)
    endt <- max(ymd_hms(idf1$pst), na.rm = TRUE)
    ref <- data.frame(pstp = seq(from = startt, to = endt, 
                                 by = '1 min'))
    ref$pst <- as.character(ref$pstp)
    
    idf1 <- merge(ref, idf1[, c("pm25", "temp", "rh", "pst")], by = "pst", all.x = TRUE)
    idf1$site <- col_name[[i]]
    return(idf1)
  }) 
  names(sensor) <- col_name[1:last]
  return(sensor)
}

##---------------------------------------------------------------------------------
## read_ref_monitor_data: function to read in data from reference nephelometers/monitors
## inputs: monitor_file_list (required)
##        monitor_name_list (required)
## List the reference monitor files to be read in in 'monitor_file_list'.
## Specify short names for the reference monitors in 'monitor_name_list'.
## The code expects the pm2.5 column to be labeled 'neph', and expects a'Date' and 'Time' column.
## Returns a list of reference monitor dataframes, indexed by names specified in 'monitor_name_list'.
## For example, df[['SEL']]
##---------------------------------------------------------------------------------

read_ref_monitor_data <- function(monitor_file_list, monitor_name_list) {
  lf <- length(monitor_file_list)
  neph <- lapply(1:lf, function(i) {
    print(paste0(i, " started......", monitor_name_list[[i]]))
    df <- read.csv(monitor_file_list[[i]], stringsAsFactors = FALSE)
    df$pst <- paste0(mdy(df$Date), " ", str_pad(df$Time, 5, pad = "0"), ":00")
    df$neph[df$neph < 0] <- NA
    df$neph[!is.na(df$neph) & df$neph > 100] <- NA
    startt <- min(ymd_hms(df$pst), na.rm = TRUE)
    endt <- max(ymd_hms(df$pst), na.rm = TRUE)
    ref <- data.frame(pstp = seq(from = startt, to = endt, 
                                 by = '1 min'))
    ref$pst <- as.character(ref$pstp)
    df <- merge(ref, df[, c("neph", "pst")], by = "pst", all.x = TRUE)
    return(df)
  })
  names(neph) <- monitor_name_list
  return(neph)
}


##---------------------------------------------------------------------------------
## create_1hour_df: function that takes in a list of dataframes with pm25/neph data
##                  and returns a dataframe with 1-hour data
## inputs: sensor_list (required)
##         time_column, default = "pst
## Takes a list of sensor dataframes (as returned from the read functions) and returns 
## a list of dataframes with all numeric columns averaged to the hour.
##---------------------------------------------------------------------------------
create_1hour_df <- function(sensor_list, time_column = "pst") {
  df_1hour <- lapply(sensor_list, function(df) {
    df$day <- as.Date(df[, time_column])
    df$hour <- hour(df[,time_column])
    df1h <- df %>% group_by(day, hour) %>% 
      summarize_if(is.numeric, mean, na.rm = TRUE) %>%
      mutate_all(.funs = funs(ifelse(is.nan(.), NA, .))) %>% 
      ungroup() %>% as.data.frame()
  })
  return(df_1hour)
}

##--------------------------------------------------------------------------------
## summary stats and graphs 
##--------------------------------------------------------------------------------
create_time_sequence <- function(sensor_list, time_column, time_interval = '1 min') {
  startt <- as.POSIXct(min(sapply(sensor_list, function(df) min(df[, time_column], na.rm = TRUE)), 
                           na.rm = TRUE), origin = "1970-01-01", tz = "UTC")
  endt <- as.POSIXct(max(sapply(sensor_list, function(df) max(df[, time_column], na.rm = TRUE)), 
                         na.rm = TRUE), origin = "1970-01-01", tz = "UTC")
  ref <- data.frame(pstp = seq(from = startt, to = endt, by = time_interval))
  ref$pst <- as.character(ref$pstp)
  colnames(ref) <- c(time_column, paste0(time_column, "_char"))
  return(ref)
}


##---------------------------------------------------------------------------------
## stats_summ: function that takes in a list of dataframes with pm25/neph data
##             and returns a single dataframe with a statistical summary of the data.
##             Summary includes count, first date, last date, min, avg, median, max   
## inputs: sensor_list (required)
##         stat_column, default = "pm25"
##---------------------------------------------------------------------------------

stats_summ <- function(sensor_list, stat_column = "pm25") {
  stat <- data.frame(
    count = sapply(sensor_list, function(df) { sum(!(is.na(df[, stat_column])))}),
    first_date = sapply(sensor_list, function(df) { as.character(min(df$day, na.rm = TRUE))}),
    last_date = sapply(sensor_list, function(df) { as.character(max(df$day, na.rm = TRUE))}),
    min = sapply(sensor_list, function(df) { round(min(df[, stat_column], na.rm = TRUE),2)}),
    avg = sapply(sensor_list, function(df) { round(mean(df[, stat_column], na.rm = TRUE),2)}),
    median = sapply(sensor_list, function(df) { round(median(df[, stat_column], na.rm = TRUE),2)}),
    max = sapply(sensor_list, function(df) { round(max(df[, stat_column], na.rm = TRUE),2)}),
    std_dev = sapply(sensor_list, function(df) { round(sd(df[, stat_column], na.rm = TRUE),2)})
  )
  return(stat)
}


##---------------------------------------------------------------------------------
## corr_summ: function that takes in a list of Purple Air dataframes and a list of
##            refernce dataframes, and returns a single dataframe with the 
##            correlation of each Purple Air with each reference monitor.
## inputs: sensor_list (required)
##         ref_monitor_list (required)
##---------------------------------------------------------------------------------

corr_summ <- function(sensor_list, ref_monitor_list) {
  corr_df <- data.frame(site = names(sensor_list))
  num_sen <- length(sensor_list)
  num_ref <- length(ref_monitor_list)
  for (r in 1:num_ref) {
    coln <- paste0("corr", names(ref_monitor_list)[r])
    corr_df[, coln] <- NA
    for(s in 1:num_sen) {
      df <- merge(sensor_list[[s]], ref_monitor_list[[r]], by = c("day", "hour"))
      corr_df[s, coln] <- round(cor(df$pm25, df$neph, use = "complete.obs"),3)
    }
  }
  return(corr_df)
}

##---------------------------------------------------------------------------------
## r2_summ: function that takes in a list of Purple Air dataframes and a list of
##            refernce dataframes, and returns a single dataframe with the 
##            Adjusted R2  of each Purple Air with each reference monitor.
##            Two R2 values are returned: one for the goodness of fit with the reference
##            monitor, and one for the goodness-of-fit including temperature and 
##            relative humidity.
## inputs: sensor_list (required)
##         ref_monitor_list (required)
##---------------------------------------------------------------------------------
r2_summ <- function(sensor_list, ref_monitor_list) {
  r2_df <- data.frame(site = names(sensor_list))
  num_sen <- length(sensor_list)
  num_ref <- length(ref_monitor_list)
  for (r in 1:num_ref) {
    coln_r2 <- paste0("r2", names(ref_monitor_list)[r])
    coln_r2th <- paste0("r2ll", names(ref_monitor_list)[r])
    r2_df[, coln_r2] <- NA
    r2_df[, coln_r2th] <- NA
    for(s in 1:num_sen) {
      df <- merge(sensor_list[[s]], ref_monitor_list[[r]], by = c("day", "hour"))
      res1 <- lm(data = df, pm25 ~ neph)
      res2 <- lm(data = df, pm25 ~ neph + temp + rh)
      r2_df[s, coln_r2] <- summary(res1)$adj.r.squared
      r2_df[s, coln_r2th] <- summary(res2)$adj.r.squared
    }
  }
  return(r2_df)
}


##---------------------------------------------------------------------------------
## create_distance_matrix: function that takes in a dataframe with the names and lat and lon
##                         for the reference monitors and a similar dataframe with the names
##                         and lat-lon of the Purple Airs used in the analysis.
##                         Returns a dataframe with distance (in km) of each sensor (rows) from
##                         each reference monitor (columns).
## inputs: ref_lat_long, default = NULL
##         sensor_lat_long, default = NULL
##---------------------------------------------------------------------------------
create_distance_matrix <- function(ref_lat_long = NULL, sensor_lat_long = NULL) {
  if (is.null(ref_lat_long)) { ref_lat_long <- ref_lat_lon }
  if (is.null(sensor_lat_long)) {
    in_files <- as.list(list.files("data/", pattern = "*Primary*"))
    col_name <- sapply(in_files, function(x) substring(x, 1, 4))
    sensor_lat_long <- data.frame( site = col_name,
                                   lat = sapply(in_files, function(f) {
                                     st_lat <- regexpr("\\(", f) + 1
                                     rem <- substring(f, st_lat)
                                     en_lat <- regexpr("-", rem) + st_lat - 2
                                     st_lon <- en_lat 
                                     en_lon <- regexpr("\\)", f) - 1
                                     lat <- trimws(substr(f, st_lat, en_lat))}),
                                   lon = sapply(in_files, function(f) {
                                     st_lat <- regexpr("\\(", f) + 1
                                     rem <- substring(f, st_lat)
                                     en_lat <- regexpr("-", rem) + st_lat - 2
                                     st_lon <- en_lat 
                                     en_lon <- regexpr("\\)", f) - 1
                                     lon <- trimws(substr(f, st_lon, en_lon)) 
                                   }), stringsAsFactors = FALSE)
  }
  sensor_lat_long$lat <- as.numeric(sensor_lat_long$lat)
  sensor_lat_long$lon <- as.numeric(sensor_lat_long$lon)
  num_sens <- nrow(sensor_lat_long)
  num_ref <- nrow(ref_lat_long)
  for (n in 1:num_ref) {
    coln <- paste0("dist", ref_lat_long$site[n])
    sensor_lat_long[, coln] <- NA
    for (i in 1:num_sens) {
      sensor_lat_long[i, coln ] <- distm(c(ref_lat_long$lon[n], ref_lat_long$lat[n]), 
                                         c(sensor_lat_long$lon[i], sensor_lat_long$lat[i]), 
                                         fun = distHaversine)
      
    }
    sensor_lat_long[,coln] <- round(sensor_lat_long[,coln]/1000, 1)
  }
  return(sensor_lat_long)
  
}
