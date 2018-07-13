library(ggplot2)
library(reshape2)
library(dplyr)

data <- read.csv("Data-Plate Reader-190817.csv", stringsAsFactors = F)
layout <- read.csv("Data-Plate Reader-190817_layout.csv", stringsAsFactors = F)

#Clean up layout
for (i in 1:nrow(layout)) {
  for (j in 2:ncol(layout)) {
    if (!(is.na(layout[i, j]) | nchar(layout[i, j]) < 1)) {
      current_val <- layout[i, j]
      cntr <- 1
    }
    layout[i, j] <- paste(current_val, LETTERS[cntr], sep = "_")
    cntr <- cntr + 1
  }
}

#Merge layout & data
layout_mlt <- melt(layout, id.vars = "X", value.name = "Contents")
layout_mlt$Well <- paste(layout_mlt$X, substr(layout_mlt$variable, 2, 
                                              nchar(as.character(layout_mlt$variable))), 
                         sep = "")

data_mlt <- melt(data, id.vars = c("Time", "Temperature"), variable.name = "Well",
                 value.name = "OD600")
data_mlt$Contents <- layout_mlt$Contents[match(as.character(data_mlt$Well), 
                                               as.character(layout_mlt$Well))]

#Clean up merged dataset
data_mlt <- data_mlt[complete.cases(data_mlt), ]
data_mlt$Time <- substr(data_mlt$Time, 1, (nchar(data_mlt$Time)-1))
data_mlt <- data_mlt[, c("Time", "Temperature", "Contents", "OD600")]

#Split contents up
split_contents <- function(input_frame) {
  my_split <- strsplit(input_frame$Contents, "_")
  input_frame$Strain <- NA
  input_frame$Stress <- NA
  input_frame$Rep <- NA
  for (i in 1:nrow(input_frame)) {
    temp <- my_split[[i]][1]
    if (grepl("\\(*\\)", temp)) {
      subsplit <- strsplit(temp, " ")
      input_frame$Strain[i] <- subsplit[[1]][1]
      input_frame$Stress[i] <- subsplit[[1]][2]
    } else {
      input_frame$Strain[i] <- temp
      input_frame$Stress[i] <- "(None)"
    }
    input_frame$Rep[i] <- my_split[[i]][length(my_split[[i]])]
  }
  input_frame$Stress <- substr(input_frame$Stress, 2, nchar(input_frame$Stress)-1)
  return(input_frame)
}

data_mlt <- split_contents(data_mlt)

#Make new variable for each indiv growth curve
data_mlt$Group <- match(paste(data_mlt$Strain, data_mlt$Stress, data_mlt$Rep),
                    unique(paste(data_mlt$Strain, data_mlt$Stress, data_mlt$Rep)))

#smooth OD data
smooth_data <- function(my_data, smooth_over, subset_by) {
  #data must be sorted sequentially before fed into function
  #my_data is a vector of the data to be smoothed
  #smooth over is how many sequential entries to average
  #the unique values of subset_by will be what is iterated over
  out_list <- rep(NA, length(my_data))
  cntr = 1
  for (my_uniq in unique(subset_by)) {
    my_sub <- subset(my_data, subset_by == my_uniq)
    out_list[cntr:(cntr+length(my_sub)-smooth_over)] <- 0
    for (i in 1:smooth_over) {
      out_list[(cntr):(cntr+length(my_sub)-smooth_over)] <-
        out_list[(cntr):(cntr+length(my_sub)-smooth_over)] + 
        my_sub[i:(length(my_sub)-smooth_over+i)]
    }  
    cntr <- cntr+length(my_sub)
  }
  out_list <- out_list/smooth_over
  return(out_list)
}

data_mlt$sm_od <- smooth_data(data_mlt$OD600, smooth_over = 10, subset_by = data_mlt$Group)

#extraction function: to get OD peak height and time
analyze_curves <- function(od_data, time_data, bandwidth = 10, return) {
  prev_max_pos <- 0
  cnt_max_pos <- bandwidth
  while (cnt_max_pos != prev_max_pos) {
    prev_max_pos <- cnt_max_pos
    search_start <- pmax(cnt_max_pos-bandwidth, 1) #start of the search window
    search_end <- pmin(search_start + 2*bandwidth, length(od_data)) #end of the search window
    cnt_max_pos <- which.max(od_data[search_start:search_end]) + search_start - 1
  }
  if (return == "max") {return(od_data[cnt_max_pos])
  } else if (return == "maxtime") {return(time_data[cnt_max_pos])}
}

#Group data by indiv growth curves
data_grp <- group_by(data_mlt[!is.na(data_mlt$sm_od), ], 
                     Group, Strain, Stress, Rep)

#Get OD peak height & time for each growth curve
data_out <- summarize(data_grp, 
                      max = analyze_curves(sm_od, Time, 
                                           bandwidth = 20, return = "max"),
                      maxtime = analyze_curves(sm_od, Time, 
                                               bandwidth = 20, return = "maxtime"))

data_mlt$Time <- as.numeric(data_mlt$Time)
data_out$maxtime <- as.numeric(data_out$maxtime)

# #Plots to visually inspect peak designation accuracy
# for (start_group in seq(from = 1, to = 96, by = 9)) {
#   my_groups <- start_group:(start_group+8)
#   print(ggplot(data = data_mlt[data_mlt$Group %in% my_groups, ],
#          aes(x = Time, y = sm_od)) + geom_line() +
#     geom_point(data = data_out[data_out$Group %in% my_groups, ],
#                aes(x = maxtime, y = max),
#                size = 3, pch = 13) +
#       facet_wrap(~Group))
# }

#Plotting tests
# for (i in seq(from = 1, to = 287, by = 20)) {
#   print(ggplot(data = data[data$Group %in% i:(i+19), ], aes(x = Time..s., y = OD600)) +
#     geom_line() + facet_wrap(~Group))
#   print(ggplot(data = data[data$Group %in% i:(i+8), ], aes(x = Time..s., y = sm_od)) +
#           geom_line() + facet_wrap(~Group))
# }

data_out$Stress <- factor(data_out$Stress,
                             levels = c("None", "0", "5", "90", "180", "270", "360"))
data_out$Strain[data_out$Strain == "Pure LB Control"] <- "LB"
data_out$Strain[data_out$Strain == "Pure Shocked LB Control"] <- "LB Shock"
data_out$Strain[data_out$Strain == "LB+PAO1 Control"] <- "PAO1"
data_out$Strain[data_out$Strain == "Shocked LB+PAO1 Control"] <- "PAO1 Shock"

data_out$Strain <- factor(data_out$Strain, 
                          levels = c("PAO1", "PAO1 Shock", "S3", "S8", "S11", "S16",
                                     "R3", "LB", "LB Shock"))
                                     
ggplot(data = data_out, aes(x = Stress, y = max)) +
  geom_point(size = 2) + facet_grid(~Strain)
