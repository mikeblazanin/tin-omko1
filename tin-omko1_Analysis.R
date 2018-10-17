library(ggplot2)

data1 <- read.csv("Data-Plate Reader-190817.csv", stringsAsFactors = F)
layout1 <- read.csv("Data-Plate Reader-190817_layout.csv", stringsAsFactors = F)
data2 <- read.csv("Data-Plate Reader-150817.csv", stringsAsFactors = F)
layout2 <- read.csv("Data-Plate Reader-150817_layout.csv", stringsAsFactors = F)
data3 <- read.csv("2018-10-12 Growth Curve.csv", stringsAsFactors = F)
layout3 <- read.csv("2018-10-9 Plate Layout.csv", stringsAsFactors = F)

layout3[layout3=="LB"] <- NA

#Clean up layout
layout_cleanup <- function(layout) {
  for (i in 1:nrow(layout)) {
    for (j in 2:ncol(layout)) {
      if(!is.na(layout[i, j])) {
        if (nchar(layout[i, j]) > 1) {
          current_val <- layout[i, j]
          cntr <- 1
        }
        layout[i, j] <- paste(current_val, LETTERS[cntr], sep = "_")
        cntr <- cntr + 1
      }
    }
  }
  return(layout)
}

cln_lay1 <- layout_cleanup(layout1)
cln_lay2 <- layout_cleanup(layout2)
cln_lay3 <- layout_cleanup(layout3)

#Merge layout & data
merge_lay_data <- function(layout, data) {
  layout_mlt <- reshape2::melt(layout, id.vars = 1, value.name = "Contents", 
                               variable.name = "column")
  layout_mlt$Well <- paste(layout_mlt[, 1], substr(layout_mlt$column, 2, 
                                                nchar(as.character(layout_mlt$column))), 
                           sep = "")
  
  data_mlt <- reshape2::melt(data, id.vars = c("Time", "Temperature"), 
                             variable.name = "Well", value.name = "OD600")
  data_mlt$Contents <- layout_mlt$Contents[match(as.character(data_mlt$Well), 
                                                 as.character(layout_mlt$Well))]
  data_mlt$format <- paste(colnames(layout)[1], "_Rep", sep = "")
  return(data_mlt)
}

#Clean up merged datasets
clean_mdata <- function(mdata) {
  mdata <- mdata[complete.cases(mdata), ]
  mdata$Time <- substr(mdata$Time, 1, (nchar(mdata$Time)-1))
  return(mdata)
}

mdata1 <- clean_mdata(merge_lay_data(cln_lay1, data1))
mdata2 <- clean_mdata(merge_lay_data(cln_lay2, data2))
mdata3 <- clean_mdata(merge_lay_data(cln_lay3, data3))

#Split contents up
split_contents <- function(input_frame) {
  colnames <- strsplit(input_frame$format[1], "_")[[1]]
  input_frame[colnames] <- do.call(rbind, strsplit(input_frame$Contents, "_"))
  return(subset(input_frame, select=-c(format, Contents)))
}

spl_data1 <- split_contents(mdata1)
spl_data2 <- split_contents(mdata2)
spl_data3 <- split_contents(mdata3)


data_mlt <- split_contents(data_mlt)

#Re-order factor levels
data_mlt$Time <- as.numeric(data_mlt$Time)
data_mlt$Stress[data_mlt$Stress == "None"] <- "0"
data_mlt$Stress <- factor(data_mlt$Stress,
                          levels = c("0", "5", "90", "180", "270", "360"))

data_mlt$Strain[data_mlt$Strain == "Pure LB Control"] <- "LB"
data_mlt$Strain[data_mlt$Strain == "Pure Shocked LB control"] <- "LB Shock"
data_mlt$Strain[data_mlt$Strain == "LB+PAO1 Control"] <- "PAO1"
data_mlt$Strain[data_mlt$Strain == "Shocked LB+PAO1 Control"] <- "PAO1 Shock"
data_mlt$Strain <- factor(data_mlt$Strain, 
                          levels = c("PAO1", "PAO1 Shock", "S3", "S8", "S11", "S16",
                                     "R3", "LB", "LB Shock"))

#Make new variable for each indiv growth curve
data_mlt$Group <- paste(data_mlt$Strain, data_mlt$Stress,
                        data_mlt$Rep)
# data_mlt$Group <- match(paste(data_mlt$Strain, data_mlt$Stress, data_mlt$Rep),
#                     unique(paste(data_mlt$Strain, data_mlt$Stress, data_mlt$Rep)))

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

#Plots to visually inspect peak designation accuracy
for (start_group in seq(from = 1, to = length(unique(data_mlt$Group)), by = 9)) {
  my_groups <- unique(data_mlt$Group)[start_group:(start_group+8)]
  print(ggplot(data = data_mlt[data_mlt$Group %in% my_groups, ],
               aes(x = Time, y = sm_od)) + geom_line() +
          facet_wrap(~Group) +
          geom_point(data = data_out[data_out$Group %in% my_groups, ],
                     aes(x = maxtime, y = max),
                     size = 3, pch = 13) +
          ylab("Smoothed OD600"))
  ggsave(filename = paste(start_group, "_growcurves.pdf", sep = ""),
         device = "pdf", width = 8, height = 8, units = "in")
}

#Plots to just look at growth curves
# for (i in seq(from = 1, to = 287, by = 20)) {
#   print(ggplot(data = data[data$Group %in% i:(i+19), ], aes(x = Time..s., y = OD600)) +
#     geom_line() + facet_wrap(~Group))
#   print(ggplot(data = data[data$Group %in% i:(i+8), ], aes(x = Time..s., y = sm_od)) +
#           geom_line() + facet_wrap(~Group))
# }

#Plots to look at summarized data
ggplot(data = data_out, aes(x = Strain, y = max, group = Stress, color = Stress)) +
  geom_point(size = 2, position = position_dodge(0.6)) +
  ylab("Max OD600 of Peak")
ggsave(filename = "gc_maxes.pdf", device = "pdf", width = 8, height = 8, units = "in")
ggplot(data = data_out, aes(x = Strain, y = maxtime, group = Stress, color = Stress)) +
  geom_point(size = 2, position = position_dodge(0.6)) +
  ylab("Time of OD600 Peak")
ggsave(filename = "gc_maxtimes.pdf", device = "pdf", width = 8, height = 8, units = "in")