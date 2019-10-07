library(ggplot2)
library(dplyr)

#Temperature Survival ----

#Read in & format data
tempr <- read.csv("Temperature-Survival.csv")
#Note temp in Celcius & Duration in minutes
colnames(tempr)[1:2] <- c("Temp", "Duration")

#Calculate percent survival
#(Note) The master stock was titered immediately before use and 200 pfu was 
# the expected Plate.count
tempr$pct_surv <- 100*tempr$Plate.count/mean(tempr$Plate.count[
  tempr$Temp == 55 & tempr$Duration == 5])

#Calculate limit of detection
tempr_limit_detec <- 100*1/mean(tempr$Plate.count[
  tempr$Temp == 55 & tempr$Duration == 5])

#Summarize
tempr <- group_by(tempr, Temp, Duration)
tempr_sum <- summarize(tempr, 
                       mean_pct_surv = mean(pct_surv))

#Flag & adjust values below detection
tempr_sum$bd <- tempr_sum$mean_pct_surv < tempr_limit_detec
tempr_sum$mean_pct_surv[tempr_sum$bd] <- tempr_limit_detec

#Make variables for plotting
tempr_sum$mean_pct_surv_lines <- tempr_sum$mean_pct_surv #this is to plot the lines
tempr_sum$ptshape <- 16

#Make pch 8 for first 0 and NA for subsequent 0's
#Stop plotting lines after the first bd point
for (temp in tempr_sum$Temp) {
  myrows <- tempr_sum$Temp == temp
  bd_rows <- (myrows & tempr_sum$bd)
  tempr_sum$ptshape[bd_rows] <- 8
  if (sum(bd_rows) > 1) {
    #identify bd rows where the previous point was also bd
    for (i in 2:sum(bd_rows)) {
      this_row <- which(bd_rows)[i]
      if (tempr_sum$bd[this_row - 1]) {
        tempr_sum$ptshape[this_row] <- NA
        tempr_sum$mean_pct_surv_lines[this_row] <- NA
      }
    }
  }
}
tempr_sum$ptshape <- as.factor(tempr_sum$ptshape)
tempr_sum$Temp <- as.factor(tempr_sum$Temp)

my_colr <- colorRampPalette(colors = c("#ffcc00", "red"))
tempr_plot <- ggplot(data = tempr_sum, aes(x = Duration, y = mean_pct_surv,
                             group = Temp, color = Temp, shape = ptshape)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_line(aes(x = Duration, y = mean_pct_surv_lines)) + 
  scale_color_manual(name = "Temperature (°C)", values = my_colr(6)) +
  scale_shape_manual(values = list("8" = 8, "16" = 16)) +
  scale_x_continuous(breaks = seq(from = 0, to = 90, by = 30),
                     limits = c(0, 90)) +
  scale_y_continuous(breaks = c(100, 10, 1),
                     labels = c("100", "10", "1"),
                     trans="log10") +
  labs(x = "Heat Shock Duration (min)",
       y = "Percent Survival (%)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_hline(yintercept = tempr_limit_detec, lty = 3, lwd = 1.15) +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15) + 
  guides(shape = FALSE) +
  NULL
tempr_plot
ggsave(filename = "temp_surv.tiff", width = 8, height = 5, units = "in")

#Statistics
#Take subset of data that has 3+ measures above limit of detection
#in treatment
tempr_sum_stats <- tempr_sum[as.numeric(as.character(tempr_sum$Temp)) <= 70, ]
tempr_sum_stats$Temp <- as.factor(tempr_sum_stats$Temp)
tempr_sum_stats$Duration <- as.numeric(tempr_sum_stats$Duration)

#Run with 55C as reference (so we can test for sig dift slope than 0)
tempr_model_55 <- lm(log10(mean_pct_surv)~Temp + Duration:Temp,
                  data = tempr_sum_stats)
anova(tempr_model_55)
summary(tempr_model_55)

#Run with 60C as reference
tempr_sum_stats$Temp <- relevel(tempr_sum_stats$Temp, "60")
tempr_model_60 <- lm(log10(mean_pct_surv)~Temp + Duration:Temp,
                 data = tempr_sum_stats)
summary(tempr_model_60)

#Run with 65C as reference
tempr_sum_stats$Temp <- relevel(tempr_sum_stats$Temp, "65")
tempr_model_65 <- lm(log10(mean_pct_surv)~Temp + Duration:Temp,
                     data = tempr_sum_stats)
summary(tempr_model_65)

#Run with 70C as reference
tempr_sum_stats$Temp <- relevel(tempr_sum_stats$Temp, "70")
tempr_model_70 <- lm(log10(mean_pct_surv)~Temp + Duration:Temp,
                     data = tempr_sum_stats)
summary(tempr_model_70)

#Results          Estimate  Std. Error  t value Pr(>|t|)    
#Temp55:Duration  0.0008093  0.0008223   0.984   0.3538  
#Temp60:Duration  0.0004729  0.0008223   0.575   0.5810    
#Temp65:Duration -0.0027540  0.0008223  -3.349   0.0101 *  
#Temp70:Duration -0.0181560  0.0008223 -22.080 1.87e-08 ***


#Temperature Duration Survival ----
temprdur <- read.csv("Temperature-Duration-Survival.csv")
#Note that all shocks were carried out at 70C, drop that column
temprdur <- subset(temprdur, select = -c(Temperature.of.shock...70))
#Redo titer calculation for safety
temprdur$Titer <- temprdur$Plate.count * 10**(temprdur$Dilution + 1)
#Rename Phage Stocks
mynames <- c("A" = "R3", "B" = "S3", "C" = "S8", "D" = "S11", "E" = "S16")
temprdur$Sample <- names(mynames)[match(temprdur$Sample, mynames)]
#Calculate pct survival relative to mean of 0 shock duration
temprdur <- group_by(temprdur, Sample, Duration.of.shock..m.)
temprdur_sum <- summarize(temprdur, titer_mean = mean(Titer),
                          titer_se = sd(Titer)/n())
temprdur$pct_surv <- 100*temprdur$Titer/temprdur_sum$titer_mean[match(temprdur$Sample,
                                                                  temprdur_sum$Sample)]
#Calculate limits of detection by assuming 1 pfu on -1 plate (so titer = 100)
#(only plot the highest limit of detection)
temprdur_limit_detec <- 100*100/temprdur_sum$titer_mean[
  temprdur_sum$Duration.of.shock..m. == 0]

temprdur <- group_by(temprdur, Sample, Duration.of.shock..m.)
temprdur_sum <- summarize(temprdur, titer_mean = mean(Titer),
                          pct_mean = mean(pct_surv))

#Make plot
temprdur_plot <- ggplot(data = temprdur_sum, aes(x = Duration.of.shock..m.,
                                y = pct_mean,
                                group = Sample,
                                color = Sample)) +
  geom_point(size = 2.5, alpha = 0.7) + geom_line() +
  scale_y_continuous(breaks = c(100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001),
                     labels = c("100", "10", "1", "0.1", "0.01", "0.001",
                                "0.0001", "0.00001"),
                     trans="log10") +
  scale_color_discrete(name = "Phage Stock") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "70°C Heat Shock Duration (min)", y = "Percent Survival (%)") +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15) +
  geom_hline(yintercept = max(temprdur_limit_detec), lty = 3, lwd = 1.15)
temprdur_plot
ggsave(filename = "temp_duration_surv.tiff", width = 8, height = 5, units = "in")

#Statistics
colnames(temprdur_sum)[2] <- "Duration"
temprdur_model <- lm(log10(pct_mean)~Sample*Duration, data = temprdur_sum)
anova(temprdur_model)


#Growth curve analysis ----

data1 <- read.csv("Data-Plate Reader-190817.csv", stringsAsFactors = F)
layout1 <- read.csv("Data-Plate Reader-190817_layout.csv", stringsAsFactors = F)
data2 <- read.csv("Data-Plate Reader-150817.csv", stringsAsFactors = F)
layout2 <- read.csv("Data-Plate Reader-150817_layout.csv", stringsAsFactors = F)
data3 <- read.csv("2018-10-12 Growth Curve.csv", stringsAsFactors = F)
layout3 <- read.csv("2018-10-9 Plate Layout.csv", stringsAsFactors = F)

layout3[layout3=="LB"] <- NA

#Clean up layout
layout_cleanup <- function(layout) {
  #This function takes a layout dataframe with ...'s where info needs to be 
  #iterated in and does so, while adding _A, _B etc for replicate wells
  #with the same contents
  #Wells labeled NA are not filled in (they are empty)
  
  #Define a matrix to track well contents we've seen before
  #So that when they come up again the replicate number can resume
  #where it left off
  vals_and_cntr <- matrix(nrow = 0, ncol = 2)
  for (i in 1:nrow(layout)) {
    for (j in 2:ncol(layout)) {
      if(!is.na(layout[i, j])) { #Non-empty well
        if (nchar(layout[i, j]) > 1) { #Non ... well
          current_val <- layout[i, j]
          if (!(current_val %in% vals_and_cntr[, 1])) {
            #This is the first time we've seen these well contents
            #So we should start numbering replicates at the beginning
            vals_and_cntr <- rbind(vals_and_cntr, c(current_val, 1))
            row <- nrow(vals_and_cntr)
          } else { 
            #this isn't the first time we've seen these contents
            #resume replicate numbering where we left off
            row <- which(vals_and_cntr[, 1] == current_val)
          }
        }
        layout[i, j] <- paste(current_val, 
                              LETTERS[as.numeric(vals_and_cntr[row, 2])], sep = "_")
        vals_and_cntr[row, 2] <- as.numeric(vals_and_cntr[row, 2]) + 1
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
  #Tidies (melts) the layout dataframe
  #Then merges the now-tidy layout & data dataframes
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
  mdata$Time <- as.numeric(substr(mdata$Time, 1, (nchar(mdata$Time)-1)))
  return(mdata)
}

mdata1 <- clean_mdata(merge_lay_data(cln_lay1, data1))
mdata2 <- clean_mdata(merge_lay_data(cln_lay2, data2))
mdata3 <- clean_mdata(merge_lay_data(cln_lay3, data3))

#Split contents up
split_contents <- function(input_frame) {
  #Splits contents at every underscore
  #Uses the column name of the first column as the format for new columns
  # where the split contents are distributed into
  colnames <- strsplit(input_frame$format[1], "_")[[1]]
  input_frame[colnames] <- do.call(rbind, strsplit(input_frame$Contents, "_"))
  return(subset(input_frame, select=-c(format)))
}

spl_data1 <- split_contents(mdata1)
spl_data2 <- split_contents(mdata2)
spl_data3 <- split_contents(mdata3)

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

spl_data1$sm_od <- smooth_data(spl_data1$OD600, smooth_over = 10,
                               subset_by = spl_data1$Contents)
spl_data2$sm_od <- smooth_data(spl_data2$OD600, smooth_over = 10,
                               subset_by = spl_data2$Contents)
spl_data3$sm_od <- smooth_data(spl_data3$OD600, smooth_over = 10,
                               subset_by = spl_data3$Contents)

#extraction function: to get OD peak height and time
analyze_curves <- function(od_data, time_data, bandwidth = 10, return) {
  #Takes vectors of the od_data and time_data
  #Bandwidth is how wide the window should be to look for a peak
  #Narrower bandwidth will get you an earlier local maxima
  #Wider bandwidth will get you a later more-global maxima
  #Designed to be run w/ group_by
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
grp_data1 <- dplyr::group_by(spl_data1[!is.na(spl_data1$sm_od), ], 
                     bacteria, phage, phageshock, bactshock, mediashock, Rep, Contents)
grp_data2 <- dplyr::group_by(spl_data2[!is.na(spl_data2$sm_od), ],
                             pctmediashocked, bacteria, phage, phageshock, Rep, Contents)
grp_data3 <- dplyr::group_by(spl_data3[!is.na(spl_data3$sm_od), ],
                             bact, phage, totalpfuinoc, Rep, Contents)
                             
#Get OD peak height & time for each growth curve
out_data1 <- dplyr::summarize(grp_data1, 
                              max = analyze_curves(sm_od, Time, bandwidth = 20, 
                                                   return = "max"),
                              maxtime = analyze_curves(sm_od, Time, bandwidth = 20, 
                                                       return = "maxtime"))
out_data2 <- dplyr::summarize(grp_data2, 
                              max = analyze_curves(sm_od, Time, 
                                                   bandwidth = 20, return = "max"),
                              maxtime = analyze_curves(sm_od, Time, 
                                                       bandwidth = 20, return = "maxtime"))
out_data3 <- dplyr::summarize(grp_data3, 
                              max = analyze_curves(sm_od, Time, 
                                                   bandwidth = 20, return = "max"),
                              maxtime = analyze_curves(sm_od, Time, 
                                                       bandwidth = 20, return = "maxtime"))

#Add information about quantity of pfu inoculated to out_data1
out_data1$pfu_inoc <- 200
out_data1$pfu_inoc[(out_data1$phage == "S3" | out_data1$phage == "S8") &
                     out_data1$phageshock == 270] <- 100
out_data1$pfu_inoc[(out_data1$phage == "S11" | out_data1$phage == "S16") &
                     out_data1$phageshock == 270] <- 50
out_data1$pfu_inoc[out_data1$phage == "R3" & out_data1$phageshock == 360] <- 50
out_data1$pfu_inoc[out_data1$phage == "NA"] <- 0

#Growth Curve Figures ----

#Plots to visually inspect peak designation accuracy
view_peaks <- function(data_mlt, data_out, plt_point = TRUE,
                       numplots = 9, lwd = 1) {
  for (start_group in seq(from = 1, to = length(unique(data_mlt$Contents)), by = numplots)) {
    my_groups <- unique(data_mlt$Contents)[start_group:(start_group+numplots-1)]
    myplot <- ggplot(data = data_mlt[data_mlt$Contents %in% my_groups, ],
                 aes(x = Time, y = sm_od)) + geom_line(lwd = lwd) +
            facet_wrap(~Contents) + ylab("Smoothed OD600")
    if (plt_point) {myplot <- myplot + 
                              geom_point(data = data_out[data_out$Contents %in% my_groups, ],
                                                  aes(x = maxtime, y = max),
                                                  size = 3, pch = 13)
    }
    print(myplot)
    # ggsave(filename = paste(start_group, "_growcurves.pdf", sep = ""),
    #        device = "pdf", width = 8, height = 8, units = "in")
  }
}

#view_peaks(grp_data1, out_data1)
#view_peaks(grp_data2, out_data2)
#view_peaks(grp_data3, out_data3)
#view_peaks(grp_data3, out_data3, plt_point = F, numplots = 1, lwd = 2)

#Plots to look at summarized data

#GC Plot 1 ####
#Create "plot" column for x axis
plot1 <- out_data1
plot1$plot <- NA
for (i in 1:nrow(plot1)) {
  if (!plot1$phage[i] == "NA") {
    plot1$plot[i] <- plot1$phageshock[i]
  } else {
    if (plot1$bacteria[i]=="NA") {plot1$plot[i] <- "- Ctrl"
    } else {
      if (plot1$mediashock[i] == 1) {
        plot1$plot[i] <- "+ Ctrl Shock"
      } else {plot1$plot[i] <- "+ Ctrl"}
    }
  }
}

#Reformat columns
plot1$phageshock <- factor(plot1$phageshock,
                           levels = c("NA", 0, 5, 90, 180, 270, 360))
plot1$plot <- factor(plot1$plot,
                        levels = c("+ Ctrl", "+ Ctrl Shock",
                                   "0", "5", "90", "180", "270",
                                   "360", "- Ctrl"))
plot1 <- plot1[-which(plot1$max < 0.1 & plot1$plot != "- Ctrl"), ]
#Rename Phage Stocks
mynames <- c("A" = "R3", "B" = "S3", "C" = "S8", "D" = "S11", "E" = "S16")
plot1$phage <- names(mynames)[match(plot1$phage, mynames)]

#Make plot
# #with all data
# ggplot(data = plot1, aes(x = plot, y = max, group = phage, color = phage)) +
#   geom_point(size = 2, position = position_dodge(0.6)) +
#   ylab("Max Bacterial Density (OD600)") + xlab("Duration of Heatshock (min)") +
#   scale_colour_discrete(name="Phage Replicate") +
#   theme_bw() + scale_y_continuous(limits = c(0, NA))
# #ggsave(filename = "gc_maxes.tiff", width = 8, height = 5, units = "in")

#With summarized data
plot1$pfu_inoc <- as.numeric(as.character(plot1$pfu_inoc))
plot1 <- group_by(plot1, plot, phage)
plot1_sum <- summarize(plot1,
                       pfu_inoc = mean(pfu_inoc),
                       maxpeak_mean = mean(max),
                       maxpeak_se = sd(max)/n())

plot1_sum$pfu_inoc <- as.factor(plot1_sum$pfu_inoc)                       
gc_plot1 <- ggplot(data = plot1_sum[plot1_sum$plot != "- Ctrl" &
                                      plot1_sum$pfu_inoc %in% c(0, 200),], 
                   aes(x = plot, y = maxpeak_mean, color = phage)) +
  geom_point(position = position_dodge(width = .4), size = 3) +
  geom_errorbar(position = position_dodge(width = .4),
                aes(ymax = maxpeak_mean + 1.96*maxpeak_se,
                    ymin = maxpeak_mean - 1.96*maxpeak_se),
                width = .7) +
  theme_bw()  +
  labs(x = "70°C Heat Shock Duration (min)", y = "Peak Bacterial Density (OD600)") +
  scale_color_discrete(name = "Phage Stock") +
  geom_hline(yintercept = plot1_sum$maxpeak_mean[plot1_sum$plot == "- Ctrl"],
             lty = 3, lwd = 1.15) +
  scale_x_discrete(labels = c("+ Ctrl", "+ Ctrl\nShock", "0", "5", 
                              "90", "180", "270")) +
  theme(panel.grid = element_blank())
  #  theme(axis.text.x = element_text(angle = 30, size = 12, hjust = 1))
gc_plot1
ggsave(filename = "gc_maxes.tiff", width = 8, height = 5, units = "in")

#Plot all the data
plot1_sum_temp <- plot1_sum
plot1_sum_temp$pfu_inoc <- factor(plot1_sum_temp$pfu_inoc,
                                  levels = c("0, 0", levels(plot1_sum_temp$pfu_inoc)))
plot1_sum_temp$pfu_inoc[plot1_sum_temp$plot == "- Ctrl"] <- "0, 0"
gc_plot1_alldata <- ggplot(data = plot1_sum_temp, 
                   aes(x = plot, y = maxpeak_mean, color = phage,
                       shape = pfu_inoc)) +
  geom_point(position = position_dodge(width = .4), size = 3) +
  geom_errorbar(position = position_dodge(width = .4),
                aes(ymax = maxpeak_mean + 1.96*maxpeak_se,
                    ymin = maxpeak_mean - 1.96*maxpeak_se),
                width = .7) +
  theme_bw()  +
  theme(panel.grid = element_blank()) +
  labs(x = "Heat Shock Duration (min)", y = "Peak Bacterial Density (OD600)") +
  scale_color_discrete(name = "Phage Stock") +
  scale_shape_manual(name = "Phage; Bacteria\nInoculated (pfu; cfu)",
                     values = c(8, 18, 15, 17, 16),
                     labels = expression("0;           0",
                                         paste("0; ", 20%.%10^6),
                                         paste("50;   ", 5%.%10^6),
                                         paste("100; ", 10%.%10^6),
                                         paste("200; ", 20%.%10^6))) +
  scale_x_discrete(labels = c("+ Ctrl", "+ Ctrl\nShock", "0", "5", 
                              "90", "180", "270", "360", "- Ctrl")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))
gc_plot1_alldata
ggsave(filename = "gc_maxes_alldata.tiff", width = 8, height = 5, units = "in")


#Statistics
gc1_stats <- plot1[plot1$plot %in%
                     c("+ Ctrl", "+ Ctrl Shock", 0, 5, 90, 180,270, 360) &
                     plot1$pfu_inoc %in% c(0, 200),]
gc1_stats$pfu_inoc <- relevel(as.factor(gc1_stats$pfu_inoc), ref = "200")
gc1_stats$plot <- relevel(as.factor(gc1_stats$plot), ref = "0")
gc1_stats$phage[is.na(gc1_stats$phage)] <- "None"
curve1_model <- lm(max~plot + phage + phage:plot,
                   data = gc1_stats)
anova(curve1_model)

#Results
##significant effect of heat treatment:
#plot        6 9.6585 1.60976 1478.546 < 2.2e-16 ***
##Significant variation between replicate stocks:
#phage       4 0.0598 0.01494   13.726 1.403e-07 ***
#plot:phage 11 0.3348 0.03043   27.952 < 2.2e-16 ***


#Contrasts of interest:
# +shock - 0
# 0-5
# 5-90
# 90-180
# 180-270

#Look for contrasts (with 0):
summary(curve1_model)

#Contrasts (with 90):
gc1_stats$plot <- relevel(as.factor(gc1_stats$plot), ref = "90")
curve1_model <- lm(max~plot + phage + phage:plot,
                   data = gc1_stats)
summary(curve1_model)

#Contrasts (with 180):
gc1_stats$plot <- relevel(as.factor(gc1_stats$plot), ref = "180")
curve1_model <- lm(max~plot + phage + phage:plot,
                   data = gc1_stats)
summary(curve1_model)

#Contrasts of interest
#(plot0) - plot+ Ctrl Shock            1.160120   0.023332  49.723  < 2e-16 ***
#(plot0) - plot5                       0.555637   0.026941  20.624  < 2e-16 ***
#(plot90) - plot5                      0.22287    0.02694   8.272 7.34e-11 ***
#(plot90) - plot180                    0.01434    0.02694   0.532 0.596859    
#(plot180) - plot270                   0.110953   0.026941   4.118 0.000146 ***



# curve1_aovmodel <- aov(max~plot + phage + phage:plot,
#                        data = plot1[plot1$plot %in%
#                                       c("+ Ctrl", "+ Ctrl Shock", 0, 5, 90, 180,
#                                         270, 360) &
#                                       plot1$pfu_inoc %in% c(0, 200),])
# summary(curve1_aovmodel) #to check that it's the same numbers
# TukeyHSD(curve1_aovmodel, "plot")

#Combined Temperature Figure ####

#Add space above & below temprdur plot
tempr_plot <- tempr_plot + theme(plot.margin = unit(c(0.1, 0.05, 0.15, 0.05), "in"))
temprdur_plot <- temprdur_plot + theme(plot.margin = unit(c(0.15, 0.05, 0.2, 0.05), "in"))
gc_plot1 <- gc_plot1 + theme(plot.margin = unit(c(0.2, 0.05, 0.2, 0.05), "in"))

comb_temp_plot <- cowplot::plot_grid(tempr_plot, temprdur_plot, gc_plot1,
                   ncol = 1, labels = c("A", "B", "C"),
                   align = "v", label_y = c(1, 1, 1.1))

cowplot::save_plot("combined_temp.tiff", comb_temp_plot,
                   base_height = 8, base_width = 5)

#GC Plot 2 ####
#Reformat & adjust column values
plot2 <- out_data2
plot2$plot <- paste(out_data2$bacteria,
                    out_data2$phage, out_data2$phageshock)
for (i in 1:nrow(plot2)) {
  if (plot2$bacteria[i] == "NA" & plot2$phage[i] == "NA") {
    plot2$plot[i] <- "- Ctrl"
  } else if (plot2$bacteria[i] != "NA" & plot2$phage[i] == "NA") {
    plot2$plot[i] <- "Bact Only"
  } else if (plot2$bacteria[i] == "NA" & plot2$phage[i] != "NA") {
    plot2$plot[i] <- "Phage Only"
  } else {plot2$plot[i] <- "Bact & Phage"}
}
plot2$phageshock <- factor(plot2$phageshock,
                           levels = c("NA", "0", "5", "360"))
plot2$pctmediashocked <- factor(plot2$pctmediashocked,
                                levels = c("0", "5", "50", "100"))
plot2$plot <- factor(plot2$plot,
                     levels = c("Bact & Phage", "Bact Only",
                                "Phage Only", "- Ctrl"))

ggplot(data = plot2, aes(x = plot, y = max, group = pctmediashocked,
                             color = pctmediashocked, shape = phageshock)) +
         geom_point(size = 2, position = position_dodge(0.6)) +
  xlab("")

#GC Plot 3 ####
out_data3$totalpfuinoc[out_data3$totalpfuinoc == "NA"] <- 0
out_data3$totalpfuinoc <- as.numeric(out_data3$totalpfuinoc)

ggplot(data = out_data3, aes(x = totalpfuinoc, y = max)) +
  geom_point(size = 2, pch = 21) + ylim(0, 1.6) +
  labs(x = "Total PFU Inoculated", y = "Max Bacterial Density (OD600)") +
  theme_bw()
ggsave(filename = "gc_maxes_varypfu.tiff", width = 8, height = 5, units = "in")

# ggplot(data = out_data3[out_data3$totalpfuinoc > 0, ], aes(x = totalpfuinoc, y = maxtime)) +
#   geom_jitter(height = 0, width = 5, pch = 21, size = 2)

#Saline ----
saline_data <- read.csv("Saline-Survival.csv", stringsAsFactors = F)
#saline_data$Saline.Concentration..M.[saline_data$Saline.Concentration..M. == "0.17 (LB)"] <- "0.17"
#saline_data$Duration.of.shock..m. <- factor(saline_data$Duration.of.shock..m.)
saline_data$Saline.Concentration..M. <- factor(saline_data$Saline.Concentration..M.)
saline_data$pct_surv <- 100*saline_data$Plate.Count/mean(
  saline_data$Plate.Count[saline_data$Saline.Concentration..M. == "0.17 (LB)" &
                            saline_data$Duration.of.shock..m. == 5])
saline_limit_detection <- 100*1/mean(
  saline_data$Plate.Count[saline_data$Saline.Concentration..M. == "0.17 (LB)" &
                            saline_data$Duration.of.shock..m. == 5])

saline_summary <- dplyr::summarize(group_by(saline_data, Saline.Concentration..M., 
                                     Duration.of.shock..m.),
                                   mean_pct_surv = mean(pct_surv))

# #Plot all data
# ggplot(data = saline_data, aes(x = Saline.Concentration..M., y = Plate.Count,
#                                color = Duration.of.shock..m., 
#                                group = Duration.of.shock..m.)) +
#   geom_point(position = position_dodge(0.4))

# #Plot summarized data, Saline Conc on x
# ggplot(data = saline_summary, aes(x = Saline.Concentration..M., 
#                                   y = mean_pct_surv,
#                                   color = Duration.of.shock..m., 
#                                   group = Duration.of.shock..m.)) +
#   geom_point() + geom_line() +
#   theme_bw() +
#   labs(x = "Saline Concentration (M)", y = "Percent Survival (%)") +
#   scale_color_discrete(name = "Duration of Shock (min)")
# ggsave(filename = "saline_byconc.tiff", width = 8, height = 5, units = "in")

#Plot summarized data, Duration on x
my_cols <- colorRampPalette(c("gray", "Navy"))
ggplot(data = saline_summary, aes(x = Duration.of.shock..m., 
                                  y = mean_pct_surv,
                                  color = Saline.Concentration..M., 
                                  group = Saline.Concentration..M.)) +
  geom_point(size = 3) + geom_line(lwd = 1.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Duration of Shock (min)", 
       y = "Percent Survival (%)") +
  scale_color_manual(name = "Saline\nConcentration (M)",
                       values = my_cols(8)) +
  scale_y_continuous(breaks = c(100, 10, 1),
                     labels = c("100", "10", "1"),
                     trans="log10") +
  scale_x_continuous(limits = c(0, 90),
                     breaks = c(0, 30, 60, 90)) +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15) +
  geom_hline(yintercept = saline_limit_detection, lty = 3, lwd = 1.15)
ggsave(filename = "saline_bydur.tiff", width = 8, height = 5, units = "in")

#Statistics
saline_summary$Duration.of.shock..m. <- as.numeric(as.character(
  saline_summary$Duration.of.shock..m.))
saline_summary$Saline.Concentration..M. <- relevel(
  saline_summary$Saline.Concentration..M., ref = "0.17 (LB)")
saline_model <- lm(log10(mean_pct_surv) ~ Saline.Concentration..M. +
                     Saline.Concentration..M.:Duration.of.shock..m.,
                   data = saline_summary)
anova(saline_model)

#Urea ----

#Read & factorize data
urea_data <- read.csv("Urea-Survival.csv")
urea_data$Urea.concentration..M. <- factor(urea_data$Urea.concentration..M.)
#urea_data$Duration.of.shock..m. <- factor(urea_data$Duration.of.shock..m.)
#Calculate percent survival
urea_data$pct_surv <- 100*urea_data$Plate.count/mean(
  urea_data$Plate.count[urea_data$Urea.concentration..M. == 0 &
    urea_data$Duration.of.shock..m. == 0])
#Calculate limit of detection in percent
urea_limit_detection <- 100*1/mean(urea_data$Plate.count[
  urea_data$Urea.concentration..M. == 0 &
    urea_data$Duration.of.shock..m. == 0])
#Summarize
urea_summary <- dplyr::summarize(group_by(urea_data, Urea.concentration..M.,
                                          Duration.of.shock..m.),
                                 mean_pct_surv = mean(pct_surv))

#Flag & adjust values below limit of detection
urea_summary$bd <- urea_summary$mean_pct_surv < urea_limit_detection
urea_summary$mean_pct_surv[urea_summary$mean_pct_surv < urea_limit_detection] <- 
  urea_limit_detection

#Make variables for plotting
urea_summary$mean_pct_surv_lines <- urea_summary$mean_pct_surv #this is to plot the lines
urea_summary$ptshape <- 16

#Make pch 8 for first 0 and NA for subsequent 0's
#Stop plotting lines after the first bd point
for (conc in urea_summary$Urea.concentration..M.) {
  myrows <- urea_summary$Urea.concentration..M. == conc
  bd_rows <- (myrows & urea_summary$bd)
  urea_summary$ptshape[bd_rows] <- 8
  if (sum(bd_rows) > 1) {
    #identify bd rows where the previous point was also bd
    for (i in 2:sum(bd_rows)) {
      this_row <- which(bd_rows)[i]
      if (urea_summary$bd[this_row - 1]) {
        urea_summary$ptshape[this_row] <- NA
        urea_summary$mean_pct_surv_lines[this_row] <- NA
      }
    }
  }
}
urea_summary$ptshape <- as.factor(urea_summary$ptshape)

#Plot summarized data, duration on X
ggplot(data = urea_summary, aes(x = Duration.of.shock..m.,
                                y = mean_pct_surv,
                                group = Urea.concentration..M.,
                                color = Urea.concentration..M.,
                                shape = ptshape)) +
  geom_point(size = 3) + 
  geom_line(lwd = 1.5, aes(x = Duration.of.shock..m., y = mean_pct_surv_lines)) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Duration of Shock (min)", y = "Percent Survival (%)") +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15) + 
  geom_hline(yintercept = urea_limit_detection, lty = 3, lwd = 1.15) +
  scale_color_manual(name = "Urea\nConcentration (M)",
                     values = my_cols(11)) +
  scale_shape_manual(values = list("8" = 8, "16" = 16)) +
  scale_y_continuous(breaks = c(100, 10, 1),
                     labels = c("100", "10", "1"),
                     trans="log10") +
  guides(shape = FALSE) #don't show shape legend
ggsave(filename = "urea_byduration.tiff", width = 8, height = 5, units = "in")

#Statistics
urea_stats <- urea_summary[as.numeric(as.character(
  urea_summary$Urea.concentration..M.)) <= 5 &
    urea_summary$bd == FALSE, ]
urea_stats$Duration.of.shock..m. <- as.numeric(as.character(
  urea_stats$Duration.of.shock..m.))

#With 0 as reference
urea_stats$Urea.concentration..M. <- relevel(urea_stats$Urea.concentration..M., "0")
urea_model_0 <- lm(log10(mean_pct_surv) ~ Urea.concentration..M. +
                   Urea.concentration..M.:Duration.of.shock..m.,
                 data = urea_stats)

#General results
anova(urea_model_0)

#                                               Df  Sum Sq Mean Sq F value    Pr(>F)    
# Urea.concentration..M.                        5 1.86958 0.37392  229.48 1.779e-11 ***
# Urea.concentration..M.:Duration.of.shock..m.  6 1.25614 0.20936  128.49 3.536e-10 ***

#Contrasts of interest:
#0M with slope 0
#1M with slope 0
#2M with slope 0
#3M with slope 0
#4M with slope 0
#5M with slope 0

#Get contrasts (all slopes are contrasted with 0 already) 
summary(urea_model_0)

#Results
# Urea.concentration..M.0:Duration.of.shock..m.  0.0013618  0.0005309   2.565  0.02478 *  
# Urea.concentration..M.1:Duration.of.shock..m. -0.0023029  0.0006328  -3.639  0.00339 ** 
# Urea.concentration..M.2:Duration.of.shock..m. -0.0019510  0.0006328  -3.083  0.00948 ** 
# Urea.concentration..M.3:Duration.of.shock..m. -0.0052271  0.0006328  -8.260 2.71e-06 ***
# Urea.concentration..M.4:Duration.of.shock..m. -0.0086706  0.0006328 -13.701 1.09e-08 ***
# Urea.concentration..M.5:Duration.of.shock..m. -0.0228416  0.0010365 -22.037 4.48e-11 ***
