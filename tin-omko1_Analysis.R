library(ggplot2)
library(dplyr)

#Temperature Survival ----
tempr <- read.csv("Temperature-Survival.csv")
#Note temp in Celcius & Duration in minutes
colnames(tempr)[1:2] <- c("Temp", "Duration")
#The master stock was titered immediately before use and 200 pfu was the expected
# Plate.count
#However, since we observe no decay at 55 & 60C, we'll use those as our
#reference values
tempr$pct_surv <- 100*tempr$Plate.count/mean(tempr$Plate.count[
  tempr$Temp %in% c(55, 60)])
tempr <- group_by(tempr, Duration, Temp)
tempr_sum <- summarize(tempr, 
                       pct_surv_mean = mean(pct_surv))
tempr_sum$Temp <- as.factor(tempr_sum$Temp)
my_colr <- colorRampPalette(colors = c("#ffcc00", "Red"))
tempr_plot <- ggplot(data = tempr_sum, aes(x = Duration, y = pct_surv_mean,
                             group = Temp, color = Temp)) +
  geom_point(size = 3, alpha = 0.7) + geom_line() + 
  scale_color_manual(name = "Temperature (°C)", values = my_colr(6)) +
  scale_x_continuous(breaks = seq(from = 0, to = 90, by = 30),
                     limits = c(0, 90)) +
  scale_y_continuous(breaks = seq(from = 0, to = 125, by = 25)) +
  labs(x = "Heat Shock Duration (min)",
       y = "Percent Survival (%)") +
  theme_bw() +
#  scale_y_continuous(trans="log10") +
  NULL
tempr_plot
ggsave(filename = "temp_surv.tiff", width = 8, height = 5, units = "in")

#Statistics
tempr_sum$Temp <- as.factor(tempr_sum$Temp)
tempr_sum$Duration <- as.numeric(tempr_sum$Duration)
tempr_model <- lm(log10(pct_surv_mean+1)~Temp + Duration:Temp, 
                  data = tempr_sum)
anova(tempr_model)
summary(tempr_model)

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
temprdur$pct_surv = 100*temprdur$Titer/temprdur_sum$titer_mean[match(temprdur$Sample,
                                                                  temprdur_sum$Sample)]
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
  labs(x = "70°C Heat Shock Duration (min)", y = "Percent Survival (%)")
temprdur_plot
ggsave(filename = "temp_duration_surv.tiff", width = 8, height = 5, units = "in")

#Statistics
colnames(temprdur_sum)[2] <- "Duration"
temprdur_model <- lm(log10(pct_mean)~Sample*Duration, data = temprdur_sum)
anova(temprdur_model)
summary(temprdur_model)


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
                           levels = c(NA, 0, 5, 90, 180, 270, 360))
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
  labs(x = "Heat Shock Duration (min)", y = "Peak Bacterial Density (OD600)") +
  scale_color_discrete(name = "Phage Stock") +
  geom_hline(yintercept = plot1_sum$maxpeak_mean[plot1_sum$plot == "- Ctrl"],
             lty = 2, lwd = 1.15) +
  scale_x_discrete(labels = c("+ Ctrl", "+ Ctrl\nShock", "0", "5", 
                              "90", "180", "270"))
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
plot1$pfu_inoc <- relevel(as.factor(plot1$pfu_inoc), ref = "200")
plot1$plot <- relevel(as.factor(plot1$plot), ref = "0")
plot1$phage[is.na(plot1$phage)] <- "None"
curve1_model <- lm(max~plot + phage + phage:plot, 
                   data = plot1[plot1$plot %in% 
                                  c("+ Ctrl", "+ Ctrl Shock", 0, 5, 90, 180, 
                                    270, 360) &
                                  plot1$pfu_inoc %in% c(0, 200),])
anova(curve1_model)
summary(curve1_model)
curve1_aovmodel <- aov(max~plot + phage + phage:plot, 
                       data = plot1[plot1$plot %in% 
                                      c("+ Ctrl", "+ Ctrl Shock", 0, 5, 90, 180, 
                                        270, 360) &
                                      plot1$pfu_inoc %in% c(0, 200),])
summary(curve1_aovmodel) #to check that it's the same numbers
TukeyHSD(curve1_aovmodel, "plot")

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
saline_data$Saline.Concentration..M.[saline_data$Saline.Concentration..M. == "0.17 (LB)"] <- "0.17"
saline_data$Duration.of.shock..m. <- factor(saline_data$Duration.of.shock..m.)
saline_data$Saline.Concentration..M. <- factor(saline_data$Saline.Concentration..M.)
saline_data$pct_surv <- 100*saline_data$Plate.Count/
  mean(saline_data$Plate.Count[saline_data$Saline.Concentration..M. == 0.17])
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
  labs(x = "Duration of Shock (min)", 
       y = "Percent Survival (%)") +
  scale_color_manual(name = "Saline\nConcentration (M)",
                       values = my_cols(8)) +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15)
ggsave(filename = "saline_bydur.tiff", width = 8, height = 5, units = "in")

#Statistics
saline_summary$Duration.of.shock..m. <- as.numeric(as.character(
  saline_summary$Duration.of.shock..m.))
saline_summary$Saline.Concentration..M. <- relevel(
  saline_summary$Saline.Concentration..M., ref = "0.17")
saline_model <- lm(log10(mean_pct_surv) ~ Saline.Concentration..M. +
                     Saline.Concentration..M.:Duration.of.shock..m.,
                   data = saline_summary)
anova(saline_model)
summary(saline_model)

#Urea ----
urea_data <- read.csv("Urea-Survival.csv")
urea_data$Urea.concentration..M. <- factor(urea_data$Urea.concentration..M.)
urea_data$Duration.of.shock..m. <- factor(urea_data$Duration.of.shock..m.)
urea_data$pct_surv <- 100*urea_data$Plate.count/mean(urea_data$Plate.count[
  urea_data$Urea.concentration..M. == 0])
urea_summary <- dplyr::summarize(group_by(urea_data, Urea.concentration..M.,
                                          Duration.of.shock..m.),
                                 mean_pct_surv = mean(pct_surv))

##Plot all data, duration on X
# ggplot(data = urea_data, aes(x = Duration.of.shock..m., y = Plate.count,
#                              group = Urea.concentration..M.,
#                              color = Urea.concentration..M.)) +
#   geom_point(position = position_dodge(0.6))

##Plot all data, conc on X
# ggplot(data = urea_data, aes(x = Urea.concentration..M.,
#                              y = Plate.count,
#                              group = Duration.of.shock..m.,
#                              color = Duration.of.shock..m.)) +
#   geom_point(position = position_dodge(0.4)) +
#   theme_bw()

#Plot summarized data, duration on X
ggplot(data = urea_summary, aes(x = Duration.of.shock..m.,
                                y = mean_pct_surv,
                                group = Urea.concentration..M.,
                                color = Urea.concentration..M.)) +
  geom_point(size = 3) + geom_line(lwd = 1.5) +
  theme_bw() + 
  labs(x = "Duration of Shock (min)", y = "Percent Survival (%)") +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15) + 
  scale_color_manual(name = "Urea\nConcentration (M)",
                     values = my_cols(11))
ggsave(filename = "urea_byduration.tiff", width = 8, height = 5, units = "in")

#Statistics
urea_summary$Duration.of.shock..m. <- as.numeric(as.character(
  urea_summary$Duration.of.shock..m.))
urea_model <- lm(log10(mean_pct_surv+1) ~ Urea.concentration..M. + 
                   Urea.concentration..M.:Duration.of.shock..m.,
                 data = urea_summary)
anova(urea_model)
summary(urea_model)
