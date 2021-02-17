library(ggplot2)
library(ggsignif)
library(dplyr)

###Needed general functions:
#Define function to fix two-tailed tests into one-tailed & apply Bonferroni correction
adjust_ttests <- function(summary_coefficients, alternative = NULL) {
  if (!(alternative %in% c("greater", "less"))) {
    warning("alternative has not been provided, or is not 'greater' or 'less'")
  }
  
  summary_coefficients <- cbind(summary_coefficients, "One-tailed Pr" = NA)
  #Fix p-values from two-tailed to one-tailed
  for (i in 1:nrow(summary_coefficients)) {
    if (alternative == "greater") { #one tailed test for higher values than expected
      if (summary_coefficients[i, 3] < 0) { #p-values will need to be halved & flipped
        summary_coefficients[i, 5] <- 1-(summary_coefficients[i, 4]/2)
      } else if (summary_coefficients[i, 3] > 0) {#p-values need to be halved
        summary_coefficients[i, 5] <- summary_coefficients[i, 4]/2
      }
    } else if (alternative == "less") { #one tailed test for lower values than expected
      if (summary_coefficients[i, 3] > 0) { #p-values will need to be halved & flipped
        summary_coefficients[i, 5] <- 1-(summary_coefficients[i, 4]/2)
      } else if (summary_coefficients[i, 3] < 0) { #p-values need to be halved
        summary_coefficients[i, 5] <- summary_coefficients[i, 4]/2
      }
    } else {} #remain a two-tailed test
  }
  return(summary_coefficients)
}

#Define function to generate letter groups for plotting of stats significance
make_cld <- function(firstmate = NULL, secondmate = NULL, sig_matrix = NULL,
                     p_vals = NULL, conf.level = 0.95) {
  #compact letter display
  #firstmate and secondmate should be matched lists of the two levels
  # which the p vals correspond to the pairwise test between
  #the order of the levels of firstmate & secondmade will be used to determine
  # lettering (with the earlier levels having lower letters)
  #where e.g. the first entry of each corresponds to the first p_vals entry
  #Otherwise, provide a matrix (sig_matrix) with TRUE/FALSE entries denoting 
  # whether each pair of levels is significantly different from one another
  
  if ((is.null(firstmate) | is.null(secondmate) | is.null(p_vals)) & 
      is.null(sig_matrix)) {
    stop("Either sig_matrix needs to be provided, or firstmate, secondmate, and p_vals must be")
  }
  
  #TODO: add check that sig_matrix is symmetric
  
  #If sig_matrix is not provided, create it
  if (is.null(sig_matrix)) {
    firstmate <- as.factor(firstmate)
    secondmate <- as.factor(secondmate)
    
    all_levels <- unique(c(as.character(levels(firstmate)),
                           as.character(levels(secondmate))))
    sig_matrix <- matrix(nrow = length(all_levels), 
                         ncol = length(all_levels),
                         FALSE)
    for (i in 1:length(firstmate)) {
      if (p_vals[i] < (1-conf.level)) {
        sig_matrix[match(firstmate[i], all_levels),
                   match(secondmate[i], all_levels)] <- TRUE
        sig_matrix[match(secondmate[i], all_levels),
                   match(firstmate[i], all_levels)] <- TRUE
      }
    }
  }
  
  #Run insert and absorb
  # Hans-Peter Piepho (2004) An Algorithm for a Letter-Based Representation of
  # All-Pairwise Comparisons, Journal of Computational and Graphical 
  # Statistics, 13:2, 456-466
  
  #Initialize letters matrix
  letters_matrix <- matrix(nrow = nrow(sig_matrix), ncol = 1,
                           data = TRUE)
  
  #Insert
  for (i in 1:nrow(sig_matrix)) {
    for (j in 1:ncol(sig_matrix)) {
      #If significant and in lower triangle of matrix
      if (sig_matrix[i, j] & lower.tri(sig_matrix)[i, j]) {
        for (k in 1:ncol(letters_matrix)) {
          #If both levels (which are sig dift) are included in this letter
          if (letters_matrix[i, k] & letters_matrix[j, k]) {
            #Duplicate the current column
            letters_matrix <- cbind(letters_matrix, letters_matrix[, k])
            #then remove one treatment from current column, and one from the
            # duplicated column
            letters_matrix[j, ncol(letters_matrix)] <- FALSE
            letters_matrix[i, k] <- FALSE
          }
        }
      }
    }
  }
  
  #Absorb
  i <- 1
  while (i <= ncol(letters_matrix)) { #this is the column we'll compare against
    j <- 1
    while (j <= ncol(letters_matrix)) {
      if (i != j) {
        #if column j TRUE's are a subset of column i
        if (all((letters_matrix[, i] | letters_matrix[, j]) == letters_matrix[, i])) {
          #get rid of column j
          letters_matrix <- letters_matrix[, -c(j)]
          if (j < i) {i <- i-1} #removing the column shifts our ith column too
          j <- j-1
        }
      }
      j <- j+1
    }
    i <- i+1
  }
  
  #TODO: implement sweeping (below is non-functional)
  #Sweep
  # i <- ncol(letters_matrix)
  # #i is the column we're considering whether it's redundant
  # # (working backwards from the end)
  # while (i > 0) {
  #   #For each pair of rows (j, k) that are TRUE in column i, we need to verify 
  #   # that there is another column where those two rows are also TRUE
  #   truerows <- which(letters_matrix[, i])
  #   redundant_pairs <- 0 #status of row i
  #   for (j in truerows) {
  #     for (k in truerows) {
  #       for (m in 1:ncol(letters_matrix)) { #check if j, k both TRUE in ea col m
  #         if (!(m != i & letters_matrix[j, m] & letters_matrix[k, m])) {
  #           redundant <- FALSE
  #         }
  #       }
  #     }
  #   }
  #   i <- i-1
  # }
  
  output <- rep(NA, nrow(letters_matrix))
  for (i in 1:length(output)) {
    output[i] <- paste(letters[1:ncol(letters_matrix)][letters_matrix[i, ]], 
                       collapse = "")
  }
  
  return(output)
}

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

#my_colr <- colorRampPalette(colors = c("#ffcc00", "red"))
my_colr <- function(n) {hcl.colors(n, palette = "lajolla")}

#my_colr <- function(n) {hcl.colors(n, palette = "lajolla")}
tempr_plot <- ggplot(data = tempr_sum[as.character(tempr_sum$Temp) < 80,], 
                     aes(x = Duration, y = mean_pct_surv,
                         group = Temp, color = Temp, shape = ptshape)) +
  geom_point(size = 3, alpha = 0.9) + 
  geom_line(aes(x = Duration, y = mean_pct_surv_lines)) + 
  scale_color_manual(name = "Temperature (°C)", 
                     values = my_colr(7)[2:7]) +
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

#Run stats
tempr_model_55 <- lm(log10(mean_pct_surv)~Temp + Duration:Temp,
                     data = tempr_sum_stats)
anova(tempr_model_55)
qqnorm(tempr_model_55$residuals)

#Get contrasts (all slopes are contrasted with 0 already) 
tempr_coeffs <- summary(tempr_model_55)$coefficients

#Take the subset we actually care about (just the slopes)
tempr_coeffs <- tempr_coeffs[5:8, ]

#Adjust t-tests to be one-tailed
tempr_coeffs <- adjust_ttests(tempr_coeffs, alternative = "less")

#Add corrected p-values
tempr_coeffs <- cbind(tempr_coeffs, 
                      "Adj p" = p.adjust(tempr_coeffs[, "One-tailed Pr"], 
                                         method = "bonferroni"))
tempr_coeffs

#For degrees of freedom
summary(tempr_model_55)

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

my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#000000")

#Make plot
temprdur_plot <- ggplot(data = temprdur_sum, aes(x = Duration.of.shock..m.,
                                y = pct_mean,
                                group = Sample,
                                fill = Sample,
                                color = Sample,
                                shape = Sample)) +
  geom_point(size = 2.5, alpha = 0.8) + geom_line() +
  scale_y_continuous(breaks = c(100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001),
                     labels = c("100", "10", "1", "0.1", "0.01", "0.001",
                                "0.0001", "0.00001"),
                     trans="log10") +
  scale_fill_manual(name = "Phage Stock", values = my_cols) +
  scale_color_manual(name = "Phage Stock", values = my_cols) +
  scale_shape_manual(name = "Phage Stock", values = c(21:25, 4)) +
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
summary(temprdur_model)

#Growth curve analysis ----

gc_data1 <- read.csv("Data-Plate Reader-190817.csv", stringsAsFactors = F)
gc_layout1 <- read.csv("Data-Plate Reader-190817_layout.csv", stringsAsFactors = F)
gc_data2 <- read.csv("Data-Plate Reader-150817.csv", stringsAsFactors = F)
gc_layout2 <- read.csv("Data-Plate Reader-150817_layout.csv", stringsAsFactors = F)
gc_data3 <- read.csv("2018-10-12 Growth Curve.csv", stringsAsFactors = F)
gc_layout3 <- read.csv("2018-10-9 Plate Layout.csv", stringsAsFactors = F)

gc_layout3[gc_layout3=="LB"] <- NA

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

gc_layout1 <- layout_cleanup(gc_layout1)
gc_layout2 <- layout_cleanup(gc_layout2)
gc_layout3 <- layout_cleanup(gc_layout3)

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

gc_data1 <- clean_mdata(merge_lay_data(gc_layout1, gc_data1))
gc_data2 <- clean_mdata(merge_lay_data(gc_layout2, gc_data2))
gc_data3 <- clean_mdata(merge_lay_data(gc_layout3, gc_data3))

#Split contents up
split_contents <- function(input_frame) {
  #Splits contents at every underscore
  #Uses the column name of the first column as the format for new columns
  # where the split contents are distributed into
  colnames <- strsplit(input_frame$format[1], "_")[[1]]
  input_frame[colnames] <- do.call(rbind, strsplit(input_frame$Contents, "_"))
  return(subset(input_frame, select=-c(format)))
}

gc_data1 <- split_contents(gc_data1)
gc_data2 <- split_contents(gc_data2)
gc_data3 <- split_contents(gc_data3)

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

gc_data1$sm_od <- smooth_data(gc_data1$OD600, smooth_over = 10,
                               subset_by = gc_data1$Contents)
gc_data2$sm_od <- smooth_data(gc_data2$OD600, smooth_over = 10,
                               subset_by = gc_data2$Contents)
gc_data3$sm_od <- smooth_data(gc_data3$OD600, smooth_over = 10,
                               subset_by = gc_data3$Contents)

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
gc_data1 <- dplyr::group_by(gc_data1[!is.na(gc_data1$sm_od), ], 
                     bacteria, phage, phageshock, bactshock, mediashock, Rep, Contents)
gc_data2 <- dplyr::group_by(gc_data2[!is.na(gc_data2$sm_od), ],
                             pctmediashocked, bacteria, phage, phageshock, Rep, Contents)
gc_data3 <- dplyr::group_by(gc_data3[!is.na(gc_data3$sm_od), ],
                             bact, phage, totalpfuinoc, Rep, Contents)
                             
#Get OD peak height & time for each growth curve
out_data1 <- dplyr::summarize(gc_data1, 
                              max = analyze_curves(sm_od, Time, bandwidth = 20, 
                                                   return = "max"),
                              maxtime = analyze_curves(sm_od, Time, bandwidth = 20, 
                                                       return = "maxtime"))
out_data2 <- dplyr::summarize(gc_data2, 
                              max = analyze_curves(sm_od, Time, 
                                                   bandwidth = 20, return = "max"),
                              maxtime = analyze_curves(sm_od, Time, 
                                                       bandwidth = 20, return = "maxtime"))
out_data3 <- dplyr::summarize(gc_data3, 
                              max = analyze_curves(sm_od, Time, 
                                                   bandwidth = 20, return = "max"),
                              maxtime = analyze_curves(sm_od, Time, 
                                                       bandwidth = 20, return = "maxtime"))

#Add other information for gc1 ----
#Add information about quantity of pfu inoculated to out_data1
out_data1$pfu_inoc <- 200
out_data1$pfu_inoc[(out_data1$phage == "S3" | out_data1$phage == "S8") &
                     out_data1$phageshock == 270] <- 100
out_data1$pfu_inoc[(out_data1$phage == "S11" | out_data1$phage == "S16") &
                     out_data1$phageshock == 270] <- 50
out_data1$pfu_inoc[out_data1$phage == "R3" & out_data1$phageshock == 360] <- 50
out_data1$pfu_inoc[out_data1$phage == "NA"] <- 0

#For gc1 add the treatment names
out_data1$plot <- NA
for (i in 1:nrow(out_data1)) {
  if (!out_data1$phage[i] == "NA") {
    out_data1$plot[i] <- out_data1$phageshock[i]
  } else {
    if (out_data1$bacteria[i]=="NA") {out_data1$plot[i] <- "- Ctrl"
    } else {
      if (out_data1$mediashock[i] == 1) {
        out_data1$plot[i] <- "+ Ctrl Shock"
      } else {out_data1$plot[i] <- "+ Ctrl"}
    }
  }
}

#Rename Phage Stocks
mynames <- c("A" = "R3", "B" = "S3", "C" = "S8", "D" = "S11", "E" = "S16")
out_data1$phage <- names(mynames)[match(out_data1$phage, mynames)]

#Summarize replicate wells ----
out_data1$pfu_inoc <- as.numeric(as.character(out_data1$pfu_inoc))
out_data1 <- group_by(out_data1, plot, phage)
sum_data1 <- summarize(out_data1,
                       pfu_inoc = mean(pfu_inoc),
                       maxpeak_mean = mean(max),
                       maxpeak_se = sd(max)/n())

#Exploratory analysis ----

#Make plot with all data
out_data1$plot <- factor(out_data1$plot,
                        levels = c("+ Ctrl", "+ Ctrl Shock", 0, 5, 90, 180,
                                   270, 360, "- Ctrl"))
ggplot(data = out_data1, aes(x = plot, y = max, group = phage, color = phage)) +
  geom_point(size = 2, position = position_dodge(0.6)) +
  ylab("Max Bacterial Density (OD600)") + xlab("Duration of Heatshock (min)") +
  scale_colour_discrete(name="Phage Replicate") +
  theme_bw() + scale_y_continuous(limits = c(0, NA))

#There's something weird going on with 180 C (S8) & D (S11)

#Look at all wells and check peak designation accuracy
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

#view_peaks(gc_data1, out_data1)
#view_peaks(gc_data2, out_data2)
#view_peaks(gc_data3, out_data3)

#There appear to be no bacteria in 180 D (S11) at all (all 3 reps)
#And there also appear to be no bacteria in 180 C (S8) replicate C

#Exclude data points which were not inoculated with bacteria
out_data1 <- out_data1[-which(out_data1$max < 0.1 & out_data1$plot != "- Ctrl"), ]

#Re-summarize
sum_data1 <- summarize(out_data1,
                       pfu_inoc = mean(pfu_inoc),
                       maxpeak_mean = mean(max),
                       maxpeak_se = sd(max)/n())

##Statistics ----

#ANOVA (with stocks unsummarized)
gc1_stats <- out_data1[out_data1$plot %in% 
                     c("+ Ctrl", "+ Ctrl Shock", 0, 5, 90, 180,270, 360) &
                     out_data1$pfu_inoc %in% c(0, 200), ]
gc1_stats$pfu_inoc <- relevel(as.factor(gc1_stats$pfu_inoc), ref = "200")
gc1_stats$plot <- relevel(as.factor(gc1_stats$plot), ref = "0")
gc1_stats$phage[is.na(gc1_stats$phage)] <- "None"
curve1_model <- lm(max~plot + phage + phage:plot,
                   data = gc1_stats)
anova(curve1_model)
#summary(curve1_model)

#Contrasts of interest:
# Each of the phage stocks at 0 w/ +Ctrl and +Shock (10x)
# +Ctrl with +Shock
# 0, 5, 90, 180, 270 - all pairwise among each other (10x), using summarized
#   stocks values

#Each phage stock at 0 w/ +Ctrl and +Shock
# with +Ctrl
curve1_coefs <- data.frame("Estimate" = vector(length = 10),
                                "Std.error" = vector(length = 10),
                                "df" = vector(length = 10),
                                "tvalue" = vector(length = 10),
                                "One.tailed.Pr" = vector(length = 10))
i <- 1
for (group1 in c("+ Ctrl", "+ Ctrl Shock")) {
  for (group2 in c("A", "B", "C", "D", "E")) {
    temp <- t.test(gc1_stats$max[gc1_stats$plot == 0 & gc1_stats$phage == group2],
                   gc1_stats$max[gc1_stats$plot == group1],
                   alternative = "less")
    curve1_coefs[i, ] <- data.frame("Estimate" = temp$estimate[1] - temp$estimate[2],
                                    "Std.error" = temp$stderr,
                                    "df" = temp$parameter,
                                    "t.value" = temp$statistic,
                                    "One.tailed.Pr" = temp$p.value)
    rownames(curve1_coefs)[i] <- paste(group1, group2, sep = "-")
    i <- i+1
  }
}

#Add corrected p-values
curve1_coefs <- cbind(curve1_coefs, 
                      "Adj p" = p.adjust(curve1_coefs[, "One.tailed.Pr"], 
                                           method = "bonferroni"))
print(round(curve1_coefs$Estimate, 2))
print(round(curve1_coefs$df, 2))
print(round(curve1_coefs$tvalue, 2))
print(round(curve1_coefs$`Adj p`, 4))

#+Ctrl with +Shock
t.test(gc1_stats$max[gc1_stats$plot == "+ Ctrl Shock"],
       gc1_stats$max[gc1_stats$plot == "+ Ctrl"],
       alternative = "two.sided")

#All pairs among 0, 5, 90, 180, 270
gc1_stats_heatsonly <- sum_data1[sum_data1$plot %in% c(0, 5, 90, 180, 270) &
                               sum_data1$pfu_inoc %in% c(0, 200), ]
curve1_model_heatshocks <- aov(maxpeak_mean~plot,
                               data = gc1_stats_heatsonly)
#just a check before Tukey that treatment is still sig w/ only a subset of the data
summary(curve1_model_heatshocks)

#Run Tukey test of all pairs
gc_heat_tukey <- TukeyHSD(curve1_model_heatshocks,
                          which = "plot")
gc_heat_tukey

#Get the treatment pairs & p-values to assign letters for plotting
my_split <- unlist(strsplit(rownames(gc_heat_tukey$plot), "-"))
treat1 <- factor(my_split[seq(from = 1, to = length(my_split), by = 2)],
                 levels = c("0", "5", "90", "180", "270"))
treat2 <- factor(my_split[seq(from = 2, to = length(my_split), by = 2)],
                 levels = c("0", "5", "90", "180", "270"))

#Get letters for plotting
gc1_groups <- data.frame(plot = levels(treat1),
                         height = NA,
                         statgroup = make_cld(firstmate = treat1, 
                                              secondmate = treat2, 
                                              p_vals = gc_heat_tukey$plot[, 4]))
for (i in 1:nrow(gc1_groups)) {
  gc1_groups$height[i] <- max(sum_data1$maxpeak_mean[
    as.character(sum_data1$plot) == as.character(gc1_groups$plot[i])])
}

#Growth Curve Publication Figures ----

#Make plot with summarized data
plot1 <- sum_data1
plot1$pfu_inoc <- as.factor(plot1$pfu_inoc)
plot1$phage[is.na(plot1$phage)] <- "NA"

#supress error bars when se < 0.008
plot1$maxpeak_se_plot <- plot1$maxpeak_se
# plot1$maxpeak_se_plot[plot1$maxpeak_se < 0.008] <- 0.0001
# 
# plot1$error_bar_n_all[plot1$maxpeak_se < 0.008] <- 0.001
# plot1$error_bar_n_subset[plot1$maxpeak_se < 0.008] <- 0.001

#Run through and calculate number of groups in each "plot" x val
# so we can adjust the error bars to the right width
plot1$error_bar_n_all <- 0
plot1$error_bar_n_subset <- 0
for (myplot in unique(plot1$plot)) {
  plot1$error_bar_n_all[which(plot1$plot == myplot)] <-
    nrow(plot1[which(plot1$plot == myplot), ])
  plot1$error_bar_n_subset[which(plot1$plot != "- Ctrl" &
                              plot1$pfu_inoc %in% c(0, 200) & 
                              plot1$plot == myplot)] <-
    nrow(plot1[which(plot1$plot != "- Ctrl" &
                       plot1$pfu_inoc %in% c(0, 200) & 
                       plot1$plot == myplot), ])
}

my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#000000")

gc_plot1 <- ggplot(data = plot1[plot1$plot != "- Ctrl" &
                                plot1$pfu_inoc %in% c(0, 200),], 
                   aes(x = plot, y = maxpeak_mean)) +
  geom_point(position = position_dodge(width = .4), size = 3,
             aes(color = phage, fill = phage, shape = phage)) +
  geom_errorbar(position = position_dodge(width = .4),
                aes(ymax = maxpeak_mean + 1.96*maxpeak_se_plot,
                    ymin = maxpeak_mean - 1.96*maxpeak_se_plot,
                    color = phage,
                    width = .7*error_bar_n_subset/5)) +
  theme_bw()  +
  ylim(NA, 1.85) +
  labs(x = "70°C Heat Shock Duration (min)", y = "Peak Bacterial Density (OD600)") +
  scale_color_manual(name = "Phage Stock", values = my_cols) +
  scale_fill_manual(name = "Phage Stock", values = my_cols) +
  scale_shape_manual(name = "Phage Stock",
                     values = c(21:25, 4)) +
  geom_hline(yintercept = plot1$maxpeak_mean[plot1$plot == "- Ctrl"],
             lty = 3, lwd = 1.15) +
  scale_x_discrete(labels = c("+ Ctrl", "+ Ctrl\nShock", "0", "5", 
                              "90", "180", "270")) +
  theme(panel.grid = element_blank()) +
  #  theme(axis.text.x = element_text(angle = 30, size = 12, hjust = 1)) +
  geom_text(data = gc1_groups,
            aes(x = plot, y = height + 0.3, label = statgroup)) +
  geom_signif(comparisons = list(c("+ Ctrl", "0")),
              tip_length = 0.01, y_position = 1.7,
              annotation = ifelse(max(curve1_coefs[1:5, "Adj p"]) < 0.001,
                                  "p<0.001",
                                  paste("p<", max(curve1_coefs[1:5, "Adj p"]),
                                        sep = ""))) +
  geom_signif(comparisons = list(c("+ Ctrl Shock", "0")),
              tip_length = 0.01, y_position = 1.5,
              annotation = ifelse(max(curve1_coefs[1:5, "Adj p"]) < 0.001,
                                  "p<0.001",
                                  paste("p<", max(curve1_coefs[1:5, "Adj p"]),
                                        sep = ""))) +
  NULL
gc_plot1
ggsave(filename = "gc_maxes.tiff", width = 8, height = 5, units = "in")

#Plot all the summarized data
plot1_sum_temp <- plot1
plot1_sum_temp$pfu_inoc <- factor(plot1_sum_temp$pfu_inoc,
                                  levels = c("0, 0", levels(plot1_sum_temp$pfu_inoc)))
plot1_sum_temp$pfu_inoc[plot1_sum_temp$plot == "- Ctrl"] <- "0, 0"
gc_plot1_alldata <- ggplot(data = plot1_sum_temp, 
                   aes(x = plot, y = maxpeak_mean, color = phage,
                       shape = pfu_inoc)) +
  geom_point(position = position_dodge(width = .4), size = 3) +
  geom_errorbar(position = position_dodge(width = .4),
                aes(ymax = maxpeak_mean + 1.96*maxpeak_se,
                    ymin = maxpeak_mean - 1.96*maxpeak_se,
                    width = .7*error_bar_n_all/5)) +
  theme_bw()  +
  theme(panel.grid = element_blank()) +
  labs(x = "70°C Heat Shock Duration (min)", y = "Peak Bacterial Density (OD600)") +
  scale_color_manual(name = "Phage Stock", values = my_cols) +
  scale_fill_manual(name = "Phage Stock", values = my_cols) +
  scale_shape_manual(name = "Phage; Bacteria\nInoculated (pfu; cfu)",
                     values = c(4, 18, 15, 17, 16),
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
# my_cols <- colorRampPalette(c("light gray", "black"))
my_cols <- function(n) {hcl.colors(n, palette = "lajolla")}

ggplot(data = saline_summary, aes(x = Duration.of.shock..m., 
                                  y = mean_pct_surv,
                                  group = Saline.Concentration..M.)) +
  geom_point(size = 3.5, pch = 21, aes(fill = Saline.Concentration..M.)) + 
  geom_line(lwd = 1.5, aes(color = Saline.Concentration..M.)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Duration of Salt Shock (min)", 
       y = "Percent Phage Survivors (%)") +
  scale_color_manual(name = "Saline\nConcentration (M)",
                       values = my_cols(8)) +
  scale_fill_manual(name = "Saline\nConcentration (M)",
                    values = my_cols(8)) +
  scale_y_continuous(breaks = c(100, 10, 1),
                     labels = c("100", "10", "1"),
                     trans="log10") +
  scale_x_continuous(limits = c(0, 90),
                     breaks = c(0, 30, 60, 90)) +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15) +
 geom_hline(yintercept = saline_limit_detection, lty = 3, lwd = 1.15) +
  NULL
ggsave(filename = "saline_bydur.tiff", width = 8, height = 5, units = "in")

ggplot(data = saline_summary, aes(x = Duration.of.shock..m., 
                                  y = mean_pct_surv,
                                  group = Saline.Concentration..M.)) +
  geom_point(size = 3.5, pch = 21, aes(fill = Saline.Concentration..M.)) + 
  geom_line(lwd = 1.5, aes(color = Saline.Concentration..M.)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Duration of Salt Shock (min)", 
       y = "Percent Phage Survivors (%)") +
  scale_color_manual(name = "Saline\nConcentration (M)",
                     values = my_cols(8)) +
  scale_fill_manual(name = "Saline\nConcentration (M)",
                    values = my_cols(8)) +
  #  guides(fill = guide_legend(override.aes = list(size=2))) +
  scale_x_continuous(limits = c(0, 90),
                     breaks = c(0, 30, 60, 90)) +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15) +
  NULL
ggsave(filename = "saline_bydur_nolimit.tiff", width = 8, height = 5, units = "in")

#Statistics
saline_summary$Duration.of.shock..m. <- as.numeric(as.character(
  saline_summary$Duration.of.shock..m.))
saline_summary$Saline.Concentration..M. <- relevel(
  saline_summary$Saline.Concentration..M., ref = "0.17 (LB)")
saline_model <- lm(log10(mean_pct_surv) ~ Saline.Concentration..M. +
                     Saline.Concentration..M.:Duration.of.shock..m.,
                   data = saline_summary)
anova(saline_model)

#Urea (Tin) ----

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
#Flag values bd
urea_data$bd <- FALSE
urea_data$bd[urea_data$Plate.count < 1] <- TRUE
#Summarize
#Here we're making any point which has any of the 3 titers below detection
# as below detection itself
urea_summary <- dplyr::summarize(group_by(urea_data, Urea.concentration..M.,
                                          Duration.of.shock..m.),
                                 mean_pct_surv = mean(pct_surv),
                                 bd = any(bd))

#Flag & adjust values below limit of detection
# urea_summary$bd <- urea_summary$mean_pct_surv < urea_limit_detection
urea_summary$mean_pct_surv[which(urea_summary$bd)] <- urea_limit_detection

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
my_cols <- function(n) {hcl.colors(n, palette = "lajolla")}
ggplot(data = urea_summary[as.numeric(as.character(
                        urea_summary$Urea.concentration..M.)) <= 6, ], 
  aes(x = Duration.of.shock..m., y = mean_pct_surv,
      group = Urea.concentration..M.)) +
  scale_shape_manual(values = list("8" = 8, "16" = 21)) +
  geom_point(size = 3, aes(fill = Urea.concentration..M.,
                           shape = ptshape)) + 
  geom_line(lwd = 1.5,
            aes(x = Duration.of.shock..m., y = mean_pct_surv_lines,
                color = Urea.concentration..M.)) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Duration of Urea Shock (min)", 
       y = "Percent Phage Survivors (%)") +
  geom_hline(yintercept = 100, lty = 2, lwd = 1.15) + 
  geom_hline(yintercept = urea_limit_detection, lty = 3, lwd = 1.15) +
  scale_color_manual(name = "Urea\nConcentration (M)",
                     values = my_cols(8)[2:8]) +
  scale_fill_manual(name = "Urea\nConcentration (M)",
                    values = my_cols(8)[2:8]) +
  scale_y_continuous(breaks = c(100, 10, 1),
                     labels = c("100", "10", "1"),
                     trans="log10") +
  scale_x_continuous(breaks = c(0, 30, 60, 90)) +
  guides(shape = FALSE, #don't show shape legend
         #override default legend shape to one that has a fill attr
         fill=guide_legend(override.aes=list(shape=21))) 
ggsave(filename = "urea_byduration.tiff", width = 8, height = 5, units = "in")

#Statistics

#Take subset of data that has 3+ observations in treatment & exclude any bd points
urea_stats <- urea_summary[as.numeric(as.character(
  urea_summary$Urea.concentration..M.)) <= 5 &
    urea_summary$bd == FALSE, ]
urea_stats$Duration.of.shock..m. <- as.numeric(as.character(
  urea_stats$Duration.of.shock..m.))

#Run stats
urea_model_0 <- lm(log10(mean_pct_surv) ~ Urea.concentration..M. +
                     Urea.concentration..M.:Duration.of.shock..m.,
                   data = urea_stats)
anova(urea_model_0)
qqnorm(urea_model_0$residuals)

#Get contrasts (all slopes are contrasted with 0 already) 
urea_coeffs <- summary(urea_model_0)$coefficients

#Take the subset we actually care about (just the slopes)
urea_coeffs <- urea_coeffs[7:12, ]

#Adjust t-tests to be one-tailed
urea_coeffs <- adjust_ttests(urea_coeffs, alternative = "less")

#Add corrected p-values
urea_coeffs <- cbind(urea_coeffs, 
                      "Adj p" = p.adjust(urea_coeffs[, "One-tailed Pr"], 
                                         method = "bonferroni"))
urea_coeffs

#For degrees of freedom:
summary(urea_model_0)

#Urea (Emma) ----

#Read & factorize data
urea_data2 <- read.csv("Urea-Survival-Emma.csv")
urea_data2$Molarity <- as.factor(urea_data2$Molarity)
#Calculate titer
urea_data2$pfu_ml <- urea_data2$Plaques/(10**(urea_data2$Dilution-1))

#Adjust values below bd
urea_data2$bd <- 0
urea_data2$bd[urea_data2$pfu_ml == 0] <- 1
urea_data2$pfu_ml[urea_data2$bd == 1] <- 100

#Calculate percent survival
urea_data2$pct_surv <- NA
for (i in 1:nrow(urea_data2)) {
  urea_data2$pct_surv[i] <- 100*urea_data2$pfu_ml[i]/mean(
    urea_data2$pfu_ml[urea_data2$Date == urea_data2$Date[i] &
                        urea_data2$Molarity == 0 &
                        urea_data2$Time == 0])
}

urea_detec_lim2 = max(urea_data2$pct_surv[urea_data2$bd == 1])

#Summarize
#Here we're making any point which has any of the 3 titers below detection
# as below detection itself
urea2_summary <- dplyr::summarize(group_by(urea_data2, 
                                           Date, Molarity, Time),
                                 mean_pfuml = mean(pfu_ml),
                                 mean_pctsurv = mean(pct_surv),
                                 bd = ifelse(any(bd == 1), 1, 0))

urea2_sum_sum <- dplyr::summarise(group_by(urea2_summary,
                                           Molarity, Time),
                                  mean_pctsurv = mean(mean_pctsurv))

#Plot summarized data, duration on X
my_cols <- function(n) {hcl.colors(n, palette = "lajolla")}
tiff("Urea_byduration2.tiff", width = 6, height = 4, units = "in", res = 300)
ggplot(data = urea2_sum_sum,
       aes(x = Time, y = mean_pctsurv, color = as.factor(Molarity))) +
  geom_line(data = urea2_summary,
            aes(x = as.numeric(Time), y = mean_pctsurv, 
                color = as.factor(Molarity),
                group = paste(Molarity, Date), lty = as.factor(bd)),
            alpha = 0.5, lwd = 0.75,
            position = position_jitter(width = 0, height = 0.07, seed = 1)) +
  geom_line(lwd = 2) +
  scale_y_continuous(breaks = c(100, 10, 1),
                     labels = c("100", "10", "1"),
                     trans="log10") +  
  scale_x_continuous(breaks = c(0, 45, 90)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)) +
  labs(x = "Duration of Urea Shock (min)",
       y = "Percent Phage Survivors (%)") +
  scale_color_manual(name = "Urea\nConcentration (M)",
                     values = my_cols(6)[2:6]) +
  geom_hline(yintercept = urea_detec_lim2, lty = 3, lwd = 1, alpha = 0.5) +
  guides(lty = FALSE) + #don't show shape legend
  #facet_grid(~Date) +
  NULL
dev.off()


#Statistics

#Take subset of data that has 3+ observations in treatment & exclude any bd points
urea2_stats <- urea2_summary[urea2_summary$bd == 0, ]

#Run stats
urea2_model_0 <- lm(log10(mean_pctsurv) ~ Molarity +
                     Molarity:Time,
                   data = urea2_stats)
anova(urea2_model_0)
qqnorm(urea2_model_0$residuals)

#Get contrasts (all slopes are contrasted with 0 already) 
urea2_coeffs <- summary(urea2_model_0)$coefficients

#Take the subset we actually care about (just the slopes)
urea2_coeffs <- urea2_coeffs[6:10, ]

#Adjust t-tests to be one-tailed
urea2_coeffs <- adjust_ttests(urea2_coeffs, alternative = "less")

#Add corrected p-values
urea2_coeffs <- cbind(urea2_coeffs, 
                     "Adj p" = p.adjust(urea2_coeffs[, "One-tailed Pr"], 
                                        method = "bonferroni"))
urea2_coeffs

#For degrees of freedom:
summary(urea2_model_0)
