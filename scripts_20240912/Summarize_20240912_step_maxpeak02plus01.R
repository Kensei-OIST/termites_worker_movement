## Comparative analysis of termite movement
## This file is for Summarize
## Get individual based parameters of movement patterns
## Or obtain species level data by pooling data wihitn species
## 12/21/2021 N Mizumoto last updated

# Setup -------------------------------------------------------------------
{
  rm(list = ls())
  setwd("/Volumes/Phosphorus/All_movement_data/Analysis_data/Japan_5sp")
  PROJHOME = normalizePath(getwd())
  source(file.path(PROJHOME, "20240912/scripts_20240912/Source_20240912.R"))
  
  # Data location
  rda.place <- file.path(PROJHOME, "data", "rda")
  
  other.place <- file.path(PROJHOME, "data/other")
  if (!dir.exists(other.place)) {dir.create(other.place)}
  
  
  # Load data
  load(file.path(rda.place, "All_Trajectory_Data.rda"))
  
  #colors
  color1 <- c("OP"="#9E2A2B", "MP"="#02C39A", "CP"="#1C77C3")
  color2 <- c("Neo_sug"="#9E2A2B", "Zoo_nev"= "#540B0E", "Hod_sjo"="#02C39A", "Cop_for"="#1C77C3", "Odo_for"="#39A9DB")
}

# Species level parameter (pool within species) -------------------------------------------------------------------

{
  Species.list = unique(data.all$Species)
  
  RUN_FILTERING <- FALSE
  if(RUN_FILTERING){
    # filter outlier steps (eg. slips)
    df.all.filtered <- data.frame()
    df.speed.all <- data.frame()
    for (i in 1:length(Species.list)){
      df.species = data.all[data.all$Species == Species.list[i],]
      print(paste("Species name:", Species.list[i]))
      Name.list <- unique(df.species$Name)
      for (i in 1:length(Name.list)){
        Name.temp = Name.list[i]
        df.temp = df.species[df.species$Name == Name.temp,]
        df.temp$speed.all <- df.temp$step / 0.2 # 0.2 sec = 1 frame, speed = mm/sec
        df.temp$speed.bl.all <- df.temp$step.bl / 0.2   # speed normalized by BL of individual
        df.temp$speed.meanBL.all <- df.temp$step.bl.sp / 0.2  # speed normalized by mean BL of species
        speed.all <- df.temp$speed.all[!is.na(df.temp$speed.all)]
        #speed.bl.all <- df.temp$speed.bl.all[!is.na(df.temp$speed.bl.all)]
        mean_speed.all <- mean(speed.all)
        sd_speed.all <- sd(speed.all)
        #mean_speed.bl.all <- mean(speed.bl.all)
        #sd_speed.bl.all <- sd(speed.bl.all)
        # Define the number of standard deviations
        #upper_bound <- mean_speed.all + 2 * sd_speed.all
        upper_bound <- (mean_speed.all + 2 * sd_speed.all)/unique(df.temp$BL)
        #upper_bound.bl <- mean_speed.bl.all + 2 * sd_speed.bl.all
        # Identify outliers (values above the upper bound)
        #outliers <- subset(df.temp, speed >= upper_bound)
        #df.temp <- subset(df.temp, speed.all < upper_bound)
        df.temp <- subset(df.temp, speed.bl.all < upper_bound)
        df.speed.all <- rbind(df.speed.all, 
                              data.frame(Species = unique(df.temp$Species),
                                         Nesting = unique(df.temp$Nesting),
                                         Name = unique(df.temp$Name),
                                         Colony = unique(df.temp$Colony),
                                         mean_speed_all = mean(df.temp$speed.all),
                                         sd_speed_all = sd(df.temp$speed.all),
                                         mean_speed_norm_all = mean(df.temp$speed.bl.all),
                                         sd_speed_norm_all = sd(df.temp$speed.bl.all)
                                         )
                              )
        df.all.filtered <- rbind(df.all.filtered, df.temp)
      }
    }
    
    df.mean.speed.all <- df.speed.all %>%
      group_by(Species, Nesting) %>%
      summarise(
        mean_mean_speed_all = mean(mean_speed_all),
        mean_sd_speed_all = mean(sd_speed_all),
        mean_mean_speed_norm_all = mean(mean_speed_norm_all),
        mean_sd_speed_norm_all = mean(sd_speed_norm_all)
      )
    write.csv(df.speed.all, file.path(other.place, "mean.speed.all.filtered.csv"), row.names=FALSE)
    write.csv(df.mean.speed.all, file.path(other.place, "mean.speed.all.filtered_spec.csv"), row.names=FALSE)
    
    save(df.all.filtered, file=file.path(rda.place, "All_Trajectory_Data_filtered.norm.rda"))
  }
  
  
  
  # correlations for  trajectory length vs BL-------------------------------------------------------------
  
  
  sink(file.path(other.place,"meanspeed_BL_cor_test.txt"))
  by(df.all.filtered, df.all.filtered$Species, function(sub_df) {
    cor(sub_df$speed.all, sub_df$BL)
  })
  by(df.all.filtered, df.all.filtered$Species, function(sub_df) {
    cor.test(sub_df$speed.all, sub_df$BL)
  })
  
  cor(df.all.filtered$speed.all, df.all.filtered$BL)
  cor.test(df.all.filtered$speed.all, df.all.filtered$BL, method = "pearson")
  sink()
  
  
  sink(file.path(other.place, "meanspeed_BL_cor_test2.txt"))
  cat("Correlation by Species:\n\n")
  by(df.all.filtered, df.all.filtered$Species, function(sub_df) {
    cor(sub_df$speed.all, sub_df$BL)
  })
  
  cat("\n\nCorrelation Test by Species:\n\n")
  by(df.all.filtered, df.all.filtered$Species, function(sub_df) {
    cor.test(sub_df$speed.all, sub_df$BL)
  })
  
  cat("\n\nOverall Correlation:\n\n")
  cat(cor(df.all.filtered$speed.all, df.all.filtered$BL), "\n")
  
  cat("\n\nOverall Correlation Test:\n\n")
  print(cor.test(df.all.filtered$speed.all, df.all.filtered$BL, method = "pearson"))

  sink()
  

  
  
  
 #  -----------------------------------------------------------------------------------------------------------
  
  
  
  
  
  
  load(file.path(rda.place, "All_Trajectory_Data_filtered.norm.rda"))
  Species.list = unique(df.all.filtered$Species)

  
  {
    # plotting displacement histograms
    pause.threshold = NULL
    displacements.dir <- file.path(PROJHOME, "plot/displacements_norm")
    if (!dir.exists(displacements.dir)) {dir.create(displacements.dir)}
    
    # find_peaks <- function(x) {
    #   peaks <- which(diff(sign(diff(x))) == -2) + 1
    #   return(peaks)
    # }
    
    for(i in 1:length(Species.list)){
      # df.temp = df.all.filtered[df.all.filtered$Species == Species.list[i],]
      # #create and save histograms
      # pdf(file.path(displacements.dir, paste0((paste(Species.list[i], "displace_filtered", sep="_")), ".pdf")),  width = 25, height = 20)
      # breaks <- seq(min(df.temp$step), max(df.temp$step)+0.001, by = 0.001)
      # h <- hist(df.temp$step, breaks = breaks, prob = TRUE,
      #          main = Species.list[i], col = color2[unname(Species.list[i])])
      # x <- h$density
      # x[1:10] <- 0
      # peaks <- find_peaks(x)
      # largest_peaks <- tail(sort(x[peaks]),10)
      # right_peak <- max(breaks[peaks[which(x[peaks] %in% largest_peaks)]])
      # #left_peak <- min(breaks[peaks[which(x[peaks] %in% largest_peaks)]])
      
      df.temp = df.all.filtered[df.all.filtered$Species == Species.list[i],]
      pdf(file.path(displacements.dir, paste0((paste(Species.list[i], "displace_filtered", sep="_")), "_step_maxpeak02plus01.pdf")),  width = 25, height = 20)
      breaks <- seq(min(df.temp$step), max(df.temp$step)+0.001, by = 0.001)
      h <- hist(df.temp$step, breaks = breaks, prob = TRUE,
                main = Species.list[i], col = color2[as.character(unname(Species.list[i]))])
      
      # if(Species.list[i] %in% c("Cop_for")) {
      #   density_est <- density(df.temp$step, bw = 0.02)
      #   min_peak_height <- 0.05
      # } else {
      #   density_est <- density(df.temp$step) #adjust bw for sensitivity
      #   min_peak_height <- 0.05
      # }
      
      density_est <- density(df.temp$step) #adjust bw for sensitivity
      min_peak_height <- 0.1
      
      peaks <- density_est$y[diff(sign(diff(density_est$y))) < 0] # Peaks are local maxima
      valleys <- density_est$y[diff(sign(diff(density_est$y))) > 0]  # Identify local minima in the negative density estimate
      print(valleys)
      print(peaks)
      peaks <- peaks[peaks > min_peak_height]
      peak_values <- density_est$x[which(density_est$y %in% peaks)]
      valley_values <- density_est$x[which(density_est$y %in% valleys)]
      print(peak_values)
      print(valley_values)
      lines(density_est, col = "blue")
      points(peak_values, peaks, col = "red", pch = 19)
      points(valley_values, valleys, col = "green", pch = 19)
      right_peak <- max(peak_values)
      left_peak <- min(peak_values)
      print(right_peak)
      print(left_peak)
      
      # target_peak_indexes <- which(peak_values > left_peak & peak_values < right_peak)
      # target_peak_values <- peak_values[target_peak_indexes]
      # target_peak_index <- which(peaks == min(peaks[target_peak_indexes]))
      # target_peak <- peak_values[target_peak_index]

      # target_valley_indexes <- which(valley_values > left_peak & valley_values < right_peak)
      # target_valley_values <- valley_values[target_valley_indexes]
      # target_valley_index <- which(valleys == min(valleys[target_valley_indexes]))
      # if(Species.list[i] %in% c("Cop_for")){
      #   target_valley <- target_valley_values[3]
      # } else {
      #   target_valley <- valley_values[target_valley_index]
      # }
      # target_valley <- valley_values[target_valley_index]
      # print(target_valley)
      # target <- target_valley - ((target_valley-left_peak)*0.2)
      
      target <- right_peak*0.2 + 0.1

      pause.threshold = c(pause.threshold, target)
      arrows(right_peak, 5, right_peak, 0.1, col=2)
      #arrows(right_peak*0.2+0.01 , 5, right_peak*0.2+0.01, 0.1, col=2)
      arrows(target, 5, target, 0.1, col=3)

      
      # if(Species.list[i] %in% c("Neo_sug")) {
      #   temp.threshold <- 0.22
      #  } else if (Species.list[i] %in% c("Hod_sjo"))  {
      #    temp.threshold <- 0.8
      #  } else if (Species.list[i] %in% c("Zoo_nev"))  {
      #    temp.threshold <- 0.43
      #  } else if (Species.list[i] %in% c("Cop_for"))  {
      #    temp.threshold <- 1.2
      #  } else if (Species.list[i] %in% c("Odo_for"))  {
      #    temp.threshold <- 1.65
      #  } else {
      #    temp.threshold <- 0.5
      #  }
      # 
      # arrows(temp.threshold, 5, temp.threshold, 0.1, col=2)
      # pause.threshold = c(pause.threshold, temp.threshold)
      
      dev.off()
      }
    
    names(pause.threshold) = Species.list
    save(pause.threshold, file=file.path(rda.place, "Pause_Threshold.rda"))
  }



# Individual level parameter -------------------------------------------------------------------
{
  Name.list = unique(df.all.filtered$Name)
  df.ind.data = NULL
  for(i in 1:length(Name.list)){
    df.temp = df.all.filtered[df.all.filtered$Name == Name.list[i],]
    Traveled.dis = sum(df.temp$step, na.rm=T)
    Traveled.dis.bl = sum(df.temp$step.bl, na.rm=T)
    Traveled.dis.bl.sp = sum(df.temp$step.bl.sp, na.rm=T)
    df = data.frame(df.temp[1,c("Genus", "Specific", "Species", "Colony", "Rep", "Name", "Nesting")],Traveled.dis, Traveled.dis.bl, Traveled.dis.bl.sp)
    df.ind.data = rbind(df.ind.data, df)
  }
  save(df.ind.data, file=file.path(rda.place, "Ind_data_norm.rda"))
}




# Creating vector for pause and move duration -------------------------------------------------
{
  Name.list = unique(df.all.filtered$Name)
  df.duration = NULL
  for(i in 1:length(Name.list)){
    df.temp = df.all.filtered[df.all.filtered$Name == Name.list[i],]
    print(paste(i, "Species: ", unique(df.temp$Species), ";  File: ",  Name.list[i]))
    pause <- df.temp$step < pause.threshold[as.character(unique(df.temp$Species))]
    move <- !pause
    sum_step_move <- sum(df.temp$step[move])
    sum_step_pause <- sum(df.temp$step[pause])
    sum_step.bl_move <- sum(df.temp$step.bl[move])
    sum_step.bl_pause <- sum(df.temp$step.bl[pause])
    
    # creating vector for pause durations and move durations
    L <- length(pause)
    count <- 0
    in_pause <- 0
    pause_label <- rep(NA,L)
    move_label <- rep(NA,L)
    for(i in 2:L){
      #print(i)
      if(pause[i]){
        if(in_pause==0){
          count <- count + 1
        }
        in_pause <- 1
        pause_label[i] <- count
      } else {
        in_pause <- 0
        move_label[i] <- count
      }
    }
    pause_sec <- 
      as.vector(tapply(pause[1:L], pause_label[1:L], sum))*(1/Analysis.FPS)
    move_sec <- 
      as.vector(tapply(move[1:L], move_label[1:L], sum))*(1/Analysis.FPS)
    
    df.duration.temp <- data.frame(
      df.temp[1,c("Species","Name", "Colony", "Nesting")],
      behav = "move",
      duration=move_sec,
      sum_step = sum_step_move,
      sum_step.bl = sum_step.bl_move
    )
    df.duration.temp <- rbind(df.duration.temp,
                              data.frame(
                                df.temp[1,c("Species","Name", "Colony", "Nesting")],
                                behav = "pause",
                                duration=pause_sec,
                                sum_step = sum_step_pause,
                                sum_step.bl = sum_step.bl_pause
                              ))
    
    df.duration <- rbind(df.duration, df.duration.temp)
  }
  save(df.duration, file=file.path(rda.place, "Duration_norm.rda"))
}


# get summed time for move and pause

{
  df.sum_mean_duration <- df.duration %>%
    group_by(Species, Name, Colony, Nesting) %>%
    summarise(
      sum_move_sec = sum(duration[behav == "move"]),
      sum_pause_sec = sum(duration[behav == "pause"]),
      mean_move_sec = mean(duration[behav == "move"]),
      mean_pause_sec = mean(duration[behav == "pause"]),
      sum_steps_move = unique(sum_step[behav == "move"]),
      sum_steps_pause = unique(sum_step[behav == "pause"]),
      sum_steps_move_norm = unique(sum_step.bl[behav == "move"]),
      sum_steps_pause_norm = unique(sum_step.bl[behav == "pause"]),
      mean_speed_move = sum_steps_move/sum_move_sec,
      mean_speed_move_norm = sum_steps_move_norm/sum_move_sec,
      .groups = "drop"
    ) %>%
    arrange(Species, Name) %>%
    group_by(Species) %>%
    as.data.frame()
  
  df.sum_mean_duration <- merge(df.sum_mean_duration, df.ind.data[, c("Name", "Traveled.dis", "Traveled.dis.bl")], by = "Name", all.x = TRUE)
  df.sum_mean_duration$mean_speed <- df.sum_mean_duration$Traveled.dis/df.sum_mean_duration$sum_move_sec
  df.sum_mean_duration$mean_speed_norm <- df.sum_mean_duration$Traveled.dis.bl/df.sum_mean_duration$sum_move_sec
  
  save(df.sum_mean_duration, file=file.path(rda.place, "Duration_Summarize_all_norm.rda"))
  write.csv(df.sum_mean_duration, file.path(other.place, "durations_summarize_all_norm_step_maxpeak02plus01.csv"), row.names=FALSE)
  
  df.mean_duration_spec <- df.sum_mean_duration %>%
    group_by(Species, Nesting) %>%
    summarise(
      mean_sum_move_sec = mean(sum_move_sec),
      mean_sum_pause_sec = mean(sum_pause_sec),
      mean_speed_move = mean(mean_speed_move),
      mean_speed_move_norm = mean(mean_speed_move_norm),
      mean_speed = mean(mean_speed),
      mean_speed_norm = mean(mean_speed_norm)
    )
  write.csv(df.mean_duration_spec, file.path(other.place, "durations_summarize_spec_norm_step_maxpeak02plus01.csv"), row.names=FALSE)
  
} 


# LMM for mean velocity; pause time vs Nesting, Species-------------------------------------------------------------

#mean speed
sink(file.path(other.place,"mean_speed_norm_LMM_step_maxpeak02plus01.txt"))

r <- lmer(mean_speed ~ Nesting + (1|Species/Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Nesting="Tukey"))
summary(multicomparison)

r <- lmer(mean_speed ~ Species + (1|Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
summary(multicomparison)

r <- lmer(mean_speed_norm ~ Nesting + (1|Species/Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Nesting="Tukey"))
summary(multicomparison)

r <- lmer(mean_speed_norm ~ Species + (1|Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
summary(multicomparison)

sink()




sink(file.path(other.place,"pause_sec_LMM_step_maxpeak02plus01.txt"))

r <- lmer(sum_pause_sec ~ Nesting + (1|Species/Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Nesting="Tukey"))
summary(multicomparison)

r <- lmer(sum_pause_sec ~ Species + (1|Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
summary(multicomparison)

sink()


sink(file.path(other.place,"mean_speed_move_norm_LMM_step_maxpeak02plus01.txt"))

r <- lmer(mean_speed_move ~ Nesting + (1|Species/Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Nesting="Tukey"))
summary(multicomparison)

r <- lmer(mean_speed_move ~ Species + (1|Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
summary(multicomparison)

r <- lmer(mean_speed_move_norm ~ Nesting + (1|Species/Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Nesting="Tukey"))
summary(multicomparison)

r <- lmer(mean_speed_move_norm ~ Species + (1|Colony), data=df.sum_mean_duration)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
summary(multicomparison)

sink()


# boxplot for summed move pause in minutes-------------------------------------------------------------

{
  Name.list = unique(df.sum_mean_duration$Name)
  df.sum_mean_duration$Nesting <- factor(df.sum_mean_duration$Nesting , levels=c("OP", "MP", "CP"))
  df.sum_mean_duration$Species <- factor(df.sum_mean_duration$Species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  
  
  df.sum_mean_duration$Species <- recode(df.sum_mean_duration$Species,
                                         "Cop_for" = "C. formosanus",
                                         "Hod_sjo" = "H. sjostedti",
                                         "Neo_sug" = "N. sugioi",
                                         "Odo_for" = "O. formosanus",
                                         "Zoo_nev" = "Z. nevadensis")
  #unique(df.outper$Species)
  
  df.sum_mean_duration$Species <- as.character(df.sum_mean_duration$Species)
  df.sum_mean_duration <- df.sum_mean_duration %>%
    mutate(Species = case_when(
      Species == "Cop_for" ~ "C. formosanus",
      Species == "Hod_sjo" ~ "H. sjostedti",
      Species == "Neo_sug" ~ "N. sugioi",
      Species == "Odo_for" ~ "O. formosanus",
      Species == "Zoo_nev" ~ "Z. nevadensis",
      TRUE ~ Species  # Keep the original value if no match
    ))
  
  
  df.sum_mean_duration$Species <- factor(df.sum_mean_duration$Species , levels=c("N. sugioi", "Z. nevadensis", "H. sjostedti", "C. formosanus", "O. formosanus"))
  
  color3 <- c("N. sugioi"="#9E2A2B", "Z. nevadensis"= "#540B0E", "H. sjostedti"="#02C39A", "C. formosanus"="#1C77C3", "O. formosanus"="#39A9DB")
  
  move_pause.place <- file.path(PROJHOME, "plot/move_pause")
  if (!dir.exists(move_pause.place)) {dir.create(move_pause.place)}
  
  
  #move sum
  move_sum.box<-
    ggplot(df.sum_mean_duration, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Species, y=sum_move_sec/60), lwd=0.5) + 
    #geom_jitter(aes(x=Species, y=sum_move_sec/60), color = "black", size = 1, alpha = 0.5, position = position_jitter(width = 0.2)) +
    labs(x = "Species", y="Total time moved (min)") +
    coord_cartesian(ylim = c(0, 60)) +
    scale_fill_manual(name = "Nesting-type",
                      values = color1) +
    theme_bw() +
    theme(
      legend.position = "none",  # Remove legend
      axis.text.x = element_text(size = 6, family = "Times", face = "italic", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6, family = "Helvetica"),
      axis.title.x = element_text(size = 7, family = "Helvetica"),
      axis.title.y = element_text(size = 7, family = "Helvetica"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    )
  print(move_sum.box)
  ggsave(file.path(move_pause.place, "move_sum_box_step_maxpeak02plus01.pdf"), move_sum.box, width = 6, height = 5.5, units = "cm", device='pdf', dpi=600)
  
  #pause sum
  pause_sum.box<-
    ggplot(df.sum_mean_duration, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Species, y=sum_pause_sec/60), lwd=0.5) + 
    labs(x = "Species", y="Total time paused (min)") +
    coord_cartesian(ylim = c(0, 60)) +
    scale_fill_manual(name = "Nesting-type",
                      values = color1) +
    theme_bw() +
    theme(
      legend.position = "none",  # Remove legend
      axis.text.x = element_text(size = 6, family = "Times", face = "italic", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6, family = "Helvetica"),
      axis.title.x = element_text(size = 7, family = "Helvetica"),
      axis.title.y = element_text(size = 7, family = "Helvetica"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    )
  print(pause_sum.box)
  ggsave(file.path(move_pause.place, "pause_sum_box_step_maxpeak02plus01.pdf"), pause_sum.box, width = 6, height = 5.5, units = "cm", device='pdf', dpi=600)

}