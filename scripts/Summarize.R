## Comparative analysis of termite movement
## This file is for Summarize
## Get individual based parameters of movement patterns
## Or obtain species level data by pooling data wihitn species
## 12/21/2021 N Mizumoto last updated

# Setup -------------------------------------------------------------------
{
  rm(list = ls())
  source("scripts_20240912/Source_20240912.R")
  
  # Load data
  load("data_fmt/All_Trajectory_Data.rda")
  
  #colors
  color1 <- c("OP"="#9E2A2B", "MP"="#02C39A", "CP"="#1C77C3")
  color2 <- c("Neo_sug"="#9E2A2B", "Zoo_nev"= "#540B0E", 
              "Hod_sjo"="#02C39A", "Cop_for"="#1C77C3", "Odo_for"="#39A9DB")
}

# Filtering
# Sometimes termites try to climb up the plastic wall and fell down.
# This creates unusual high speed movements(especially in N. sugioi). 
# We filtered out these unusual movements by detecting them using (larger than 2sd)
{
  Species.list = unique(data.all$Species)
  
  df.all.filtered <- data.frame()
  df.speed.all <- data.frame()
  for (i_s in 1:length(Species.list)){
    df.species = data.all[data.all$Species == Species.list[i_s],]
    print(paste("Species name:", Species.list[i_s]))
    Name.list <- unique(df.species$Name)
    for (i_n in 1:length(Name.list)){
      Name.temp = Name.list[i_n]
      df.temp = df.species[df.species$Name == Name.temp,]
      
      df.temp$speed.all <- df.temp$step / 0.2 # 0.2 sec = 1 frame, speed = mm/sec
      df.temp$speed.bl.all <- df.temp$step.bl / 0.2   # speed normalized by BL of individual
      df.temp$speed.meanBL.all <- df.temp$step.bl.sp / 0.2  # speed normalized by mean BL of species
      
      mean_speed.all <- mean(df.temp$speed.all[!is.na(df.temp$speed.all)], na.rm = T)
      sd_speed.all <- sd(df.temp$speed.all[!is.na(df.temp$speed.all)], na.rm = T)
      
      upper_bound <- (mean_speed.all + 2 * sd_speed.all) / unique(df.temp$BL)
  
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
  write.csv(df.speed.all, "output/mean.speed.all.filtered.csv", row.names=FALSE)
  write.csv(df.mean.speed.all, "output/mean.speed.all.filtered_spec.csv", row.names=FALSE)
  
  save(df.all.filtered, file="data_fmt/All_Trajectory_Data_filtered.norm.rda")
}
# ---------------------------------------------------------------------------- #

# Obtain move/pause thresholds -------------------------------------------------
{
  load("data_fmt/All_Trajectory_Data_filtered.norm.rda")
  Species.list = unique(df.all.filtered$Species)
  
  # plotting displacement histograms
  pause.threshold = NULL
  displacements.dir <- "output/displacements_norm"
  if (!dir.exists(displacements.dir)) {dir.create(displacements.dir)}
  
  for(i_s in 1:length(Species.list)){
    df.temp = df.all.filtered[df.all.filtered$Species == Species.list[i_s],]
    
    # calculate
    breaks <- seq(min(df.temp$step), max(df.temp$step)+0.001, by = 0.001)
    h <- hist(df.temp$step, breaks = breaks, prob = TRUE,
              main = Species.list[i_s], col = color2[as.character(unname(Species.list[i_s]))])
    
    density_est <- density(df.temp$step) #adjust bw for sensitivity
    min_peak_height <- 0.1
    peaks <- density_est$y[diff(sign(diff(density_est$y))) < 0] # Peaks are local maxima
    valleys <- density_est$y[diff(sign(diff(density_est$y))) > 0]  # Identify local minima in the negative density estimate

    peaks <- peaks[peaks > min_peak_height]
    peak_values <- density_est$x[which(density_est$y %in% peaks)]
    valley_values <- density_est$x[which(density_est$y %in% valleys)]
    
    right_peak <- max(peak_values)
    left_peak <- min(peak_values)

    target <- right_peak*0.2 + 0.1
    
    pause.threshold = c(pause.threshold, target)
    
    # plot
    pdf(file.path(displacements.dir, 
                  paste0((paste(Species.list[i_s], "displace_filtered", sep="_")),
                         "_step_maxpeak02plus01.pdf")),  width = 25, height = 20)
    hist(df.temp$step, breaks = breaks, prob = TRUE,
         main = Species.list[i_s], col = color2[as.character(unname(Species.list[i_s]))])
    lines(density_est, col = "blue")
    points(peak_values, peaks, col = "red", pch = 19)
    points(valley_values, valleys, col = "green", pch = 19)
    arrows(right_peak, 5, right_peak, 0.1, col=2)
    #arrows(right_peak*0.2+0.01 , 5, right_peak*0.2+0.01, 0.1, col=2)
    arrows(target, 5, target, 0.1, col=3)
    dev.off()
  }
  
  names(pause.threshold) = Species.list
  save(pause.threshold, file="data_fmt/Pause_Threshold.rda")
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
    df = data.frame(df.temp[1,c("Genus", "Specific", "Species",
                                "Colony", "Rep", "Name", "Nesting")],
                    Traveled.dis, Traveled.dis.bl, Traveled.dis.bl.sp)
    df.ind.data = rbind(df.ind.data, df)
  }
  save(df.ind.data, file="data_fmt/Ind_data_norm.rda")
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
  save(df.duration, file="data_fmt/Duration_norm.rda")
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
  
  df.sum_mean_duration <- merge(df.sum_mean_duration, 
                                df.ind.data[, c("Name", "Traveled.dis", "Traveled.dis.bl")],
                                by = "Name", all.x = TRUE)
  df.sum_mean_duration$mean_speed <- df.sum_mean_duration$Traveled.dis/df.sum_mean_duration$sum_move_sec
  df.sum_mean_duration$mean_speed_norm <- df.sum_mean_duration$Traveled.dis.bl/df.sum_mean_duration$sum_move_sec
  
  save(df.sum_mean_duration, file="data_fmt/Duration_Summarize_all_norm.rda")
  write.csv(df.sum_mean_duration, "output/durations_summarize_all_norm_step_maxpeak02plus01.csv", row.names=FALSE)
  
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
  write.csv(df.mean_duration_spec, "output/durations_summarize_spec_norm_step_maxpeak02plus01.csv", row.names=FALSE)
  
} 

# LMM for mean velocity; pause time vs Nesting, Species ------------------------
{
  # mean speed
  sink(file.path(other.place,"mean_speed_norm_LMM_step_maxpeak02plus01.txt"))
  
  df.sum_mean_duration$Forager <- df.sum_mean_duration$Nesting != "OP"
  
  # speed in mm/sec while moving
  r <- lmer(mean_speed ~ Forager + (1|Species/Colony), data=df.sum_mean_duration)
  Anova(r)
  multicomparison<-glht(r,linfct=mcp(Nesting="Tukey"))
  summary(multicomparison)
  
  r <- lmer(mean_speed ~ Species + (1|Colony), data=df.sum_mean_duration)
  Anova(r)
  multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
  summary(multicomparison)
  
  # speed in bl/sec while moving
  r <- lmer(mean_speed_norm ~ Forager + (1|Species/Colony), data=df.sum_mean_duration)
  Anova(r)
  
  r <- lmer(mean_speed_norm ~ Species + (1|Colony), data=df.sum_mean_duration)
  Anova(r)
  multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
  summary(multicomparison)
  
  sink()
}

# LMM for total pausing time vs Nesting, Species ------------------------
{
  sink(file.path(other.place,"pause_sec_LMM_step_maxpeak02plus01.txt"))
  
  df.sum_mean_duration$pause_prop <- 
    df.sum_mean_duration$sum_pause_sec / (df.sum_mean_duration$sum_move_sec + df.sum_mean_duration$sum_pause_sec)
  
  y1 = df.sum_mean_duration$pause_prop + 0.01
  df.sum_mean_duration$pause_prop_logit = log((y1) / (1-(y1)))
  
  r <- lmer(pause_prop_logit ~ Forager + (1|Species/Colony), data=df.sum_mean_duration)
  Anova(r)
  
  r <- lmer(pause_prop_logit ~ Species + (1|Colony), data=df.sum_mean_duration)
  Anova(r)
  multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
  summary(multicomparison)
  
  sink()
}

## not checked L302-330 (9/20/2024)
##TODO: confirm what are mean_speed and mean_speed_move

# LMM for total pausing time vs Nesting, Species ------------------------
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
  
  move_pause.place <- "output/move_pause"
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