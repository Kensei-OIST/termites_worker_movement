## Comparative analysis of termite movement
## This file is for Summarize
## Get individual based parameters of movement patterns
## Or obtain species level data by pooling data wihitn species
## 12/21/2021 N Mizumoto last updated

# Setup -------------------------------------------------------------------
{
  rm(list = ls())
  PROJHOME = normalizePath(getwd())
  source(file.path(PROJHOME, "scripts/Source.R"))
  
  # Data location
  rda.place <- file.path(PROJHOME, "data", "rda")
  
  # Load data
  load(file.path(rda.place, "All_Trajectory_Data.rda"))
  
  #colors
  color1 <- c("OP"="#9E2A2B", "MP"="#02C39A", "CP"="#1C77C3")
  color2 <- c("Neo_sug"="#9E2A2B", "Zoo_nev"= "#540B0E", "Hod_sjo"="#02C39A", "Cop_for"="#1C77C3", "Odo_for"="#39A9DB")
}

# Species level parameter (pool within species) -------------------------------------------------------------------

{
  Species.list = unique(data.all$Species)
  pause.threshold = NULL
  
  # T if plotting displacement histograms
  displacements.dir <- file.path(PROJHOME, "plot/displacements")
  if (!dir.exists(displacements.dir)) {dir.create(displacements.dir)}
  
  for(i in 1:length(Species.list)){
    df.temp = data.all[data.all$Species == Species.list[i],]
  
    #create and save histograms
    png(file.path(displacements.dir, paste0((paste(Species.list[i], "displace", sep="_")), ".png")),  width = 250, height = 200, units='mm',res = 600)
    
    h <- hist(df.temp$step.bl[df.temp$step.bl<2], breaks = seq(0,2,0.01), prob= T,
              main=Species.list[i], xlim=c(0,2), col = color2[unname(Species.list[i])])
    x <- h$density
    x[1:10] <- 0
    local.maximum <- h$breaks[ which(diff( sign( diff(x) ) ) == -2)+1 ]
    local.maximum <- h$breaks[which(x == max(x[which(diff( sign( diff(x) ) ) == -2)+1]))]
    arrows(local.maximum, 5, local.maximum, 1, col=1)
    arrows(local.maximum*0.2, 5, local.maximum*0.2, 1, col=1)
    pause.threshold = c(pause.threshold, local.maximum*0.2 + 0.1)
    dev.off()
    }
  
  names(pause.threshold) = Species.list
  save(pause.threshold, file=file.path(rda.place, "Pause_Threshold.rda"))
}



# Individual level parameter -------------------------------------------------------------------
{
  Name.list = unique(data.all$Name)
  df.ind.data = NULL
  for(i in 1:length(Name.list)){
    df.temp = data.all[data.all$Name == Name.list[i],]
    Traveled.dis = sum(df.temp$step, na.rm=T)
    df = data.frame(df.temp[1,c("Genus", "Specific", "Species", "Colony", "Rep", "Name", "Nesting")], Traveled.dis)
    df.ind.data = rbind(df.ind.data, df)
  }
  save(df.ind.data, file=file.path(rda.place, "Ind_data.rda"))
}



# Creating vector for pause and move duration -------------------------------------------------
{
  Name.list = unique(data.all$Name)
  df.duration = NULL
  for(i in 1:length(Name.list)){
    df.temp = data.all[data.all$Name == Name.list[i],]
    
    pause <- df.temp$step.bl < pause.threshold[df.temp$Species[i]]
    move <- !pause
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
      df.temp[1,c("Species","Name")],
      behav = "move",
      duration=move_sec
    )
    df.duration.temp <- rbind(df.duration.temp,
                              data.frame(
                                df.temp[1,c("Species","Name")],
                                behav = "pause",
                                duration=pause_sec
                              )
    )
    
    df.duration <- rbind(df.duration, df.duration.temp)
  }
  save(df.duration, file=file.path(rda.place, "Duration.rda"))
}


# get summed time for move and pause
{
  other.dir <- file.path(PROJHOME, "plot/other")
  if (!dir.exists(other.dir)) {dir.create(other.dir)}
  
  mean_duration<- df.duration %>% group_by(Species, Name) %>% 
    summarise(sum_move_sec = sum(duration[behav=="move"]),
              sum_pause_sec= sum(duration[behav=="pause"]),
              .groups = "drop") %>% 
    group_by(Species) %>% 
    summarise(mean_move_min = mean(sum_move_sec)/60,
              mean_pause_min= mean(sum_pause_sec)/60,
              .groups = "drop") %>% 
    as.data.frame()
  print(mean_duration)
  write.csv(mean_duration, file.path(other.dir, "mean_durations.csv"), row.names=FALSE)
}  
