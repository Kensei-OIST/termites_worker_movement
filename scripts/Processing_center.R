## Comparative analysis of termite movement
## This file is for processing row data
## Load .csv files, preprocess for movement analysis, and save as .rda file
## 12/16/2023 K Kikuchi last updated

# Setup -------------------------------------------------------------------
{
  rm(list = ls())
  source("scripts/Source.R")
  
  # Data location
  rda.place <- file.path("data", "rda")
  raw.place <- file.path("data", "raw")
  plot.place <- file.path("plot")
  if (!dir.exists(plot.place)) {dir.create(plot.place)}
  
  # Nesting information
  Species.info <- data.frame(
    Genus = c("Zoo", "Neo", "Hod", "Cop", "Odo"),
    Nesting = c("OP", "OP", "MP", "CP", "CP")
  )
  
  #Species list
  Species.list <- list.files(raw.place)
    
  #colors
  color1 <- c("OP"="#9E2A2B", "MP"="#02C39A", "CP"="#1C77C3")
  color2 <- c("Neo_sug"="#9E2A2B", "Zoo_nev"= "#540B0E", "Hod_sjo"="#02C39A", "Cop_for"="#1C77C3", "Odo_for"="#39A9DB")
}



# Check trajectories -------------------------------------------------------------------
{
  Plot = F
  Dataframe = T
  Traj_Dataframe = T
  data.all <- NULL
  Species.list <- list.files(raw.place)
  #Species.list <- c("Cop_for", "Hod_sjo", "Neo_sug" ,"Odo_for", "Zoo_nev")
 
  Colonyplot = F
  BLTraj_Plot = F
  Traj_boxPlot = F
  df.trajlen <- NULL
  
  for(i in 1:length(Species.list)){
    #get scale data (pixels per plate diameter (=14.5 mm) and BL) from file "res.csv" created by videoScaleBL.py
    scaledata= list.files(file.path(raw.place, Species.list[i]), full.names = TRUE, pattern = "res.csv")
    scaledf<-  data.frame(fread(scaledata))
    pixelpermm = scaledf[[c("scale")]]/Dish.Size
    names(pixelpermm) = scaledf[[2]]
    Body.Length = scaledf[[c("bodyLength0")]]/pixelpermm
    names(Body.Length) = scaledf[[2]]
    
    #get center from body length data
    center_x = scaledf[[c("center_x")]]
    names(center_x) = scaledf[[2]]
    center_y = scaledf[[c("center_y")]]
    names(center_y) = scaledf[[2]]

    rawdata = list.files(file.path(raw.place, Species.list[i]), full.names = TRUE, pattern = "position.csv")
    dataname = list.files(file.path(raw.place, Species.list[i]), full.names = FALSE, pattern = "position.csv")
  
    if(length(rawdata) == 0){ next; }
    
    for(v in 1:length(dataname)){
      # Video information
      Genus = sub("_\\w*", "", Species.list[i])
      Specific = sub("[A-Za-z0-9]*_", "", Species.list[i])
      Species = paste(Genus, Specific, sep="_")
      Colony = sub("\\-[0-9]*-position.csv", "", dataname[v])
      Rep = sub(".*-(.*)-.*", "\\1", dataname[v])
      
      Name = paste(Genus, Specific, Colony, Rep, sep="_")
      Name.scale = paste(Colony, Rep, sep="-")
      Nesting = Species.info[Species.info$Genus==Genus,"Nesting"]
      BL= unname(Body.Length[Name.scale])
      
      print(paste("Species: ", i, "/", length(Species.list), ";  File: ",  v, "/", length(dataname), "->", Name))
      
      d <- data.frame(fread(rawdata[v], header=T))
      
      # extract first 60 minutes
      d <- d[1:(Original.FPS * 60 * 60), ]

      # Down sampling
      d <- d[seq(1,dim(d)[1],Original.FPS/Analysis.FPS), ]

      # scaling from measurements
      d_scaled <- d
      #check NA
      which(is.na(d_scaled))
      #Original.FPS = 25 FPS, one step = 5 frames, 5 steps = 1 sec, one step = 1/5 sec
      d_scaled[,1] <- d_scaled[,1] / Original.FPS
      x = d_scaled[,2]; y = d_scaled[,3]

      #get center coordinates from BL data
      center.x= center_x[Name.scale]
      center.y= center_y[Name.scale]
      
      #subtract centerx, y 
      x = (x-center.x)/pixelpermm[Name.scale]
      y = (y-center.y)/pixelpermm[Name.scale]
      d_scaled[,2] = x; d_scaled[,3] = y
      
      #Plot trajectories (high res by species)
      if(Plot){
        rawtraj_dir= file.path(PROJHOME, "plot/raw-trajectories", Species.list[i])
        if (!dir.exists(rawtraj_dir)) {dir.create(rawtraj_dir)}
        rawtrajplot <- ggplot() +
          geom_path(data=d_scaled, aes(x0, y0, color=Species), show.legend=FALSE) +
          coord_fixed() +
          ggtitle(paste(Name)) +
          scale_color_manual(name = "Species",
                             values = color2) +
          theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
                plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                panel.border = element_rect(color="black", fill=NA),
                plot.title = element_text(hjust = 0.5, face="bold", size=18)
          ) +
          labs(x = "x", y = "y")
        ggsave(filename=paste0(Name, ".png"), path=rawtraj_dir, width = 5, height = 5, device='png', dpi=600)
      }
      
      
      #create dataframe
      if(Dataframe){
        df.temp <- data.frame(
          Genus, Specific, Species, Colony, Rep, Name, Nesting, BL,
          time = d_scaled[, 1],
          x = d_scaled[, 2],
          y = d_scaled[, 3],
          step = c(NA, sqrt(diff(x)^2 + diff(y)^2)),
          step.bl = c(NA, sqrt(diff(x)^2 + diff(y)^2))/BL,
          angle = c(NA, angle_cal(x, y, length(x)), NA),
          centerdis = sqrt(x^2+y^2),
          trajlen.bl= sum(c(NA, sqrt(diff(x)^2 + diff(y)^2))/BL, na.rm=TRUE)
        )
        data.all <- rbind(data.all, df.temp)
        } 
      
      if(Traj_Dataframe){
      # Create dataframe for BL and trajectory lengths
      df.trajlen.temp <- data.frame(
        sample= Name.scale,
        Name= Name,
        species = Species,
        BL,
        Nesting = Nesting,
        trajlen= sum(df.temp$step, na.rm=TRUE),
        trajlen.bl= sum(df.temp$step.bl, na.rm=TRUE))
      df.trajlen<- rbind(df.trajlen,df.trajlen.temp )
      }
      
    }
  }
}  
  
  # remove data based on trajectory--------------------------------------------------------
{
  removedata <- c()
  data.all <-subset(data.all, !(Name %in% removedata))
  df.traj <- subset(df.trajlen, !(Name %in% removedata))
  
  
  save(data.all, file=file.path(rda.place, "All_Trajectory_Data.rda"))
  save(data.trajl, file=file.path(rda.place, "Species_Trajectory_Data.rda"))
}
  
  # plot BL vs. trajectory length-------------------------------------------------------------
{
  trajBL_dir <- file.path(PROJHOME, "plot/BL-trajectories")
  if (!dir.exists(trajBL_dir)) {dir.create(trajBL_dir)}
  
  #std.bl, by nesting
  if(BLTraj_Plot){
    trajplot<- ggplot(df.traj, aes(x=BL, y=trajlen.bl, col=Nesting)) +
      geom_point() +
      ggtitle("BL vs. Std. Trajectory_length") +
      labs(x = "BL (mm)", y="Std. Trajectory Length (mm)") +
      scale_color_manual(name = "Nesting",
                         values = color1) +
      theme_bw()  + 
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(trajplot)
    ggsave(filename="BL-stdtrajL_nestcol.png", path=trajBL_dir, width = 5, height = 5, device='png',  dpi = 600)
  }
  
  
  #non-std by nesting
  if(BLTraj_Plot){
    #png(file.path(PROJHOME, "plot/BL-trajectories/))
    trajplot<- ggplot(df.traj, aes(x=BL, y=trajlen, col=Nesting)) +
      geom_point() +
      #geom_text_repel(aes(label = sample), size = 3, show.legend=FALSE) +
      ggtitle("BL vs. Trajectory_length") +
      labs(x = "BL (mm)", y="Trajectory Length (mm)") +
      scale_color_manual(name = "Nesting",
                         values = color1) +
      theme_bw()  + 
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(trajplot)
    ggsave(filename="BL_trajL_nestcol.png", path=trajBL_dir, width = 5, height = 5, device='png',  dpi = 600)
    #dev.off()
  }
  
  #std.bl, by species
  if(BLTraj_Plot){
    trajplot<- ggplot(df.traj, aes(x=BL, y=trajlen.bl, col=species)) +
      geom_point() +
      ggtitle("BL vs. Std. Trajectory_length") +
      labs(x = "BL (mm)", y="Std. Trajectory Length (mm)") +
      scale_color_manual(name = "Nesting",
                         values = color2) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(trajplot)
    ggsave(filename="BL-stdtrajL_speccol.png", path=trajBL_dir, width = 5, height = 5, device='png', dpi=600)
  }
  

  #non-std, by species
  if(BLTraj_Plot){
    trajplot<- ggplot(df.traj, aes(x=BL, y=trajlen, col=species)) +
      geom_point() +
      #geom_text_repel(aes(label = sample), size = 3, show.legend=FALSE) +
      ggtitle("BL vs. Trajectory_length") +
      labs(x = "BL (mm)", y="Trajectory Length (mm)") +
      scale_color_manual(name = "Nesting",
                         values = color2) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(trajplot)
    ggsave(filename="BL_trajL_speccol.png", path=trajBL_dir, width = 5, height = 5, device='png', dpi=600)
  }
  
}
  
# plot trajectory length box plots-------------------------------------------------------------
{
  df.traj$Nesting <- factor(df.traj$Nesting , levels=c("OP", "MP", "CP"))
  df.traj$species <- factor(df.traj$species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  df.traj$species <- recode(df.traj$species,
                                   "Cop_for" = "C. formosanus",
                                   "Hod_sjo" = "H. sjostedti",
                                   "Neo_sug" = "N. sugioi",
                                   "Odo_for" = "O. formosanus",
                                   "Zoo_nev" = "Z. nevadensis")
  trajbox_dir <- file.path(PROJHOME, "plot/trajectory_boxplots")
  if (!dir.exists(trajbox_dir)) {dir.create(trajbox_dir)}
  
  if(Traj_boxPlot){
    nestplot<- ggplot(df.traj, aes(x=Nesting, y=trajlen.bl, col=Nesting)) +
      geom_boxplot(outlier.shape = NA, lwd=0.5) +
      geom_point(size=0.3) +
      ggtitle("Nesting vs. Std. Trajectory_length") +
      labs(x = "Nesting", y="Std. Trajectory Length") +
      scale_color_manual(name = "Nesting-type",
                         values = color1) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(nestplot)
    ggsave(filename="nest-stdtrajL_box_nestcol.png", path=trajbox_dir, width = 5, height = 5, device='png',  dpi = 600)
  }
  
  if(Traj_boxPlot){
    nestplot<- ggplot(df.traj, aes(x = Nesting, y=trajlen.bl)) +
      geom_boxplot(aes(fill=species), outlier.shape = NA, lwd=0.5) +
      geom_point(position=position_jitterdodge(), size=0.5, aes(fill=species, alpha=0.1), show.legend = FALSE) +
      ggtitle("Nesting vs. Std. Trajectory_length") +
      labs(x = "Nesting", y="Std. Trajectory Length") +
      scale_fill_manual(name = "species", values=color2) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(nestplot)
    ggsave(filename="nest-stdtrajL_box_speccol.png", path=trajbox_dir, width = 6, height = 5, device='png',  dpi = 600)
    }
  
  if(Traj_boxPlot){
    nestplot<- ggplot(df.traj, aes(x=species, y=trajlen.bl)) +
      geom_boxplot(aes(fill=Nesting), outlier.shape = NA, lwd=0.5) +
      geom_point(position=position_jitterdodge(jitter.width =0.15), size=0.5, aes(alpha =0.1), show.legend = FALSE)+
      ggtitle("Species vs. Std. Trajectory_length") +
      labs(x = "Species", y="Std. Trajectory Length") +
      scale_fill_manual(name = "Nesting-type",
                         values = color1) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(nestplot)
    ggsave(filename="spec-stdtrajL_box_nestcol.png", path=trajbox_dir, width = 6, height = 5, device='png',  dpi = 600)
  }
  

