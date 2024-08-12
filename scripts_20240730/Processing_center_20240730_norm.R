## Comparative analysis of termite movement
## This file is for processing row data
## Load .csv files, preprocess for movement analysis, and save as .rda file
## 12/16/2023 K Kikuchi last updated

# Setup -------------------------------------------------------------------
{
  rm(list = ls())
  PROJHOME = normalizePath(getwd())
  source(file.path(PROJHOME, "scripts_20240401/Source_20240401.R"))
  
  # Data location
  rda.place <- file.path(PROJHOME, "data", "rda")
  raw.place <- file.path(PROJHOME, "data", "raw")
  plot.place <- file.path(PROJHOME, "plot")
  if (!dir.exists(plot.place)) {dir.create(plot.place)}
  other.place <- file.path(PROJHOME, "data", "other")
  if (!dir.exists(other.place)) {dir.create(other.place)}
  
  
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
  df.traj <- NULL
  
  for(i in 1:length(Species.list)){
    #get scale data (pixels per plate diameter (=14.5 mm) and BL) from file "res.csv" created by videoScaleBL.py
    scaledata= list.files(file.path(raw.place, Species.list[i]), full.names = TRUE, pattern = "res.csv")
    scaledf<-  data.frame(fread(scaledata))
    pixelpermm = scaledf[[c("scale")]]/Dish.Size
    names(pixelpermm) = scaledf[[2]]
    Body.Length = scaledf[[c("bodyLength0")]]/pixelpermm
    names(Body.Length) = scaledf[[2]]
    Body.Length.sp = mean(Body.Length)
    
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
      BL.sp= Body.Length.sp
      
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
      #Original.FPS = 25 FPS, one step = 5 frames, 5 steps = 1 sec, one step = 1/5 sec = 0.2 sec
      d_scaled[,1] <- d_scaled[,1] / Original.FPS
      x = d_scaled[,2]; y = d_scaled[,3]

      #get center coordinates from BL data
      center.x= center_x[Name.scale]
      center.y= center_y[Name.scale]
      
      #subtract centerx, y and convert to mm
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
          Genus, Specific, Species, Colony, Rep, Name, Nesting, BL, BL.sp,
          time = d_scaled[, 1],
          x = d_scaled[, 2],
          y = d_scaled[, 3],
          step = c(NA, sqrt(diff(x)^2 + diff(y)^2)),
          step.bl = c(NA, sqrt(diff(x)^2 + diff(y)^2))/BL,
          step.bl.sp = c(NA, sqrt(diff(x)^2 + diff(y)^2))/BL.sp,
          angle = c(NA, angle_cal(x, y, length(x)), NA),
          centerdis = sqrt(x^2+y^2),
          trajlen = sum(c(NA, sqrt(diff(x)^2 + diff(y)^2)), na.rm=TRUE),
          trajlen.bl= sum(c(NA, sqrt(diff(x)^2 + diff(y)^2))/BL, na.rm=TRUE)
        )
        data.all <- rbind(data.all, df.temp)
        } 
      
      if(Traj_Dataframe){
      # Create dataframe for BL and trajectory lengths
      df.traj.temp <- data.frame(
        sample= Name.scale,
        Name= Name,
        Species = Species,
        Colony = Colony,
        Rep = Rep,
        BL,
        Nesting = Nesting,
        trajlen= sum(df.temp$step, na.rm=TRUE),
        trajlen.bl= sum(df.temp$step.bl, na.rm=TRUE))
      df.traj<- rbind(df.traj,df.traj.temp )
      }
      
    }
  }
}  

save(df.traj,file=file.path(rda.place, "Species_Trajectory_Data.rda"))

  # remove data based on trajectory--------------------------------------------------------
{
  removedata <- c()
  data.all <-subset(data.all, !(Name %in% removedata))
  df.traj <- subset(df.traj, !(Name %in% removedata))
  
  
  save(data.all, file=file.path(rda.place, "All_Trajectory_Data.rda"))
  save(df.traj, file=file.path(rda.place, "Species_Trajectory_Data.rda"))
}




load(file.path(rda.place, "All_Trajectory_Data.rda"))
load(file.path(rda.place, "Species_Trajectory_Data.rda"))


# correlations for  trajectory length vs BL-------------------------------------------------------------

trajBL_dir <- file.path("plot/BL-trajectories")
if (!dir.exists(trajBL_dir)) {dir.create(trajBL_dir)}


sink(file.path(trajBL_dir,"trajlen_BL_cor_test.text"))
by(df.traj, df.traj$Species, function(sub_df) {
  cor(sub_df$trajlen, sub_df$BL)
})
by(df.traj, df.traj$Species, function(sub_df) {
  cor.test(sub_df$trajlen, sub_df$BL)
})

cor(df.traj$trajlen, df.traj$BL)
cor.test(df.traj$trajlen, df.traj$BL, method = "pearson")
sink()


#correlations
library("ggpubr")
library(tidyverse)
library(broom)

sink(file.path(trajBL_dir,"trajlen_BL_cor.text"))

#colinearity
cor(df.traj$trajlen, df.traj$BL, method=c("pearson"))
cor.test(df.traj$trajlen, df.traj$BL, method=c("pearson"))

#Building a regression model
reg_model <- lm(trajlen ~ BL, data = df.traj)
reg_model
reg_model.diag.metrics <- augment(reg_model)

reg_model_within_species <- by(df.traj, df.traj$Species, function(sub_df) {
  lm(trajlen ~ BL, data = sub_df)
})
reg_model_within_species
lapply(reg_model_within_species, summary)
reg_model_within_species.diag.metrics <- lapply(reg_model_within_species, function(model) {
  augment(model)
})

reg_model_between_species <- lm(trajlen ~ BL * Species, data = df.traj)
reg_model_between_species
lapply(reg_model_between_species, summary)
reg_model_between_species.diag.metrics <- augment(reg_model_between_species)

sink()


#plot residuals error 
res_error_plot <- ggplot(reg_model.diag.metrics, aes(BL, trajlen)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = BL, yend = .fitted), color = "red", size = 0.3)
print(res_error_plot)
ggsave(filename="trajlen_BL_res_error.pdf", path=trajBL_dir, width = 7, height = 6, device='pdf')

pdf(file = file.path(trajBL_dir, "trajlen_BL_reg-diag.pdf"),  width = 7, height = 6)
par(mfrow = c(2, 2))
plot(reg_model)
dev.off()


res_error_plots <- list()
for (species in names(reg_model_within_species.diag.metrics)) {
  augmented_data <- reg_model_within_species.diag.metrics[[species]]
  res_error_plot <- ggplot(augmented_data, aes(BL, trajlen)) +
    geom_point() +
    stat_smooth(method = lm, se = FALSE) +
    geom_segment(aes(xend = BL, yend = .fitted), color = "red", size = 0.3) +
    labs(title = paste("Residual Error Plot for", species),
         x = "Body Length (BL)", y = "Residuals") +
    theme_minimal()
  res_error_plots[[species]] <- res_error_plot
}
for (species in names(res_error_plots)) {
  ggsave(filename = paste0(species, "_res_error_plot.pdf"),
         plot = res_error_plots[[species]],
         path = trajBL_dir, width = 7, height = 6, device = 'pdf')
}

for (species in names(reg_model_within_species)) {
  pdf(file = file.path(trajBL_dir, paste0("trajlen_BL_reg-diag_", species, ".pdf")),  width = 7, height = 6)
  par(mfrow = c(2, 2))
  model <- reg_model_within_species[[species]]
  plot(model)
  dev.off()
}


res_error_plot <- ggplot(reg_model_between_species.diag.metrics, aes(BL, trajlen)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = BL, yend = .fitted), color = "red", size = 0.3)
print(res_error_plot)
ggsave(filename="trajlen_BL_res_error_between_spec.pdf", path=trajBL_dir, width = 7, height = 6, device='pdf')

pdf(file = file.path(trajBL_dir, "trajlen_BL_reg-diag_between_species.pdf"),  width = 7, height = 6)
par(mfrow = c(2, 2))
plot(reg_model_between_species)
dev.off()


# #within species chieck for no linearity  traveleddistance - BL
# 
# library(ggfortify)
# autoplot(model)
# model.diag.metrics <- model.diag.metrics %>%
#   mutate(index = 1:nrow(model.diag.metrics)) %>%
#   select(index, everything(), -.se.fit, -.sigma)
# 
# plot(model, 1)
# plot(model, 3)
# model2 <- lm(log(trajlen) ~ BL, data = df.traj)
# plot(model2, 3)
# plot(model, 2)


# LMM for  trajectory length vs Nesting-------------------------------------------------------------

# colony <- NULL
# for(i in 1:dim(df.traj)[1]){
#   colony <- c(colony, str_split(df.traj$Name[i], "_")[[1]][3])
# }
# df.traj$colony <- colony
#r <- lmer(trajlen.bl ~ Nesting + (1|Species/colony), data=df.traj)


sink(file.path(other.place,"trajlen_LMM.text"))

r <- lmer(trajlen ~ Nesting + (1|Species/Colony), data=df.traj)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Nesting="Tukey"))
summary(multicomparison)

r <- lmer(trajlen.bl ~ Nesting + (1|Species/Colony), data=df.traj)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Nesting="Tukey"))
summary(multicomparison)

r <- lmer(trajlen ~ Species + (1|Colony), data=df.traj)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
summary(multicomparison)


r <- lmer(trajlen.bl ~ Species + (1|Colony), data=df.traj)
Anova(r)
multicomparison<-glht(r,linfct=mcp(Species="Tukey"))
summary(multicomparison)


sink()








# plot BL vs. trajectory length-------------------------------------------------------------
{
  trajBL_dir <- file.path("plot/BL-trajectories")
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
    ggsave(filename="BL-stdtrajL_nestcol.pdf", path=trajBL_dir, width = 5, height = 5, device='pdf',  dpi = 600)
  }
  
  
  #non-std BL, by nesting
  if(BLTraj_Plot){
    trajplot<- ggplot(df.traj, aes(x=BL, y=trajlen, col=Nesting)) +
      geom_point() +
      ggtitle("BL vs. Trajectory_length") +
      labs(x = "BL (mm)", y="Trajectory Length (mm)") +
      scale_color_manual(name = "Nesting",
                         values = color1) +
      theme_bw()  + 
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(trajplot)
    ggsave(filename="BL_trajL_nestcol.png", path=trajBL_dir, width = 5, height = 5, device='png',  dpi = 600)
    ggsave(filename="BL_trajL_nestcol.pdf", path=trajBL_dir, width = 5, height = 5, device='pdf',  dpi = 600)
  }
  
  
  #std.bl, by Species, PDF
  if(TrajBL_Plot){
    trajplot <- ggplot(df.traj, aes(x=BL, y=trajlen.bl, col=Species)) +
      geom_point() +
      labs(x = "BL (mm)", y="Trajectory Length (mm)") +
      scale_color_manual(name = "Nesting", values = color2) +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 8, family = "Helvetica"),
        axis.text.y = element_text(size = 8, family = "Helvetica"),
        axis.title.x = element_text(size = 10, family = "Helvetica"),
        axis.title.y = element_text(size = 10, family = "Helvetica"),
        legend.position = "right",  # Remove legend
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm")
      )
    print(trajplot)
    ggsave(filename="BL-stdtrajL_speccol.pdf", path=trajBL_dir, width = 8, height = 7, device='pdf', dpi=600)
  }
  
  

  #non-std, by Species
  if(BLTraj_Plot){
    trajplot<- ggplot(df.traj, aes(x=BL, y=trajlen, col=Species)) +
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
  
  #non-std, by Species, PDF
  if(BLTraj_Plot){
    trajplot <- ggplot(df.traj, aes(x=BL, y=trajlen, col=Species)) +
      geom_point() +
      labs(x = "BL (mm)", y="Trajectory Length (mm)") +
      scale_color_manual(name = "Nesting",
                         values = color2) +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 8, family = "Helvetica"),
        axis.text.y = element_text(size = 8, family = "Helvetica"),
        axis.title.x = element_text(size = 10, family = "Helvetica"),
        axis.title.y = element_text(size = 10, family = "Helvetica"),
        legend.position = "none",  # Remove legend
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm")
      )
    print(trajplot)
    ggsave(filename="BL-trajL_speccol.pdf", path=trajBL_dir, width = 8, height = 7, device='pdf', dpi=600)
  }
  
  
  
  
}
  
# plot trajectory length box plots-------------------------------------------------------------
{
  df.traj$Nesting <- factor(df.traj$Nesting , levels=c("OP", "MP", "CP"))
  df.traj$Species <- factor(df.traj$Species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  df.traj$Species <- recode(df.traj$Species,
                                   "Cop_for" = "C. formosanus",
                                   "Hod_sjo" = "H. sjostedti",
                                   "Neo_sug" = "N. sugioi",
                                   "Odo_for" = "O. formosanus",
                                   "Zoo_nev" = "Z. nevadensis")
  
  df.traj<- df.traj %>%
    mutate(Species = case_when(
      Species == "Cop_for" ~ "C. formosanus",
      Species == "Hod_sjo" ~ "H. sjostedti",
      Species == "Neo_sug" ~ "N. sugioi",
      Species == "Odo_for" ~ "O. formosanus",
      Species == "Zoo_nev" ~ "Z. nevadensis",
      TRUE ~ Species  # Keep the original value if no match
    ))
  

  
  color3 <- c("N. sugioi"="#9E2A2B", "Z. nevadensis"= "#540B0E", "H. sjostedti"="#02C39A", "C. formosanus"="#1C77C3", "O. formosanus"="#39A9DB")
  
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
    ggsave(filename="nest-stdtrajL_box_nestcol.pdf", path=trajbox_dir, width = 5, height = 5, device='pdf',  dpi = 600)
  }
  
  if(Traj_boxPlot){
    nestplot<- ggplot(df.traj, aes(x = Nesting, y=trajlen.bl)) +
      geom_boxplot(aes(fill=Species), outlier.shape = NA, lwd=0.5) +
      geom_point(position=position_jitterdodge(), size=0.5, aes(fill=Species, alpha=0.1), show.legend = FALSE) +
      ggtitle("Nesting vs. Std. Trajectory_length") +
      labs(x = "Nesting", y="Std. Trajectory Length") +
      scale_fill_manual(name = "Species", values=color3) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    print(nestplot)
    ggsave(filename="nest-stdtrajL_box_speccol.pdf", path=trajbox_dir, width = 5, height = 5, device='pdf',  dpi = 600)
    }
  
  if(TrajNest_Plot){
    nestplot<- ggplot(df.traj, aes(x=Species, y=trajlen.bl)) +
      geom_boxplot(aes(fill=Nesting), outlier.shape = NA, lwd=0.5) +
      geom_point(position=position_jitterdodge(jitter.width =0.15), size=0.5, aes(alpha =0.1), show.legend = FALSE) +
      labs(x = NULL, y="Std. Trajectory Length") +
      scale_fill_manual(name = "Nesting-type",
                        values = color1) +
      theme_bw() +
      theme(legend.position = "none",  # Remove legend
            axis.text.x = element_text(size = 8, family = "Times", face = "italic"),
            axis.text.y = element_text(size = 8, family = "Helvetica"),
            axis.title.x = element_text(size = 10, family = "Helvetica"),
            axis.title.y = element_text(size = 10, family = "Helvetica"),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
    print(nestplot)
    ggsave(filename="nest-std_trajL_specbox-1-2.pdf", path=trajbox_dir, width = 11, height = 8, units = "cm", device='pdf', dpi=600)
  }
  
  
}

