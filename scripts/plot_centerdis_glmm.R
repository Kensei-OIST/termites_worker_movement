## Comparative analysis of termite movement
## This file is for Summarize
## Get individual based parameters of movement patterns
## Or obtain species level data by pooling data wihitn species

{
  rm(list = ls())
  source("scripts/Source.R")
  
  #plot location
  centerdis.dir <- "output/center_distance"
  if (!dir.exists(centerdis.dir)) {dir.create(centerdis.dir)}
  
  # Load data
  load("data_fmt/All_Trajectory_Data_filtered.norm.rda")
  
  #colors
  color1 <- c("OP"="#9E2A2B", "MP"="#02C39A", "CP"="#1C77C3")
  color2 <- c("Neo_sug"="#9E2A2B", "Zoo_nev"= "#540B0E", "Hod_sjo"="#02C39A", "Cop_for"="#1C77C3", "Odo_for"="#39A9DB")
  
  centerBar_Plot = F
  centerViolin_Plot = F
  centerHist_Plot = F

}

# Species level parameter (pool within species) -------------------------------------------------------------------

{
  Species.list = unique(df.all.filtered$Species)
  Name.list = unique(df.all.filtered$Name)
  max.dis=Dish.Size/2
}
  

{  
  df.inout.all <- NULL
  for(i in 1:length(Name.list)){
    Name = Name.list[i]
    df.Name = df.all.filtered[df.all.filtered$Name == Name,]
    Species = unique(df.Name$Species)
    Nesting = unique(df.Name$Nesting)
    Colony = unique(df.Name$Colony)
    under.dis = sum(df.Name$centerdis<max.dis/sqrt(2))*0.2/60  #make inside and outside areas equal
    over.dis = sum(df.Name$centerdis>=max.dis/sqrt(2))*0.2/60  # multiply counts by sec per step (0.2 sec), divide to min
    
    df.inout.all.temp <- data.frame(Species, Nesting, Name, Colony, inner=under.dis, outer=over.dis)
    df.inout.all <-rbind(df.inout.all, df.inout.all.temp)
  }
  
  df.inout.all$Species <- str_replace_all(df.inout.all$Species, "(\\w{3})\\w*_(\\w{3})\\w*", "\\1_\\2")
  
  df.inout <- melt(as.data.table(df.inout.all), id = c("Species","Nesting", "Name", "Colony")) 
  df.inout$Nesting <- factor(df.inout$Nesting , levels=c("OP", "MP", "CP"))
  df.inout$Species <- factor(df.inout$Species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  df.inout$variable <- factor(df.inout$variable , levels=c("inner","outer" ))
  
  df.inout.mean <-aggregate(cbind(inner, outer) ~ Species + Nesting, data =df.inout.all, FUN=mean)
  df.inout.mean$Nesting <- factor(df.inout.mean$Nesting , levels=c("OP", "MP", "CP"))
  df.inout.mean$Species <- factor(df.inout.mean$Species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  df.inout.sd <-aggregate(cbind(inner, outer) ~ Species + Nesting, data =df.inout.all, FUN=sd)
  df.inout.agg <- merge(df.inout.mean, df.inout.sd, by = c("Species", "Nesting"))
  
  # Rename the columns
  colnames(df.inout.agg) <- c("Species", "Nesting", "inner.mean", "outer.mean", "inner.sd", "outer.sd")
  write.csv(df.inout.agg, "output/inout.csv", row.names=FALSE)
}


{
  df.outper <- df.inout.all[,1:(ncol(df.inout.all)-2)]
  df.outper[,c("outer_min")] <- df.inout.all[,c("outer")]
  df.outper[,c("outer_sec")] <- df.inout.all[,c("outer")]*60
  df.outper[,c("outper")] <- df.inout.all[,c("outer")]/60
  df.outper$Nesting <- factor(df.outper$Nesting , levels=c("OP", "MP", "CP"))
  df.outper$Species <- factor(df.outper$Species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  summary(df.outper$outper)
  write.csv(df.outper, "output/out_sum.csv", row.names=FALSE)
  save(df.outper,file= "data_fmt/Outter_Sum.rda")
  
  df.inper <- df.inout.all[,1:(ncol(df.inout.all)-2)]
  df.inper[,c("inner_min")] <- df.inout.all[,c("inner")]
  df.inper[,c("inner_sec")] <- df.inout.all[,c("inner")]*60
  df.inper[,c("inper")] <- df.inout.all[,c("inner")]/60
  df.inper$Nesting <- factor(df.inper$Nesting , levels=c("OP", "MP", "CP"))
  df.inper$Species <- factor(df.inper$Species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  summary(df.inper$inper)
}

save(df.inper, file = "data_fmt/df_wallfollowing.rda")



# LMM -------------------------------------------------------------
load("data_fmt/df_wallfollowing.rda")
df.inper$Forager <- df.inper$Nesting != "OP"
#Anova
{ # inner time ~ species/nesting
  # logit transformation
  y1 = df.inper$inper + 0.01
  df.inper$inper.logit = log((y1) / (1-(y1)))

  #create model
  sink("output/inper.nesting_Anova.txt")
  r1 <- lmer(inper.logit ~ Forager + (1|Species/Colony),  data = df.inper)
  summary(r1)
  Anova(r1, type=2) 
  multicomparison<-glht(r1,linfct=mcp(Nesting="Tukey"))
  summary(multicomparison)
  
  sink()
  
  sink(file.path(other.place, "inper.species_Anova.txt"))
  r2 <- lmer(inper.logit ~ Species + (1|Colony),  data = df.inper)
  summary(r2)
  Anova(r2, type=2) 
  multicomparison<-glht(r2,linfct=mcp(Species="Tukey"))
  summary(multicomparison)
  sink()
  
}

# outer time ~ species/nesting
{ 
  # logit transformation
  y = df.outper$outper + 0.0001
  df.outper$outper.logit = log((y) / (1-(y)))
  
  #create model
  sink(file.path(other.place, "outper.nesting_Anova.txt"))
  r1 <- lmer(outper.logit ~ Nesting + (1|Species/Colony),  data = df.outper)
  summary(r1)
  Anova(r1, type=2)  
  sink()
  
  sink(file.path(other.place, "outper.species_Anova.txt"))
  r2 <- lmer(outper.logit ~ Species + (1|Colony),  data = df.outper)
  summary(r2)
  Anova(r2, type=2)  
  sink()
  
}


#  outer time ~ species/nesting
{
  y = df.outper$outer_min + 0.0001
  df.outper$outer_min.logit = log((y) / (1-(y)))
  
  #create model
  sink(file.path(other.place, "outer_min.nesting_Anova.txt"))
  r1 <- lmer(outer_min ~ Nesting + (1|Species/Colony),  data = df.outper)
  summary(r1)
  Anova(r1, type=2)  
  sink()
  
  sink(file.path(other.place, "outer_min.species_Anova.txt"))
  r2 <- lmer(outer_min ~ Species + (1|Colony),  data = df.outper)
  summary(r2)
  Anova(r2, type=2)  
  sink()
  # 
  # #create model
  # sink(file.path(other.place, "outer_min.nesting_Anova.txt"))
  # r1 <- lmer(outer_min.logit ~ Nesting + (1|Species/Colony),  data = df.outper)
  # summary(r1)
  # Anova(r1, type=2)  
  # sink()
  # 
  # sink(file.path(other.place, "outer_min.species_Anova.txt"))
  # r2 <- lmer(outer_min.logit ~ Species + (1|Colony),  data = df.outper)
  # summary(r2)
  # Anova(r2, type=2)  
  # sink()

}




  
#boxplots by outper~nesting-------------------------------------------------------------

{
  outper_nest.box<-
    ggplot(df.outper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Nesting, y=outper))+ 
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
  print(outper_nest.box)
  ggsave(file.path(centerdis.dir, "outper_nest.box.png"), outper_nest.box, dpi = 200)
  
  outper_nestspec.box<-
    ggplot(df.outper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Species, y=outper))+ 
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
  print(outper_nestspec.box)
  ggsave(file.path(centerdis.dir, "outper_nestspec.box.png"), outper_nestspec.box, dpi = 200)
}


{
  #PDF
  df.outper$Species <- recode(df.outper$Species,
                              "Cop_for" = "C. formosanus",
                              "Hod_sjo" = "H. sjostedti",
                              "Neo_sug" = "N. sugioi",
                              "Odo_for" = "O. formosanus",
                              "Zoo_nev" = "Z. nevadensis")
  outper_nestspec.box<-
    ggplot(df.outper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Species, y=outper), lwd=0.5) + 
    labs(x = "Species", y="Time spent in outer region (%)") +
    coord_cartesian(ylim = c(0, max(df.outper$outper))) +
    scale_fill_manual(name = "Nesting-type",
                      values = color1) +
    theme_bw() +
    theme(
      legend.position = "none",  # Remove legend
      axis.text.x = element_text(size = 8, family = "Times", face = "italic"),
      axis.text.y = element_text(size = 8, family = "Helvetica"),
      axis.title.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    )
  print(outper_nestspec.box)
  ggsave(file.path(PROJHOME, "plot/center_distance/outper_nestspec_box.pdf"), outper_nestspec.box, width = 10, height = 7, units = "cm", device='pdf', dpi=600)
}


#boxplots by inper~nesting
{
  inper_nest.box<-
    ggplot(df.inper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Nesting, y=inper))+ 
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
  print(inper_nest.box)
  ggsave(file.path(centerdis.dir,"inper_nest.box.png"), inper_nest.box, dpi = 200)
  
  inper_nestspec.box<-
    ggplot(df.inper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Species, y=inper))+ 
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
  print(inper_nestspec.box)
  ggsave(file.path(centerdis.dir, "inper_nestspec.box.png"), inper_nestspec.box, dpi = 200)
}

