## Comparative analysis of termite movement
## This file is for Summarize
## Get individual based parameters of movement patterns
## Or obtain species level data by pooling data wihitn species

{
  rm(list = ls())
  PROJHOME = normalizePath(getwd())
  source(file.path(PROJHOME, "scripts/Source.R"))

  # Data location
  rda.place <- file.path(PROJHOME, "data", "rda")
  
  #plot location
  centerdis.dir <- file.path(PROJHOME, "plot/center_distance")
  if (!dir.exists(centerdis.dir)) {dir.create(centerdis.dir)}
  
  
  #other output locations
  other.dir <- file.path(PROJHOME, "plot/other")
  if (!dir.exists(other.dir)) {dir.create(other.dir)}

  # Load data
  load(file.path(rda.place, "All_Trajectory_Data.rda"))
  
  #colors
  color1 <- c("OP"="#9E2A2B", "MP"="#02C39A", "CP"="#1C77C3")
  color2 <- c("Neo_sug"="#9E2A2B", "Zoo_nev"= "#540B0E", "Hod_sjo"="#02C39A", "Cop_for"="#1C77C3", "Odo_for"="#39A9DB")
  
  centerBar_Plot = F
  centerViolin_Plot = F
  centerHist_Plot = F

}

# Species level parameter (pool within species) -------------------------------------------------------------------

{
  Species.list = unique(data.all$Species)
  Name.list = unique(data.all$Name)
  max.dis=Dish.Size/2
}
  

{  
  df.inout.all <- NULL
  for(i in 1:length(Name.list)){
    Name = Name.list[i]
    df.Name = data.all[data.all$Name == Name,]
    Species = unique(df.Name$Species)
    Nesting = unique(df.Name$Nesting)
    under.dis = sum(df.Name$centerdis<max.dis/sqrt(2))*0.2/60  #make inside and outside areas equal
    over.dis = sum(df.Name$centerdis>=max.dis/sqrt(2))*0.2/60  # multiply counts by sec per step (0.2 sec), divide to min
    
    df.inout.all.temp <- data.frame(Species, Nesting, Name, inner=under.dis, outer=over.dis)
    df.inout.all <-rbind(df.inout.all, df.inout.all.temp)
  }
  
  
  df.inout <- melt(as.data.table(df.inout.all), id = c("Species","Nesting", "Name")) 
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
  write.csv(df.inout.agg, file.path(other.dir, "inout.csv"), row.names=FALSE)
}


{
  df.outper <- df.inout.all[,1:(ncol(df.inout.all)-2)]
  df.outper[,c("outer_min")] <- df.inout.all[,c("outer")]
  df.outper[,c("outer")] <- df.inout.all[,c("outer")]/60
  df.outper$Nesting <- factor(df.outper$Nesting , levels=c("OP", "MP", "CP"))
  df.outper$Species <- factor(df.outper$Species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  summary(df.outper$outer)
  
  {
    outper.his<-
    ggplot(df.outper, aes(outer, fill=Nesting)) + 
    geom_histogram(bins = 10) + 
    facet_wrap(~Species)+
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
    print(outper.his)
    ggsave(file.path(centerdis.dir, "outper.his.png"), outper.his, dpi = 200)
  }
  
  df.inper <- df.inout.all[,1:(ncol(df.inout.all)-2)]
  df.inper[,c("inner")] <- df.inout.all[,c("inner")]/60
  df.inper$Nesting <- factor(df.inper$Nesting , levels=c("OP", "MP", "CP"))
  df.inper$Species <- factor(df.inper$Species , levels=c("Neo_sug", "Zoo_nev", "Hod_sjo", "Cop_for", "Odo_for"))
  summary(df.inper$inner)
  
  {
    inper.his<-
      ggplot(df.inper, aes(inner, fill=Nesting)) + 
      geom_histogram(bins = 10) + 
      facet_wrap(~Species)+
      scale_fill_manual(name = "Nesting-type",
                        values = color1) 
    print(inper.his)
    ggsave(file.path(centerdis.dir, "inper.his.png"), inper.his, dpi = 200)
  }
}



{ # add colony row
  Colony <- NULL
  for(row_i in 1:dim(df.inper)[1]){
    Colony = c(Colony, str_split(df.inper$Name, pattern = "_")[[row_i]][3])
  }
  df.inper$Colony = Colony
  
  # logit transformation
  y = df.inper$inner + 0.01
  df.inper$inner.logit = log((y) / (1-(y)))
  
  #create model
  sink(file.path(other.dir, "inper.nesting_Anova.txt"))
  r1 <- lmer(inner.logit ~ Nesting + (1|Species/Colony),  data = df.inper)
  summary(r1)
  Anova(r1, type=2)  
  sink()
  }
  
#boxplots by outper~nesting
{
  outper_nest.box<-
    ggplot(df.outper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Nesting, y=outer))+ 
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
  print(outper_nest.box)
  ggsave(file.path(centerdis.dir, "outper_nest.box.png"), outper_nest.box, dpi = 200)
  
  outper_nestspec.box<-
    ggplot(df.outper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Species, y=outer))+ 
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
  print(outper_nestspec.box)
  ggsave(file.path(centerdis.dir, "outper_nestspec.box.png"), outper_nestspec.box, dpi = 200)
}


#boxplots by inper~nesting
{
  inper_nest.box<-
    ggplot(df.inper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Nesting, y=inner))+ 
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
  print(inper_nest.box)
  ggsave(file.path(centerdis.dir,"inper_nest.box.png"), inper_nest.box, dpi = 200)
  
  inper_nestspec.box<-
    ggplot(df.inper, aes(fill=Nesting)) + 
    geom_boxplot(aes(x=Species, y=inner))+ 
    scale_fill_manual(name = "Nesting-type",
                      values = color1) 
  print(inper_nestspec.box)
  ggsave(file.path(centerdis.dir, "inper_nestspec.box.png"), inper_nestspec.box, dpi = 200)
}

