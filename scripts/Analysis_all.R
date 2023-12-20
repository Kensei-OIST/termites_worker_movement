## Comparative analysis of termite movement
## This file is for Analysis (e.g., statistical analysis)
## 12/18/2023 K K

# Setup -------------------------------------------------------------------
{
  rm(list = ls())
  source("scripts/Source.R")
  
  # Data location
  rda.place <- "data/rda"
  
  #plot location
  movement.dir <- "plot/move_pause"
  if (!dir.exists(movement.dir)) {dir.create(movement.dir)}
  
  #other output locations
  other.dir <- "plot/other"
  if (!dir.exists(other.dir)) {dir.create(other.dir)}
  
  # Load data
  load(file.path(rda.place, "Pause_Threshold.rda"))
  load(file.path(rda.place, "Ind_data.rda"))
  load(file.path(rda.place, "Duration.rda"))
  
  color1 <- c("OP"="#9E2A2B", "MP"="#02C39A", "CP"="#1C77C3")
  color2 <- c("Neo_sug"="#9E2A2B", "Zoo_nev"= "#540B0E", "Hod_sjo"="#02C39A", "Cop_for"="#1C77C3", "Odo_for"="#39A9DB")
}

#functions -------------------------------------------------------------------
## Truncated power-law
# r: 0-1
TP <- function(r,myu,xmin,xmax){
  return( ( xmax^(1-myu) - (1-r)*(xmax^(1-myu)-xmin^(1-myu) ) )^(1/(1-myu)) )
}

## Stretched exponention
# r (0-1)
SE <- function(r, lambda, beta, xmin){
  return( (xmin^beta - 1/lambda*log(1-r)  )^(1/beta))
}

##### Log Likelihood fuction #####
## Truncated Power-Law
# param = myu
TP_LLF <- function(param, data, xmin, xmax){
  length(data)*(log(param-1)-log(xmin^(1-param)-xmax^(1-param))) - param*sum(log(data))
}

## Stretched exponential
# (param[1] = Beta, param[2]= lambda)
SE_LLF <- function(param, data, xmin){
  length(data)*log(param[1])+length(data)*log(param[2])+(param[1]-1)*sum(log(data))-param[2]*sum(data^param[1]-xmin^param[1])
}

## Truncated Power-Law (for binned dataset)
# param = myu
TP_bin_LLF <- function(param, data, xmin, xmax){
  j = 1:(max(data)/0.2)
  dj <- j
  for(i in j){
    dj[i] = sum(round(data*5) == i)
  }
  return(-length(data)*log(xmin^(1-param)-xmax^(1-param)) +sum(dj*log( (xmin+(j-1)*0.2)^(1-param) - (xmin+j*0.2)^(1-param) )))
}

## Stretched exponential (for binned dataset)
# (param[1] = Beta, param[2]= lambda)
SE_bin_LLF <- function(param, data, xmin){
  j = 1:(max(data)/0.2)
  dj <- j
  for(i in 1:length(j)){
    dj[i] = sum(round(data*5) == i)
  }
  return( length(data)*param[2]*xmin^param[1] + sum( dj*log( exp(-param[2]*(xmin+(j-1)*0.2)^param[1]) - exp(-param[2]*(xmin+j*0.2)^param[1]) ) ) )
}

#-------------------------------------------------------------------

# Comparison of traveling distance ------
{
  sink(file.path(other.dir, "Traveled.dis_Anova.txt"))
  r <- lmer(Traveled.dis ~ Nesting + (1|Species/Colony), data=df.ind.data)
  Anova(r)
  sink()
}

# Fitting ------
Plot_behav = F
Plot_behav.all = F

{
  df.behav.all <-NULL
  df.fitting.res <- NULL
  Species.List <- unique(df.duration$Species)
  for(species in 1:length(Species.List)){
    for(behavior in 1:2){
      Species = Species.List[species]
      Behav = c("move","pause")[behavior]
      print(paste(Species, Behav))
      StepData <- df.duration[df.duration$Species == Species &
                                df.duration$behav == Behav , "duration"]
      
      # Obtain info
      min_sec <- min(StepData)
      max_sec <- max(StepData) # Max x set as Max value
      Step_lab <- sort(unique(StepData))
      Step_count <- as.vector(table(StepData))
      Step_ratio <- Step_count / sum(Step_count)
      Step_cum <- rep(0, length(Step_ratio))
      for(i in 1:length(Step_ratio)){
        Step_cum[i] <- sum(Step_ratio[i:length(Step_ratio)])
      }
      
      # Fit with truncated Power-law (pre-binned)
      idea.x <- seq(0,0.99999,0.00001)
      #param_TP <- optimize(TP_LLF, interval=c(1,10), data=StepData, xmin=min_sec, xmax=max_sec, maximum=T)
      param_TP <- optimize(TP_bin_LLF, interval=c(0,10), data=StepData, xmin=min_sec, xmax=max_sec, maximum=T)
      
      y1.est_TP <- TP(idea.x, param_TP$maximum[1], min_sec, max_sec)
      AIC1.TP <- -2*TP_LLF(param_TP$maximum[1], StepData, min_sec, max_sec)+2*1
      
      # fit with stretched exponential (pre-binned)
      param_SE <- optim(c(0.1,0.1), SE_bin_LLF, data=StepData, xmin=min_sec, control = list(fnscale = -1), method="Nelder-Mead")
      #param_SE <- optim(c(0.2,0.2), SE_LLF, data=StepData, xmin=min_sec, control = list(fnscale = -1), method="Nelder-Mead")
      y2.est_SE <- SE(idea.x, param_SE$par[2], param_SE$par[1], min_sec)
      AIC2.SE <- -2*SE_LLF(param_SE$par, StepData, min_sec)+2*2
      
      ## Selection
      AIC <- c(AIC1.TP, AIC2.SE)
      delta <- AIC - min(AIC)
      w1.TP_exp <- exp(-delta[1]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
      w2.SE_exp <- exp(-delta[2]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
      
      if(w1.TP_exp > w2.SE_exp){
        judge <- "TP"
      } else { judge <- "SE" }
      
      ## result plot
      if(Plot_behav){
        png(file.path(PROJHOME, "plot/move_pause", paste0((paste(Species, Behav, sep="_")), ".png")))
        plot(Step_lab,Step_cum, log="xy",axes=F,
             ylab="inverse cumulative frequency", xlab="time", ylim=c(0.00001,1), xlim=c(min_sec, max_sec),
             main=paste(Species, Behav))
        axis(1, c(0.2,0.5,1,2,5,10,20,50,100,200,500))
        axis(2, c(0.001,0.01,0.1,1))
        lines(y1.est_TP, rev(idea.x), col="red", lwd = 1.5)
        lines(y2.est_SE, rev(idea.x), col="blue", lwd = 1.5)
        dev.off()
      }
      
      
      df.all_temp <- data.frame(Species, Behav, Step_lab, Step_cum)
      df.behav.all<- rbind(df.behav.all, df.all_temp)
      
      df <- data.frame(
        Species, Behav, min_sec, max_sec, judge, w1.TP_exp, w2.SE_exp, NumStep = length(StepData),
        TP_myu = param_TP$maximum, SE_beta = param_SE$par[1], SE_lambda = param_SE$par[2]
      )
      df.fitting.res <- rbind(df.fitting.res, df)
      
    }
  }
}



save(df.fitting.res, file=file.path(rda.place, "fitting_res.rda"))
save(df.behav.all, file=file.path(rda.place, "fitting_res_all.rda"))

#if loading from saved file
#load(file.path(rda.place, "fitting_res_all.rda"))

#plot all move and pause
if(Plot_behav.all){
  df.move.all <- df.behav.all[df.behav.all$Behav == "move",]
  df.stop.all <- df.behav.all[df.behav.all$Behav == "pause" ,]
  
  # move_all plot
  png(file.path(PROJHOME, "plot/move_pause/move_all.png"),width = 240, height = 300, units='mm',res=600)
  move_all.p <- ggplot(df.move.all, aes(x=Step_lab,y=Step_cum, col=unname(Species))) +
    geom_point() +
    ggtitle("Move all") +
    #labs(x = "Nesting", y="Std. Trajectory Length") +
    scale_x_continuous(name="time",trans='log10',  breaks=c(0.2,0.5,1,2,5,10,20,50,100,200,500)) +
    scale_y_continuous(name="inverse cumulative frequency", trans='log10', breaks=c(0.001,0.01,0.1,1)) +
    theme_bw() +
    scale_color_manual(name="Species", values = color2)
  
  print(move_all.p)
  dev.off()
  
  # pause_all plot
  png(file.path(PROJHOME, "plot/move_pause/pause_all.png"), width = 240, height = 300, units='mm',res=600)
  stop_all.p <- ggplot(df.stop.all) +
    geom_point(aes(x=Step_lab,y=Step_cum, col=unname(Species))) +
    ggtitle("Stop all") +
    scale_x_continuous(name="time",trans='log10',  breaks=c(0.2,0.5,1,2,5,10,20,50,100,200,500)) +
    scale_y_continuous(name="inverse cumulative frequency", trans='log10', breaks=c(0.001,0.01,0.1,1)) +
    theme_bw()+
    scale_color_manual(name="Species", values = color2)
  print(stop_all.p)
  dev.off()
}

