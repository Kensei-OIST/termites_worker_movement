## Comparative analysis of termite movement
## This file is for dependencies
## 03/02/2021 last updated

# Packages ---------------------
{
  library(data.table)
  library(stringr)
  
  library(viridis)
  library(ggplot2)
  library(ggrepel)
  require(grid)
  require(gridExtra)
  library(ggridges)
  
  library(exactRankTests)
  library(CircStats)
  
  library(survival)
  library("survminer")
  require(coxme)
  
  library(lme4)
  library(car)
  library(multcomp)
  
  require(phytools)
  library(stringr)
  library(scales)
  
  library(extrafont)
  #font_import(pattern="PT")
  library(hrbrthemes)
  library(rstatix)
  
  library(tidyr)
  library(dplyr)
  loadfonts()
}


# parameters ---------------------
{
  today <- Sys.Date()
  Original.FPS = 25
  Analysis.FPS = 5
  #Dish.Size = 150
  # dish size is inner diameter measured in body lengths python program
  Dish.Size = 145
}

# Functions ---------------------
{
  # Standard Error
  se  <-  function(x){
    y  <-  x[!is.na(x)]  #  remove  the  missing  values
    sqrt(var(as.vector(y))/length(y))
  }
  
  angle_cal <- function(X, Y, Length){
    Ax <- (X[3:Length-1] - X[3:Length-2])
    Bx <- (X[3:Length] - X[3:Length-1])
    Ay <- (Y[3:Length-1] - Y[3:Length-2])
    By <- (Y[3:Length] - Y[3:Length-1])
    hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
    cos <- round((Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5,14)
    return(acos(cos)*hugo)
  }
}
