# Installation packages
install.packages("biomod2", dep = T)
install.packages("colorRamps", dep = T)
install.packages("dismo", dep = T)
install.packages("dplyr", dep = T)
install.packages("maps", dep = T)
install.packages("maptools", dep = T)
install.packages("rgdal", dep = T)
install.packages("plotKML", dep = T)
install_github("envirometrix/plotKML")
install.packages("raster", dep = T)
install.packages("usdm")
install.packages("foreach", dep = T)
install.packages("doParallel", dep = T)
install.packages("virtualspecies", dep = T)
install.packages("filesstrings")
install.packages("geodata", dep = T)
install.packages("gdalUtils", dep = T)
devtools::install_github("gearslaboratory/gdalUtils")
install.packages("tidyverse")
install.packages("ggcorrplot")

# Loading libraries

library(biomod2)
library(colorRamps)
library(dismo)
library(dplyr)
library(maps)
library(maptools)
library(plotKML)
library(raster)
library(rgdal)
library(RStoolbox)
library(usdm)
library(foreach)
library(doParallel)
library(virtualspecies)
library(filesstrings)
library(geodata)
library(gdalUtils)
library(tools)
library(tidyverse)
library(ggcorrplot)

# Raster file calling - VIF ----

setwd("D:/Envdata/Multicorr_var")
MultiVar <- list.files(path = ".", pattern = "grd")

TempName <- c("bio13", "bio14", "bio15", "bio2", "bio3", "bio4", "gsp", "gst", "kg0", "kg2", "kg3",
              "kg4", "kg5", "npp", "bdod", "cfvo", "ocd", "phh2o", "sand", "silt", "soc", "elev", "rough", "srad", "wavg")

for (i in (1:length(MultiVar))){
  AA <- MultiVar[i]
  BB <- substr(AA, 1, nchar(AA)-4)
  temp <- raster::stack(MultiVar[i]) # for multiband image
  names(temp) <- TempName
  assign(BB, temp)
}

setwd("D:/Envdata/")  

# Loading species occurrence data

crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

args <- list.files("Species/", pattern = "*.csv$",full.names = TRUE) 

args <- sort(args)

for(i in (1:length(args))){
  inputDataFile <- args[i]
  outputFolder <- inputDataFile %>%
    basename %>%
    file_path_sans_ext
  outputFolder
  
  outputFolder1 <- gsub("_", ".", outputFolder) 
  setwd("D:/Envdata/Process")
  
  if (!dir.exists(outputFolder)) {
    dir.create(outputFolder, recursive = TRUE)
  }
  
  setwd("D:/Envdata/")
  occsData <- readr::read_csv(inputDataFile)
  if (dim(occsData)[1] > 500) {occsData <- occsData[sample(1:dim(occsData)[1], 500, replace=FALSE),] }
  occsData <- select(occsData, "Longitude", "Latitude", "sp")
  
  names(occsData)[names(occsData)== outputFolder] <- outputFolder
  
  names(occsData)[3]<-outputFolder 
  sp::coordinates(occsData) <- c("Longitude", "Latitude")
  sp::proj4string(occsData) <- crs.wgs84
  
  
  # 9- Data formatting ####
  occurrences<-data.frame(occsData)
  names(occurrences)[names(occurrences)=="outputFolder"] <- outputFolder
  occurrences<-dplyr::select(occurrences, c("Longitude", "Latitude",outputFolder))
  names(occurrences)
  
  # 10- Modeling ####
  DataSpecies<-occurrences
  myRespName <- outputFolder
  DataSpecies[,myRespName]
  
  #myResp <- as.numeric(DataSpecies[,myRespName])
  DataSpecies[,myRespName] <- 1
  myResp <- DataSpecies[,myRespName]
  
  myExpl <- bio_present
  myRespCoord <- DataSpecies[c("Longitude", "Latitude")]
  
  # Checking !!!! 
  
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName,
                                       filter.raster = TRUE,
                                       PA.nb.rep = 1,
                                       PA.nb.absences = 1000,
                                       PA.strategy = 'random')
  
  plot(myBiomodData)
  str(myBiomodData)
  
  setwd("D:/Envdata/Process/")
  
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('ANN', 'GBM', 'GAM', 'CTA', 'RF','MARS'), #SRE, FDA, MAXENT
    bm.options = myBiomodOption,
    CV.perc = 0.7, 
    CV.nb.rep = 3, 
    metric.eval = c('KAPPA','TSS','ROC'),
    scale.models = FALSE,
    modeling.id = paste(outputFolder)
  )
  
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  write.csv(myBiomodModelEval, file = file.path(outputFolder, "myBiomodModelEval.csv"),
            row.names = FALSE)
  
  getwd()
  myBiomodModelEval
  
  # 11- Projection ####
  myBiomodProj <- BIOMOD_Projection(
    bm.mod = myBiomodModelOut,
    new.env = myExpl,
    proj.name = outputFolder,
    models.chosen = 'all',
    metric.binary = "TSS",
    compress = 'xz',
    build.clamping.mask = TRUE,
    output.format = '.grd')
  
  
  # 12- Ensemble ####
  myBiomodEM <- BIOMOD_EnsembleModeling(
    bm.mod = myBiomodModelOut,
    models.chosen = 'all',
    em.by='all',
    em.algo = c('EMcv', 'EMci', 'EMwmean'),
    metric.select = c('ROC'),
    prob.ci.alpha = 0.05,
    metric.select.thresh = c(0.7),
    var.import = 1)
  
  myVarImportEM<-data.frame(get_variables_importance(myBiomodEM))
  myVarImportEM<-myVarImportEM[5]
  write.csv(myVarImportEM, file = file.path(outputFolder, "myVarImportEM.csv"),
            row.names = T)
  
  #---------------------------- developing from here ------------------------------------
  
  myVarImportEM<-data.frame(get_variables_importance(myBiomodEM))
  myVarImportEM
  
  myVarImportEM<-myVarImportEM[c("full.name", "filtered.by", "algo", "expl.var", "var.imp")]
  write.csv(myVarImportEM, file = file.path(outputFolder, "myVarImportEM.csv"),
            row.names = T)
  
  
  # 13- Ensemble validation ####
  myBiomodEMEval<-get_evaluations(myBiomodEM)
  write.csv(myBiomodEMEval, file = file.path(outputFolder, "myBiomodEMEval.csv"),
            row.names = FALSE)
  
  # 14- Projection under current climatic conditions ####
  myBiomodEM_proj <-BIOMOD_EnsembleForecasting(bm.em  = myBiomodEM,
                                               bm.proj = myBiomodProj,
                                               models.chosen = 'all',
                                               proj.name = outputFolder,
                                               metric.binary = "TSS",
                                               output.format = '.grd')

    outpath<-file.path(outputFolder1,paste0("proj_",outputFolder))
  currentPred <- raster::stack(file.path(outputFolder1,paste("proj_",outputFolder,"/proj_",outputFolder,"_",outputFolder1,"_ensemble_TSSbin.grd",sep="")))
  raster::writeRaster(currentPred,
              file.path(outpath,paste0(outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  
  currentPred_c <- raster::stack(file.path(outputFolder1,paste("proj_",outputFolder,"/proj_",outputFolder,"_",outputFolder1,"_ensemble.grd",sep="")))
  raster::writeRaster(currentPred_c,
              file.path(outpath,paste0(outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  ClampMask <- raster::stack(file.path(outputFolder1,paste("proj_",outputFolder,"/proj_",outputFolder,"_ClampingMask.grd",sep="")))
  raster::writeRaster(ClampMask,file.path(outpath,paste0(outputFolder,"_ClampingMask.tif")),overwrite= TRUE)
  
  getwd()
  
  # 15- Projections under current climatic conditions ####
  myBiomodEM_1140_126 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_1140_126,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_1140_126",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_1140_126/proj_SSP_1140_126_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1,paste0("proj_SSP_1140_126/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_1140_126/proj_SSP_1140_126_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1,paste0("proj_SSP_1140_126/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_1140_126)
  
  myBiomodEM_1140_370 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_1140_370,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_1140_370",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_1140_370/proj_SSP_1140_370_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1,paste0("proj_SSP_1140_370/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_1140_370/proj_SSP_1140_370_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1,paste0("proj_SSP_1140_370/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_1140_370)
  
  
  myBiomodEM_1140_585 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_1140_585,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_1140_585",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_1140_585/proj_SSP_1140_585_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1,paste0("proj_SSP_1140_585/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_1140_585/proj_SSP_1140_585_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1,paste0("proj_SSP_1140_585/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_1140_585)
  
  myBiomodEM_4170_126 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_4170_126,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_4170_126",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_4170_126/proj_SSP_4170_126_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1,paste0("proj_SSP_4170_126/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_4170_126/proj_SSP_4170_126_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1,paste0("proj_SSP_4170_126/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_4170_126)
  
  myBiomodEM_4170_370 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_4170_370,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_4170_370",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_4170_370/proj_SSP_4170_370_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1,paste0("proj_SSP_4170_370/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_4170_370/proj_SSP_4170_370_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1, paste0("proj_SSP_4170_370/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_4170_370)
  
  
  myBiomodEM_4170_585 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_4170_585,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_4170_585",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_4170_585/proj_SSP_4170_585_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1,paste0("proj_SSP_4170_585/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_4170_585/proj_SSP_4170_585_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1,paste0("proj_SSP_4170_585/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_4170_585)
  
  myBiomodEM_7100_126 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_7100_126,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_7100_126",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_7100_126/proj_SSP_7100_126_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1,paste0("proj_SSP_7100_126/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_7100_126/proj_SSP_7100_126_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1,paste0("proj_SSP_7100_126/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_7100_126)
  
  
  myBiomodEM_7100_370 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_7100_370,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_7100_370",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_7100_370/proj_SSP_7100_370_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1, paste0("proj_SSP_7100_370/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_7100_370/proj_SSP_7100_370_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1, paste0("proj_SSP_7100_370/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_7100_370)
  
  
  myBiomodEM_7100_585 <-BIOMOD_EnsembleForecasting( bm.em  = myBiomodEM,
                                                    new.env = bio_7100_585,
                                                    models.chosen = 'all',
                                                    proj.name = "SSP_7100_585",
                                                    metric.binary = "TSS",
                                                    output.format = '.grd')
  
  FutureProj <- raster::stack(file.path(outputFolder1, paste("proj_SSP_7100_585/proj_SSP_7100_585_",outputFolder1,"_ensemble_TSSbin.grd", sep="")))
  raster::writeRaster(FutureProj,
              file.path(outputFolder1, paste0("proj_SSP_7100_585/", outputFolder, "_TSSbin.tif")),
              overwrite= TRUE)
  FutureProj_c <- raster::stack(file.path(outputFolder1, paste("proj_SSP_7100_585/proj_SSP_7100_585_",outputFolder1,"_ensemble.grd", sep="")))
  raster::writeRaster(FutureProj_c,
              file.path(outputFolder1, paste0("proj_SSP_7100_585/", outputFolder, "_c.tif")),
              overwrite= TRUE)
  
  rm(FutureProj, FutureProj_c, myBiomodEM_7100_585)
}

