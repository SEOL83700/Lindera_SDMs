packages <- c(
  "viridis",
  "raster",
  "nhSDM",
  "terra",
  "tidyverse",
  "ggplot2",
  "sf",
  "magrittr",
  "khroma",
  "wesanderson",
  "dplyr"
)

install.packages(setdiff(packages, installed.packages()[,"Package"]))

lapply(packages, library, character.only = TRUE)

library(viridis)
library(raster)
library(nhSDM)
library(terra)
library(tidyverse)
library(ggplot2)
library(sf)
library(magrittr)
library(khroma)
library(wesanderson)
library(dplyr)


setwd("D:/Envdata/")

csv.path <- "D:/Envdata/Backup_species/"
original.path <- "D:/Envdata/Backup_process/"
biomod.path <- "D:/Envdata/Biomod/"

# Make dir : Operate these code for only first time ----

dir.create(paste0(biomod.path, "Sum")) 
dir.create(paste0(biomod.path, "Stats"))
dir.create(paste0(biomod.path, "Stack"))
dir.create(paste0(biomod.path, "Sum_normal"))
dir.create(paste0(biomod.path, "DSR"))
dir.create(paste0(biomod.path, "Plot"))
dir.create(paste0(biomod.path, "Plot/Total"))
dir.create(paste0(biomod.path, "Plot/Species"))
dir.create(paste0(biomod.path, "Plot/Species/Latitude"))
dir.create(paste0(biomod.path, "Plot/Species/Longitude"))
dir.create(paste0(biomod.path, "Plot/Species/Elevation"))
dir.create(paste0(biomod.path, "Plot/Species/Binary"))
dir.create(paste0(biomod.path, "Plot/Species/RangeDifference_km2"))
dir.create(paste0(biomod.path, "Plot/Species/RangeDifference_Percent"))
dir.create(paste0(biomod.path, "Plot/DSR"))
dir.create(paste0(biomod.path, "Stack/csv"))

codelist <- list.files(csv.path)

codes <- NULL

for (i in (1:length(codelist))){
  AA <- codelist[i]
  BB <- substr(AA, 1, nchar(AA)-4)
  codes <- c(codes, BB)
}

codes_dot <- gsub("_", ".", codes)
codes_empty <- gsub("_", " ", codes)
Scenario <- c("SSP_Current", "SSP_1140_126", "SSP_1140_370", "SSP_1140_585","SSP_4170_126", "SSP_4170_370", "SSP_4170_585",
             "SSP_7100_126", "SSP_7100_370", "SSP_7100_585")

# Rearrangement of result file ----

 for(i in 1:length(Scenario)){
   dir.create(paste0(biomod.path, Scenario[i]))
 }
 
 for(k in 1:length(codes_dot)){
   bin.path <- list.files(paste0(original.path, codes_dot[k], "/proj_", codes[k]), full.names = T, pattern = "*_TSSbin.tif")
   file.copy(bin.path, to = paste0(biomod.path, Scenario[1],"/", codes_dot[k], "_bin.tif"))
   for(i in 2:length(Scenario)){
     bin.path <- sort(list.files(paste0(original.path, codes_dot[k], "/proj_", Scenario[i]), full.names = T, pattern = "*_TSSbin.tif"))
     file.copy(bin.path, to = paste0(biomod.path, Scenario[i],"/", codes_dot[k], "_bin.tif"), overwrite =  T)    
   }
 }

 bin.path <- list.files(paste0(original.path, codes_dot[1], "/proj_", codes[1]), full.names = T, pattern = "*_TSSbin.tif")
 bin.path <- sort(list.files(paste0("D:/Envdata/Backup_process/", codes_dot[k], "/proj_", Scenario[i]), full.names = T, pattern = "*_TSSbin.tif"))
 str(codes)
 
 
Total <- 3240351 # total Non-na cell count
DeToKm <- 0.0681247 # Area(km2) of Korean Peninsula / Total cell count (for calculating of area(km2))
df1 <- data.frame(stringsAsFactors = T)

length(list)

for(i in 1:length(Scenario)){
  list <- sort(list.files(paste0(biomod.path, Scenario[i]), full.names = T, pattern = "*_bin.tif"))
  temp <- stack(x = list)
  names(temp) <- codes_empty
  assign(paste0("stack_", Scenario[i]), temp)
  writeRaster(temp,filename = paste0(biomod.path, "Stack/", Scenario[i], "_stack"), overwrite= TRUE)
  rs1 <- calc(temp, sum)
  rs2 <- rs1 / length(list)
  writeRaster(rs1, filename = paste0(biomod.path, "Sum/", Scenario[i], "_sum"), overwrite= TRUE)
  writeRaster(rs2, filename = paste0(biomod.path, "Sum_normal/", Scenario[i], "_sum_normal"), overwrite= TRUE)
  df <- data.frame(stringsAsFactors = T)
  for(j in 1:length(codes)){
    sp <- temp[[j]]
    pre <- freq(sp, value=1)
    abs <- freq(sp, value=0)
    Total_km <- Total * DeToKm
    pre_km <- pre * DeToKm
    abs_km <- abs * DeToKm
    tmp_tb <- cbind(codes_empty[j], Total, pre, abs, Total_km, pre_km, abs_km, Scenario[i])
    df <- rbind(df, tmp_tb)
  }
  names(df) <- c("sp", "total", "pre", "abs", "total_km2", "pre_km2", "abs_km2", "Scenario")
  assign(paste0("df_", Scenario[i]), df)
  df1 <- rbind(df1, df)
  assign("df_TotalScenario", df1)
}
write.csv(x = df_TotalScenario, file = paste0(biomod.path, "Stats/TotalAreaPerSpecies.csv"))

MultiVar <- list.files(path = paste0(biomod.path, "Sum_normal/"), pattern = "grd")

for (i in (1:length(MultiVar))){
  AA <- MultiVar[i]
  BB <- substr(AA, 1, nchar(AA)-4)
  temp <- stack(paste0(biomod.path, "Sum_normal/", MultiVar[i])) # for multiband image
  names(temp) <- BB
  assign(BB, temp)
}

plot(SSP_1140_126_sum_normal)

# Calculating of range change index(RCI) ----

Sn <- c("SSP_1140_126", "SSP_1140_370", "SSP_1140_585","SSP_4170_126", "SSP_4170_370", "SSP_4170_585",
        "SSP_7100_126", "SSP_7100_370", "SSP_7100_585", "SSP_Current")

ls_stack <- sort(list.files(paste0(biomod.path, "Stack"), full.names = T, pattern = "*.grd"))
Current <- stack(ls_stack[length(ls_stack)])
df1 <- data.frame(stringsAsFactors = T)
for(i in 1:(length(ls_stack)-1)){
  Future <- stack(ls_stack[i])
  df <- data.frame(stringsAsFactors = T)
  for(j in 1:length(codes)){
    tmp_c <- Current[[j]]
    tmp_f <- Future[[j]]
    cur <- freq(tmp_c, value = 1)
    tmp_chn <- tmp_f - tmp_c
    tmp_chn_col <- freq(tmp_chn, value = 1)
    tmp_chn_ext <- freq(tmp_chn, value = -1)
    rci <- round((tmp_chn_col - tmp_chn_ext) / cur, 1)
    rci_percent <- round(rci * 100, 1)
    tmp_tb <- cbind(codes_empty[j], cur, tmp_chn_col, tmp_chn_ext, rci, rci_percent, Sn[i]) 
    df <- rbind(df, tmp_tb)
  }
  names(df) <- c("sp", "pres", "col", "ext", "rci", "rci_percent", "Scenario")
  assign(paste0("df_", Sn[i]), df)
  df1 <- rbind(df1, df)
  assign("df_rci_TotalScenario", df1)
}
write.csv(x = df_rci_TotalScenario, file = paste0(biomod.path, "Stats/RciPerSpecies.csv"))

# Calculation of SR difference----
Sn <- c("SSP_1140_126", "SSP_1140_370", "SSP_1140_585","SSP_4170_126", "SSP_4170_370", "SSP_4170_585",
        "SSP_7100_126", "SSP_7100_370", "SSP_7100_585", "SSP_Current")

ls_sum_normal <- sort(list.files(paste0(biomod.path, "Sum_normal"), full.names = T, pattern = "*.grd"))
ls_sum_normal

Current <- stack(ls_sum_normal[length(ls_sum_normal)])

for(i in 1:(length(ls_sum_normal)-1)){
  Future <- stack(ls_sum_normal[i])
  DSR <- Future - Current
  assign(paste0("DSR_", Sn[i]), DSR)
  writeRaster(DSR, filename = paste0(biomod.path, "DSR/", Sn[i], "_DSR"), overwrite= TRUE)
}

MultiVar <- list.files(path = paste0(biomod.path, "DSR/"), pattern = "grd")

for (i in (1:length(MultiVar))){
  AA <- MultiVar[i]
  BB <- substr(AA, 1, nchar(AA)-4)
  temp <- stack(paste0(biomod.path, "DSR/", MultiVar[i])) # for multiband image
  names(temp) <- BB
  assign(BB, temp)
}

# plot(DSR, col = inferno(256), zlim = c(-1, 1))


# Extract of environment range for each species----

Scenario <- c("SSP_Current", "SSP_1140_126", "SSP_1140_370", "SSP_1140_585","SSP_4170_126", "SSP_4170_370", "SSP_4170_585",
              "SSP_7100_126", "SSP_7100_370", "SSP_7100_585")

# Calling environment variable
setwd("D:/Envdata/Multicorr_var")
MultiVar <- list.files(path = ".", pattern = "grd")

TempName <- c("bio13", "bio14", "bio15", "bio2", "bio3", "bio4", "gsp", "gst", "kg0", "kg2", "kg3",
              "kg4", "kg5", "npp", "bdod", "cfvo", "ocd", "phh2o", "sand", "silt", "soc", "elev", "rough", "srad", "wavg")

for (i in (1:length(MultiVar))){
  AA <- MultiVar[i]
  BB <- substr(AA, 1, nchar(AA)-4)
  temp <- stack(MultiVar[i]) # for multiband image
  names(temp) <- TempName
  assign(BB, temp)
}
setwd("D:/Envdata/")  

Scntemp <- c("bio_present", "bio_1140_126", "bio_1140_370", "bio_1140_585", "bio_4170_126", "bio_4170_370", "bio_4170_585", "bio_7100_126", "bio_7100_370", "bio_7100_585")

df1 <- data.frame(stringsAsFactors = T)

for(i in 1:length(Scenario)){
  list <- sort(list.files(paste0(biomod.path, "/", Scenario[i]), full.names = T, pattern = "*_bin.tif"))
  temp <- stack(x = list)
  names(temp) <- codes_empty
  env <- get(Scntemp[i])
  df <- data.frame(stringsAsFactors = T)
  for(j in 1:length(codes)){
    ssp <- temp[[j]]
    tpts <- rasterToPoints(ssp, spatial = T)
    tpts@data <- data.frame(tpts@data, long = coordinates(tpts)[,1], lat = coordinates(tpts)[,2])
    sp <- as.data.frame(tpts@data)
    sp <- sp[sp[,1] == 1,c(2,3)]
    env1 <- raster::extract(env, sp) 
    total <- cbind(sp, env1)
    aa <- total %>% select(long, lat, elev) %>% 
      summarise(long_mean = mean(long, na.rm = T), long_sd = sd(long, na.rm = T), long_low = min(long, na.rm = T), long_upp = max(long, na.rm = T), long_se=long_sd/fnxmsqrt(179)
                lat_mean = mean(lat, na.rm = T), lat_sd = sd(lat, na.rm = T), lat_low = min(lat, na.rm = T), lat_upp = max(lat, na.rm = T),
                ele_mean = mean(elev, na.rm = T), ele_sd = sd(elev, na.rm = T), ele_low = min(elev, na.rm = T), ele_upp = max(elev, na.rm = T))
    aa1 <- cbind(codes_empty[j], Scenario[i], aa)
    names(aa1)[1:2] <- c("sp", "Scenario")
    df <- rbind(df, aa1)
  }
  assign(paste0("env_", Scenario[i]), df)
  df1 <- rbind(df1, df)
  assign("df_TotalScenario", df1)
}
write.csv(x = df_TotalScenario, file = paste0(biomod.path, "Stats/EnvironmentVariable.csv"))

# Extraction of descriptive statistics of environment for total species---- 
setwd("D:/Envdata/Biomod/Stats")
hh <- read.csv("EnvironmentVariable.csv")

colnames(df1)
df2 <- hh %>% group_by(Scenario) %>%
  summarise(long_mean = mean(long_mean, na.rm = T), long_sd = mean(long_sd, na.rm = T), lat_mean = mean(lat_mean, na.rm = T), lat_sd = mean(lat_sd, na.rm = T),
            ele_mean = mean(ele_mean, na.rm = T), ele_sd = mean(ele_sd, na.rm = T))


df3 <- df2[1:3, ] %>% mutate(Period = "2011-2040")
df4 <- df2[4:6, ] %>% mutate(Period = "2041-2070")
df5 <- df2[7:9, ] %>% mutate(Period = "2071-2100")
df6 <- df2[10, ] %>% mutate(Period = "Current")
df7 <- rbind(df3,df4,df5,df6)
df2 <- df7

df3 <- df2[c(1,4,7), ] %>% mutate(SSP = "SSP1-2.6")
df4 <- df2[c(2,5,8), ] %>% mutate(SSP = "SSP3-7.0")
df5 <- df2[c(3,6,9), ] %>% mutate(SSP = "SSP5-8.5")
df6 <- df2[10, ] %>% mutate(SSP = "SSP1-2.6")
df7 <- df2[10, ] %>% mutate(SSP = "SSP3-7.0")
df8 <- df2[10, ] %>% mutate(SSP = "SSP5-8.5")
df7 <- rbind(df3,df4,df5,df6,df7,df8)
df2 <- df7

rm(df3, df4, df5, df6, df7, df8)


nCol <- 4
myCol <- rev(viridis(n = nCol, option = "B"))
myCol <- inferno(n = nCol)
myCol <- c("#287D8EFF", "#F8A815", "coral3")

df2$Period <- factor(df2$Period, levels = c("Current", "2011-2040", "2041-2070", "2071-2100"))

# Plotting for total species----
b1 <- ggplot(data = df2) +
  geom_point(aes(x = Period, y = long_mean, color = SSP, size = 1.5, alpha = 0.5)) + 
  geom_line(aes(x = Period, y = long_mean, group = SSP, color = SSP)) +
  geom_errorbar(aes(x = Period, ymin = long_mean - long_sd, ymax = long_mean + long_sd, group = SSP, color = SSP)) + # errorbar(SD)
  scale_color_manual(values = myCol) + 
  xlab("Period") + ylab("Longitude(decimal degree)") + 
  labs(color = "Scenario") +
  # theme(legend.box.background = element_rect(), legend.box.margin = margin(2, 2, 2, 2), legend.text = element_text(size = 8, face = "bold")) +
  guides(alpha = FALSE, size = FALSE) +
  theme_classic()

b1
ggsave(filename = "Longitude change.jpeg", path = paste0(biomod.path, "Plot/Total"), dpi = 300, units = "cm", width = 12, height = 10)
unlink("Longitude change.jpeg")
while (!is.null(dev.list()))  dev.off()

b1 <- ggplot(data = df2) +
  geom_point(aes(x = Period, y = lat_mean, color = SSP, size = 1.5, alpha = 0.5)) + 
  geom_line(aes(x = Period, y = lat_mean, group = SSP, color = SSP)) +
  geom_errorbar(aes(x = Period, ymin = lat_mean - lat_sd, ymax = lat_mean + lat_sd, group = SSP, color = SSP)) + # errorbar(SD)
  scale_color_manual(values = myCol) + 
  xlab("Period") + ylab("Latitude(decimal degree)") + 
  labs(color = "Scenario") +
  # theme(legend.box.background = element_rect(), legend.box.margin = margin(2, 2, 2, 2), legend.text = element_text(size = 8, face = "bold")) +
  guides(alpha = FALSE, size = FALSE) +
  theme_classic()
b1

ggsave(filename = "Latitude change.jpeg", path = paste0(biomod.path, "Plot/Total"), dpi = 300, units = "cm", width = 12, height = 10)
unlink("Latitude change.jpeg")
while (!is.null(dev.list()))  dev.off()

b1 <- ggplot(data = df2) +
  geom_point(aes(x = Period, y = ele_mean, color = SSP, size = 1.5, alpha = 0.5)) + 
  geom_line(aes(x = Period, y = ele_mean, group = SSP, color = SSP)) +
  geom_errorbar(aes(x = Period, ymin = ele_mean - ele_sd, ymax = ele_mean + ele_sd, group = SSP, color = SSP)) + # errorbar(SD)
  scale_color_manual(values = myCol) + 
  xlab("Period") + ylab("Elevation(m)") + 
  labs(color = "Scenario") +
  # theme(legend.box.background = element_rect(), legend.box.margin = margin(2, 2, 2, 2), legend.text = element_text(size = 8, face = "bold")) +
  guides(alpha = FALSE, size = FALSE) +
  theme_classic()
b1

ggsave(filename = "Elevation change.jpeg", path = paste0(biomod.path, "Plot/Total"), dpi = 300, units = "cm", width = 12, height = 10)
unlink("Elevation change.jpeg")
while (!is.null(dev.list()))  dev.off()

# Plotting for each species----

aa <- unique(df_TotalScenario$Scenario)  
aa
df3 <- df_TotalScenario %>% filter(Scenario == "SSP_Current") %>% mutate(Period = "Period 1")
df4 <- df_TotalScenario %>% filter(Scenario %in% aa[2:4]) %>% mutate(Period = "Period 2")
df5 <- df_TotalScenario %>% filter(Scenario %in% aa[5:7]) %>% mutate(Period = "Period 3")
df6 <- df_TotalScenario %>% filter(Scenario %in% aa[8:10]) %>% mutate(Period = "Period 4")

df_sp <- rbind(df3, df4, df5, df6)

aa[c(2,5,8)]
df3 <- df_sp %>% filter(Scenario == "SSP_Current") %>% mutate(SSP = "SSP126")
df4 <- df_sp %>% filter(Scenario == "SSP_Current") %>% mutate(SSP = "SSP370")
df5 <- df_sp %>% filter(Scenario == "SSP_Current") %>% mutate(SSP = "SSP585")
df6 <- df_sp %>% filter(Scenario %in% aa[c(2,5,8)]) %>% mutate(SSP = "SSP126")
df7 <- df_sp %>% filter(Scenario %in% aa[c(3,6,9)]) %>% mutate(SSP = "SSP370")
df8 <- df_sp %>% filter(Scenario %in% aa[c(4,7,10)]) %>% mutate(SSP = "SSP585")

df_sp <- rbind(df3, df4, df5, df6, df7, df8)

rm(df3, df4, df5, df6, df7, df8)

for(i in 1:length(codes_empty)){
  df_tmp <- df_sp %>% filter(sp == codes_empty[i])
  b1 <- ggplot(data = df_tmp) +
    geom_point(aes(x = Period, y = long_mean, color = SSP, size = 1.5, alpha = 0.5)) + 
    geom_line(aes(x = Period, y = long_mean, group = SSP, color = SSP)) +
    geom_errorbar(aes(x = Period, ymin = long_mean - long_sd, ymax = long_mean + long_sd, group = SSP, color = SSP)) + # errorbar(SD)
    scale_color_manual(values = myCol) + 
    xlab("Period") + ylab("Longitude(decimal degree)") + 
    labs(color = "Scenario") +
    # theme(legend.box.background = element_rect(), legend.box.margin = margin(2, 2, 2, 2), legend.text = element_text(size = 8, face = "bold")) +
    guides(alpha = FALSE, size = FALSE) +
    theme_classic()
  
  ggsave(filename = paste0(codes_empty[i], "_Longitude change.jpeg"), path = paste0(biomod.path, "Plot/Species/Longitude"), dpi = 300, units = "cm", width = 10, height = 10)
  unlink(paste0(codes_empty[i], "_Longitude change.jpeg"))
  while (!is.null(dev.list()))  dev.off()
  
  b1 <- ggplot(data = df_tmp) +
    geom_point(aes(x = Period, y = lat_mean, color = SSP, size = 1.5, alpha = 0.5)) + 
    geom_line(aes(x = Period, y = lat_mean, group = SSP, color = SSP)) +
    geom_errorbar(aes(x = Period, ymin = lat_mean - lat_sd, ymax = lat_mean + lat_sd, group = SSP, color = SSP)) + # errorbar(SD)
    scale_color_manual(values = myCol) + 
    xlab("Period") + ylab("Latitude(decimal degree)") + 
    labs(color = "Scenario") +
    # theme(legend.box.background = element_rect(), legend.box.margin = margin(2, 2, 2, 2), legend.text = element_text(size = 8, face = "bold")) +
    guides(alpha = FALSE, size = FALSE) +
    theme_classic()
  
  ggsave(filename = paste0(codes_empty[i], "_Latitude change.jpeg"), path = paste0(biomod.path, "Plot/Species/Latitude"), dpi = 300, units = "cm", width = 10, height = 10)
  unlink(paste0(codes_empty[i], "_Latitude change.jpeg"))
  while (!is.null(dev.list()))  dev.off()
  
  b1 <- ggplot(data = df_tmp) +
    geom_point(aes(x = Period, y = ele_mean, color = SSP, size = 1.5, alpha = 0.5)) + 
    geom_line(aes(x = Period, y = ele_mean, group = SSP, color = SSP)) +
    geom_errorbar(aes(x = Period, ymin = ele_mean - ele_sd, ymax = ele_mean + ele_sd, group = SSP, color = SSP)) + # errorbar(SD)
    scale_color_manual(values = myCol) + 
    xlab("Period") + ylab("Elevation(m)") + 
    labs(color = "Scenario") +
    # theme(legend.box.background = element_rect(), legend.box.margin = margin(2, 2, 2, 2), legend.text = element_text(size = 8, face = "bold")) +
    guides(alpha = FALSE, size = FALSE) +
    theme_classic()
  
  ggsave(filename = paste0(codes_empty[i], "_Elevation change.jpeg"), path = paste0(biomod.path, "Plot/Species/Elevation"), dpi = 300, units = "cm", width = 10, height = 10)
  unlink(paste0(codes_empty[i], "_Elevation change.jpeg"))
  while (!is.null(dev.list()))  dev.off()
  
}


# Plotting of binary map for each species

ls_dsr <- c("DSR_SSP_1140_126","DSR_SSP_1140_370","DSR_SSP_1140_585","DSR_SSP_4170_126","DSR_SSP_4170_370","DSR_SSP_4170_585","DSR_SSP_7100_126","DSR_SSP_7100_370","DSR_SSP_7100_585")
dsrcode <- c("SSP126: 11-40","SSP370: 11-40","SSP585: 11-40","SSP126: 41-70","SSP370: 41-70","SSP585: 41-70","SSP126: 71-00","SSP370: 71-00","SSP585: 71-00")

dsr_df <- data.frame(stringsAsFactors = T)
for(i in 1:length(ls_dsr)){
  dsr <- get(ls_dsr[i])
  df_dsr <- as.data.frame(dsr, xy = T)
  df_na <- na.omit(df_dsr)
  df_na$Scenario <- dsrcode[i]
  dsr_df <- rbind(dsr_df, df_na)
}

fill_min <- min(dsr_df$layer)
fill_max <- max(dsr_df$layer)
plot.dsr <- ggplot(data = dsr_df) +
  geom_raster(aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(limits = c(fill_min, fill_max)) +
  # scale_fill_viridis(limits = c(-1, 1), option = "inferno") +# color option : "inferno", "magma", "turbo"  
  theme_void() +
  theme(legend.position = "bottom") +
  coord_equal() +
  labs(fill = "Difference of Species Range") +
  facet_wrap(Scenario~.)

ggsave(filename = "Difference of Species Range.jpeg", path = paste0(biomod.path, "Plot/DSR"), dpi = 300, units = "cm", width = 16, height = 25)
unlink("Difference of Species Range.jpeg")
#?str(test.df)

Scenario <- c("SSP_Current", "SSP_1140_126", "SSP_1140_370", "SSP_1140_585","SSP_4170_126", "SSP_4170_370", "SSP_4170_585",
              "SSP_7100_126", "SSP_7100_370", "SSP_7100_585")
dsrcode_c <- c("Current", "SSP126: 11-40","SSP370: 11-40","SSP585: 11-40","SSP126: 41-70","SSP370: 41-70","SSP585: 41-70","SSP126: 71-00","SSP370: 71-00","SSP585: 71-00")


for(i in 1:length(codes)){
  bin.df1 <- data.frame()
  for(j in 1:length(Scenario)){
    env <- get(paste0("stack_", Scenario[j]))
    env <- env[[i]]
    env.df <- as.data.frame(env, xy = T)
    names(env.df) <- c("x", "y", "value")
    na.df <- na.omit(env.df)
    na.df$sp <- codes_empty[i]
    na.df$Scenario <- dsrcode_c[j]
    na.df <- na.df %>% mutate(preab = ifelse(value == 1, "Presence", "Absence"))
    bin.df1 <- rbind(bin.df1, na.df)
  }
  for(k in 1:length(dsrcode_c)){
    df.c <- bin.df1 %>% filter(Scenario == dsrcode_c[k])
    plot.df <- ggplot(data = df.c) +
      geom_raster(aes(x = x, y = y, fill = preab, alpha = 0.1)) +
      theme_void() +
      theme(legend.position = "bottom") +
      coord_equal() +
      labs(fill = "Distribution") +
      scale_fill_manual(values = rev(wes_palette(n = 2, name = "Darjeeling1", type = "discrete")))
    ggsave(filename = paste0(codes_empty[i], "_", Scenario[k], "_Distribution.jpeg"), path = paste0(biomod.path, "Plot/Species/Binary"), dpi = 300, units = "cm", width = 16, height = 25)
    unlink(paste0(codes_empty[i], "_", Scenario[k], "_Distribution.jpeg"))
  }
}


# Pre-processing to make graphs of range difference for each species ----
#names(df_range)
df_rci <- read.csv(paste0(biomod.path, "Stats/RciPerSpecies_re.csv"), stringsAsFactors = T)
df_range <- read.csv(paste0(biomod.path, "Stats/TotalAreaPerSpecies.csv"), stringsAsFactors = T)
df_range <- df_range %>% mutate(preratio = pre_km2 / total_km2 * 100)
tmp1 <- df_range %>% select(sp, pre_km2, Scenario) %>% pivot_wider(names_from = "Scenario", values_from = "pre_km2")
tmp1
tmp2 <- tmp1 %>% select(sp, SSP_Current, SSP_1140_126) %>% mutate(diff = SSP_1140_126 - SSP_Current, Scenario = "SSP_1140_126")
tmp3 <- tmp1 %>% select(sp, SSP_1140_126, SSP_4170_126) %>% mutate(diff = SSP_4170_126 - SSP_1140_126, Scenario = "SSP_4170_126")

tmp_126 <- tmp1 %>% select(sp, SSP_Current, SSP_1140_126, SSP_4170_126, SSP_7100_126) %>% mutate(SSP = "SSP126")
tmp_370 <- tmp1 %>% select(sp, SSP_Current, SSP_1140_370, SSP_4170_370, SSP_7100_370) %>% mutate(SSP = "SSP370")
tmp_585 <- tmp1 %>% select(sp, SSP_Current, SSP_1140_585, SSP_4170_585, SSP_7100_585) %>% mutate(SSP = "SSP585")


tmp_test <- tmp_126 %>% select(sp, SSP_Current,  SSP_1140_126, SSP) %>% mutate(diff = SSP_1140_126 - SSP_Current, Period = "P1-P2")
tmp_test1 <- tmp_126 %>% select(sp, SSP_1140_126, SSP_4170_126, SSP) %>% mutate(diff = SSP_4170_126 - SSP_1140_126, Period = "P2-P3")
tmp_test2 <- tmp_126 %>% select(sp, SSP_4170_126, SSP_7100_126, SSP) %>% mutate(diff = SSP_7100_126 - SSP_4170_126, Period = "P3-P4")
tmp_test_diff <- tmp_test %>% select(sp, diff, Period, SSP)
tmp_test1_diff <- tmp_test1 %>% select(sp, diff, Period, SSP)
tmp_test2_diff <- tmp_test2 %>% select(sp, diff, Period, SSP)

tmp_diff_126 <- rbind(tmp_test_diff, tmp_test1_diff, tmp_test2_diff)

tmp_test <- tmp_370 %>% select(sp, SSP_Current,  SSP_1140_370, SSP) %>% mutate(diff = SSP_1140_370 - SSP_Current, Period = "P1-P2")
tmp_test1 <- tmp_370 %>% select(sp, SSP_1140_370, SSP_4170_370, SSP) %>% mutate(diff = SSP_4170_370 - SSP_1140_370, Period = "P2-P3")
tmp_test2 <- tmp_370 %>% select(sp, SSP_4170_370, SSP_7100_370, SSP) %>% mutate(diff = SSP_7100_370 - SSP_4170_370, Period = "P3-P4")
tmp_test_diff <- tmp_test %>% select(sp, diff, Period, SSP)
tmp_test1_diff <- tmp_test1 %>% select(sp, diff, Period, SSP)
tmp_test2_diff <- tmp_test2 %>% select(sp, diff, Period, SSP)

tmp_diff_370 <- rbind(tmp_test_diff, tmp_test1_diff, tmp_test2_diff)

tmp_test <- tmp_585 %>% select(sp, SSP_Current,  SSP_1140_585, SSP) %>% mutate(diff = SSP_1140_585 - SSP_Current, Period = "P1-P2")
tmp_test1 <- tmp_585 %>% select(sp, SSP_1140_585, SSP_4170_585, SSP) %>% mutate(diff = SSP_4170_585 - SSP_1140_585, Period = "P2-P3")
tmp_test2 <- tmp_585 %>% select(sp, SSP_4170_585, SSP_7100_585, SSP) %>% mutate(diff = SSP_7100_585 - SSP_4170_585, Period = "P3-P4")
tmp_test_diff <- tmp_test %>% select(sp, diff, Period, SSP)
tmp_test1_diff <- tmp_test1 %>% select(sp, diff, Period, SSP)
tmp_test2_diff <- tmp_test2 %>% select(sp, diff, Period, SSP)

tmp_diff_585 <- rbind(tmp_test_diff, tmp_test1_diff, tmp_test2_diff)

tmp_df_range <- rbind(tmp_diff_126, tmp_diff_370, tmp_diff_585)

rm(tmp1, tmp_126, tmp_370, tmp_585, tmp_test, tmp_test1, tmp_test2, tmp_test_diff, tmp_test1_diff, tmp_test2_diff, tmp_diff_126, tmp_diff_370, tmp_diff_585)

tmp1 <- df_range %>% select(sp, preratio, Scenario) %>% pivot_wider(names_from = "Scenario", values_from = "preratio")

tmp_126 <- tmp1 %>% select(sp, SSP_Current, SSP_1140_126, SSP_4170_126, SSP_7100_126) %>% mutate(SSP = "SSP126")
tmp_370 <- tmp1 %>% select(sp, SSP_Current, SSP_1140_370, SSP_4170_370, SSP_7100_370) %>% mutate(SSP = "SSP370")
tmp_585 <- tmp1 %>% select(sp, SSP_Current, SSP_1140_585, SSP_4170_585, SSP_7100_585) %>% mutate(SSP = "SSP585")


tmp_test <- tmp_126 %>% select(sp, SSP_Current,  SSP_1140_126, SSP) %>% mutate(diff = SSP_1140_126 - SSP_Current, Period = "P1-P2")
tmp_test1 <- tmp_126 %>% select(sp, SSP_1140_126, SSP_4170_126, SSP) %>% mutate(diff = SSP_4170_126 - SSP_1140_126, Period = "P2-P3")
tmp_test2 <- tmp_126 %>% select(sp, SSP_4170_126, SSP_7100_126, SSP) %>% mutate(diff = SSP_7100_126 - SSP_4170_126, Period = "P3-P4")
tmp_test_diff <- tmp_test %>% select(sp, diff, Period, SSP)
tmp_test1_diff <- tmp_test1 %>% select(sp, diff, Period, SSP)
tmp_test2_diff <- tmp_test2 %>% select(sp, diff, Period, SSP)

tmp_diff_126 <- rbind(tmp_test_diff, tmp_test1_diff, tmp_test2_diff)

tmp_test <- tmp_370 %>% select(sp, SSP_Current,  SSP_1140_370, SSP) %>% mutate(diff = SSP_1140_370 - SSP_Current, Period = "P1-P2")
tmp_test1 <- tmp_370 %>% select(sp, SSP_1140_370, SSP_4170_370, SSP) %>% mutate(diff = SSP_4170_370 - SSP_1140_370, Period = "P2-P3")
tmp_test2 <- tmp_370 %>% select(sp, SSP_4170_370, SSP_7100_370, SSP) %>% mutate(diff = SSP_7100_370 - SSP_4170_370, Period = "P3-P4")
tmp_test_diff <- tmp_test %>% select(sp, diff, Period, SSP)
tmp_test1_diff <- tmp_test1 %>% select(sp, diff, Period, SSP)
tmp_test2_diff <- tmp_test2 %>% select(sp, diff, Period, SSP)

tmp_diff_370 <- rbind(tmp_test_diff, tmp_test1_diff, tmp_test2_diff)

tmp_test <- tmp_585 %>% select(sp, SSP_Current,  SSP_1140_585, SSP) %>% mutate(diff = SSP_1140_585 - SSP_Current, Period = "P1-P2")
tmp_test1 <- tmp_585 %>% select(sp, SSP_1140_585, SSP_4170_585, SSP) %>% mutate(diff = SSP_4170_585 - SSP_1140_585, Period = "P2-P3")
tmp_test2 <- tmp_585 %>% select(sp, SSP_4170_585, SSP_7100_585, SSP) %>% mutate(diff = SSP_7100_585 - SSP_4170_585, Period = "P3-P4")
tmp_test_diff <- tmp_test %>% select(sp, diff, Period, SSP)
tmp_test1_diff <- tmp_test1 %>% select(sp, diff, Period, SSP)
tmp_test2_diff <- tmp_test2 %>% select(sp, diff, Period, SSP)

tmp_diff_585 <- rbind(tmp_test_diff, tmp_test1_diff, tmp_test2_diff)

tmp_df_ratio <- rbind(tmp_diff_126, tmp_diff_370, tmp_diff_585)

rm(tmp1, tmp_126, tmp_370, tmp_585, tmp_test, tmp_test1, tmp_test2, tmp_test_diff, tmp_test1_diff, tmp_test2_diff, tmp_diff_126, tmp_diff_370, tmp_diff_585)

names(tmp_df_range) <- c("sp", "range_diff", "Period", "SSP")
names(tmp_df_ratio) <- c("sp", "ratio_diff", "Period", "SSP")
df_diff <- cbind(tmp_df_range, tmp_df_ratio)
names(df_diff)
df_diff <- df_diff[ ,c(1,2,6,7,8)]

df_range_summ <- df_diff %>%
  group_by(Period, SSP) %>%
  summarise(
    range_diff_mean = mean(range_diff, na.rm = TRUE),
    range_diff_sd = sd(range_diff, na.rm = TRUE),
    .groups = 'drop'
  )

df_ratio_summ <- df_diff %>%
  group_by(Period, SSP) %>%
  summarise(
    ratio_diff_mean = mean(ratio_diff, na.rm = TRUE),
    ratio_diff_sd = sd(ratio_diff, na.rm = TRUE),
    .groups = 'drop'
  )

# Join the two summary dataframes
df_summ <- left_join(df_range_summ, df_ratio_summ, by = c("Period", "SSP"))
print(df_summ)

names(df_summ) <- c("Period", "SSP", "range_diff_mean", "range_diff_sd", "ratio_diff_mean", "ratio_diff_sd")
write.csv(x = df_summ, file = paste0(biomod.path, "Stats/Total_DifferencePerPeriod.csv"))

b1 <- ggplot(data = df_summ) +
  geom_point(aes(x = Period, y = range_diff_mean, color = SSP, size = 1.5, alpha = 0.5)) + 
  geom_line(aes(x = Period, y = range_diff_mean, group = SSP, color = SSP)) +
  geom_errorbar(aes(x = Period, ymin = range_diff_mean - range_diff_sd, ymax = range_diff_mean + range_diff_sd, group = SSP, color = SSP)) + # errorbar(SD)
  scale_color_manual(
    values = myCol,
    breaks = c("SSP126", "SSP370", "SSP585"),  
    labels = c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")  
  ) + 
  xlab("Period") + ylab("range difference (square km)") + 
  labs(color = "Scenario") +
  guides(alpha = FALSE, size = FALSE) +
  theme_classic()

b1 <- b1 + scale_x_discrete(labels = c("P1-P2" = "(2011-2040)-Current", "P2-P3" = "(2041-2070)-(2011-2040)", "P3-P4" = "(2071-2100)-(2041-2070)"))


ggsave(filename = "Range change_km2.jpeg", path = paste0(biomod.path, "Plot/Total"), dpi = 300, units = "cm", width = 20, height = 10)
unlink("Range change_km2.jpeg")
while (!is.null(dev.list()))  dev.off()

b1 <- ggplot(data = df_summ) +
  geom_point(aes(x = Period, y = ratio_diff_mean, color = SSP, size = 1.5, alpha = 0.5)) + 
  geom_line(aes(x = Period, y = ratio_diff_mean, group = SSP, color = SSP)) +
  geom_errorbar(aes(x = Period, ymin = ratio_diff_mean - ratio_diff_sd, ymax = ratio_diff_mean + ratio_diff_sd, group = SSP, color = SSP)) + # errorbar(SD)
  scale_color_manual(
    values = myCol,
    breaks = c("SSP126", "SSP370", "SSP585"),  
    labels = c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5")  
  ) + 
  xlab("Period") + ylab("range difference(%)") + 
  labs(color = "Scenario") +
  guides(alpha = FALSE, size = FALSE) +
  theme_classic()

b1 <- b1 + scale_x_discrete(labels = c("P1-P2" = "(2011-2040)-Current", "P2-P3" = "(2041-2070)-(2011-2040)", "P3-P4" = "(2071-2100)-(2041-2070)"))

b1

ggsave(filename = "Range change_Percent.jpeg", path = paste0(biomod.path, "Plot/Total"), dpi = 300, units = "cm", width = 20, height = 10)
unlink("Range change_Percent.jpeg")
while (!is.null(dev.list()))  dev.off()


# plotting graphs of range difference for each species ----

for(i in 1:length(codes_empty)){
  df_tmp <- df_diff %>% filter(sp == codes_empty[i])
  b1 <- ggplot(data = df_tmp) +
    geom_point(aes(x = Period, y = range_diff, color = SSP, size = 1.5, alpha = 0.5)) + 
    geom_line(aes(x = Period, y = range_diff, group = SSP, color = SSP)) +
    scale_color_manual(values = myCol) + 
    xlab("Period") + ylab("range difference(square km)") + 
    labs(color = "Scenario") +
    guides(alpha = FALSE, size = FALSE) +
    theme_classic()
  
  ggsave(filename = paste0(codes_empty[i], "_Range difference_km2.jpeg"), path = paste0(biomod.path, "Plot/Species/RangeDifference_km2"), dpi = 300, units = "cm", width = 10, height = 10)
  unlink(paste0(codes_empty[i], "_Range difference_km2.jpeg"))
  while (!is.null(dev.list()))  dev.off()
  
  b1 <- ggplot(data = df_tmp) +
    geom_point(aes(x = Period, y = ratio_diff, color = SSP, size = 1.5, alpha = 0.5)) + 
    geom_line(aes(x = Period, y = ratio_diff, group = SSP, color = SSP)) +
    scale_color_manual(values = myCol) + 
    xlab("Period") + ylab("Ratio difference(%)") + 
    labs(color = "Scenario") +
    guides(alpha = FALSE, size = FALSE) +
    theme_classic()
  
  ggsave(filename = paste0(codes_empty[i], "_Range difference_percent.jpeg"), path = paste0(biomod.path, "Plot/Species/RangeDifference_Percent"), dpi = 300, units = "cm", width = 10, height = 10)
  unlink(paste0(codes_empty[i], "_Range difference_percent.jpeg"))
  while (!is.null(dev.list()))  dev.off()
  
}

rm(b1)

df_rci_tmp1 <- df_rci %>% filter(Scenario %in% c("SSP_1140_126", "SSP_1140_370", "SSP_1140_585")) %>% mutate(Period = "P1-P2")
df_rci_tmp2 <- df_rci %>% filter(Scenario %in% c("SSP_4170_126", "SSP_4170_370", "SSP_4170_585")) %>% mutate(Period = "P1-P3")
df_rci_tmp3 <- df_rci %>% filter(Scenario %in% c("SSP_7100_126", "SSP_7100_370", "SSP_7100_585")) %>% mutate(Period = "P1-P4")
df_rci_tmp <- rbind(df_rci_tmp1,df_rci_tmp2,df_rci_tmp3)
df_rci <- df_rci_tmp

df_rci_tmp1 <- df_rci %>% filter(Scenario %in% c("SSP_1140_126", "SSP_4170_126", "SSP_7100_126")) %>% mutate(SSP = "SSP126")
df_rci_tmp2 <- df_rci %>% filter(Scenario %in% c("SSP_1140_370", "SSP_4170_370", "SSP_7100_370")) %>% mutate(SSP = "SSP370")
df_rci_tmp3 <- df_rci %>% filter(Scenario %in% c("SSP_1140_585", "SSP_4170_585", "SSP_7100_585")) %>% mutate(SSP = "SSP585")
df_rci_tmp <- rbind(df_rci_tmp1,df_rci_tmp2,df_rci_tmp3)

df_rci <- df_rci_tmp
rm(df_rci_tmp1, df_rci_tmp2, df_rci_tmp3, df_rci_tmp)
?sort
df_rci <- df_rci %>% arrange(rci_percent, sp)
df_rci$change <- ifelse(df_rci$rci_percent >= 0, "gain", "loss")

df_rci_126 <- df_rci %>% filter(SSP == "SSP126")
df_rci_370 <- df_rci %>% filter(SSP == "SSP370")
df_rci_585 <- df_rci %>% filter(SSP == "SSP585")

b1 <- ggplot(data = df_rci_126) +
  geom_bar(aes(x = reorder(sp, rci_percent, decreasing = T), y = rci_percent, fill = change), stat = "identity", position = "dodge" ) + 
  labs(fill = "Habitat") + 
  xlab("Species") + ylab("Range change index(RCI)") + 
  scale_fill_manual(values = c("aquamarine4", "firebrick2")) + 
  theme_classic() + 
  coord_flip() +
  facet_wrap(Period~.)  

ggsave(filename = "RCI_SSP126.jpeg", path = paste0(biomod.path, "Plot/Total"), dpi = 300, units = "cm", width = 20, height = 20)
unlink("RCI_SSP126.jpeg")
while (!is.null(dev.list()))  dev.off()


b1 <- ggplot(data = df_rci_370) +
  geom_bar(aes(x = reorder(sp, rci_percent, decreasing = T), y = rci_percent, fill = change), stat = "identity", position = "dodge" ) + 
  labs(fill = "Habitat") + 
  xlab("Species") + ylab("Range change index(RCI)") + 
  scale_fill_manual(values = c("aquamarine4", "firebrick2")) + 
  theme_classic() + 
  coord_flip() +
  facet_wrap(Period~.)  

ggsave(filename = "RCI_SSP370.jpeg", path = paste0(biomod.path, "Plot/Total"), dpi = 300, units = "cm", width = 20, height = 20)
unlink("RCI_SSP370.jpeg")
while (!is.null(dev.list()))  dev.off()

b1

b1 <- ggplot(data = df_rci_585) +
  geom_bar(aes(x = reorder(sp, rci_percent, decreasing = T), y = rci_percent, fill = change), stat = "identity", position = "dodge" ) + 
  labs(fill = "Habitat") + 
  xlab("Species") + ylab("Range change index(RCI)") + 
  scale_fill_manual(values = c("aquamarine4", "firebrick2")) +
  theme_classic() + 
  coord_flip() +
  facet_wrap(Period~.)  


ggsave(filename = "RCI_SSP585.jpeg", path = paste0(biomod.path, "Plot/Total"), dpi = 300, units = "cm", width = 20, height = 20)
unlink("RCI_SSP585.jpeg")
while (!is.null(dev.list()))  dev.off()

