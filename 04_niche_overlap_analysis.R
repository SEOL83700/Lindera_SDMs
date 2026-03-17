library(readr)
library(dplyr)
library(tidyr)
library(purrr)

setwd("D:/Envdata/Envdata_Lindera/Biomod_Lindera/Stats")

a1 <- read_csv("EnvironmentVariable.csv")
a2 <- read_csv("RciPerSpecies.csv")
a3 <- read_csv("Total_DifferencePerPeriod.csv")
a4 <- read_csv("TotalAreaPerSpecies.csv")

str(a2)

rci <- a2 %>% filter(Scenario %in% c("SSP_7100_126", "SSP_7100_370", "SSP_7100_585")) %>% 
  dplyr::select(sp, rci_percent, Scenario) %>%
  group_by(sp) %>% summarise(mean = mean(rci_percent))

rci1 <- rci %>% mutate(A_A = ifelse(mean < 0, 1, 0)) %>%
  mutate(A_B = ifelse(mean >= 0, 1, 0)) %>%
  mutate(A_A_1 = ifelse(mean < 0 & mean >= -30, 1, 0)) %>%
  mutate(A_A_2 = ifelse(mean < -30 & mean > -100, 1, 0)) %>%
  mutate(A_A_3 = ifelse(mean == -100, 1, 0)) %>%
  mutate(A_B_1 = ifelse(mean < 30 & mean >= 0, 1, 0)) %>%
  mutate(A_B_2 = ifelse(mean < 100 & mean >= 30, 1, 0)) %>%
  mutate(A_B_3 = ifelse(mean >= 100, 1, 0)) %>%
  mutate(Area_group = ifelse(A_A_1 == 1, "Dec.0-30", 
                             ifelse(A_A_2 == 1, "Dec.30-100", 
                                    ifelse(A_A_3 == 1, "Extinction", 
                                           ifelse(A_B_1 == 1, "Inc.0-30", 
                                                  ifelse(A_B_2 == 1, "Inc.30-100", "Inc.100>=")))))) %>%
  rename(Area_mean = "mean") %>%
  mutate(Area = ifelse(A_A == 1, "Loss", ifelse(A_B == 1, "Gain", "Loss"))) %>%
  mutate(Area = ifelse(Area_group == "Extinction", "Extinction", Area))


la_fu <- a1 %>% filter(Scenario %in% c("SSP_7100_126", "SSP_7100_370", "SSP_7100_585")) %>% dplyr::select(sp, long_mean) %>%
  group_by(sp) %>% summarise(mean = mean(long_mean)) %>% rename(la_mean_f = "mean")

la_c <- a1 %>% filter(Scenario == "SSP_Current") %>% dplyr::select(sp, long_mean) %>% rename(la_mean_c = "long_mean")

lon <- left_join(la_c, la_fu, by = "sp")

lon1 <- lon %>% mutate(lon_diff = la_mean_f - la_mean_c) %>%
  mutate(B_A = ifelse(lon_diff < 0, 1, 0)) %>%
  mutate(B_B = ifelse(lon_diff >= 0, 1, 0)) %>%
  mutate(B_A_1 = ifelse(lon_diff < 0 & lon_diff > -1, 1, 0)) %>%
  mutate(B_A_2 = ifelse(lon_diff <= -1, 1, 0)) %>%
  mutate(B_A_3 = ifelse(is.na(lon_diff), 1, 0)) %>%
  mutate(B_B_1 = ifelse(lon_diff < 1 & lon_diff >= 0, 1, 0)) %>%
  mutate(B_B_2 = ifelse(lon_diff >= 1, 1, 0)) %>%
  mutate(Lon_group = ifelse(B_A_1 == 1, "Dec.0-1", 
                            ifelse(B_A_2 == 1, "Dec.1>=", 
                                   ifelse(B_A_3 == 1, "Extiction", # Not operate, why?
                                          ifelse(B_B_1 == 1, "Inc.0-1", "Inc.1>="))))) %>%
  mutate(Lon_group = replace_na(Lon_group, "Extinction")) %>%
  mutate(Longitude = ifelse(B_A == 1, "Loss", ifelse(B_B == 1, "Gain", NA))) %>%
  mutate(Longitude = replace_na(Longitude, "Loss"))  

elv_fu <- a1 %>% filter(Scenario %in% c("SSP_7100_126", "SSP_7100_370", "SSP_7100_585")) %>% dplyr::select(sp, ele_mean) %>%
  group_by(sp) %>% summarise(mean = mean(ele_mean)) %>% rename(ele_mean_f = "mean")

elv_c <- a1 %>% filter(Scenario == "SSP_Current") %>% dplyr::select(sp, ele_mean) %>% rename(ele_mean_c = "ele_mean")

elv <- left_join(elv_c, elv_fu, by = "sp")

elv1 <- elv %>% mutate(ele_diff = ele_mean_f - ele_mean_c) %>%
  mutate(C_A = ifelse(ele_diff < 0, 1, 0)) %>%
  mutate(C_B = ifelse(ele_diff >= 0, 1, 0)) %>%
  mutate(C_A_1 = ifelse(ele_diff < 0 & ele_diff > -100, 1, 0)) %>%
  mutate(C_A_2 = ifelse(ele_diff <= -100, 1, 0)) %>%
  mutate(C_A_3 = ifelse(is.na(ele_diff), 1, 0)) %>%
  mutate(C_B_1 = ifelse(ele_diff < 100 & ele_diff >= 0, 1, 0)) %>%
  mutate(C_B_2 = ifelse(ele_diff >= 1, 1, 0)) %>%
  mutate(Elv_group = ifelse(C_A_1 == 1, "Dec.0-100", 
                            ifelse(C_A_2 == 1, "Dec.100>=", 
                                   ifelse(C_A_3 == 1, "Extiction", # Not operate, why?
                                          ifelse(C_B_1 == 1, "Inc.0-100", "Inc.100>="))))) %>%
  mutate(Elv_group = replace_na(Elv_group, "Extinction")) %>%
  mutate(Elevation = ifelse(C_A == 1, "Loss", ifelse(C_B == 1, "Gain", NA))) %>%
  mutate(Elevation = replace_na(Elevation, "Loss"))  

dflist <- list(rci1, lon1, elv1)

Total <- dflist %>% reduce(inner_join, by = "sp")

Group <- Total %>% select(sp, Area, Longitude, Elevation)
Anno <- Group[,-1]

names(Anno) <- c("Area", "Longitude", "Elevation")


original.path <- "D:/Envdata/Envdata_Lindera/Biomod_Lindera/"

niche.path <- "D:/Envdata/Thesis/Envdata_Quercus/Biomod_Quercus/Niche/"

Scenario <- c("SSP_current", "SSP_7100_126", "SSP_7100_370", "SSP_7100_585")

#All_Species <- as.vector(Group[Group$Area %in% c("Loss", "Gain", "Extinction"), 1])
All_Species <- Group %>%
  filter(Area %in% c("Loss", "Gain", "Extinction")) %>%
  pull(1)

Group_Gain <- as.vector(as.data.frame(Group[Group$Area == "Gain", 1])[,1])
Group_Loss <- as.vector(as.data.frame(Group[Group$Area == "Loss" | Group$Area == "Extinction", 1])[,1])

Gain <- paste0(gsub(" ", ".", Group_Gain), "_bin")
Loss <- paste0(gsub(" ", ".", Group_Loss), "_bin")

for (i in 1:length(Scenario)) {
  dir.create(paste0(niche.path, Scenario[i], "/Gain/"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(niche.path, Scenario[i], "/Loss/"), recursive = TRUE, showWarnings = FALSE)
}

# file copy (Gain)
for (i in 1:length(Scenario)) {
  for (k in 1:length(Gain)) {
    from <- paste0(original.path, Scenario[i], "/", Gain[k], ".tif")
    to <- paste0(niche.path, Scenario[i], "/Gain/", Gain[k], ".tif")
    cat("Copying from:", from, "to", to, "\n")
    # 파일 복사
    if (file.exists(from)) {
      file.copy(from, to, overwrite = TRUE)
    } else {
      warning(paste("File not found:", from))
    }
  }
}

# file copy (Loss)
for (i in 1:length(Scenario)) {
  for (k in 1:length(Loss)) {
    from <- paste0(original.path, Scenario[i], "/", Loss[k], ".tif")
    to <- paste0(niche.path, Scenario[i], "/Loss/", Loss[k], ".tif")
    cat("Copying from:", from, "to", to, "\n")
    # 파일 복사
    if (file.exists(from)) {
      file.copy(from, to, overwrite = TRUE)
    } else {
      warning(paste("File not found:", from))
    }
  }
}


###### Niche Overlap calculation----
library(tidyverse)
library(raster)
library(rgdal)
library(humboldt)
library(ENMeval)
library(dismo)

# For Gain
n.groups <- length(Group_Gain)
D <- matrix(nrow = n.groups, ncol = n.groups)
g.codenames <- data.frame(Group_Gain)
names(g.codenames) <- "sp"
rownames(D) <- colnames(D) <- g.codenames$sp

for (i in 1:length(Scenario)) {
  files <- list.files(path = paste0(niche.path, Scenario[i], "/Gain/"), pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) {
    cat("No TIF files found for scenario:", Scenario[i], "\n")
    next
  } else {
    cat("Found TIF files for scenario:", Scenario[i], " - ", length(files), " files\n")
  }
  z <- stack(files)
  for (m in 2:n.groups) {
    for (n in 1:(m - 1)) {
      x1 <- z[[m]]
      x2 <- z[[n]]
      result <- try(nicheOverlap(x1, x2, "D", mask = TRUE), silent = TRUE)
      if (inherits(result, "try-error")) {
        cat("Error in nicheOverlap for files", m, "and", n, "\n")
      } else {
        D[m, n] <- result
      }
    }
  }
  write.csv(x = D, file = paste0(niche.path, Scenario[i], "_Gain.csv"))
}


# For Loss
n.groups <- length(Group_Loss)
D <- matrix(nrow = n.groups, ncol = n.groups)
g.codenames <- data.frame(Group_Loss)
names(g.codenames) <- "sp"
rownames(D) <- colnames(D) <- g.codenames$sp

for (i in 1:length(Scenario)) {
  files <- list.files(path = paste0(niche.path, Scenario[i], "/Loss/"), pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) {
    cat("No TIF files found for scenario:", Scenario[i], "\n")
    next
  } else {
    cat("Found TIF files for scenario:", Scenario[i], " - ", length(files), " files\n")
  }
  z <- stack(files)
  for (m in 2:n.groups) {
    for (n in 1:(m - 1)) {
      x1 <- z[[m]]
      x2 <- z[[n]]
      result <- try(nicheOverlap(x1, x2, "D", mask = TRUE), silent = TRUE)
      if (inherits(result, "try-error")) {
        cat("Error in nicheOverlap for files", m, "and", n, "\n")
      } else {
        D[m, n] <- result
      }
    }
  }
  write.csv(x = D, file = paste0(niche.path, Scenario[i], "_Loss.csv"))
}

# For All Species

#Scenario <- c("SSP_7100_126", "SSP_7100_370", "SSP_7100_585")

n.groups <- length(All_Species)
D <- matrix(nrow = n.groups, ncol = n.groups)
g.codenames <- data.frame(All_Species)
names(g.codenames) <- "sp"
rownames(D) <- colnames(D) <- g.codenames$sp

for (i in 1:length(Scenario)) {
  files <- list.files(path = paste0(original.path, Scenario[i], "/"), pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) {
    cat("No TIF files found for scenario:", Scenario[i], "\n")
    next
  } else {
    cat("Found TIF files for scenario:", Scenario[i], " - ", length(files), " files\n")
  }
  z <- stack(files)
  for (m in 2:n.groups) {
    for (n in 1:(m - 1)) {
      x1 <- z[[m]]
      x2 <- z[[n]]
      result <- try(nicheOverlap(x1, x2, "D", mask = TRUE), silent = TRUE)
      if (inherits(result, "try-error")) {
        cat("Error in nicheOverlap for files", m, "and", n, "\n")
      } else {
        D[m, n] <- result
      }
    }
  }
  write.csv(x = D, file = paste0(niche.path, Scenario[i], ".csv"))
}


#####visualization_heatmap-----

library(ComplexHeatmap)
library(pheatmap)
library(tibble)
library(grid)
library(ggplot2)
library(gghalves)
library(ggpubr)
library(rstatix)


# Reshaping of correlation data-----
test <- D

aa <- test %>% as.data.frame() %>% rownames_to_column(var = "var1") %>%
  gather(key = "var2", value = "cor", -var1)
ab <- test %>% as.data.frame() %>% rownames_to_column(var = "var2") %>%
  gather(key = "var1", value = "cor", -var2) %>% rename(var2 = var1, var1 = var2)
ac <- data.frame(var1 = g.codenames$sp, var2 = g.codenames$sp, cor = 1)

tmp <- bind_rows(aa, ab, ac) %>%
  group_by(var1, var2) %>%
  summarise(cor = mean(cor, na.rm = TRUE), .groups = "drop") %>%  # 또는 first(cor)
  tidyr::pivot_wider(names_from = var2, values_from = cor) %>%
  tibble::column_to_rownames("var1") %>%
  as.matrix()


# Heatmap 
dir.create("D:/Envdata/Envdata_Lindera/Biomod_Lindera/Niche/NicheOverlap", recursive = TRUE)
plot.path <- "D:/Envdata/Envdata_Lindera/Biomod_Lindera/Niche/NicheOverlap"

annot_colors <- list(Area = c(Gain = "red", Loss = "blue", Extinction = "white"),
                     Longitude = c(Gain = "red", Loss = "blue", Extinction = "white"), 
                     Elevation = c(Gain = "red", Loss = "blue", Extinction = "white"))

setwd("D:/Envdata/Envdata_Lindera/Biomod_Lindera/Niche/")

Scenario <- c("SSP_current", "SSP_7100_126", "SSP_7100_370", "SSP_7100_585")


# All scenario-----

for (i in seq_along(Scenario)) {
  message(">> ", Scenario[i])
  
  test <- read.csv(file = paste0(Scenario[i], ".csv"))
  names(test) <- gsub(x = names(test), pattern = "\\.", replacement = " ")
  test <- test %>% tibble::column_to_rownames("X")
  
  aa <- test %>% cor_gather() %>%
    dplyr::rename(var1 = "var1", var2 = "var2", cor = "cor")
  
  ab <- test %>% cor_gather() %>%
    dplyr::select(var2, var1, cor) %>%
    dplyr::rename(var2 = "var1", var1 = "var2", cor = "cor")
  
  ac <- tibble::tibble(var1 = g.codenames$sp, var2 = g.codenames$sp, cor = 1)
  
  tmp <- dplyr::bind_rows(aa, ab, ac) %>%
    dplyr::group_by(var1, var2) %>%
    dplyr::summarise(cor = mean(cor, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = var2, values_from = cor) %>%
    tibble::column_to_rownames("var1") %>%
    as.matrix()
  
  Anno_loss <- Group %>%
    dplyr::filter(sp %in% Group_Loss) %>%
    dplyr::select(sp, Area, Longitude, Elevation) %>%
    tibble::column_to_rownames(var = "sp")
  
  newnames_row <- lapply(rownames(tmp), function(x) bquote(italic(.(x))))
  newnames_col <- lapply(colnames(tmp), function(x) bquote(italic(.(x))))
  
  png(file = paste0(plot.path, "/", Scenario[i], ".png"),
      res = 300, height = 8, width = 12, units = "cm")
  
  pheatmap::pheatmap(tmp,
                     annotation_row = Anno_loss,
                     annotation_names_col = FALSE,
                     fontsize = 6,
                     annotation_colors = annot_colors,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     display_numbers = FALSE,
                     labels_row = as.expression(newnames_row),
                     labels_col = as.expression(newnames_col)
  )
  dev.off()
}


###### Niche overlap gradient graph-----
getwd()

setwd("D:/Envdata/Envdata_Lindera/Biomod_Lindera/Niche/")


tmp_total <- data.frame()

# Read the first CSV file and clean column names
tmp1 <- read.csv("SSP_7100_585.csv")
names(tmp1) <- gsub(x = names(tmp1), pattern = "\\.", replacement = " ")
test <- tmp1 %>% tibble::column_to_rownames("X")
str(g.codenames)

# Compute correlations and reshape data
aa <- test %>% cor_gather() %>% rename(var1 = "var1", var2 = "var2", cor = "cor")
ab <- test %>% cor_gather() %>% select(var2, var1, cor) %>% rename(var2 = "var1", var1 = "var2", cor = "cor")
ac <- data.frame(g.codenames, g.codenames, 1) %>% rename(var1 = "sp", var2 = "sp.1", cor = "X1")
tmp <- rbind(aa, ab, ac) %>% arrange(var2, var1)

# Merge with group data and prepare final data frame
tmp2 <- tmp %>% 
  left_join(Group, by = c("var1" = "sp")) %>%
  left_join(Group, by = c("var2" = "sp")) %>% 
  select(var1, var2, cor, Area.x, Area.y) %>%
  mutate(Group = ifelse(Area.x == "Gain" & Area.y == "Gain", "Gain",
                        ifelse(Area.x == "Loss" & Area.y == "Loss", "Loss", "G-L"))) %>%
  mutate(class = "SSP_7100_585", code = paste0(var1, "+", var2))

str(tmp2)

# Append to the total data frame 1 
tmp_total <- rbind(tmp_total, tmp2)

# Read the second CSV file and clean column names
tmp1 <- read.csv("SSP_current.csv")
names(tmp1) <- gsub(x = names(tmp1), pattern = "\\.", replacement = " ")
test <- tmp1 %>% tibble::column_to_rownames("X")
str(g.codenames)

# Compute correlations and reshape data
aa <- test %>% cor_gather() %>% rename(var1 = "var1", var2 = "var2", cor = "cor")
ab <- test %>% cor_gather() %>% select(var2, var1, cor) %>% rename(var2 = "var1", var1 = "var2", cor = "cor")
ac <- data.frame(g.codenames, g.codenames, 1) %>% rename(var1 = "sp", var2 = "sp.1", cor ="X1")
tmp <- rbind(aa, ab, ac) %>% arrange(var2, var1)

# Merge with group data and prepare final data frame
tmp2 <- tmp %>% 
  left_join(Group, by = c("var1" = "sp")) %>% 
  left_join(Group, by = c("var2" = "sp")) %>%
  select("var1", "var2", "cor", "Area.x", "Area.y") %>%
  mutate(Group = ifelse(Area.x == "Gain" & Area.y == "Gain", "Gain", 
                        ifelse(Area.x == "Loss" & Area.y == "Loss", "Loss", "G-L"))) %>%
  mutate(class = "SSP_Current", code = paste(var1, "+", var2))

# Append to the total data frame 2 
tmp_total <- rbind(tmp_total, tmp2)

# Check final data
print(head(tmp_total))
print(table(tmp_total$class, tmp_total$Group))

# Remove duplicate data
tmp_total <- tmp_total %>% distinct()

# Check data again !
print(table(tmp_total$class, tmp_total$Group))

# Filter and prepare data for plotting
data2 <- tmp_total %>% filter(class == "SSP_Current" | class == "SSP_7100_585")
data2$class <- factor(data2$class, levels = c("SSP_Current", "SSP_7100_585"), order = TRUE)
data2$Group <- as.factor(data2$Group)
data2$code <- as.factor(data2$code)

# Adding a jittered numeric value for x-axis plotting
y <- ifelse(data2$class == "SSP_Current", 1, ifelse(data2$class == "SSP_7100_585", 2, 3))
data2$classnum <- as.factor(y)
data2$y <- jitter(y, amount = .05)
data3 <- data2 %>% filter(cor != 1)

# Check data again !
print(head(data3))
print(table(data3$class, data3$Group))

# filtering data
data_current <- data3 %>% filter(class == "SSP_Current")
data_ssp <- data3 %>% filter(class == "SSP_7100_585")

# Create the plot
library(gghalves)
library(ggpubr)
library(OPTS)
library(ggplot2)

tiff(file = "Gradient.tiff", res = 300, height = 1400, width = 800)

ggplot(data=data3, aes(y=cor, fill=classnum)) + 
  geom_half_boxplot(data=data_current, aes(x=classnum, y=cor), width=0.6, side = "l") +
  geom_half_boxplot(data=data_ssp, aes(x=classnum, y=cor), width=0.6, side = "r") +
  geom_point(aes(x=classnum, y=cor, color=classnum), alpha = 0.1, size=1) +
  geom_line(aes(x=classnum, y=cor, group=interaction(var1, var2)), color='lightgray', alpha = 0.3) +
  ylab(expression("Schoener's " * italic(D))) +
  scale_x_discrete(labels=c("1"="Current", "2"="SSP5-8.5(2071-2100)")) +
  scale_color_manual(values = c("coral3", "#287D8EFF")) +
  scale_fill_manual(values = c("coral3", "#287D8EFF")) +
  facet_grid(~Group) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank())

dev.off()