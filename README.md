---
title: "Lidar Processing"
author: "Anish Bhattarai"
date: "`r Sys.Date()`"
output: 
  html_document:
    theme: cerulean
    toc: yes
    toc_float:
      collapsed: true
---

# Lidar-Processing

Here in this workflow we will:

1.  Read and write .las and .laz files and render customized point-cloud display

2.  Process point clouds including point classification into ground and non-ground points,

3.  Develop a digital terrain models (DTM),

4.  normalization of cloud points based on DTM

5.  Develop a digital surface models (DSM)

6.  Rasterize the point clouds

7.  Get the canopy height of the plot level (h_max (100 percentile), h_mean, h_min and h_median and percentiles (80, 85, 90, 95).

8.  Do a ANOVA and correlate with the hand collected data (ground truth)

# Installing lidR packages

```{r}
# install.packages("lidR")
```

# Loading all required libraries

```{r}
library(lidR) # package to import and wrangle lidar file (laz or las)
library(dplyr) # for wrangling the normal dataframes
library(tidyverse) # for wrangling and filtering
library(ggplot2) # for data visualization
library(sf) # to manipulate the simple feature geometry objects
library(stars) # to manipulate the raster objects
library(raster) # to manipulate raster object
library(terra) 
library(janitor)
library(ggpmisc) # to plot and get statistical parameters in the plot
library(ggthemes) # to change the themes
library(car) # to perform ANOVA
library(lme4) # to perform linear mixed effect models
library(emmeans) # to perform the emmeans after running the model
library(multcomp) # to perform pairwise comparision
library(multcompView)
library(geometry)
```

# Loading the las file

```{r}
jpc_las <- readLAS("../data/all_laz_file/mv_9_160.copc.laz",
                   select = "xyzrnc") # here we are reading the laz file and selecting only X,Y,Z coordinates, return number and number of return.

summary(jpc_las@data)

```

# Clipping the las file

The larger the area bigger is the file size and it will take a lot of processing time, so by clipping with the area of interest we can reduce the area and subsequently the file size to speed up the processing.

```{r}
clip <- read_sf("../data/clip_bbox.geojson") %>% 
  st_transform(crs = st_crs(jpc_las))

jpc_las_c <- clip_roi(jpc_las, clip)

writeLAS(jpc_las_c, "../output/clipped_mv_9_160_smallfile.laz")
```

# Classification of the point clouds

The algorithm that is used is the progressive morphological filter (pmf). It has 2 arguments, window size (which is in meters), which looks at the cloud points within those window size and threshold (th), which looks into the verticle distance and classify them as ground (also in meters).

There are other algorithm for the classification like cloth simulation filter(csf), which is most widely used. Pmf is generally used when the vegetation is too dense while csf is used when the vegetation is less dense.

We can **skip** this classification process because, the classification of point clouds can be performed in DJI terra app, which is faster and produces a RGB cloud point in the app, which can be visually assessed.

```{r}
# las <- classify_ground(jpc_las, pmf(ws = 0.1, th = 0.06)) #ws = window size (pixel size in meters), th = hreshold heights above the parameterized ground surface to be considered a ground return.
# 
# las@data
```

# DTM from the ground points

Digital Terrain Models (DTM) sometimes called Digital Elevation Models (DEM) is a topographic model of the bare Earth that can be manipulated by computer. Using triangular irregular network (tin), where it looks into the ground points and creates multiple triangles joining them and what ever falls in the plane of these are used to make a DTM.

```{r}
dtm <- rasterize_terrain(jpc_las_c, res = 0.1, algorithm = tin()) # for this to work we must first classify the ground and non-ground points, the better we can quality better DTM we can get
```

# Normalization

As the point clouds are initially measured as mean sea level (MSL) coming out from the drone, we have to normalize such that we make the ground points as 0 and the point clouds above them are now turned into above ground level (AGL). As we already created the dtm we can now use that dtm as reference to normalize the height where the dtm values are now 0, as it is the ground surface. This normalization now gives us the a digital surface model (DSM), which means, if there is a bare ground it represents the ground, but incase there is a tree,crops, buildings etc. above it represents the height of those structure.

For normalization there are other algorithm like, knnidw and kriging. But these 2 methods takes a lot of processing time. And based on alot of paper, tin algorithm has a balance between processing time and preciseness, so it is most widely used.

```{r}
normalized <- normalize_height(jpc_las_c, algorithm = tin(), dtm = dtm) 

normalized@data
```

# Filtering the point clouds

While DTM is created using TIN there might be some points that falls below the triangular structure created, so we want to make sure we get rid of those noise and we also want to make sure to get rid of any higher structures that does not represents our interest.

```{r}
filt_norm <- filter_poi(normalized, Z > 0, Z < 2) # upper threshold can be changed aaccording the plant height during the growing season

filt_norm@data %>% 
  group_by(Classification) %>% 
  summarise(avg_z = mean(Z))
```

# Top view (takes time)

```{r}
# ggplot()+
#   geom_point(data = filt_norm@data, aes(x = X, y = Y, color = factor(Classification)), size = 0.1)+
#   coord_equal()+
#   scale_color_manual(values = c("forestgreen",
#                                 "burlywood"))
```

# Sectional view

```{r, warning=F, message=F}
# ggplot()+
#   geom_point(data = filt_norm@data, aes(x = X, y = Z, color = factor(Classification)), size = 0.01)+
#   coord_equal(ratio = 10)+ # might need to change this according to the requirement.
#   scale_color_manual(values = c("forestgreen",
#                                 "burlywood"))
```

# Desnity plot

```{r}
ggplot()+
  geom_density(data = filt_norm@data, aes(x = Z))+
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0,1, 0.1))
```

# CHM

Getting only the cloud points that are plants so we can rasterize the canopy height.Once the surface heights are removed we get the canopy height measurement.

```{r}
chm <- filter_poi(filt_norm, Z > 0.1, Classification == 1) # changes can be made for the z height according to the requirement. I put 0.06 here as I am considering

chm
```

```{r}
ggplot()+
  geom_density(data = chm@data, aes(x = Z))+
  scale_x_continuous(limits = c(0, 2), 
                     breaks = seq(0,2, 0.1))
```

# Canopy Rasterization

Point to raster (p2r) method is most commonly used for the canopy rasterization. p2r algorithms are conceptually simple, consisting of establishing a grid at a user defined resolution and attributing the elevation of the highest point to each pixel.

```{r}
chm_rast <- rasterize_canopy(chm, res = 0.1, algorithm = p2r(subcircle = 0.1))

chm_stars <- chm_rast %>% 
  st_as_stars()

write_stars(chm_stars, "../output/mv_default.tif",
            delete_dsm = T)
```

```{r}
ggplot()+
  geom_stars(data = chm_stars)
```

# Plot file

```{r}
plots_mv <- read_sf("../data/ten_plants_bbox.geojson") %>% 
  arrange(plot) %>% 
  st_transform(crs = st_crs(chm_stars))
```

# Intersection

```{r}
plot_df <- st_intersection(plots_mv, st_as_sf(chm_stars)) %>% 
  arrange(plot)
```

# Height parameters lidar

```{r}
lidar_calc <- plot_df %>% 
  group_by(plot) %>% 
  summarise(h_mean = mean(Z)*100,
            h_max = max(Z)*100,
            h_min = min(Z)*100,
            h_med = median(Z)*100) %>% 
  mutate(plot = as.integer(plot))

write_csv(lidar_calc, "../output/lidar_calc_default.csv")
```

# height parameters hand

```{r}
manual <- readxl::read_xlsx("../data/hand collected data/data/1_hand data mv jun 21.xlsx") %>% 
  janitor::clean_names() %>% 
  dplyr::select(-sample, -date)

manual_w <- manual %>% 
  rename(plots = plot) %>% 
  group_by(plots) %>% 
  summarise(ht_mean_man = mean(height_cm),
            ht_max_man = max(height_cm),
            ht_med_man = median(height_cm))
```

# lidar and hand merged df

```{r}
both <- left_join(lidar_calc, manual_w, by = "plots") %>% 
  mutate(h_mean = ceiling(h_mean),
         ht_mean_man = floor(ht_mean_man))
```

# regression

```{r}
both %>% 
  ggplot(aes(x = ht_mean_man, y = h_mean))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq()+
  geom_abline(slope = 1, linetype = "dashed", color = "red")
```

# Anova

## Importing info file

```{r}
info_mv <- read_csv("../data/info_midville.csv") %>% 
  dplyr::select(plots = plot_ids, blocks, treatment)
```

```{r}
manual_info <- manual %>% 
  left_join(info_mv, by = c("plot" = "plots"))

both_wo_geom <- both %>% 
  st_drop_geometry() %>% 
  left_join(info_mv, by = "plots")
```

# Model manual

```{r}
model_man <- lmer(height_cm ~ treatment + (1|blocks/plot), data = manual_info)

Anova(model_man, type = 3)
```

# Model lidar

```{r}
model_lidar <- lmer(h_mean ~ treatment + (1|blocks), data = both_wo_geom)

Anova(model_lidar, type = 3)
```