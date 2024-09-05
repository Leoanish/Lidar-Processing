<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/GiorgosXou/Random-stuff/main/Programming/StackOverflow/Answers/70200610_11465149/w.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/GiorgosXou/Random-stuff/main/Programming/StackOverflow/Answers/70200610_11465149/b.png">
  <img alt="Shows a black logo in light color mode and a white one in dark color mode." src="https://user-images.githubusercontent.com/25423296/163456779-a8556205-d0a5-45e2-ac17-42d089e3c3f8.png">
</picture>

---
title: "Lidar Processing"
date: 2024/09/04
output: rmdformats::readthedown
---

```{css}
.badCode {
background-color: red;
}
```

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
library(tidymodels)
```

# Loading the las file

```{r}
las <- readLAS("../data/clipped_mv_9_160_smallfile.laz",
                   select = "xyzrnc") # here we are reading the laz file and selecting only X,Y,Z coordinates, return number and number of return.

summary(las@data)

```

# Clipping the las file

The larger the area bigger is the file size and it will take a lot of processing time, so by clipping with the area of interest we can reduce the area and subsequently the file size to speed up the processing.

```{r}
clip <- read_sf("../data/clip_bbox.geojson") %>% 
  st_transform(crs = st_crs(las))

las_c <- clip_roi(las, clip)
```

# Classification of the point clouds

The algorithm that is used is the progressive morphological filter (pmf). It has 2 arguments, window size (which is in meters), which looks at the cloud points within those window size and threshold (th), which looks into the verticle distance and classify them as ground (also in meters).

There are other algorithm for the classification like cloth simulation filter(csf), which is most widely used. Pmf is generally used when the vegetation is too dense while csf is used when the vegetation is less dense.

We can **skip** this classification process because, the classification of point clouds can be performed in DJI terra app, which is faster and produces a RGB cloud point in the app, which can be visually assessed.

```{r}
# las <- classify_ground(jpc_las, pmf(ws = 0.1, th = 0.06)) #ws = window size (pixel size in meters), th = threshold heights above the parameterized ground surface to be considered a ground return.
# 
# las@data
```

# DTM from the ground points

Digital Terrain Models (DTM) sometimes called Digital Elevation Models (DEM) is a topographic model of the bare Earth that can be manipulated by computer. Using triangular irregular network (tin), where it looks into the ground points and creates multiple triangles joining them and what ever falls in the plane of these are used to make a DTM.

```{r}
dtm <- rasterize_terrain(las_c, res = 0.1, algorithm = tin()) # for this to work we must first classify the ground and non-ground points, the better we can quality better DTM we can get
```

# Normalization

As the point clouds are initially measured as mean sea level (MSL) coming out from the drone, we have to normalize such that we make the ground points as 0 and the point clouds above them are now turned into above ground level (AGL). As we already created the dtm we can now use that dtm as reference to normalize the height where the dtm values are now 0, as it is the ground surface. This normalization now gives us the a digital surface model (DSM), which means, if there is a bare ground it represents the ground, but incase there is a tree,crops, buildings etc. above it represents the height of those structure.

For normalization there are other algorithm like, knnidw and kriging. But these 2 methods takes a lot of processing time. And based on alot of paper, tin algorithm has a balance between processing time and preciseness, so it is most widely used.

```{r}
normalized <- normalize_height(las_c, algorithm = tin(), dtm = dtm) 

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
ggplot()+
  geom_point(data = filt_norm@data, aes(x = X, y = Y, color = factor(Classification)), size = 0.1)+
  coord_equal()+
  scale_color_manual(values = c("forestgreen",
                                "burlywood"))
```

# Sectional view

```{r, warning=F, message=F}
ggplot()+
  geom_point(data = filt_norm@data, aes(x = X, y = Z, color = factor(Classification)), size = 0.01)+
  coord_equal(ratio = 10)+ # might need to change this according to the requirement.
  scale_color_manual(values = c("forestgreen",
                                "burlywood"))
```

# Desnity plot

```{r}
ggplot()+
  geom_density(data = filt_norm@data, aes(x = Z))+
  scale_x_continuous(limits = c(0, 1.5), 
                     breaks = seq(0,1.5, 0.1))
```

# CHM

Getting only the cloud points that are plants so we can rasterize the canopy height.Once the surface heights are removed we get the canopy height measurement.

```{r}
chm <- filter_poi(filt_norm, Z > 0.1, Classification == 1) # changes can be made for the z height according to the requirement. I put 0.1 here as I am considering small bumps, weeds and dead plants as noise.

chm
```

```{r}
ggplot()+
  geom_density(data = chm@data, aes(x = Z))+
  scale_x_continuous(limits = c(0, 1.5), 
                     breaks = seq(0,1.5, 0.1))
```

# Canopy Rasterization

Point to raster (p2r) method is most commonly used for the canopy rasterization. p2r algorithms are conceptually simple, consisting of establishing a grid at a user defined resolution and attributing the elevation of the highest point to each pixel.

```{r}
chm_rast <- rasterize_canopy(chm, res = 0.6, algorithm = p2r(subcircle = 0.1))

chm_stars <- chm_rast %>% 
  st_as_stars()

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
  arrange(plot) %>% 
  filter(!plot == 307)
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

write_csv(lidar_calc, "../output/lidar_calc_9_160.csv")
```

# height parameters hand

```{r}
manual <- readxl::read_xlsx("../data/5_hand data mv aug 26.xlsx") %>% 
  janitor::clean_names() %>% 
  dplyr::select(-sample, -date)

manual_w <- manual %>% 
  group_by(plot) %>% 
  summarise(ht_mean_man = mean(height_cm),
            ht_max_man = max(height_cm),
            ht_med_man = median(height_cm))
```

# lidar and hand merged df

```{r}
both <- left_join(lidar_calc, manual_w, by = "plot")
  # mutate(h_mean = ceiling(h_mean),
  #        ht_mean_man = floor(ht_mean_man))

summary(lm(ht_mean_man ~ h_mean, data = both))
```

# regression

```{r}
both %>% 
  ggplot(aes(x = ht_mean_man, y = h_mean))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label(labels = c("eq", "R2")))+
  geom_abline(slope = 1, linetype = "dashed", color = "red")+
  coord_obs_pred()
```

```{r}
metrica::scatter_plot(obs = both$ht_mean_man, pred = both$h_mean, print_metrics = T, metrics_list = c('R2','RMSE', "MAE"), position_eq = c(x = 0.7*90, y = 1.25*86), position_metrics = c(x = 63, y = 1.05*98))+
                       labs(title = "9_160")
```


# Anova

## Importing info file

```{r}
chm_stars_sp <- chm_stars %>% 
  as("SpatRaster")
```

```{r}
w <- matrix(1,3,3)

filled <- terra::focal(chm_stars_sp, w, fun = max, na.rm = T) %>% 
  st_as_stars()

ggplot()+
  geom_stars(data = filled, aes(fill = focal_max))
```

```{r}
plot_df_sp <- st_intersection(plots_mv, st_as_sf(filled)) %>% 
  arrange(plot) %>% 
  filter(!plot == 307)
```

# Height parameters lidar

```{r}
lidar_calc_sp <- plot_df_sp %>% 
  group_by(plot) %>% 
  summarise(h_mean = mean(focal_max)*100,
            h_max = max(focal_max)*100,
            h_min = min(focal_max)*100,
            h_med = median(focal_max)*100) %>% 
  mutate(plot = as.integer(plot))

```

```{r}
both_sp <- left_join(lidar_calc_sp, manual_w, by = "plot")
```

```{r}
both_sp %>% 
  ggplot(aes(x = ht_mean_man, y = h_mean))+
  geom_point()+
  stat_poly_line()+
  stat_poly_eq(use_label(labels = c("eq", "R2")))+
  geom_abline(slope = 1, linetype = "dashed", color = "red")+
  coord_obs_pred()
```

