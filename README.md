# geobuffer: Geodesic buffer using metric radius

# Overview

R package that allows the possibility of creating **geodesic buffers** when the radius is given in metric units. A geodesic buffer is not affected by the distortions introduced by projected coordinate systems.

In order to use `rgeos::gBuffer()` with a metric radius, one has to project the coordinates. Projecting and then applying `gBuffer` means actually producing **Euclidean** buffers as opposed to **Geodesic** ones. Euclidean buffers are affected by distortions.

<p float="left">
  <img src="https://i.stack.imgur.com/nr2bP.jpg" width="200" />
  <img src="https://i.stack.imgur.com/5wyFG.jpg" width="200" /> 
</p>
<sup>Left – Euclidian buffers, affected by distortions; Right – geodesic buffers, no distortions.</sup>

The `geobuffer` package avoids this problem, producing directly geodesic buffers with a given metric radius by wrapping around the `geosphere::destPoint()`.

The idea for the code was first expressed on gis.stackexchange - [Euclidean and Geodesic Buffering in R](https://gis.stackexchange.com/questions/250389/euclidean-and-geodesic-buffering-in-r). A related question was addressed on Stack Overflow [here](https://stackoverflow.com/questions/25411251/buffer-geospatial-points-in-r-with-gbuffer). Relevant is also the ESRI article [Understanding Geodesic Buffering](https://www.esri.com/news/arcuser/0111/geodesic.html).

Also is worth checking the package [dggridR](https://github.com/r-barnes/dggridR/), where the "messy distortions" are tackled when your analysis involves some spatial binning. Their README page explains very well the distortions problems and their consequences.

# Installation

You can install `geobuffer` from GitHub with:

``` r
# install.packages("devtools") # if you do not have devtools, then install it
devtools::install_github("valentinitnelav/geobuffer")
```
# Examples

## Example 1 – simple geodesic buffer from given point

Currently there is only one main function in the package, ` geobuffer_pts()`, which produces geodesic buffers around given points.

``` r
library(geobuffer)
library(mapview)
library(sf)
#> Linking to GEOS 3.6.1, GDAL 2.2.3, PROJ 4.9.3

bucharest_500km <- geobuffer_pts(xy = data.frame(lon = 26.101390,
                                                 lat = 44.427764),
                                 dist_m = 500*10^3,
                                 output = "sf")
bucharest_500km
#> Simple feature collection with 1 feature and 0 fields
#> geometry type:  POLYGON
#> dimension:      XY
#> bbox:           xmin: 19.83443 ymin: 39.92637 xmax: 32.36835 ymax: 48.9256
#> epsg (SRID):    4326
#> proj4string:    +proj=longlat +datum=WGS84 +no_defs
#>                         geometry
#> 1 POLYGON ((26.10139 48.9256,...

mapView(as(bucharest_500km, "Spatial"), alpha.regions = 0.2)
```
<img src="https://i.imgur.com/wRuhqYF.png" width="400" />

## Example 2 – Multiple points; data.frame for `ggplot2` and allow other buffer shapes

The function ` geobuffer_pts` is vectorized, accepting as input multiple points (be it `SpatialPoints`, `SpatialPointsDataFrame`, `sf` points, or two columns `matrix`, `data.frame` or `data.table`).

By adjusting the `step_dg` argument, one can obtain hexagons, triangles or circle-like shapes.

For convenience, the returned output can also be a data.frame that can be used directly for plotting with `ggplot2::ggplot()`.

``` r
library(geobuffer)
library(ggplot2)

buffers <- geobuffer_pts(xy = data.frame(lon = c(26.101390, 25.6112233),
                                         lat = c(44.427764, 45.6523994)),
                         dist_m = c(50*10^3, 30*10^3),
                         step_dg = 60,
                         output = "data.frame")
str(buffers)
#> 'data.frame':    14 obs. of  3 variables:
#>  $ lon: num  26.1 26.6 26.6 26.1 25.6 ...
#>  $ lat: num  44.9 44.7 44.2 44 44.2 ...
#>  $ id : int  1 1 1 1 1 1 1 2 2 2 ...
ggplot(data = buffers,
       aes(x = lon, y = lat, group = id)) +
  geom_polygon(aes(fill = as.factor(id))) +
  coord_fixed()
```
<img src="https://i.imgur.com/BVrXNsq.png" width="500" />

Also can plot the buffers with `tmap`.

``` r
library(geobuffer)
library(tmap)
data("World")

buffers <- geobuffer_pts(xy = data.frame(lon = c(26.101390, 25.6112233),
                                         lat = c(44.427764, 45.6523994)),
                         dist_m = c(50*10^3, 30*10^3),
                         step_dg = 60)

tm_shape(World[World$iso_a3 == "ROU",]) +
  tm_polygons() +
  tm_shape(buffers) +
  tm_polygons(col = "red")
```
<img src="https://i.imgur.com/Ic0SpJC.png" width="500" />

<sup>Created on 2019-02-11 by the [reprex package](https://reprex.tidyverse.org) (v0.2.1)</sup>

https://doi.org/10.5281/zenodo.7127818

# How to cite the package?

I just uploaded this packge on Zenodo after almost 4 years :)
Would be nice to cite it if you make use of it.

You can cite the package as:

> Valentin Ștefan (2022). R package for constructing geodesic buffers using metric radiuses. https://doi.org/10.5281/zenodo.7127818

# To do

- [ ] Make use of the faster alternative proposed [here](https://gis.stackexchange.com/a/251873/62753);
- [ ] Fix the buffer display when crossing the antemeridian;
- [ ] Geodesic buffer for lines and polygons;
