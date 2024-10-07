libs <- c(
  "tidyr", "dplyr", "terra",
  "sf", "exactextractr",
  "rgeoboundaries", "ggplot2",
  "plotly", "grid", "ggtern", "gganimate",
  "magick", "elevatr", "leaflet"
)

installed_libraries <- libs %in% rownames(
  installed.packages()
)

if (any(installed_libraries == F)) {
  install.packages(libs[!installed_libraries])
}

invisible(
  lapply(
    libs, library,
    character.only = T
  )
)
getwd()  # Check the current working directory
setwd("~/Uni/2024S/Bachelor/Spatio-Temporal-Analysis-Of-Sustainability-Data/Central Europe Analysis") 
getwd()



load_vrt_layers <- function(years) {
  virtual_layers <- list()
  
  for (i in seq_along(years)) {  
    year <- years[i]
    vrt_file <- paste0("central_europe_vrt", year, ".vrt")
    
    if (file.exists(vrt_file)) {
      virtual_layers[[i]] <- terra::rast(vrt_file)
    } else {
      warning(paste("VRT file for year", year, "not found."))
    }
  }
  
  return(virtual_layers)
}

years <- 2018:2023
virtual_layers <- load_vrt_layers(years)

terra::plot(virtual_layers[[1]])

zug <- rbind(gb_adm1("Switzerland")) %>%
  filter(shapeName == "Zug")

crs <- "EPSG:4326"
lc_list <- list()
for (index in seq_along(virtual_layers)) {
  year <- years[index]
  rasters <- virtual_layers[[index]]
  
  for (i in seq_along(rasters)) {
    raster <- rasters[[i]]
    
    state <- zug %>%
      sf::st_transform(crs = terra::crs(raster))
    
    land_cover <- raster %>%
      crop(vect(state), snap = "in") %>%
      mask(vect(state)) %>%
      terra::project(crs)
    
    lc_list[[index]] <-land_cover
  }
}

terra::plot(lc_list[[1]])

process_layer <- function(layer) {
  terra::ifel(layer == 2L,
              1,
              terra::ifel(layer == 7L,
                          2,
                          terra::ifel(layer == 5L,
                                      3,
                                      terra::ifel(layer == 1L, 4,
                                                  terra::ifel(layer == 11L, 5, NA)))))
}

land_cover <- lapply(lc_list, process_layer)

images <- list()

for (i in seq_along(land_cover)) {
  year <- years[i]
  layer <- land_cover[[i]]
  
  plot(st_geometry(zug),
       main = paste("Land Cover in Zug, Switzerland - Year:", year))
  plot(layer, col = c("#228B22", "#8B0000", "#f5ce42", "#9dcff2", "#f8fade"), add = TRUE, legend = FALSE)
  legend("topleft", legend = c("Tree location", "Built Areas", "Crops Location", "Water", "Rangeland"), 
         pch = 20, xpd=NA, cex = 0.6, bg="white", col= c("#228B22", "#8B0000", "#f5ce42", "#9dcff2", "#f8fade"))
  
  img_file <- paste0("zug_cover", year, ".png")
  dev.copy(png, img_file)
  dev.off()
  
  images[[i]] <- image_read(img_file)
}

animation <- image_animate(image_join(images), fps = 2)
image_write(animation, "zug_landcover_animation.gif")


somogy <- rbind(gb_adm1("Hungary")) %>%
  filter(shapeName == "Somogy")

crs <- "EPSG:4326"
lc_list <- list()
for (index in seq_along(virtual_layers)) {
  year <- years[index]
  rasters <- virtual_layers[[index]]
  
  for (i in seq_along(rasters)) {
    raster <- rasters[[i]]
    
    state <- somogy %>%
      sf::st_transform(crs = terra::crs(raster))
    
    land_cover <- raster %>%
      crop(vect(state), snap = "in") %>%
      mask(vect(state)) %>%
      terra::project(crs)
    
    lc_list[[index]] <-land_cover
  }
}

terra::plot(lc_list[[1]])

process_layer <- function(layer) {
  terra::ifel(layer == 2L,
              1,
              terra::ifel(layer == 7L,
                          2,
                          terra::ifel(layer == 5L,
                                      3,
                                      terra::ifel(layer == 1L, 4,
                                                  terra::ifel(layer == 11L, 5, NA)))))
}

land_cover <- lapply(lc_list, process_layer)

images <- list()

for (i in seq_along(land_cover)) {
  year <- years[i]
  layer <- land_cover[[i]]
  
  plot(st_geometry(somogy),
       main = paste("Land Cover in Somogy, Hungary - Year:", year))
  plot(layer, col = c("#228B22", "#8B0000", "#f5ce42", "#9dcff2", "#f8fade"), add = TRUE, legend = FALSE)
  legend("topleft", legend = c("Tree location", "Built Areas", "Crops Location", "Water", "Rangeland"), 
         pch = 20, xpd=NA, cex = 0.6, bg="white", col= c("#228B22", "#8B0000", "#f5ce42", "#9dcff2", "#f8fade"))
  
  img_file <- paste0("somogy_cover", year, ".png")
  dev.copy(png, img_file)
  dev.off()
  
  images[[i]] <- image_read(img_file)
}

animation <- image_animate(image_join(images), fps = 2)
image_write(animation, "somogy_landcover_animation.gif")
