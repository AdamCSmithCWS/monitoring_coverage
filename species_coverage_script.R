

# setup -------------------------------------------------------------------
library(ebirdst)
library(tidyverse)
library(bbsBayes2)
library(terra)
library(sf)
library(patchwork)

hex_grid <- readRDS("data/hexagonal_grid_west_hemisphere_120km.rds") %>%
  rename(strata_name = hex_name)
west_hemi <- st_union(hex_grid)
eBird_crs <- readRDS("data/eBird_crs.rds")
west_crop <- west_hemi %>%
  st_buffer(.,10000) %>%
  st_transform(.,crs = eBird_crs) %>%
  st_make_valid()
species_list <- readRDS("data/species_list.rds")
bbs_strata <- load_map("prov_state") %>%
  st_transform(.,crs = st_crs(hex_grid))
re_run <- FALSE
min_p_area <- 0.25 # minimum proportion of a grid cell with species range
# grid cells with < min_p_area are not considered part of the species' range

# species loop ------------------------------------------------------------

pdf("figures/coverage_summary.pdf",
    width = 11,
    height = 8.5)
for(i in rev(1:nrow(species_list))){

  species <- as.character(species_list[i,"english"])
  aou <- as.integer(species_list[i,"aou"])
  if(grepl("(all forms)",species)){
    speciesn <- gsub(" (all forms)",replacement = "",
                     species,fixed = TRUE)
  }else{
    speciesn <- species
  }
down <- try(ebirdst_download_status(speciesn,
                        download_ranges = TRUE,
                        download_abundance = FALSE),
            silent = TRUE)

if(class(down) == "try-error"){
  species_list[i,"eBird_range_data"] <- "unavailable"
  next
  }
species_list[i,"eBird_range_data"] <- "downloaded"

abd_seasonal_range <- load_ranges(species = speciesn, #metric = "mean",
                                  resolution = "27km")  #3km high resolution


season_sel <- "breeding"

if(season_sel %in% abd_seasonal_range$season){
abd_range <- abd_seasonal_range %>%
  filter(season == season_sel)
}else{
  season_sel <- "resident"
  abd_range <- abd_seasonal_range %>%
    filter(season == season_sel)
}



abd_range <- abd_range %>%
  st_make_valid() %>%
  st_crop(.,west_crop) %>%
  st_make_valid() %>%
  st_transform(.,crs = st_crs(hex_grid))

if(file.exists(paste0("output/",
                      aou,"_coverage_surface.rds")) &
   !re_run){

  coverage_alt <- readRDS(paste0("output/",
                                 aou,"_coverage_surface.rds"))

}else{

hex_range <- st_intersection(abd_range,hex_grid) %>%
  st_make_valid()

hex_range <- hex_range %>%
  mutate(area_km2 = as.numeric(st_area(hex_range)/1e6))

strata_w_eBird_range <- hex_range %>%
  select(strata_name) %>%
  mutate(ebird = TRUE)


# Survey data -------------------------------------------------------------


bbs_range_strat <- stratify(species,by = "hexagons",
                            strata_custom = hex_grid,
                            return_omitted = TRUE)

bbs_range <- prepare_data(bbs_range_strat,
                          min_n_routes = 1,
                          min_max_route_years = 1)

strata_w_bbs <- bbs_range$meta_strata %>%
  select(strata_name) %>%
  mutate(bbs = TRUE)



# combine surveyed and range ----------------------------------------------


strata_w_either <- full_join(strata_w_bbs,
                             strata_w_eBird_range,
                             by = "strata_name") %>%
  mutate(ebird = ifelse(is.na(ebird),FALSE,ebird),
         bbs = ifelse(is.na(bbs),FALSE,bbs),
         survey = ifelse(ebird,"ebird","bbs"),
         survey = ifelse(ebird & bbs,"both",survey))



# coverage map and calculation --------------------------------------------

coverage <- hex_grid %>%
  inner_join(strata_w_either,
             by = "strata_name") %>%
  mutate(covered = ifelse(survey %in% c("both","bbs"),
                          TRUE,FALSE)) %>%
  st_make_valid()



# area of range according to eBird
full_hex_w_range <- hex_grid %>%
  rename(area_km2_stratum = area_km2) %>%
  filter(strata_name %in% strata_w_either$strata_name) %>%
  st_make_valid()


full_hex_in_range <- st_intersection(full_hex_w_range,
                                     abd_range) %>%
  select(strata_name,area_km2_stratum) %>%
  st_make_valid()

p_area_in_range <- full_hex_in_range %>%
  mutate(area_km2_in_range = as.numeric(st_area(full_hex_in_range)/1e6),
         range = "inside",
         p_area_in_range = area_km2_in_range/area_km2_stratum) %>%
  st_drop_geometry() %>%
  select(strata_name,range,p_area_in_range)

drop_outside_range <- p_area_in_range %>%
  filter(p_area_in_range < min_p_area)

coverage_alt <- coverage %>%
  filter(survey %in% c("both","bbs") |
           !(survey == "ebird" & strata_name %in% drop_outside_range$strata_name))

}
coverage_sum <- coverage_alt %>%
  st_drop_geometry() %>%
  ungroup() %>%
  group_by(covered) %>%
  summarise(area_km2 = sum(area_km2))
covered <- coverage_sum %>%
  filter(covered) %>%
  select(area_km2) %>%
  unlist()
uncovered <- coverage_sum %>%
  filter(!covered) %>%
  select(area_km2) %>%
  unlist()
cov_overall <- signif(covered/sum(covered,uncovered),3)
cov_cap <- paste("Based on land area of cells covered vs uncovered.",
                 "Assuming that any grid cell with monitoring data \n
                 is completely covered and only grid cells with >",
                 (min_p_area)*100,"% overlap with the species' range
                 map are included in the uncovered range")
tit <- paste(species,"=",cov_overall*100,"%")
stit <- paste(season_sel,"season coverage of BBS surveys with species' data")
xlimits <- st_bbox(coverage_alt)[c("xmin","xmax")]
ylimits <- st_bbox(coverage_alt)[c("ymin","ymax")]

alt_coverage <- ggplot()+
  geom_sf(data = west_hemi)+
  geom_sf(data = bbs_strata)+
  geom_sf(data = coverage_alt,
          aes(fill = covered))+
  geom_sf(data = abd_range,
          fill = NA,
          colour = "red")+
  coord_sf(xlim = xlimits,
           ylim = ylimits)+
  labs(title = tit,
       subtitle = stit,
       caption = cov_cap)+
  scale_fill_viridis_d()+
  theme_bw()
print(alt_coverage)

saveRDS(coverage_alt,file = paste0("output/",
                                   aou,"_coverage_surface.rds"))
}


dev.off()

