

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
min_nyears <- 2 # minimum number of years of species observations on at least one BBS route
# within a grid cell to consider that cell "covered"
selected_species <- TRUE

if(selected_species){

species_sel <- c("American Robin",
                 "European Starling",
                 "American Kestrel",
                 "Nelson's Sparrow",
                 "Western Grebe",
                 "Curve-billed Thrasher",
                 "Caspian Tern",
                 "Lesser Yellowlegs",
                 "Peregrine Falcon",
                 "Connecticut Warbler")

species_list2 <- species_list %>%
  filter(english %in% species_sel)


}


pdf_title <- ifelse(selected_species,
                    "coverage_selection.pdf",
                    "coverage_summary.pdf")


# species loop ------------------------------------------------------------

pdf(paste0("figures/",pdf_title),
    width = 11,
    height = 8.5)
for(i in rev(1:nrow(species_list2))){

  species <- as.character(species_list2[i,"english"])
  aou <- as.integer(species_list2[i,"aou"])
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
  species_list2[i,"eBird_range_data"] <- "unavailable"
  next
  }
species_list2[i,"eBird_range_data"] <- "downloaded"

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
                          min_max_route_years = min_nyears)

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

canada <- bbs_strata %>%
  filter(country == "Canada") %>%
  select(country)
coverage_sum_can <- coverage_alt %>%
  st_join(.,canada,
          join = st_intersects,
          left = FALSE) %>%
  st_drop_geometry() %>%
  ungroup() %>%
  group_by(covered) %>%
  summarise(area_km2 = sum(area_km2))

covered_can <- coverage_sum_can %>%
  filter(covered) %>%
  select(area_km2) %>%
  unlist()
uncovered_can <- coverage_sum_can %>%
  filter(!covered) %>%
  select(area_km2) %>%
  unlist()

US <- bbs_strata %>%
  filter(country == "United States of America") %>%
  select(country)
coverage_sum_us <- coverage_alt %>%
  st_join(.,US,
          join = st_intersects,
          left = FALSE) %>%
  st_drop_geometry() %>%
  ungroup() %>%
  group_by(covered) %>%
  summarise(area_km2 = sum(area_km2))
covered_us <- coverage_sum_us %>%
  filter(covered) %>%
  select(area_km2) %>%
  unlist()
uncovered_us <- coverage_sum_us %>%
  filter(!covered) %>%
  select(area_km2) %>%
  unlist()
cov_overall <- signif(covered/sum(covered,uncovered),3)
cov_overall_us <- signif(covered_us/sum(covered_us,uncovered_us),3)
cov_overall_can <- signif(covered_can/sum(covered_can,uncovered_can),3)
cov_cap <- paste("Based on western hemisphere land area of 12,000km^2
hexagonal grid cells. Assuming that any grid cell
                with monitoring data is completely covered.
                Range is both the land area of all grid cells with >",
                 (min_p_area)*100,"% overlap
                 with the species' range layer from eBird, and all grid cells with
                 monitoring data for > ", min_nyears," years.")
tit <- paste(species,"=",cov_overall*100,"% overall")
stit <- paste0(str_to_title(season_sel)," season coverage of BBS surveys with \n at least 2 years of species observations \n",
              cov_overall_us*100,"% coverage of range in in USA \n",cov_overall_can*100,"% coverage of range in Canada")
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

if(!re_run){
saveRDS(coverage_alt,file = paste0("output/",
                                   aou,"_coverage_surface.rds"))
}
}


dev.off()

