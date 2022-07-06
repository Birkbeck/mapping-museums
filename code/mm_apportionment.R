
source("code/mm_setup.R")

outfolder = 'tmp'

# gen hex 10km grid 
popgrid = load_uk_dataset(paste0('uk_census_2011/uk_population_5kmgrid-2011-sdf.rds'))
hexgrid = load_uk_dataset(paste0('uk_hex_grids/uk_hex_grid_2011-10km.rds'))
hexgridpop = FUN_apportion_geographies(popgrid, hexgrid, 'pop2011')

summary(hexgridpop$pop2011); sum(hexgridpop$pop2011)
summary(popgrid$pop2011); sum(popgrid$pop2011)
saveRDS(hexgridpop, paste0(outfolder,'uk_pop_hex_grid-10km.rds'))

# gen hex 5km grid
popgrid = load_uk_dataset(paste0('uk_census_2011/uk_population_2kmgrid-2011-sdf.rds'))
hexgrid = load_uk_dataset(paste0('uk_hex_grids/uk_hex_grid_2011-5km.rds'))
hexgridpop = FUN_apportion_geographies(popgrid, hexgrid, 'pop2011')

summary(hexgridpop$pop2011); sum(hexgridpop$pop2011)
summary(popgrid$pop2011); sum(popgrid$pop2011)
saveRDS(hexgridpop, paste0(outfolder,'uk_pop_hex_grid-5km.rds'))

#write_geojson(spTransform(hexgridpop,ll_crs),'../tmp/new_geog.geojson','new_geo')
#write_geojson(spTransform(popgrid,ll_crs),'../tmp/old_geog.geojson','old_geo')

rm(popgrid,hexgrid,hexgridpop)