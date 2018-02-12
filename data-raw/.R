library(raster)

urbanUS <- readAll(raster('data-raw/urbanUS.tif'))
save(urbanUS, file='data/urbanUS.rda', compress='bzip2')
