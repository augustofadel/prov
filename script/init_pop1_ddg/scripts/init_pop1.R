# k-means with regenaration [Aloise and Araujo, 2015]

# library(magga)
# 
# dat <- readRDS('~/!filesync/instancias/sdc/DS1-200DATA.rds')
# aggr <- 5
# 
# kmeans_ma(dat, aggr) %>% table()


cl <- makeCluster(ncores)
clusterExport(cl, list('kmeans_ma', 'dat', 'aggr_levels'))

for(data_file in calib_files) {
  # cat(data_file, '\n')
  dat <- readRDS(file.path(data_path, data_file))
  dissim <- dist(dat)
  pop <- lapply(as.list(aggr_levels), function(x) {
    parSapply(cl, 1:pop_size, function(y) kmeans_ma(dat, x, dissim = dissim))
  })
  saveRDS(pop, file.path(main_path, 'results', paste0(data_file, '_pop_kmeans.rds')))
}

stopCluster(cl)
