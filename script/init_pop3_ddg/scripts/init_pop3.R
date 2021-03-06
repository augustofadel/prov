# Degree-constrained minimum spanning tree [Ning et al., 2008]

dissim <- dist(dat) %>% as.matrix
n_obj <- nrow(dat)

cl <- makeCluster(ncores)
clusterExport(cl, list('init_population', 'diameter', 'dat', 'n_obj', 'pop_size', 'aggr_levels'))

for(data_file in calib_files) {
  # cat(data_file, '\n')
  dat <- readRDS(file.path(data_path, data_file))
  dissim <- dist(dat)
  pop <- lapply(
    as.list(aggr_levels),
    function(x) init_population(
      k = x,
      n_obj = n_obj,
      p = pop_size,
      input = init_mst(dat, dissim, k = x),
      cl = cl
    )
  )
  saveRDS(pop, file.path(main_path, 'results', paste0(data_file, '_pop_mst.rds')))
}

stopCluster(cl)
