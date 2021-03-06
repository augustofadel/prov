# Traveling Salesman Problem
# farthest insertion algorithm [Rosenkrantz et al., 1977]

n_obj <- nrow(dat)

cl <- makeCluster(ncores)
clusterExport(cl, list('init_population', 'dat', 'n_obj', 'pop_size', 'aggr_levels'))

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
      input = init_tsp(dat, dissim),
      cl = cl
    )
  )
  saveRDS(pop, file.path(main_path, 'results', paste0(data_file, '_pop_tsp.rds')))
}

stopCluster(cl)


# Proporcao de solucoes tsp identicas
# tot <- 1000
# teste <- sapply(1:tot, function(x) init_tsp(dat))
# result <- apply(teste, 2, function(x) {
#   apply(teste, 2, function(y) {identical(x, y)})
# })
# sum(apply(result, 1, sum) - 1) / (2 * tot^2)

