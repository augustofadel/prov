result <- readRDS(file.path(main_path, 'results', 'exp_results.rds'))

cl <- makeCluster(ncores)
clusterExport(cl, 'result')
unexec <- parApply(
  cl, 
  result[, c('t', 'IL1', 'IL2', 'DLD')], 
  1, 
  function(x) {any(is.na(x))}
)
stopCluster(cl)

data_file <- ''

for (i in 1:3) {
  
  param <- result[i,]
  cat('Execution', i, 'from', nrow(result), '\n')
  print(param)
  
  if (data_file != param$dataset) {
    data_file <- param$dataset
    dat <- readRDS(file.path(data_path, paste0(data_file, '.rds'))) %>% 
      as.data.table
    dissim <- dist(dat)
  }
  
  set.seed(param$replica)
  pop <- 
    file.path(main_path, 'results', paste0(data_file, '_pop_', param$init,'.rds')) %>% 
    readRDS %>% 
    lapply(function(x) {x[, sort(sample(1:pop_size, param$pop_size, replace = F))]})
  
  tmp <- try({
    system.time({
      sol <- magga(
        dataset = dat,
        dissim = dissim,
        n_agreg = param$aggr,
        metricas = param$fobj,
        pop = pop,
        init_method = param$init,
        pk = param$pop_size,
        tot_ger = param$generations,
        pe = param$elite,
        pm = param$mutant,
        pr = param$survival,
        nuc = ncores,
        verbose = F
      )
    })
  })
  
  if (class(tmp) == 'try-error') {
    result[i,'err'] <- T
  } else {
    result[i, 't'] <- tmp[[3]]
    aux <- fit(dat, sol[[1]][, 1], c('IL1', 'IL2', 'DLD'))
    result[i, names(aux)] <- aux
  }
  
}