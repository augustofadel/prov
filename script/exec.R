# install.packages('gtools')

# https://github.com/End-to-end-provenance/RDataTracker

# library(devtools)
# httr::set_config(httr::config(ssl_verifypeer = F))
# install_github('End-to-end-provenance/RDataTracker')
# install_github('augustofadel/magga')
# httr::set_config(httr::config(ssl_verifypeer = T))

library(RDataTracker)
library(profvis)
library(parallel)
library(dplyr)
library(stringr)
library(TSP)
library(data.table)
library(magga)



# ddg.run(
#   r.script.path = NULL,   # the full path to the file containing the R script that is being executed.
#   ddgdir = NULL,          # the directory where the DDG should be saved.
#   overwrite = TRUE,       # Defaults is TRUE, if FALSE, adds timestamp to ddg directory to prevent overwriting
#   f = NULL,               # A function to run. Data provenance is collected within the function.
#   enable.console = TRUE,  # If TRUE, any commands executed in the console, either by typing, copying and pasting, or selecting and running, will result in a procedure node created in the provenance graph, with data nodes created for each variable assigned and data flow edges for variables used and set.
#   annotate.inside = TRUE, # specifies whether automatic annotation of functions and control constructs should be enabled.
#   first.loop = 1,         # The number of the first iteration to be annotated in a for, while, or repeat loop.
#   max.loops = 1,          # The maximum number of times that a for, while, or repeat loop will be annotated. If max.loops is -1, there is no limit. If max.loops = 0, no loops will be annotated.
#   max.snapshot.size = 10, # The maximum size for objects that should be output in snapshot files. If max.snapshot.size is -1, there is no limit. If max.snapshot.size is 0, snapshot nodes are created but no snapshot files are saved.
#   debug = FALSE,          # If TRUE, enable script debugging.
#   save.debug = FALSE,     # If TRUE, save debug files to debug directory.
#   display = FALSE         # If TRUE, display the DDG when the R script completes.
# )
# 
# 
# 
# # granularidade
# ddg.get.detail()
# ddg.set.detail(n = 3) # n = [0, 3]
# 
# ddg.annotate.on(fnames = NULL)
# ddg.annotate.off(fnames = NULL)
# # fnames = A list of one or more function names.
# 
# 
# 
# # breakpoints
# # set a breakpoint at the specified line number of the specified script
# ddg.set.breakpoint(script.name, line.num)
# # display all current breakpoints
# ddg.list.breakpoints()
# # clear all current breakponts
# ddg.clear.breakpoints()




# parameters --------------------------------------------------------------

# main_path <- '~/!filesync/mestrado/18.1_e-science/trabalho'
# data_path <- '~/!filesync/instancias/sdc'
# main_path <- '~/Desktop/filesync/mestrado/18.1_e-science/trabalho'
# data_path <- '~/Desktop/filesync/instancias/sdc'
main_path <- '~/sdc/prov'
data_path <- '~/instancias/sdc'
# files <- list.files(data_path, '.rds')
calib_files <- c(
  'DS1-200DATA.rds',
  'DS2-Uniform400.rds',
  'DS3-censo_Rio de Janeiro.rds',
  'DS4-Tarragona.rds'
)
# calib_files <- 'DS1-200DATA.rds'
aggr_levels <- c(3, 4, 5)
pop_size <- 1000
ncores <- 22



# initial population ------------------------------------------------------

# for (script in c('init_pop1', 'init_pop2', 'init_pop3')) {
  # p <- profvis({
    # source(file.path(main_path, 'script', paste0(script, '.R')))
  # })
  # htmlwidgets::saveWidget(p, file.path(main_path, 'results', paste0(script, '.html')))
  #browseURL('init_pop1.html')
  
  # p <- profvis({
  #   ddg.run(
  #     r.script.path = file.path(main_path, 'script', paste0(script, '.R'))
  #   )
  # })
  # htmlwidgets::saveWidget(p, file.path(main_path, 'results', paste0(script, '_ddg.html')))
# }




# optimization ------------------------------------------------------------

# result <- data.frame(
#   expand.grid(
#     generations = c(50, 100, 200),
#     init = c('kmeans', 'tsp', 'mst'),
#     fobj = c('IL1', 'IL2', 'DLD'),
#     pop_size = c(50, 100, 200),
#     elite = c(0.1, 0.2, 0.3),
#     mutant = c(0.1, 0.2, 0.3),
#     survival = c(0.7, 0.8, 0.9),
#     aggr = aggr_levels,
#     dataset = str_replace(calib_files, '.rds', ''),
#     replica = 1:10,
#     stringsAsFactors = F
#   ),
#   gen_time = NA,
#   total_time = NA,
#   total_time_prov = NA,
#   IL1 = NA,
#   IL2 = NA,
#   DLD = NA,
#   err = F
# )
# saveRDS(result, file.path(main_path, 'results', 'exp_results.rds'))


result <- readRDS(file.path(main_path, 'results', 'exp_results.rds'))

cl <- makeCluster(ncores)
clusterExport(cl, 'result')
unexec <- parApply(
  cl, 
  result[, c('gen_time', 'IL1', 'IL2', 'DLD')], 
  1, 
  function(x) {any(is.na(x))}
  )
stopCluster(cl)

data_file <- ''
max_gen <- max(result$generations)

for (i in 1:nrow(result)) {
  
  if (!unexec[i] | result$err[i] | result$generations[i] != max_gen)
    next
  
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
        save_progress = T,
        verbose = F
      )
    })
  })
  
  if (class(tmp) == 'try-error') {
    result[(i-2):i,'err'] <- T
  } else {
    result[(i-2):i, 'total_time'] <- tmp[[3]]
    for (n_gen in 0:(length(unique(result$generations)) - 1)) {
      # result[i - n_gen, 'gen_time'] <- sum(sol$progress[[1]]$t[1:result[i - n_gen, 'generations']])
      result[i - n_gen, 'gen_time'] <- 
        sol$progress[[1]]$t[result[i - n_gen, 'generations']] - sol$progress[[1]]$t[1]
      aux <- fit(dat, 
                 sol$progress[[1]]$best[, result[i - n_gen, 'generations']], 
                 c('IL1', 'IL2', 'DLD'))
      result[i - n_gen, names(aux)] <- aux
    }
    # aux <- fit(dat, sol$final_pop[[1]][, 1], c('IL1', 'IL2', 'DLD'))
    # result[i, names(aux)] <- aux
  }
  
  saveRDS(result, file.path(main_path, 'results', 'exp_results.rds'))
  
}






# ddg.run(
#   file.path(main_path, 'script', 'teste.R'),
#   display = T
# )

# ddg.display()
