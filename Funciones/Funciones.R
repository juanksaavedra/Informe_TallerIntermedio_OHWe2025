# Funciones -----
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Autor: Juan Carlos Saavedra-Nievas
# email: juank.saavedra@gmail.com
# Fecha: 15-10-2025
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Importar datos netCDF en formato RasterBrick-----
# Función para importar datos netCDF en el tiempo y/o profundidad
# Generando una salida con el objeto RasterBrick y una tabla con los datos
brick_mod <- function(.x, dt.time = TRUE, dt.depth = FALSE, level = 1){
  rb.list <- list()
  case.dt <- case_when(
    dt.time & dt.depth ~ 'all',
    dt.time & !dt.depth ~ 'time',
    TRUE ~ 'depth'
  )
  
  switch (case.dt,
          'all' = {
            rb.list$dph <- brick(path_sst, lvar=4, level=1)
            n_depht <- nlayers(rb.list$dph)
            nm_depth <- getZ(rb.list$dph)
            
            rb.list$tmp[[as.character(nm_depth[1])]] <- brick(path_sst, level=1)
            aux.tmp <- tb_rb(rb.list$tmp[[1]], .time = TRUE) %>% 
              mutate(Depth = nm_depth[1])
            
            for(i in 2:n_depht){
              rb.list$tmp[[as.character(nm_depth[i])]] <- brick(path_sst, level=i)
              aux.tmp <- aux.tmp %>% 
                bind_rows(
                  tb_rb(rb.list$tmp[[i]], .time = TRUE) %>% 
                    mutate(Depth = nm_depth[i])
                )
            }
            aux.tmp <- aux.tmp %>% 
              relocate(Depth, .after = Time)
            salidas <-list(brick.tmp = rb.list$tmp, dt.td = aux.tmp)
            return(salidas)
          },
          'time' = {
            rb.list$dph <- brick(path_sst, lvar=4, level=1)
            n_depht <- nlayers(rb.list$dph)
            nm_depth <- getZ(rb.list$dph)
            
            rb.list$tmp[[as.character(nm_depth[level])]] <- brick(path_sst, level=level)
            aux.tmp <- tb_rb(rb.list$tmp[[1]], .time = TRUE) %>% 
              mutate(Depth = nm_depth[level]) %>% 
              relocate(Depth, .after = Time)
            salidas <-list(brick.tmp = rb.list$tmp, dt.tmp = aux.tmp)
            return(salidas)
            
          },
          'depth' = {
            rb.list$tmp <- brick(path_sst, lvar=3, level=1)
            n_time <- nlayers(rb.list$tmp)
            nm_time <- getZ(rb.list$tmp)
            
            rb.list$dph[[nm_time[level]]] <- brick(path_sst, lvar=4, level=level)
            aux.tmp <- tb_rb(rb.list$dph[[1]], .time = FALSE) %>% 
              mutate(Time = nm_time[level]) %>% 
              relocate(Time, .before = Depth)
            
            salidas <-list(brick.dpht = rb.list$dph, dt.dpht = aux.tmp)
            return(salidas)
            
          }
  )
  invisible()
}

# Genera tabla de datos de objeto RasterBrick-----
tb_rb <- function(.rb, .time=TRUE){
  ncap <- nlayers(.rb)
  x <- coordinates(.rb)[,1] 
  y <- coordinates(.rb)[,2] 
  tb_rb <- tibble(
    # z = ifelse(.time, as.Date(getZ(.rb[[1]])), getZ(.rb[[1]])),
    z = getZ(.rb[[1]]),
    x = x, 
    y = y, 
    sst = getValues(.rb[[1]]))
  for(i in 2:ncap){
    tb_rb <- tb_rb %>% 
      bind_rows(
        tibble(
          # z = ifelse(.time, as.Date(getZ(.rb[[i]])), getZ(.rb[[i]])),
          z = getZ(.rb[[i]]),
          x = x, 
          y = y, 
          sst = getValues(.rb[[i]]))
      )
  }
  if(.time){
    tb_rb <- tb_rb %>%
      mutate(
        Time = as.Date(z),
        tipo = 'RasterBrick') %>% 
      select(Time, x, y, sst, tipo) %>% 
      arrange(x, y)
  }else{
    tb_rb <- tb_rb %>%
      mutate(
        Depth = z,
        tipo = 'RasterBrick') %>% 
      select(Depth, x, y, sst, tipo) %>% 
      arrange(x, y)
  }
  return(tb_rb)
}

# Modifica función formatXTP ----
# Función modificada para darle formato tibble a datos extraídos con extractPts()
# asociada a datos importados desde el paquete satin
formatXTP <- function(pts, plot.error = FALSE){
  error <- pts$d
  if (plot.error == TRUE){
    x11(); hist(error)
  }
  npts <- nrow(pts)
  ncols <- ncol(pts)
  att <- attributes(pts)
  pts <- as_tibble(pts) %>% 
    pivot_longer(cols = starts_with('p'), names_to = 'tp', values_to = 'sst')
  p.var <- pts %>% distinct(tp) %>% pull(tp)
  td <- tibble(expand_grid(Depth = att$attribs$depth, Time = as.Date(att$attribs$period$tmStart)), tp=p.var) %>% 
    relocate(Time, .before = 1)
  ans <- td %>% 
    full_join(pts, join_by(tp)) %>% 
    arrange(Time, Depth, id)
  ans
}

# Modifica función print.satin de paquete satin -----
# Modificación de función que imprime resumen de datos netCDF importados
# Incopora los nombres de las dimensiones impresas Lat, Lon, Time y Depth
print.satin_mod <- function (x, ...)
  {
  X <- list()
  X[["class"]] <- class(x)
  X[["attribs"]] <- x@attribs
  X[["dims"]] <- dim(x@data)
  
  # Modificación incorporada
  names(X$dims) <- c("Lat", "Lon", "Time", "Depth")  
  
  vn <- x@attribs$name
  rng.lon <- range(x@lon)
  rng.lat <- range(x@lat)
  rng.data <- c(range(as.vector(x@data), na.rm = TRUE))
  rng.per <- c(format(min(x@period$tmStart), "%Y-%m-%d"), 
               format(max(x@period$tmEnd), "%Y-%m-%d"))
  ans <- data.frame(lon = rng.lon, lat = rng.lat, rng.data, 
                    period = rng.per)
  row.names(ans) <- c("min", "max")
  names(ans)[3] <- vn
  if (length(x@depth) > 0) {
    rng.dep <- range(x@depth)
    ans$depth <- rng.dep
  }
  X[["ans"]] <- ans
  cat(paste("Object of class ", X[["class"]], "\n", sep = ""))
  cat("\n", "Title:", unlist(X[["attribs"]])[1], "\n", "Long name:", 
      unlist(X[["attribs"]])[2], "\n", "Name:", unlist(X[["attribs"]])[3], 
      "\n", "Units:", unlist(X[["attribs"]])[4], "\n", "Temporal range:", 
      unlist(X[["attribs"]])[5], "\n", "Spatial resolution:", 
      unlist(X[["attribs"]])[6], "\n")
  # cat("\nData dimensions:\n", X[["dims"]], "\n")
  cat("\nData dimensions:\n")
  print(X[["dims"]])
  cat("\nData ranges:\n")
  print(X[["ans"]])
  invisible(X)
  }

# Funcion resumen rasterbrick modificado ----
res.rb <- function(.rb){
  ndepht <- .rb[[2]] %>% reframe(n = n_distinct(Depth)) %>% pull(n)
  ntimes <- .rb[[2]] %>% reframe(n = n_distinct(Time)) %>% pull(n)
  
  if(ndepht==1 | ntimes==1){
    if(ndepht==1){
      .rb[[2]] %>% 
        reframe(across(c(Time,x,y), ~ n_distinct(.x), .names ="{.col}_{'n'}"),
                Depth = unique(Depth),
                min = min(sst, na.rm = TRUE),
                max = max(sst, na.rm = TRUE),
                mean = mean(sst, na.rm = TRUE),
                median = median(sst, na.rm = TRUE),
                q25 = unname(quantile(sst, prob=0.25, na.rm=TRUE)),
                q75 = unname(quantile(sst, prob=0.75, na.rm=TRUE)),
                sd = sd(sst, na.rm = TRUE),
                riq = q75 - q25) %>% 
        relocate(Depth, .before = 1)
      
    }else{
      .rb[[2]] %>% 
        reframe(across(c(Depth,x,y), ~ n_distinct(.x), .names ="{.col}_{'n'}"),
                Time = unique(Time),
                min = min(sst, na.rm = TRUE),
                max = max(sst, na.rm = TRUE),
                mean = mean(sst, na.rm = TRUE),
                median = median(sst, na.rm = TRUE),
                q25 = unname(quantile(sst, prob=0.25, na.rm=TRUE)),
                q75 = unname(quantile(sst, prob=0.75, na.rm=TRUE)),
                sd = sd(sst, na.rm = TRUE),
                riq = q75 - q25) %>% 
        relocate(Time, .before = 1)
    }
  }
  else{
    .rb[[2]] %>% 
      reframe(across(c(Time,Depth,x,y), ~ n_distinct(.x), .names ="{.col}_{'n'}"),
              min = min(sst, na.rm = TRUE),
              max = max(sst, na.rm = TRUE),
              mean = mean(sst, na.rm = TRUE),
              median = median(sst, na.rm = TRUE),
              q25 = unname(quantile(sst, prob=0.25, na.rm=TRUE)),
              q75 = unname(quantile(sst, prob=0.75, na.rm=TRUE)),
              sd = sd(sst, na.rm = TRUE),
              riq = q75 - q25)
  }
  
}
