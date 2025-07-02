
# Suponemos que M_UM es la matriz original y contiene las columnas con las especies y una columna 'año'.

# Número de simulaciones
n_sim <- 2000
matriz_z<- NULL
# Inicializamos un listado vacío donde se guardarán los resultados


# Loop de simulaciones
for (sim in 1:n_sim) {
  
  # Paso 1: Dividir la matriz según el año
  años_unicos <- unique(M_UM$año)
  M_UM_aux <- NULL
  
  for (año in años_unicos) {
    
    # Filtramos la matriz por cada año
    M_año <- M_UM[M_UM$año == año, ]
    
    # Paso 2: Aleatorizar las especies presentes manteniendo el número de especies presentes
    especies_presentes <- which(M_año == 1, arr.ind = TRUE)  # Encuentra las especies presentes
    
    # Aleatorizamos las especies presentes
    set.seed(1234)  # Fijamos la semilla para reproducibilidad
    especies_aleatorias <- sample(especies_presentes[, 2])  # Aleatorizamos solo las especies (columna 2)
    
    # Reconstruimos la matriz con las especies aleatorias
    M_año_aleatorio <- M_año
    M_año_aleatorio[especies_presentes] <- 0  # Ponemos todo a cero
    
    # Aquí se asignan correctamente las especies aleatorizadas en las posiciones correspondientes
    for (i in 1:length(especies_aleatorias)) {
      fila <- especies_presentes[i, 1]
      col <- especies_aleatorias[i]
      M_año_aleatorio[fila, col] <- 1
    }
    # Paso 3: Pegamos las matrices
    M_UM_aux <- rbind(M_UM_aux, M_año_aleatorio)
  }
  
  
  M_año_aux <- M_UM_aux %>%
    group_by(año) %>%
    summarise(across(starts_with("Acmella_"):last_col(), sum), .groups = "drop")
  
  M_año_aux<- as.matrix(M_año_aux)
  # Filtra las columnas donde la suma de presencia es 0
  presence_sum_columns <- apply(as.matrix(M_año_aux[, -1]), 2, sum)  # Suma por columna (especies)
  M_año_filtrado <- M_año_aux[, !(presence_sum_columns == 0)]  # Elimina las columnas con presencia 0
  species_common <- intersect(colnames(M_año_filtrado), rownames(traits))
  traits_aux<- traits[species_common,]
  if (!all(colnames(M_año_aux[, -1]) == rownames(traits_aux))) {
    stop("Los nombres de especies entre la matriz de comunidad y los traits no coinciden.")
  }

    D_año <- daisy(traits_aux, metric = "gower", stand = TRUE)
    equit_año_aux <- dbFD(x = D_año, a = as.matrix(M_año_filtrado[, -1]), corr = "cailliez", m = 10)
  #}
  temp_z<- NULL
  for(years in 1:19){
  temp_z.t<- cbind(Año = c(2004+years),
                 FEve = equit_año_aux$FEve[years],
                 FRic = equit_año_aux$FRic[years],
                 FDiv = equit_año_aux$FEve[years],
                 FDis = equit_año_aux$FRic[years],
                 FRao = equit_año_aux$RaoQ[years])
  temp_z<- rbind(temp_z, temp_z.t)}
  
  matriz_z<- rbind(matriz_z, temp_z)
  
  # Eliminamos la matriz auxiliar para el siguiente ciclo
  rm(M_UM_aux)
  print("Vamos por: ", sim, " de: ", n_sim)
}
  
resultados_año <- data.frame(
  año = M_año$año,
  FRic = equit_año$FRic,
  FEve = equit_año$FEve,
  FDis = equit_año$FDis,
  FDiv = equit_año$FDiv,
  RaoQ = equit_año$RaoQ)

Modelo_nulo<- NULL

for(years in unique( matriz_z$Año)){
  t.matriz_z<- as.data.frame(matriz_z[which((as.data.frame(matriz_z)$Año)==(2004+years)),])
  sub<- resultados_año[1,]
  value_z <- c(
    z_FRic = (sub$FRic - mean(t.matriz_z$FRic, na.rm = TRUE)) / sd(t.matriz_z$FRic, na.rm = TRUE),
    z_FEve = (sub$FEve - mean(t.matriz_z$FEve, na.rm = TRUE)) / sd(t.matriz_z$FEve, na.rm = TRUE),
    z_FDis = (sub$FDis - mean(t.matriz_z$FDis, na.rm = TRUE)) / sd(t.matriz_z$FDis, na.rm = TRUE),
    z_FDiv = (sub$FDiv - mean(t.matriz_z$FDiv, na.rm = TRUE)) / sd(t.matriz_z$FDiv, na.rm = TRUE),
    z_FRao = (sub$FRao - mean(t.matriz_z$FRao, na.rm = TRUE)) / sd(t.matriz_z$FRao, na.rm = TRUE)
  )
  
  resultados_año_z<- cbind(sub, t(value_z))
  Modelo_nulo<- rbind(Modelo_nulo, resultados_año_z)
}
  
