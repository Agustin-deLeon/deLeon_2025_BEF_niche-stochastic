########################################################
########################################################
################  de León et al 2025  ##################
########################################################
########################################################

########################################################
################# CATS per year #########################
########################################################

# Set directory
setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")

### Import databases 
library(readxl)
br<-read_excel("br.xlsx")                       # br is the matrix of species occurrences
br<-matrix(unlist(br),byrow = F, nrow=nrow(br), dimnames = list(rownames(br),colnames(br)))       # transform to matrix

bm.nuevo=NULL
for(i in 2005:2024){
  as.matrix(br[which(br[,1]==i),]) -> bm.temp
  if(i==2006){bm.temp[which(bm.temp[,2]==6),]-> bm.temp}
  if(i==2007){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2008){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2009){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  if(i==2010){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  rbind(bm.nuevo, bm.temp)-> bm.nuevo
}
as.matrix(bm.nuevo)-> br


traits<- read.csv("traits_numerico.csv")
traits[-c(15,35,119),-1]-> traits #when i import it it adds an extra col, i eliminate it

# Edit the format of trait variables (factors)
names.traits<-traits[,1]
traits<-matrix(as.numeric(unlist(traits)),byrow = F, 
               nrow=nrow(traits), dimnames = list(rownames(traits),
                                                  colnames(traits)))
rownames(traits) <- as.vector(unlist(names.traits))
traits <- traits[,-1]
traits<-as.data.frame(traits)
sapply(traits, class) ### ok, numeric
### pass categorical variables to factors
# colnames(traits)
cols <- c(1,3, 4, 5, 7, 11, 12, 13, 14, 15, 16,17) ## factors
cols.ord <- c(2, 6, 8, 18) ## ordered factors :  "Seed size"  "Leaf size"  "Sculthorpe" "Raunkiaer" 
colnames(traits[,cols.ord]); 
colnames(traits[,cols]) ## check
for(i in cols.ord){
  traits[,i] <- ordered(traits[,i])
} 
for(i in cols){
  traits[,i] <- factor(traits[,i])
} 

######################### PCoA ############################
library(ape)
library(cluster)
library(scales)
library(vegan)

### Gower distances
### computing distances among species based on their functional traits
spp_plants<-rownames(traits)
# spp_plants_viejo <- rownames(traits.pcoa)
# rownames(traits)<-spp_plants #### HAY QUE HACER LISTA CORTA DE ESTOS NOMBRES Y TEXTO PARA PONER EN LA FIG
rownames(traits) <- 1:120
D <- daisy(traits, metric = "gower", stand = TRUE) ### species x species distances (100 x 100)
# D1 <- daisy(traits, metric = "gower")
# identical(D, D1); rm(D1) ## ok
PCoA<- pcoa(D, rn = rownames(traits))
traits.pcoa.ejes<-PCoA$vectors
pcoa_axes <- as.data.frame(PCoA$vectors[, 1:5])  # Ejes 1 a 5
colnames(pcoa_axes) <- c("Axis1", "Axis2", "Axis3", "Axis4", "Axis5")
rownames(pcoa_axes)<- names.traits
proporciones <- PCoA$values$Relative_eig
porcentaje_acumulado <- cumsum(proporciones) * 100
print(proporciones)
print(porcentaje_acumulado)


library(glmmTMB)
library(MuMIn)
library(performance)
out.final=NULL
M<- br[,-c(19,39,121)] #Delate spp with 0 occurrences in selected sampling periods
apply(br[,c(19,39,121)],2,sum) #Check 0 anundance: Cirsium_vulgare, Enydra_sessilis, Stellaria_media 
years<-sort(unique(M[,1]))                                       # Ponds observed 
metacomm<-apply(M[,5:ncol(M)],2,sum)
#  id.spp<-match(colnames(M[,5:ncol(M)]),rownames(Traits))
#  Traits<-Traits[id.spp,]
#  metacomm<-metacomm[ii.observed]           
Traits<- traits.pcoa.ejes[,1:5] # reduce pool to species observed on this year
#  Traits<-Traits[ii.observed,]                                     # reduce traits matrix to observed species
it= 200 # lo ideal, lo estoy corriendo con 10 porque demora jeje
for(i in years){ # }, .combine=rbind)%dopar%  {
  charco<-as.numeric(paste(M[which(M[,1]==i),3]))                    # Create a vector indicating survey times
  metacomm<-apply(M[,5:ncol(M)],2,sum)                             # Metacommunity abundances: uso el global, algo a considerar
  metacomm.U<-rep(sum(metacomm)/length(metacomm),length(metacomm)) # Uniform expectation
  print(i)
  out<-NULL
  sigma.obs.nt<-NULL                                               # Null objects for deviance decomposition...
  sigma.obs.ut<-NULL                                               # 
  sigma.null.np<-NULL                                              # 
  sigma.obs.up<-NULL                                               # 
  m<-M[which(M[,1]==i),]                                         # pond matrix
  ss<-charco                                 # sampling times
  m.temp<-NULL                                                   # temporal out
  # if more than a sample unit was observed
  for(j in 1:nrow(m)){                                         # for each sample unit
    cat("row ",c(j, "of year",i),"\n")    
    y<-m[j,5:ncol(m)]
    m.temp.2<-cbind(metacomm,ss[j], Traits, metacomm.U,y)
    rbind(m.temp, m.temp.2)-> m.temp
  }
  #m.temp[1:2214,]-> m.temp
  m.temp <- na.omit(m.temp) 
  colnames(m.temp)[2]<-"charco.id"
  colnames(m.temp)[ncol(m.temp)]<-"y"
  colnames(m.temp)[3:7]<- c("Axis1" ,    "Axis2"    , "Axis3"    , "Axis4"   ,  "Axis5")
  as.data.frame(m.temp)-> m.temp
  # Variables defined in the space bcecause glmmTMB can not use "data=", SI DEJA, solo que tiene que ser un data.frame      
  rm(Axis1, Axis2, Axis3, Axis4, metacomm, metacomm.U)
  attach(m.temp)
  
  print("ok?")
  ######### FIRST MODEL
  # neutral mass effect and observed traits 
  M1.n.t<-glmmTMB(y~Axis1+Axis2+Axis3+Axis4+Axis5+
                    I(Axis1^2)+I(Axis2^2)+I(Axis3^2)+I(Axis4^2)+I(Axis5^2)+(1|charco.id),
                  offset=log(metacomm),
                  data=m.temp, family=binomial)
  print("ok?")
  
  M1.n.t_r2<- r2(M1.n.t)
  r2.n.t<-as.numeric(M1.n.t_r2$R2_marginal)
  # if(adjustr2.to.random==T)if(is.na(M1.n.t_r2$R2_conditional)==F)r2.n.t<-r2.n.t/(1-(M1.n.t_r2$R2_conditional-r2.n.t))
  sel.coeff<-as.vector(unlist(fixef(M1.n.t)$cond))
  print("model1 ok")
  
  ##### SECOND MODEL
  # neutral mass effect and observed traits 
  M1.u.t<-glmmTMB(y~Axis1+Axis2+Axis3+Axis4+Axis5+
                    I(Axis1^2)+I(Axis2^2)+I(Axis3^2)+I(Axis4^2)+I(Axis5^2)+(1|charco.id),
                  offset=log(metacomm.U),
                  data=m.temp, family=binomial)
  sigma.obs.ut<-c(sigma.obs.ut,M1.u.t$sigma) #desviación estandar? para qué la quiere
  #M0.u.t<-glmmTMB(y~1 +(1|year.month),
  #                offset=log(metacomm.U),
  #                data=NULL, family=binomial)
  M1.u.t_r2<-r2(M1.u.t)
  r2.u.t<-as.numeric(M1.u.t_r2$R2_marginal)
  #  if(adjustr2.to.random==T)if(is.na(M1.u.t_r2$R2_conditional)==F)r2.u.t<-r2.u.t/(1-(M1.u.t_r2$R2_conditional-r2.u.t))
  print("model2 ok")
  
  r2.n.p<-NULL #no se pierden si se borran ya acá? ok no, estos son p de permutados (los next)
  r2.u.p<-NULL
  ############
  #####. Randomized traits
  for(rand in 1:it){ #cuantas veces aleatorizo los traits creo que sería
    print(rand)
    y.p<-m.temp[sample(nrow(m.temp)),ncol(m.temp)]
    m.temp.p<-cbind(m.temp[,-ncol(m.temp)],y.p)
    
    #######THIRD MODEL        
    # neutral mass effect and permuted traits         
    M1.n.p<-glmmTMB(y.p~Axis1+Axis2+Axis3+Axis4+Axis5+
                      I(Axis1^2)+I(Axis2^2)+I(Axis3^2)+I(Axis4^2)+I(Axis5^2)+(1|charco.id),
                    offset=log(metacomm),
                    data=NULL, family=binomial)
    M1.n.p_r2<-r2(M1.n.p)
    #r.squaredGLMM(M1.n.p)[1,]-> M1.n.p_r2
    r2.n.p.temp<-as.numeric(M1.n.p_r2$R2_marginal)
    #r2.n.p.temp<-as.numeric(M1.n.p_r2[1])
    # if(adjustr2.to.random==T) if(is.na(M1.n.p_r2$R2_conditional)==F)r2.n.p.temp<-r2.n.p.temp/(1-(M1.n.p_r2$R2_conditional-r2.n.p.temp))
    #   if(adjustr2.to.random==T) if(is.na(M1.n.p_r2[2])==F)r2.n.p.temp<-r2.n.p.temp/(1-(M1.n.p_r2[2]-r2.n.p.temp))
    ######### FOURTH MODEL
    # Uniform metacomm and permuted traits      
    M1.u.p<-glmmTMB(y.p~Axis1+Axis2+Axis3+Axis4+Axis5+
                      I(Axis1^2)+I(Axis2^2)+I(Axis3^2)+I(Axis4^2)+I(Axis5^2)+(1|charco.id),
                    offset=log(metacomm.U),
                    data=NULL, family=binomial)
    M1.u.p_r2<-r2(M1.u.p)
    #r.squaredGLMM(M1.u.p)[1,]-> M1.u.p_r2
    r2.u.p.temp<-as.numeric(M1.u.p_r2$R2_marginal)
    #r2.u.p.temp<-as.numeric(M1.u.p_r2[1])
    #if(adjustr2.to.random==T)if(is.na(M1.u.p_r2$R2_conditional)==F)r2.u.p.temp<-r2.u.p.temp/(1-(M1.u.p_r2$R2_conditional-r2.u.p.temp))
    #   if(adjustr2.to.random==T) if(is.na(M1.u.p_r2[2])==F)r2.u.p.temp<-r2.u.p.temp/(1-(M1.u.p_r2[2]-r2.u.p.temp))
    r2.n.p<-c(r2.n.p,r2.n.p.temp)
    r2.u.p<-c(r2.u.p,r2.u.p.temp)
  }
  r2.n.p<-mean(r2.n.p) #ok, me quedo con la media de todas las veces que lo corrí
  r2.u.p<-mean(r2.u.p)
  
  # Following Shipley 2014 R2 values are adjusted
  r2.u.t<-max(r2.u.t,r2.u.p)  # R2 with observed traits could not be inferior to R2 with randomized traits
  r2.n.p<-max(r2.n.p,r2.u.p)  # R2 with observed metacomm could not be inferior to R2 with Uniform metacomm
  r2.n.t<-max(r2.n.t,r2.u.t)   # metacommunity prior and observed trait adjusted. It should not be less than R2.u.t.
  out<-rbind(out,c(as.numeric(i),r2.n.t,r2.u.t,r2.n.p,r2.u.p, as.vector(sel.coeff)))
  
  
  colnames(out)<-c("year","r2.n.t","r2.u.t","r2.n.p","r2.u.p", "intercept","axis1", "axis2", "axis3", "axis4", "axis5", "axis1_2", "axis2_2", "axis3_2", "axis4_2", "axis5_2") #VER bien que son estos valores
  out-> out.temp
  out.final<- rbind(out.final, out.temp)
  print(date())
}

as.data.frame(out.final)-> nn
nn.1<-cbind(nn[,1:5],(1-(nn[,2]))/(1-nn[,5]),nn[,(6:ncol(nn))])
colnames(nn.1)[6]<-"R2.unexp"        # Estimates DRIFT corrected by expected drift (R2) from randomized trait 
head(nn.1)
nn.1-> CATS_salida



####################################
# Armar las matrices a utilizar
####################################

setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
library(openxlsx)
library(glmmTMB)
library(quantreg)  
library(ggplot2) 
library(car)
library(visreg)
library(MuMIn) 
library(ggpubr)
library(gridExtra)
library(caret)
library(broom.mixed)
library(ggeffects)
library(DHARMa)
library(bestglm)
read.xlsx("bm.xlsx")-> bm
read.xlsx("br.xlsx")-> br2

read.csv("lluvias.csv")-> lluvias
bm[,7:10]<- c(rep(1, nrow(bm)))
#for(i in unique(lluvias[,1])){
#  bm[which(bm[,1]==i),7:10]<- lluvias[which(lluvias[,1]==i),2:5]
#}
#colnames(bm)
colnames(bm)<-  c("año","mes","charco","um","riq","biom","lluvia_muestreo","lluvia_anual","temp_muestreo", "temp_anual")  

read.xlsx("Environmental.xlsx")-> Environment
bm[,11:14]<- c(rep(1, nrow(bm)))
for (i in unique(Environment[,1])) {
  print(i)
  bm_subset <- bm[bm[,3] == i, 11:14]
  env_subset <- Environment[Environment[,1] == i, c(1:18)]
  bm[bm[,3] == i, 11:28]<- env_subset
}
colnames(bm)
colnames(bm)<-  c("año","mes","charco","um","riq","biom","lluvia_muestreo","lluvia_anual","temp_muestreo", "temp_anual",
                  "Pond.id"  ,      "X"           ,    "Y"       ,        "DM"   ,           "ddmm"     ,      
                  "Shape"     ,      "Islands"    ,     "log.Area"     ,   "log.Volumen"    , "Mean.Depth"  ,   
                  "Sd.Depth"    ,    "CV.Depth"   ,     "Hydroperiod"   ,  "Degree"      ,    "log.Betweenness",
                  "Closenness"   ,   "grado.percol" ,   "clos.percol"   )  
#Keep sampling events of interest 
bm.nuevo=NULL
for(i in 2005:2023){
  as.matrix(bm[which(bm[,1]==i),]) -> bm.temp
  if(i==2006){bm.temp[which(bm.temp[,2]==6),]-> bm.temp}
  if(i==2007){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2008){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2009){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  if(i==2010){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  rbind(bm.nuevo, bm.temp)-> bm.nuevo}

as.data.frame(bm.nuevo)-> bm.nuevo
ambientales<- read.xlsx("Ambientales_resumen.xlsx")
ambientales[,c(1:3, 9:16)]-> ambientales 
amb.nuevo=NULL
for(i in 2005:2023){
  as.matrix(ambientales[which(ambientales[,1]==i),]) -> amb.temp
  if(i==2006){amb.temp[which(amb.temp[,2]==6),]-> amb.temp}
  if(i==2007){amb.temp[which(amb.temp[,2]==8),]-> amb.temp}
  if(i==2008){amb.temp[which(amb.temp[,2]==8),]-> amb.temp}
  if(i==2009){amb.temp[which(amb.temp[,2]==7),]-> amb.temp}
  if(i==2010){amb.temp[which(amb.temp[,2]==7),]-> amb.temp}
  rbind(amb.nuevo, amb.temp)-> amb.nuevo
  amb.temp=NULL}
rm(amb.temp, ambientales)
as.data.frame(amb.nuevo)-> amb.nuevo

bm.nuevo[,29:36]<- c(rep(1, nrow(bm.nuevo)))
bm.act=NULL
for (i in unique(amb.nuevo[,1])) {
  print(i)
  bm.temp<- bm.nuevo[which(bm.nuevo[,1]==i),] #divido para cada año
  amb.temp<- amb.nuevo[which(amb.nuevo[,1] == i),]
  for(ii in unique(bm.temp[,3])){ #ii es id de cada charco de ese año
    print(ii)
    bm_subset <- bm.temp[which(bm.temp[,3] == ii), 29:36]
    env_subset <- amb.temp[which(amb.temp[,3] == ii), c(4:11)]
    if(length(which(amb.temp[,3] == ii))>0){
      bm.temp[which(bm.temp[,3] == ii), 29:36]<- env_subset } 
    else{ bm.temp[which(bm.temp[,3] == ii), 29:36]<- NA }} 
  rbind(bm.act, bm.temp)-> bm.act }
colnames(bm.act)
colnames(bm.act)<-  c("año","mes","charco","um","riq","biom","lluvia_muestreo","lluvia_anual","temp_muestreo", "temp_anual",
                      "Pond.id"  ,      "X"           ,    "Y"       ,        "DM"   ,           "ddmm"     ,      
                      "Shape"     ,      "Islands"    ,     "log.Area"     ,   "log.Volumen"    , "Mean.Depth"  ,   
                      "Sd.Depth"    ,    "CV.Depth"   ,     "Hydroperiod"   ,  "Degree"      ,    "log.Betweenness",
                      "Closenness"   ,   "grado.percol" ,   "clos.percol" , 
                      "cortes"   ,  "prof.media", "sd.prof" ,   "CVprof"   ,  "area"  ,
                      "log.area" ,  "vol",   "log.vol"   )  
bm.nuevo.repuesto<- bm.nuevo
bm.nuevo<- bm.act
bm.nuevo.temp<- bm.nuevo 

bm.st.total=NULL
for(i in 2005:2023){
  as.matrix(bm[which(bm[,1]==i),]) -> bm.temp
  if(i==2006){bm.temp[which(bm.temp[,2]==6),]-> bm.temp}
  if(i==2007){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2008){bm.temp[which(bm.temp[,2]==8),]-> bm.temp}
  if(i==2009){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  if(i==2010){bm.temp[which(bm.temp[,2]==7),]-> bm.temp}
  
  for(ii in unique(bm.temp[,3])){
    sum(bm.temp[which(bm.temp[,3]==ii),6])-> biomasa
    max(bm.temp[which(bm.temp[,3]==ii),4])-> um
    (biomasa/um)-> biomasa.st
    c(bm.temp[max(which(bm.temp[,3]==ii)),1:3], um , biomasa, biomasa.st)-> bm.st
    rbind(bm.st.total, bm.st)-> bm.st.total
  }
}

br.muestreos =NULL
for(i in 2005:2023){
  as.matrix(br2[which(br2[,1]==i),]) -> br2.temp
  if(i==2006){br2.temp[which(br2.temp[,2]==6),]-> br2.temp}
  if(i==2007){br2.temp[which(br2.temp[,2]==8),]-> br2.temp}
  if(i==2008){br2.temp[which(br2.temp[,2]==8),]-> br2.temp}
  if(i==2009){br2.temp[which(br2.temp[,2]==7),]-> br2.temp}
  if(i==2010){br2.temp[which(br2.temp[,2]==7),]-> br2.temp}
  rbind(br.muestreos, br2.temp)-> br.muestreos
}

charcos.evaluar <- c("1", "2", "3", "4", "5", "6", "7", "8", 
                     "10", "11", "13", "14", "15", "16", "17", "12", 
                     "21", "24", "25", "26", "27", "29", "30", "32", 
                     "33", "38", "40", "41", "42", "43", "44", "45", 
                     "47", "48", "49", "50", "51", "55", "56", "666", 
                     "91", "10022", "21022", "134", "137")

#I replace the missing data of the environmental variables with the average of that puddle for the total of the samples taken.
bm.nuevo$cortes <- ifelse(is.na(bm.nuevo$cortes), bm.nuevo$Islands, bm.nuevo$cortes)
bm.nuevo$prof.media <- ifelse(is.na(bm.nuevo$prof.media), bm.nuevo$Mean.Depth, bm.nuevo$prof.media)
bm.nuevo$CVprof <- ifelse(is.na(bm.nuevo$CVprof), bm.nuevo$CV.Depth, bm.nuevo$CVprof)
bm.nuevo$log.area <- ifelse(is.na(bm.nuevo$log.area), bm.nuevo$log.Area, bm.nuevo$log.area)
bm.nuevo$sd.prof <- ifelse(is.na(bm.nuevo$sd.prof), bm.nuevo$Sd.Depth, bm.nuevo$sd.prof)


profundidad<- NULL
for(ii in 2005:2023){
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:16,29,34,30:32,23:25)]
  colnames(df.um.ch) <- c("um.año", "ch.id", "um.biom", "um.rich", "lluvia_muestreo", "lluvia_anual", "temp_muestreo", 
                          "temp_anual", "DM", "ddmm", "Shape", "Islands", "log.Area", "Mean.Depth", 
                          "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree", "log.Betweenness")
  df.um.ch <- na.omit(df.um.ch)
  mean(df.um.ch$Mean.Depth)-> media_prof
  sd(df.um.ch$Mean.Depth)-> sd_prof
  (media_prof)/(sd_prof)-> media_cv
  cbind(ii, media_prof, media_cv)-> prof.temp
  rbind(profundidad, prof.temp)-> profundidad
}


####################################
# What is each one?:
# bm is the matrix with biomasses, bm new is selecting the sampling months I want to retain
# br2 is the matrix of species presence by um, br.sampling is selecting the months I want to retain
# puddles.evaluate has the most consistent puddles over time, which I can use to study them
# depth has the average depth and cv of the average depth for each year
####################################

#NAs in each one?
n.datos=NULL
for(ii in 2005:2023){
  bm.temp <- bm.nuevo[which(bm.nuevo[,1] == ii),]
  nrow(bm.temp)-> datos
  nrow(na.omit(bm.temp))-> datos.na
  cbind(ii, datos, datos.na) -> datos.rest
  n.datos<- rbind(n.datos, datos.rest)
}
print(n.datos)

####################################
####################################
####   PARTE 0: Neutrality vs   ####
####      Richness & Biomass    ####
####################################
####################################
setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]

as.data.frame(nn.1)->nn.1
nn.1$R2.unexp-> R2.unexp
par(mfrow=c(1,1), mar=c(5,5,2,2))
ggplot(nn.1, aes(x = year, y = R2.unexp)) +
  geom_point(color = "darkgreen", alpha = 0.6, size = 5) +
  geom_line(color = "darkgreen", alpha = 0.6, size = 1) +   
  theme_pubr() +
  labs(y = "Neutrality", x = " ") +
  theme(
    axis.title.y = element_text(size = 35),  
    axis.text.y = element_text(size = 20),   
    axis.text.x = element_text(size = 20)    
  )

library(dplyr)
###############


read.csv("lluvias.csv")->lluvias
lluvias[-c(1:3),-1]->lluvias
prof_media <- (aggregate(prof.media ~ año, data = bm.nuevo, mean))$prof.media
prof_CV <- (aggregate(CVprof ~ año, data = bm.nuevo, mean))$CVprof
cbind(año = nn.1[,-2], lluvias, prof_media, prof_CV, 
      R2.unexp=nn.1[,2])-> nn.1

##################
## Standarizing richness with iNEXT
#################
as.data.frame(br.muestreos)-> br.muestreos_df
# Cargar paquetes
library(dplyr)
library(iNEXT)

M<- as.data.frame(br.muestreos_df[,-c(2:4)])
M_año <- M %>%
  group_by(año) %>%
  summarise(across(starts_with("Acmella_"):last_col(), sum), .groups = "drop")

comunidades_por_año <- split(M_año[, -1], M_año$año) %>%
  lapply(as.numeric)  

riqueza_año <- iNEXT(comunidades_por_año, q = 0, datatype = "abundance")
riqueza_resultados <- riqueza_año$AsyEst %>%
  filter(Diversity == "Species richness") %>%
  select(Assemblage, Estimator) %>%
  rename(año = Assemblage, Riqueza = Estimator)

riq_neutral <- cbind(nn.1, riqueza = riqueza_resultados$Riqueza)
attach(riq_neutral)
bestglm(as.data.frame(na.omit(riq_neutral[,-c(1)])), family = gaussian, IC = "AIC", nvmax=2)$BestModels
#The one with only neutrality does not differ by more than 2 from the one with other covariates (AIC)
summary(lm(round(riqueza)~R2.unexp))
p<-coefficients(lm(round(riqueza)~R2.unexp))
par(mfrow=c(1,1), mar=c(5,8,1,1))
plot(round(riqueza) ~ R2.unexp, 
     ylab = "Richness", 
     xlab = "Neutrality", 
     bty = "l", 
     pch = 19, 
     col = "burlywood3",
     cex.lab = 2,   
     cex.axis = 1.3)  
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
text(x=0.12,y=93, "F(2,17): 4.3\nP: 0.05\nR2: 0.16", cex=1, pos=4)

#########################
# Neutrality vs biomass
#########################

um_por_año <- M %>%
  group_by(año) %>%
  summarise(n_UM = n())
bm.nuevo_df<- as.data.frame(bm.nuevo[,1:6])
na.omit(bm.nuevo_df)-> bm.nuevo_df
biomasa_total_anual <- bm.nuevo_df %>%
  group_by(año) %>%
  summarise(BiomasaTotal = sum(biom))
biomasa_estandarizada <- biomasa_total_anual %>%
  left_join(um_por_año, by = "año") %>%
  mutate(Biomasa_Estandarizada = BiomasaTotal / n_UM)
biomasa_estandarizada<- as.matrix(biomasa_estandarizada[,4])
biom_neutral <- cbind(nn.1, biomasa = biomasa_estandarizada)
colnames(biom_neutral)[9]<- "biom"
bestglm(as.data.frame(na.omit(biom_neutral[,-c(1)])), family = gaussian, IC = "AIC", nvmax=2)$BestModels
#The one with the lowest AIC is in fact the one with neutrality and annual rainfall
attach(biom_neutral)
summary(lm(biom~R2.unexp+lluvia_anual))
p<-coefficients(lm(biom~R2.unexp+lluvia_anual))
par(mfrow=c(1,2), mar=c(5,8,1,1))
Biom.neutral<-biom-p[3]*lluvia_anual+mean(p[3]*lluvia_anual)
plot(Biom.neutral ~ R2.unexp, ylab = "Biomass", xlab = "Neutrality", 
     bty = "l", pch = 19, col = "burlywood3",cex.lab = 1.8, cex.axis = 1.3)  # Agranda los números de los eje
curve((p[1]+mean(p[3]*lluvia_anual)+p[2]*x), add=T, col="red",lwd=3)
text(x=0.5,y=7, "F(2,16): 9.8\nP: 0.002\nR2: 0.5", cex=1, pos=4)
Biom.precipitation<-biom-p[2]*R2.unexp+mean(p[2]*R2.unexp)
plot(Biom.precipitation~lluvia_anual, xlab="Anual rain", bty="l", pch=19, col="navy", ylab = " ",
     cex.lab = 1.8, cex.axis = 1.3)
points(Biom.precipitation~lluvia_anual, cex=.6, pch=19, col="cornflowerblue")
curve((p[1]+p[3]*x+mean(p[2]*R2.unexp)), add=T, col="darkgreen",lwd=3)







####################################
####################################
###### PARTE 1:  BEF slopes ########
####################################
####################################

####################################
#Lineal
par(mfrow=c(5,4), mar=c(5,5,4,3))
pendientes.lineal <- list()
pendientes.charcos<- list()
for(ii in 2005:2023){
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]  
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:26)]
  colnames(df.um.ch) <- c("um.año", "ch.id", "um.biom", "um.rich", "lluvia_muestreo", "lluvia_anual", "temp_muestreo", "temp_anual",
                          "DM", "ddmm", "Shape", "Islands", "log.Area", "log.Volumen", "Mean.Depth", "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree", "log.Betweenness", "Closenness")
  df.um.ch <- na.omit(df.um.ch)
  df.um.ch$ch.id <- as.factor(df.um.ch$ch.id)
  titulo <- paste(ii)
  m1 <- glmmTMB(um.biom ~ um.rich + (um.rich | ch.id),  
                family = gaussian(), data = df.um.ch)
  fixef <- fixef(m1)$cond
  ranef_intercepts <- ranef(m1)$cond$ch.id[, "(Intercept)"]
  ranef_slopes <- ranef(m1)$cond$ch.id[, "um.rich"]
  pendientes.charcos[[as.character(ii)]] <- ranef_slopes
  unique_ids <- levels(df.um.ch$ch.id)
  colors <- hcl.colors(length(unique_ids), palette = "Pastel1")
  names(colors) <- unique_ids
  plot(df.um.ch$um.rich, df.um.ch$um.biom,  
       col = colors[as.factor(df.um.ch$ch.id)], pch = 19,  
       xlab = "Richness", ylab = "Biomass",  
       main = titulo,
       cex.lab = 2.5,  
       cex.axis = 2, 
       cex.main = 1.7)
  for (i in 1:length(ranef_intercepts)) {
    current_data <- df.um.ch[df.um.ch$ch.id == unique_ids[i],]
    intercept <- fixef["(Intercept)"] + ranef_intercepts[i]
    slope <- fixef["um.rich"] + ranef_slopes[i]
    x_min <- min(current_data$um.rich)
    x_max <- max(current_data$um.rich)
    curve(intercept + slope * x,  
          from = x_min, to = x_max, 
          add = TRUE, col = colors[unique_ids[i]], lwd = 1, lty = 2)
  }
  
  general_intercept <- fixef["(Intercept)"]
  general_slope <- fixef["um.rich"]
  curve(general_intercept + general_slope * x,  
        from = min(df.um.ch$um.rich), to = max(df.um.ch$um.rich),
        add = TRUE, col = "black", lwd = 3)
  
  resumen <- summary(m1)
  sigma.temp <- resumen$sigma^2
  stddev_intercept <- attr(resumen$varcor$cond$ch.id, "stddev")["(Intercept)"]
  stddev_umrich <- attr(resumen$varcor$cond$ch.id, "stddev")["um.rich"]
  
  pend.temp <- cbind(t(as.matrix(fixef)), sigma.temp, r.squaredGLMM(m1)[1], r.squaredGLMM(m1)[2], stddev_intercept, stddev_umrich, resumen$coefficients$cond["um.rich", "Pr(>|z|)"])
  colnames(pend.temp) <- c("intercepto", "pendiente", "sigma", "R2 m", "R2 c", "stddev int", "stddev pend", "pvalor")
  rownames(pend.temp) <- ii
  pendientes.lineal[[as.character(ii)]] <- pend.temp
}
pendientes.lineal <- do.call(rbind, pendientes.lineal)



pendientes_df <- do.call(rbind, lapply(names(pendientes.charcos), function(año) {
  data.frame(año = as.numeric(año), pendiente = pendientes.charcos[[año]])
}))

pendientes_df$year <- as.numeric(pendientes_df$año)

setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
read.csv("lluvias.csv")->lluvias
lluvias[-c(1:3),-1]->lluvias
prof_media <- (aggregate(prof.media ~ año, data = bm.nuevo, mean))$prof.media
prof_CV <- (aggregate(CVprof ~ año, data = bm.nuevo, mean))$CVprof

cbind(nn.1, lluvias, prof_media=profundidad[,2], prof_CV=profundidad[,2],BEF=pendientes.lineal[,2])-> nn.1
as.data.frame(nn.1)-> nn.1

df_completo <- merge(pendientes_df, nn.1, by = "year")
df_completo[,-c(1)]-> df_completo
df_completo <- df_completo[, c(1, 3, 4:ncol(df_completo),2)]

par(mfrow=c(1,2))
bestglm(na.omit(df_completo[,-c(1,10)]), family = gaussian, IC = "AIC")$BestModels

df_reducido <- df_completo[, !(names(df_completo) %in% c("lluvia_muestreo", "prof_CV"))]  # Ejemplo, ajusta según tu análisis

##### ESTE ES EL QUE TIENE TODAS LAS PENDIENTES
bestglm(na.omit(df_reducido[,-c(1,7)]), family = gaussian, IC = "AIC")$BestModels


setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
read.csv("lluvias.csv")->lluvias
lluvias[-c(1:3),-1]->lluvias
prof_media <- (aggregate(prof.media ~ año, data = bm.nuevo, mean))$prof.media
prof_CV <- (aggregate(CVprof ~ año, data = bm.nuevo, mean))$CVprof

cbind(nn.1, lluvias, prof_media=profundidad[,2], prof_CV=profundidad[,2],BEF=pendientes.lineal[,2])-> nn.1
as.data.frame(nn.1)-> nn.1

m1<-(lm(BEF~R2.unexp+lluvia_anual, data=nn.1))
summary(m1)
r.squaredGLMM(m1)
par(mfrow = c(1, 2), mar = c(5, 4, 2, 1), oma = c(0, 0, 3, 0))
visreg(m1, "R2.unexp", 
       xlab = "Neutrality", 
       ylab = "BEF", 
       main = " ", 
       line.par = list(col = "navy", lwd = 2), 
       fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
visreg(m1, "lluvia_anual", 
       xlab = "Annual rain", 
       ylab = "BEF", 
       main = " ", 
       line.par = list(col = "navy", lwd = 2), 
       fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
mtext("p-valor: 0.004\nR2: 0.5", side = 3, outer = TRUE, line = 0, cex = 1)


#####################
# Representation
attach(nn.1)
summary(lm(BEF~R2.unexp + lluvia_anual))
p<-coefficients(lm(BEF~R2.unexp + lluvia_anual))
par(mfrow=c(1,2), mar=c(5,5,1,1))
BEF.neutral<-BEF-p[3]*lluvia_anual+mean(p[3]*lluvia_anual)
plot(BEF.neutral~R2.unexp, xlab="Neutrality", bty="l", pch=19, col="burlywood3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "BEF")
curve((p[1]+mean(p[3]*lluvia_anual)+p[2]*x), add=T, col="red",lwd=3)

BEF.precipitation<-BEF-p[2]*R2.unexp+mean(p[2]*R2.unexp)
plot(BEF.precipitation~lluvia_anual, xlab="Anual rain", bty="l", pch=19, col="navy",
     cex.lab = 1.8, cex.axis = 1.3, ylab = " ")
points(BEF.precipitation~lluvia_anual, cex=.6, pch=19, col="cornflowerblue")
curve((p[1]+p[3]*x+mean(p[2]*R2.unexp)), add=T, col="darkgreen",lwd=3)
text(x=55,y=1, "F(2,16): 8.4\nP: 0.003\nR2: 0.45", cex=1, pos=4)


library(tidyverse)

max_bef <- max(nn.1$BEF, na.rm = TRUE)
min_bef <- min(nn.1$BEF, na.rm = TRUE)
max_r2 <- max(nn.1$R2.unexp, na.rm = TRUE)
min_r2 <- min(nn.1$R2.unexp, na.rm = TRUE)

nn.2 <- nn.1 %>%
  mutate(BEF_scaled = min_r2 + (BEF - min_bef) * (max_r2 - min_r2) / (max_bef - min_bef))

ggplot(nn.2, aes(x = year)) +
  # Neutrality con eje izquierdo (sin escalar)
  geom_point(aes(y = R2.unexp), color = "darkgreen", alpha = 0.6, size = 5) +
  geom_line(aes(y = R2.unexp), color = "darkgreen", alpha = 0.6, size = 1) +
  # BEF con eje derecho (escalado)
  geom_point(aes(y = BEF_scaled), color = "darkblue", alpha = 0.6, size = 5) +
  geom_line(aes(y = BEF_scaled), color = "darkblue", alpha = 0.6, size = 1) +
  scale_y_continuous(
    name = "Neutrality",  # Neutrality a la izquierda
    sec.axis = sec_axis(~ min_bef + (. - min_r2) * (max_bef - min_bef) / (max_r2 - min_r2), 
                        name = "BEF")  # BEF a la derecha
  ) +
  theme_pubr() +
  labs(x = "Year") +
  theme(
    axis.title.y = element_text(size = 25, color = "darkgreen"),  # Neutrality en verde a la izquierda
    axis.text.y = element_text(size = 20, color = "darkgreen"),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y.right = element_text(size = 25, color = "navy"),  # BEF en azul a la derecha
    axis.text.y.right = element_text(size = 20, color = "navy")
  )





#############
# time lag BEF~Drift(t-1)

library(tidyverse)

# Desplazar R2.unexp un año hacia adelante
nn.3 <- nn.1 %>%
  arrange(year) %>%  # Asegurar orden temporal
  mutate(R2.unexp_lag = lag(R2.unexp)) %>%  # Desplazar un año
  filter(!is.na(R2.unexp_lag))  # Eliminar el primer año que no tiene lag

# Normalizar BEF para que comparta escala con R2.unexp_lag
max_bef <- max(nn.3$BEF)
min_bef <- min(nn.3$BEF)
max_r2 <- max(nn.3$R2.unexp_lag)
min_r2 <- min(nn.3$R2.unexp_lag)

nn.3 <- nn.3 %>%
  mutate(BEF_scaled = min_r2 + (R2.unexp_lag - min_r2) * (max_bef - min_bef) / (max_r2 - min_r2))

ggplot(nn.3, aes(x = year)) +
  # BEF con eje izquierdo (sin cambios)
  geom_point(aes(y = BEF), color = "darkblue", alpha = 0.6, size = 5) +
  geom_line(aes(y = BEF), color = "darkblue", alpha = 0.6, size = 1) +
  # R2.unexp con eje derecho (laggeado y escalado)
  geom_point(aes(y = BEF_scaled), color = "darkgreen", alpha = 0.6, size = 5) +
  geom_line(aes(y = BEF_scaled), color = "darkgreen", alpha = 0.6, size = 1) +
  scale_y_continuous(
    name = "BEF", 
    sec.axis = sec_axis(~ min_r2 + (. - min_bef) * (max_r2 - min_r2) / (max_bef - min_bef), 
                        name = "Neutrality (t-1)")
  ) +
  theme_pubr() +
  labs(x = " ") +
  theme(
    axis.title.y = element_text(size = 35, color = "darkblue"),
    axis.text.y = element_text(size = 20, color = "darkblue"),
    axis.text.x = element_text(size = 20),
    axis.title.y.right = element_text(size = 35, color = "darkgreen"),
    axis.text.y.right = element_text(size = 20, color = "darkgreen")
  )






####################################
#Pendientes evaluando linealidad
par(mfrow=c(4,5), mar=c(4,4,4,3))
pendientes=NULL
source("~/Desktop/Facultad/Maestría/BEF/Scripts/Years_evaluando.R") #1000 lineas como para pegarlo acá
setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
read.csv("lluvias.csv")->lluvias
lluvias[-c(1:3),-1]->lluvias
prof_media <- (aggregate(prof.media ~ año, data = bm.nuevo, mean))$prof.media
prof_CV <- (aggregate(CVprof ~ año, data = bm.nuevo, mean))$CVprof

cbind(nn.1, lluvias, prof_media=profundidad[,2], prof_CV=profundidad[,3],BEF=pendientes[,2])-> nn.1
as.data.frame(nn.1)-> nn.1
na.omit(nn.1)-> nn.1
par(mfrow=c(1,2))
bestglm((nn.1), family = gaussian, IC = "AIC")$BestModels
m1<-(lm(BEF~R2.unexp, data=nn.1))
summary(m1)
r.squaredGLMM(m1)
par(mfrow = c(1, 1), mar = c(5, 4, 2, 1), oma = c(0, 0, 3, 0))
visreg(m1, "R2.unexp", 
       xlab = "Neutrality", 
       ylab = "BEF", 
       main = " ", 
       line.par = list(col = "navy", lwd = 2), 
       fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
mtext("p-valor: 0.023\nR2: 0.27", side = 3, outer = TRUE, line = 0, cex = 1)


####################################
####################################
# Pendientes con covariables ambientales 
find.x <- function(m, y_en, x_en, KK = 3) {
  out=NULL
  require(bestglm)
  require(mgcv)
  X <- m[, x_en]
  y <- m[, y_en]
  Xy <- as.data.frame(cbind(X, y))
  # Modelo de regresión
  model <- bestglm(Xy, family = gaussian, IC = "AIC")
  model$BestModels-> modelos
  names(model$BestModel$coefficients)[-1]-> vars_best
  a.t <- NULL
  for (i in 1:nrow(modelos)) {
    mod.temp <- modelos[i, modelos[i, ] == TRUE]
    vars_temp <- colnames(modelos)[which(modelos[i, ] == TRUE)]
    formula <- as.formula(paste("y", "~", paste(vars_temp, collapse = "+")))
    modelo <- lm(formula, data = Xy)
    aa <- summary(modelo)$coefficients
    if ("um.rich" %in% rownames(aa)) {
      a <- aa["um.rich", 1]  # Coeficiente estimado de "um.rich"
    } else {
      a <- NA  # Si "um.rich" no está, asignar NA
    }
    a.t <- rbind(a.t, cbind(row = i, a))
  }
  print(a.t)
  cbind(modelos, a=a.t[,2])-> modelos
  modelos[which(modelos$Criterion==median(modelos$Criterion)),]-> temp
  vars_modelo <- colnames(temp)[which(temp[1, ] == TRUE)]
  out[[1]]<-(vars_modelo);
  out[[2]]<-(modelos)
  out[[3]]<- vars_best
  return(out)
}

library(visreg)
#source("~/Desktop/Facultad/Maestría/Tesis/Seleccion_mejor_modelo.R")
#source("~/Desktop/Facultad/Maestría/Tesis/Seleccion_modelo_year.R")
pendientes.covariables <- NULL
par(mfrow = c(4, 5), mar = c(4, 4, 4, 2))
for (ii in 2005:2023) {
  pendientes.covariables.temp <- as.data.frame(matrix(NA, ncol = 18, nrow = 1))
  colnames(pendientes.covariables.temp) <- c("Intercepto", "um.rich", "p.valor", "R2", "Islands", "log.Area",
                                             "Mean.Depth", "CV.Depth", "Hydroperiod", "Degree", "log.Betweenness", 
                                             "um.rich_Islands", "um.rich_log.Area", 
                                             "um.rich_Mean.Depth", "um.rich_CV.Depth", "um.rich_Hydroperiod", 
                                             "um.rich_Degree", "um.rich_log.Betweenness")
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  #df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:18, 20:25)]
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:16,29,34,30:32,23:25)]
  colnames(df.um.ch) <- c("um.año", "ch.id", "um.biom", "um.rich", "lluvia_muestreo", "lluvia_anual", "temp_muestreo", 
                          "temp_anual", "DM", "ddmm", "Shape", "Islands", "log.Area", "Mean.Depth", 
                          "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree", "log.Betweenness")
  df.um.ch <- na.omit(df.um.ch)
  df.um.ch$ch.id <- as.factor(df.um.ch$ch.id)
  df.um.ch[, c(4, 12:19)] <- scale(df.um.ch[, c(4, 12:19)])
  df.um.ch[, 20:27] <- df.um.ch[, 4] * df.um.ch[, c(12:19)]
  colnames(df.um.ch)[20:27] <- paste0("um.rich_", colnames(df.um.ch)[c(12:19)])
  
  titulo <- paste("Year:", ii)
  variables <- find.x(m = df.um.ch, y_en = 3, x_en = c(4, 12:27))
  
  if (ii == 2021) {
    variables_modelo2 <- variables[[1]]
    formula <- as.formula(paste("um.biom", "~", paste(variables_modelo2, collapse = "+")))
    modelo <- lm(formula, data = df.um.ch)
    residuos <- residuals(modelo)
    residuales.2021 <- cbind(residuos, df.um.ch)
    modelo.residual.2021 <- lm(residuos ~ um.rich, data = residuales.2021)
    visreg(modelo.residual.2021, "um.rich", xlab = "Richness", ylab = " ", main = "Year: 2021*",
           line.par = list(col = "darkgreen", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
    summary(modelo.residual.2021)
    mtext("Estimado: 0.01, p-valor: 0.93", side = 3, line = 0, cex = 0.8, col = "black")
    pendientes.covariables.temp$um.rich <- c(0.01)
    pendientes.covariables.temp$p.valor <- c(0.93)
    coeficientes <- summary(modelo)$coefficients
    
  } else {
    if (length(variables[[1]]) == 0) {
      plot.new()
      title(main = titulo)
      print("El vector está vacío")
    } else {
      variables_modelo2 <- variables[[1]]
      variables_modelo <- variables[[3]]
      formula <- as.formula(paste("um.biom", "~", paste(variables_modelo2, collapse = "+")))
      modelo <- lm(formula, data = df.um.ch)
      coeficientes <- summary(modelo)$coefficients
      
      if ("um.rich" %in% variables_modelo) {
        p_valor <- coeficientes["um.rich", 4]
        estimado <- coeficientes["um.rich", 1]
        visreg(modelo, "um.rich", xlab = "Richness", ylab = " ", main = titulo,
               line.par = list(col = "darkgreen", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
        texto <- ifelse(p_valor > 0.001, sprintf("p-valor: %.3f, Estimado: %.3f", p_valor, estimado),
                        sprintf("p-valor < 0.001, Estimado: %.3f", estimado))
        mtext(texto, side = 3, line = 0, cex = 0.8, col = "black")
        pendientes.covariables.temp$um.rich <- estimado
        pendientes.covariables.temp$p.valor <- p_valor
      } else {
        plot.new()
        title(main = titulo)
      }
    }
  }
  pendientes.covariables.temp$Intercepto <- coeficientes["(Intercept)", 1]
  covariables <- c("Islands", "log.Area", "log.Volumen", "Mean.Depth", "CV.Depth", "Hydroperiod", 
                   "Degree", "log.Betweenness", "Closenness", "um.rich_Islands", "um.rich_log.Area", 
                   "um.rich_log.Volumen", "um.rich_Mean.Depth", "um.rich_CV.Depth", "um.rich_Hydroperiod", 
                   "um.rich_Degree", "um.rich_log.Betweenness", "um.rich_Closenness")
  for (var in covariables) {
    if (var %in% variables_modelo2) pendientes.covariables.temp[[var]] <- coeficientes[var, 1]
  }
  pendientes.covariables.temp$R2 <- r.squaredGLMM(modelo)[1]
  pendientes.covariables <- rbind(pendientes.covariables, pendientes.covariables.temp)
}

############################################################
############################################################

setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
read.csv("lluvias.csv")->lluvias
lluvias[-c(1:3),-1]->lluvias
cbind(nn.1, lluvias, prof_media=profundidad[,2], prof_CV=profundidad[,3],
      BEF=pendientes.covariables[,2])-> nn.1
as.data.frame(nn.1)-> nn.1
par(mfrow=c(1,2))
bestglm(na.omit(nn.1[,-1]), family = gaussian, IC = "AIC", nvmax=2)$BestModels
m1<-(lm(BEF~R2.unexp+lluvia_anual, data=nn.1))
summary(m1)
r.squaredGLMM(m1)
par(mfrow = c(1, 2), mar = c(5, 4, 2, 1), oma = c(0, 0, 3, 0))
visreg(m1, "R2.unexp", 
       xlab = "Neutrality",  
       ylab = "BEF", 
       main = " ", 
       line.par = list(col = "navy", lwd = 2), 
       fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
visreg(m1, "lluvia_anual", 
       xlab = "Annual rain", 
       ylab = "BEF", 
       main = " ", 
       line.par = list(col = "navy", lwd = 2), 
       fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
mtext("p-valor: 0.009\nR2: 0.44", side = 3, outer = TRUE, line = 0, cex = 1)


###################
attach(nn.1)
summary(lm(BEF~R2.unexp + lluvia_anual))
p<-coefficients(lm(BEF~R2.unexp + lluvia_anual))
par(mfrow=c(1,2), mar=c(5,5,1,1))
BEF.neutral<-BEF-p[3]*lluvia_anual+mean(p[3]*lluvia_anual)
plot(BEF.neutral~R2.unexp, xlab="Neutrality", bty="l", pch=19, col="burlywood3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "BEF")
curve((p[1]+mean(p[3]*lluvia_anual)+p[2]*x), add=T, col="red",lwd=3)

BEF.precipitation<-BEF-p[2]*R2.unexp+mean(p[2]*R2.unexp)
plot(BEF.precipitation~lluvia_anual, xlab="Anual rain", bty="l", pch=19, col="navy",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "BEF")
points(BEF.precipitation~lluvia_anual, cex=.6, pch=19, col="cornflowerblue")
curve((p[1]+p[3]*x+mean(p[2]*R2.unexp)), add=T, col="darkgreen",lwd=3)
text(x=50,y=1.5, "F(2,16): 6.4\nP: 0.009\nR2: 0.38", cex=1, pos=4)





library(ggplot2)
library(ggbeeswarm)
pendientes_filtrado <- pendientes.covariables[, c(5:18)]
# Transformar la matriz en formato largo
pendientes_largo <- reshape2::melt(pendientes_filtrado, na.rm = TRUE)
# Contar las observaciones por cada variable
counts <- pendientes_largo %>%
  group_by(variable) %>%
  summarise(n = n())
# Calcular el valor máximo de 'value' para cada variable
max_values <- pendientes_largo %>%
  group_by(variable) %>%
  summarise(max_value = max(value, na.rm = TRUE))
# Crear el gráfico
ggplot(pendientes_largo, aes(x = variable, y = value, fill = variable)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.5) +
  # Gráfico de beeswarm
  geom_quasirandom(aes(color = variable), alpha = 0.7) +
  # Agregar el promedio de cada variable
  stat_summary(fun = "mean", geom = "point", shape = 18, color = "red", size = 3) +
  # Agregar el número de datos sobre cada variable
  geom_text(data = counts, aes(x = variable, y = -10, 
                               label = paste("n =", n)),color = "black", size = 4) +
  labs(title = " ", x = " ",  y = "Valores estimados") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none", color = "none")
#Alternativa en la cual el punto es pintado en función de la neutralidad del año
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(dplyr)
pendientes_filtrado <- pendientes.covariables[, c(5:18)]
pendientes_largo <- reshape2::melt(pendientes_filtrado, na.rm = TRUE)
pendientes_largo$neutralidad <- NA
pendientes_largo$neutralidad[!is.na(pendientes_largo$value)] <- rep(R2.unexp, length.out = sum(!is.na(pendientes_largo$value)))
counts <- pendientes_largo %>%
  group_by(variable) %>%
  summarise(n = n())
ggplot(pendientes_largo, aes(x = variable, y = value)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.5, aes(fill = variable)) +
  geom_quasirandom(aes(color = neutralidad), alpha = 0.7, size=3) +
  stat_summary(fun = "mean", geom = "point", shape = 18, color = "red", size = 3) +
  geom_text(data = counts, aes(x = variable, y = -10, 
                               label = paste("n =", n)), color = "black", size = 4) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.8) + 
  labs(title = " ", x = " ",  y = "Estimated values") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_color_gradient(low = "lightblue", high = "navy") +
  labs(shape = "Neutrality", color = "Neutrality")+
  theme(legend.position = "right") 


####################################
# Pendientes con covariables ambientales (quantil)
par(mfrow=c(4,5), mar=c(4,4,4,3))
pendientes.quantil=NULL
#Funcion creada para hacer algo similar al bestglm pero con regresiones quantiles
#install.packages("combinat")
library(combinat)
seleccion_mejores_modelos_cuantil <- function(data, respuesta, para_cuadratica, tau = 0.5, top_n = 5) {
  variables <- setdiff(names(data), respuesta)
  variables <- c(variables, paste0("I(", para_cuadratica, "^2)"))  # Agregar la versión cuadrática
  # Generar todas las combinaciones posibles de las variables
  combinaciones <- unlist(lapply(1:length(variables), function(m) {
    combn(variables, m, simplify = FALSE)
  }), recursive = FALSE)
  # Almacenar los modelos y sus AICs
  modelos_seleccionados <- list()
  aics <- numeric()
  for (combo in combinaciones) {
    # Crear la fórmula del modelo
    formula_text <- paste(respuesta, "~", paste(combo, collapse = "+"))
    formula <- tryCatch(as.formula(formula_text), error = function(e) NULL)
    if (is.null(formula)) next
    # Ajustar el modelo cuantílico
    modelo <- tryCatch(rq(formula, data = data, tau = tau), error = function(e) NULL)
    if (is.null(modelo)) next
    criterio <- AIC(modelo)
    modelos_seleccionados[[length(modelos_seleccionados) + 1]] <- list(
      variables = combo, aic = criterio)
    aics <- c(aics, criterio)
  }
  if (length(aics) > 0) {
    top_indices <- order(aics)[1:min(top_n, length(aics))]
    mejores_modelos <- modelos_seleccionados[top_indices]
    # Construir la matriz de salida
    salida <- data.frame(matrix(FALSE, nrow = top_n, ncol = length(variables) + 1))
    colnames(salida) <- c(variables, "AIC")
    for (i in seq_along(mejores_modelos)) {
      modelo <- mejores_modelos[[i]]
      salida[i, modelo$variables] <- TRUE
      salida[i, "AIC"] <- mejores_modelos[[i]]$aic
    }
    return(salida)
  } else { return(NULL)}
}

for(ii in 2005:2023){
  pendientes.covariables.temp<- matrix(NA, ncol=14, nrow=1)
  colnames(pendientes.covariables.temp)<- c("Intercepto", "um.rich", "p.valor", "R2","I(um.rich^2)","Islands","log.Area","log.Volumen",
                                            "Mean.Depth","CV.Depth","Hydroperiod","Degree",         
                                            "log.Betweenness",
                                            "Closenness")
  as.data.frame(pendientes.covariables.temp)-> pendientes.covariables.temp
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  #df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:18, 20:25)]
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:16,29,34,30:32,23:25)]
  colnames(df.um.ch) <- c("um.año", "ch.id", "um.biom", "um.rich", "lluvia_muestreo", "lluvia_anual", "temp_muestreo", 
                          "temp_anual", "DM", "ddmm", "Shape", "Islands", "log.Area", "Mean.Depth", 
                          "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree", "log.Betweenness")
  df.um.ch <- na.omit(df.um.ch)
  df.um.ch$ch.id <- as.factor(df.um.ch$ch.id)
  df.um.ch[, c(4, 12:19)] <- scale(df.um.ch[, c(4, 12:19)])
  df.um.ch[, 20:27] <- df.um.ch[, 4] * df.um.ch[, c(12:19)]
  colnames(df.um.ch)[20:27] <- paste0("um.rich_", colnames(df.um.ch)[c(12:19)])
  
  titulo <- paste("Year:", ii)
  salida <- seleccion_mejores_modelos_cuantil(data= df.um.ch[,-c(1,2,5:11,14,16,19,20)], respuesta = "um.biom", tau=0.8, para_cuadratica = "um.rich")
  print(salida)
  a<- which(salida$AIC==median(salida$AIC))
  colnames(salida)[which(salida[a, ] == TRUE)]-> variables
  if (length(variables) == 0) {
    plot.new()  # Inicia un nuevo gráfico vacío
    title(main = titulo)  # Agrega solo el título
    # Acción si el vector está vacío
    print("El vector está vacío")
  } else{
    formula <- as.formula(paste("um.biom", "~", paste(variables, collapse = "+")))
    model<- rq(formula, data = df.um.ch, tau = 0.8)
    resumen <- summary(model, se = "boot")
    p_valores <- as.data.frame(t(as.matrix(round(resumen$coefficients[, 4],5))))
    coeficientes <- as.data.frame(t(summary(model)$coefficients[,1]))
    if ("um.rich" %in% variables) {
      p_valor <- p_valores$um.rich
      estimado <- coeficientes$um.rich
      visreg(model, "um.rich",
             xlab = "Richness", 
             ylab = " ", 
             main = titulo,
             line.par = list(col = "darkgreen", lwd = 2),    # Personaliza línea
             fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))  # Personaliza relleno
      if (p_valor > 0.001) {
        mtext(sprintf("p-valor: %.3f, Estimado: %.3f", round(p_valor, 3), round(estimado, 3)),
              side = 3, line = 0, cex = 0.8, col = "black")
      } else {
        mtext(sprintf("p-valor < 0.001, Estimado: %.3f", round(estimado, 3)),
              side = 3, line = 0, cex = 0.8, col = "black")
      }
      pendientes.covariables.temp$um.rich<-estimado
      pendientes.covariables.temp$p.valor<- p_valor
    } else {
      plot.new()  # Inicia un nuevo gráfico vacío
      title(main = titulo)  # Agrega solo el título
    }
  }
  pendientes.covariables.temp$Intercepto<-coeficientes$`(Intercept)`
  if(("I(um.rich^2)" %in% variables)){pendientes.covariables.temp$`I(um.rich^2)`<-coeficientes$`I(um.rich^2)`} else{ }
  covariables <- c("Islands", "log.Area", "log.Volumen", "Mean.Depth", "CV.Depth", "Hydroperiod", 
                   "Degree", "log.Betweenness", "Closenness", "um.rich_Islands", "um.rich_log.Area", 
                   "um.rich_log.Volumen", "um.rich_Mean.Depth", "um.rich_CV.Depth", "um.rich_Hydroperiod", 
                   "um.rich_Degree", "um.rich_log.Betweenness", "um.rich_Closenness")
  for (var in covariables) {
    if (var %in% variables_modelo2) pendientes.covariables.temp[[var]] <- coeficientes[var, 1]
  }
  rho <- sum(model$residuals * ifelse(model$residuals >= 0, 0.5, 1 - 0.5))
  rho_nulo <- sum(rq(df.um.ch$um.rich ~ 1, tau = 0.8)$residuals * ifelse(rq(df.um.ch$um.rich ~ 1, tau = 0.8)$residuals >= 0, 0.5, 1 - 0.5))
  pseudo_r2 <- 1 - (rho / rho_nulo)
  pseudo_r2 -> pendientes.covariables.temp$R2
  pendientes.quantil<- rbind(pendientes.quantil, pendientes.covariables.temp)
}







setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
read.csv("lluvias.csv")->lluvias
lluvias[-c(1:3),-1]->lluvias
cbind(nn.1, lluvias, prof_media=profundidad[,2], prof_CV=profundidad[,2],BEF=pendientes.quantil[,2])-> nn.1
as.data.frame(nn.1)-> nn.1
na.omit(nn.1)-> nn.1
par(mfrow=c(1,2))
bestglm(nn.1, family = gaussian, IC = "AIC", nvmax=2)$BestModels
m1<-(lm(BEF~R2.unexp+lluvia_anual, data=nn.1))
summary(m1)
r.squaredGLMM(m1)
par(mfrow = c(1, 2), mar = c(5, 4, 2, 1), oma = c(0, 0, 3, 0))
visreg(m1, "R2.unexp", 
       xlab = "Neutrality", 
       ylab = "BEF", 
       main = " ", 
       line.par = list(col = "navy", lwd = 2), 
       fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
visreg(m1, "lluvia_anual", 
       xlab = "Annual rain", 
       ylab = "BEF", 
       main = " ", 
       line.par = list(col = "navy", lwd = 2), 
       fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
mtext("p-valor: 0.06\nR2: 0.33", side = 3, outer = TRUE, line = 0, cex = 1)







###########################################
###########################################
# PARTE 3: Intentando explicar Neutralidad 
###########################################
###########################################


####################################################
# EXPLICANDO NEUTRALIDAD
####################################################

R2.unexp<- nn.1$R2.unexp
pendientes.lineal[,2]-> BEF
year<- c(2005:2023)

#Auto-correlación temporal?
install.packages("lmtest")
library(lmtest)
# Ajustar un modelo de regresión (ejemplo)
modelo <- lm(R2.unexp ~ year, data = nn.1)
# Realizar la prueba de Durbin-Watson
dwtest(modelo)

library(lmtest)
# Ajustar un modelo de regresión (ejemplo)
modelo <- lm(BEF ~ year, data = nn.1)
# Realizar la prueba de Durbin-Watson
dwtest(modelo)


read.csv("lluvias.csv")-> lluvias
lluvias[-c(1:2),]-> lluvias
setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
attach(lluvias)
delta_lluvia_anual <- c(diff(lluvia_anual, lag=1))
delta_neutralidad <- c(NA, diff(nn.1$R2.unexp, lag=1))
delta_temp_anual <- c(diff(temp_anual, lag=1))
delta_lluvia_muestreo <- c(diff(lluvia_muestreo, lag=1))
delta_temp_muestreo <- c(diff(temp_muestreo, lag=1))
neutralidad_previo<- c(NA, nn.1$R2.unexp[-(nrow(nn.1))])
prof_media <- profundidad[,2]
prof_CV <- profundidad[,3]
df_model <- data.frame(
  año = 2005:2023,  # Asegúrate de ajustar los años según tus datos
  R2_unexp = nn.1$R2.unexp,
  delta_lluvia_anual = delta_lluvia_anual,
  delta_neutralidad = delta_neutralidad,
  delta_temp_anual = delta_temp_anual, 
  delta_lluvia_muestreo = delta_lluvia_muestreo,
  delta_temp_muestreo = delta_temp_muestreo,
  neutralidad_previo= neutralidad_previo,
  prof_media = prof_media,
  prof_CV = prof_CV)
library(ape)
library(cluster)
library(scales)
library(vegan)

### Gower distances
op1<-dist(df_model[,-1])
op2<-daisy(df_model[,-1], stand = TRUE) ### species x species distances (100 x 100)
#mantel(op1,op2) significativo
PCoA_neutralidad<- pcoa(op1)
neutralidad.pcoa.ejes<-PCoA_neutralidad$vectors
eigenvalues.neutralidad <- PCoA_neutralidad$values$Eigenvalues
explained_variance <- eigenvalues.neutralidad / sum(eigenvalues.neutralidad) * 100
print(explained_variance)
df_model$PCOA1<- neutralidad.pcoa.ejes[,1]
df_model$PCOA2<- neutralidad.pcoa.ejes[,2]



##########################
# This is a function for data visualization and exploration of potential relationships
# a function for data exploration
# this function generate a matrix for the bestglm packake, adding quadratic variables
# in addition a GAM model is fitted and plotted for each potential independent variable
# m i sthe matrix
# y_en identify the column with the dependent variable
# x_en identify the column(s) with independent variable(s)
# vmax is the maximum number of independent variables to consider
# KK is the parameter for GAM 
find.x<-function(m, y_en, x_en, vmax=2, KK=3){
  require(bestglm)
  X<-m[,x_en]
  y<-m[,y_en]
  for(i in x_en){
    # library(AED)
    library(mgcv)
    x<-m[,i]
    plot(y~x,main=colnames(m)[i],type="p", col="navy", pch=19)
    M0 <- gam(y~s(x, fx=FALSE, k=KK, bs="cr"))
    M0pred <- predict(M0, se = TRUE, type = "response")
    I1 <- order(x)
    lines(x[I1], M0pred$fit[I1], lty=1, col="red", lwd=2)
    lines(x[I1], M0pred$fit[I1]+2*M0pred$se[I1],lty=2, col="red", lwd=2)
    lines(x[I1], M0pred$fit[I1]-2*M0pred$se[I1],lty=2, col="red", lwd=2)
  }
  X<-cbind(X,X^2)
  colnames(X)<-c(colnames(m)[x_en], paste(colnames(m)[x_en],2,sep="_"))
  Xy<-as.data.frame(cbind(X,y))
  bestglm(Xy,family=gaussian, nvmax=vmax, IC="AIC", TopModels = 10)->bb
  return(bb)
}
par(mfrow=c(1,2));find.x(m=na.omit(df_model), y_en=4, x_en=c(3,5:12), vmax=2, KK=3)
par(mfrow=c(1,2));find.x(m=na.omit(df_model), y_en=2, x_en=c(3,5:12), vmax=2, KK=3)


m1 <- lm(R2.unexp ~ delta_lluvia_muestreo + I(delta_lluvia_muestreo^2), data = df_model)
summary(m1)
r.squaredGLMM(m1)
par(mfrow=c(1,1))
visreg(m1, "delta_lluvia_muestreo",
       xlab = "Survey rain variation", 
       ylab = "Neutrality", 
       main = "p-valor: 0.085\nR2: 0.24",
       line.par = list(col = "navy", lwd = 2), 
       fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(R2.unexp~delta_lluvia_muestreo, data= df_model,
     xlab = "∆ Sampling-period rainfall", ylab = "Neutrality",
     pch = 19, col = "blue", bty = "l", cex.lab = 2.2, cex.axis= 1.4)
curve(coef(m1)[1] + coef(m1)[2] * x + coef(m1)[3]*x^2, 
      from = min(df_model$delta_lluvia_muestreo), to = max(df_model$delta_lluvia_muestreo), 
      add = TRUE, col = "red", lwd = 2)
text(x=0-75,y=0.65, "F(2,16): 2.8\nP: 0.09\nR2: 0.17", cex=1.3, pos=4)




###########################################
###########################################
# PARTE 4: Intentando explicar mecanismos BEF detras de la Neutralidad 
###########################################
###########################################

###################################

# Librerías necesarias
library(dplyr)
library(cluster)
library(FD)

M<- as.data.frame(br.muestreos)
br.muestreos_df<- as.data.frame(br.muestreos)
M<-M[,-c(which((apply(M,2,sum))==0))]

rownames(traits)<- colnames(M)[-c(1:4)]

# Paso 1: Agrupar y sumar por Charco
M_charco <- M %>%
  group_by(año, mes, Charco) %>%
  summarise(across(starts_with("Acmella_"):last_col(), sum), .groups = "drop")
if (!all(colnames(M_charco[, -c(1:3)]) == rownames(traits))) {
  stop("Los nombres de especies entre la matriz de comunidad y los traits no coinciden.")
}
D_charco <- daisy(traits, metric = "gower", stand = TRUE)
equit_charco <- dbFD(x = D_charco, a = M_charco[, -c(1:3)], corr = "cailliez", m = 10)

# Paso 1: Agrupar y sumar por Marco y Año
M_UM <- br.muestreos_df[,-(which(apply(br.muestreos_df,2,sum)==0))]
D_UM <- daisy(traits, metric = "gower", stand = TRUE)
equit_UM <- dbFD(x = D_UM, a = M_UM[-c(which((apply(M_UM[, -c(1:4)],1,sum))==0)), -c(1:4)], corr = "cailliez", m = 10)


M_año <- M %>%
  group_by(año) %>%
  summarise(across(starts_with("Acmella_"):last_col(), sum), .groups = "drop")
# Paso 2: Comprobar la consistencia entre las matrices para Año
if (!all(colnames(M_año[, -1]) == rownames(traits))) {
  stop("Los nombres de especies entre la matriz de comunidad y los traits no coinciden.")
}
D_año <- daisy(traits, metric = "gower", stand = TRUE)
equit_año <- dbFD(x = D_año, a = M_año[, -c(1)], corr = "cailliez", m = 10)


# Paso 1: Crear un data frame con las métricas por Marco y Año
M_UM<- M_UM[-(which((apply(M_UM[, -c(1:4)],1,sum))==0)),]
resultados_UM <- data.frame(
  Marco = M_UM$Marco,
  año = M_UM$año,
  FRic = equit_UM$FRic,
  FEve = equit_UM$FEve,
  FDis = equit_UM$FDis,
  FDiv = equit_UM$FDiv,
  RaoQ = equit_UM$RaoQ)

# Paso 2: Calcular el promedio de las métricas por año para Marco
promedios_por_año_UM <- resultados_UM %>%
  group_by(año) %>%
  summarise(
    FRic_prom = mean(FRic, na.rm = TRUE),
    FEve_prom = mean(FEve, na.rm = TRUE),
    FDis_prom = mean(FDis, na.rm = TRUE),
    FDiv_prom = mean(FDiv, na.rm = TRUE),
    Rao_prom= mean(RaoQ, na.rm=TRUE)
  )


# Paso 1: Crear un data frame con las métricas por charco y año
resultados_charcos <- data.frame(
  charco = M_charco$Charco,
  año = M_charco$año,
  FRic = equit_charco$FRic,
  FEve = equit_charco$FEve,
  FDis = equit_charco$FDis,
  FDiv = equit_charco$FDiv,
  RaoQ = equit_charco$RaoQ)

# Paso 2: Calcular el promedio de las métricas por año
promedios_por_año <- resultados_charcos %>%
  group_by(año) %>%
  summarise(
    FRic_prom = mean(FRic, na.rm = TRUE),
    FEve_prom = mean(FEve, na.rm = TRUE),
    FDis_prom = mean(FDis, na.rm = TRUE),
    FDiv_prom = mean(FDiv, na.rm = TRUE),
    Rao_prom= mean(RaoQ, na.rm=TRUE)
  )

# Paso 1: Crear un data frame con las métricas por Año
resultados_año <- data.frame(
  año = M_año$año,
  FRic = equit_año$FRic,
  FEve = equit_año$FEve,
  FDis = equit_año$FDis,
  FDiv = equit_año$FDiv,
  RaoQ = equit_año$RaoQ)

# Paso 2: Calcular el promedio de las métricas por año
promedios_por_año_año <- resultados_año %>%
  group_by(año) %>%
  summarise(
    FRic_prom = mean(FRic, na.rm = TRUE),
    FEve_prom = mean(FEve, na.rm = TRUE),
    FDis_prom = mean(FDis, na.rm = TRUE),
    FDiv_prom = mean(FDiv, na.rm = TRUE),
    Rao_prom = mean(RaoQ, na.rm = TRUE)
  )

# Visualizar los resultados
promedios_por_año_UM
promedios_por_año
promedios_por_año_año


library(ggplot2)
library(tidyr)
library(dplyr)
df_time<- as.data.frame(cbind(year=c(2005:2023), neutralidad= R2.unexp,
                              fric_ch = promedios_por_año$FRic_prom, 
                              fdis_ch= promedios_por_año$FDis_prom,
                              feve_ch= promedios_por_año$FEve_prom,
                              fdiv_ch= promedios_por_año$FDiv_prom))

# Convertir el df a formato largo
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor")

# Graficar todas las variables en función de year
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor", color = "Variable")


# Convertir a formato largo y escalar por el máximo
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor") %>%
  group_by(Variable) %>%
  mutate(Valor = Valor / max(Valor, na.rm = TRUE)) # Normalización

# Graficar
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor (Escalado)", color = "Variable")


# Convertir a formato largo y escalar entre 0 y 1
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor") %>%
  group_by(Variable) %>%
  mutate(Valor = (Valor - min(Valor, na.rm = TRUE)) / 
           (max(Valor, na.rm = TRUE) - min(Valor, na.rm = TRUE)))

# Graficar
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor Normalizado (0-1)", color = "Variable")






###################################
library(ggplot2)
library(tidyr)
library(dplyr)
df_time<- as.data.frame(cbind(year=c(2005:2023), neutralidad= R2.unexp,
                              fric_um = promedios_por_año_UM$FRic_prom, 
                              fdis_um= promedios_por_año_UM$FDis_prom,
                              feve_um= promedios_por_año_UM$FEve_prom,
                              fdiv_um= promedios_por_año_UM$FDiv_prom))

# Convertir el df a formato largo
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor")

# Graficar todas las variables en función de year
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor", color = "Variable")


# Convertir a formato largo y escalar por el máximo
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor") %>%
  group_by(Variable) %>%
  mutate(Valor = Valor / max(Valor, na.rm = TRUE)) # Normalización

# Graficar
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor (Escalado)", color = "Variable")


# Convertir a formato largo y escalar entre 0 y 1
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor") %>%
  group_by(Variable) %>%
  mutate(Valor = (Valor - min(Valor, na.rm = TRUE)) / 
           (max(Valor, na.rm = TRUE) - min(Valor, na.rm = TRUE)))

# Graficar
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor Normalizado (0-1)", color = "Variable")





######################################################
df_time<- as.data.frame(cbind(year=c(2005:2023), neutralidad= R2.unexp,
                              fric_metacomm = promedios_por_año_año$FRic_prom, 
                              fdis_metacomm= promedios_por_año_año$FDis_prom,
                              feve_metacomm= promedios_por_año_año$FEve_prom,
                              fdiv_metacomm= promedios_por_año_año$FDiv_prom))

# Convertir el df a formato largo
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor")

# Graficar todas las variables en función de year
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor", color = "Variable")


# Convertir a formato largo y escalar por el máximo
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor") %>%
  group_by(Variable) %>%
  mutate(Valor = Valor / max(Valor, na.rm = TRUE)) # Normalización

# Graficar
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor (Escalado)", color = "Variable")


# Convertir a formato largo y escalar entre 0 y 1
df_long <- df_time %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Valor") %>%
  group_by(Variable) %>%
  mutate(Valor = (Valor - min(Valor, na.rm = TRUE)) / 
           (max(Valor, na.rm = TRUE) - min(Valor, na.rm = TRUE)))

# Graficar
ggplot(df_long, aes(x = year, y = Valor, color = Variable)) +
  geom_line() + 
  geom_point() + 
  theme_pubr() +
  labs(x = "Año", y = "Valor Normalizado (0-1)", color = "Variable")



print(promedios_por_año)
setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
as.data.frame(nn.1)->nn.1
cbind(promedios_por_año, nn.1[,2], BEF)-> nn.1

library(MuMIn)
library(visreg)
par(mfrow=c(1,2))
lm(BEF~ FEve_prom, data=nn.1)-> m1
summary(m1); r.squaredGLMM(m1)
visreg(m1, "FEve_prom", xlab = "FEve", ylab = "BEF", 
       main = "p-valor: 0.93\nR2: 0.0004", line.par = list(col = "navy", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
lm(R2.unexp~ FEve_prom, data=nn.1)-> m1
summary(m1); r.squaredGLMM(m1)
visreg(m1, "FEve_prom", xlab = "FEve", ylab = "Neutrality", 
       main = "p-valor: 0.38\nR2: 0.04", line.par = list(col = "navy", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))

par(mfrow=c(1,2))
lm(BEF~ FRic_prom, data=nn.1)-> m1
summary(m1); r.squaredGLMM(m1)
visreg(m1, "FRic_prom", xlab = "FRic", ylab = "BEF", 
       main = "p-valor: 0.01\nR2: 0.35", line.par = list(col = "navy", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
lm(R2.unexp~ FRic_prom, data=nn.1)-> m1
summary(m1); r.squaredGLMM(m1)
visreg(m1, "FRic_prom", xlab = "FRic", ylab = "Neutrality", 
       main = "p-valor< 0.001\nR2: 0.69", line.par = list(col = "navy", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))

par(mfrow=c(1,2))
lm(BEF~ FDis_prom, data=nn.1)-> m1
summary(m1); r.squaredGLMM(m1)
visreg(m1, "FDis_prom", xlab = "FDis", ylab = "BEF", 
       main = "p-valor: 0.9\nR2: 0.0008", line.par = list(col = "navy", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
lm(R2.unexp~ FDis_prom, data=nn.1)-> m1
summary(m1); r.squaredGLMM(m1)
visreg(m1, "FDis_prom", xlab = "FDis", ylab = "Neutrality", 
       main = "p-valor< 0.001\nR2: 0.6", line.par = list(col = "navy", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))

par(mfrow=c(1,1))
lm(BEF~ Rao_prom, data=nn.1)-> m1
summary(m1); r.squaredGLMM(m1)
visreg(m1, "Rao_prom", xlab = "Indice Rao", ylab = "BEF", 
       main = "p-valor: 0.6\nR2: 0.013", line.par = list(col = "navy", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))
lm(R2.unexp~ Rao_prom, data=nn.1)-> m1
summary(m1); r.squaredGLMM(m1)
visreg(m1, "Rao_prom", xlab = "Indice Rao", ylab = "Neutrality", 
       main = "p-valor< 0.001\nR2: 0.64", line.par = list(col = "navy", lwd = 2), fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)))


#Chequear correlacion con riqueza a ver que onda (i.e de la saldia de equit.all)
cbind(riqueza_prom_comunidades, promedios_por_año[,-1])-> prueba.cor
prueba.cor<- as.data.frame(prueba.cor)
library(corrplot)
par(mfrow=c(1,1))
my_colors <- colorRampPalette(c("navy", "white", "firebrick"))(200)
corrplot(cor(prueba.cor),
         method = "number", 
         type = "upper", 
         col = my_colors, 
         tl.col = "black",    # Color de las etiquetas
         number.cex = 0.8)  


###########################
# Neutralidad y BEF vs metricas DF
############

# A nivel UM (10% de las UM sólo fueron vistas por debajo del )
setwd("/Users/agustindeleon/Desktop/Facultad/Maestría/Tesis/Sintesis")
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
read.csv("lluvias.csv")->lluvias
lluvias[-c(1:3),-1]->lluvias
prof_media <- (aggregate(prof.media ~ año, data = bm.nuevo, mean))$prof.media
prof_CV <- (aggregate(CVprof ~ año, data = bm.nuevo, mean))$CVprof

cbind(nn.1, lluvias, prof_media=profundidad[,2], prof_CV=profundidad[,2],BEF=pendientes.lineal[,2])-> nn.1
as.data.frame(nn.1)-> nn.1

colnames(nn.1)[1:2]<- c("año", "R2.unexp")

resultados_año_completo <- merge(resultados_año, nn.1[, c("año", "R2.unexp", "BEF")], by = "año", all.x = TRUE)

resultados_UM_completo <- merge(resultados_UM, nn.1[, c("año", "R2.unexp", "BEF")], by = "año", all.x = TRUE)

resultados_charcos_completo <- merge(resultados_charcos, nn.1[, c("año", "R2.unexp", "BEF")], by = "año", all.x = TRUE)



# a nivel um
par(mfrow=c(1,5), mar=c(5,5,1,1))
mRic<-(lm(data= resultados_UM_completo, FRic~R2.unexp))
summary(mRic)
p<- coefficients(mRic)
plot(FRic~R2.unexp, data=resultados_UM_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FRic (UM)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mEve<-(lm(data= resultados_UM_completo, FEve~R2.unexp))
summary(mEve)
p<- coefficients(mEve)
plot(FEve~R2.unexp, data=resultados_UM_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FEve (UM)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mDis<-(lm(data= resultados_UM_completo, FDis~R2.unexp))
summary(mDis)
p<- coefficients(mDis)
plot(FDis~R2.unexp, data=resultados_UM_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDis (UM)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mDiv<-(lm(data= resultados_UM_completo, FDiv~R2.unexp))
summary(mDiv)
p<- coefficients(mDiv)
plot(FDiv~R2.unexp, data=resultados_UM_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDiv (UM)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mRao<-(lm(data= resultados_UM_completo, RaoQ~R2.unexp))
summary(mRao)
p<- coefficients(mRao)
plot(RaoQ~R2.unexp, data=resultados_UM_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "Rao (UM)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)


par(mfrow=c(1,5), mar=c(5,5,1,1))
mRic<-(lm(data= resultados_UM_completo, FRic~BEF))
summary(mRic)
p<- coefficients(mRic)
plot(FRic~BEF, data=resultados_UM_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FRic (UM)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="red",lwd=3)
mEve<-(lm(data= resultados_UM_completo, FEve~BEF))
summary(mEve)
p<- coefficients(mEve)
plot(FEve~BEF, data=resultados_UM_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FEve (UM)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="red",lwd=3)
mDis<-(lm(data= resultados_UM_completo, FDis~BEF))
summary(mDis)
p<- coefficients(mDis)
plot(FDis~BEF, data=resultados_UM_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDis (UM)", main= "p-value: 0.08")
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
mDiv<-(lm(data= resultados_UM_completo, FDiv~BEF))
summary(mDiv)
p<- coefficients(mDiv)
plot(FDiv~BEF, data=resultados_UM_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDiv (UM)", main= "p-value: 0.1")
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
mRao<-(lm(data= resultados_UM_completo, RaoQ~BEF))
summary(mRao)
p<- coefficients(mRao)
plot(RaoQ~BEF, data=resultados_UM_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "Rao (UM)", main= "p-value: 0.001")
curve((p[1]+p[2]*x), add=T, col="red",lwd=3)




# a nivel charco
par(mfrow=c(1,5), mar=c(5,5,1,1))
mRic<-(lm(data= resultados_charcos_completo, FRic~R2.unexp))
summary(mRic)
p<- coefficients(mRic)
plot(FRic~R2.unexp, data=resultados_charcos_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FRic (ponds)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mEve<-(lm(data= resultados_charcos_completo, FEve~R2.unexp))
summary(mEve)
p<- coefficients(mEve)
plot(FEve~R2.unexp, data=resultados_charcos_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FEve (ponds)", main= "p-value: 0.21")
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
mDis<-(lm(data= resultados_charcos_completo, FDis~R2.unexp))
summary(mDis)
p<- coefficients(mDis)
plot(FDis~R2.unexp, data=resultados_charcos_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDis (ponds)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mDiv<-(lm(data= resultados_charcos_completo, FDiv~R2.unexp))
summary(mDiv)
p<- coefficients(mDiv)
plot(FDiv~R2.unexp, data=resultados_charcos_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDiv (ponds)", main= "p-value: 0.003")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mRao<-(lm(data= resultados_charcos_completo, RaoQ~R2.unexp))
summary(mRao)
p<- coefficients(mRao)
plot(RaoQ~R2.unexp, data=resultados_charcos_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "Rao (ponds)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)


par(mfrow=c(1,5), mar=c(5,5,1,1))
mRic<-(lm(data= resultados_charcos_completo, FRic~BEF))
summary(mRic)
p<- coefficients(mRic)
plot(FRic~BEF, data=resultados_charcos_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FRic (ponds)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="red",lwd=3)
mEve<-(lm(data= resultados_charcos_completo, FEve~BEF))
summary(mEve)
p<- coefficients(mEve)
plot(FEve~BEF, data=resultados_charcos_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FEve (ponds)", main= "p-value: 0.95")
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
mDis<-(lm(data= resultados_charcos_completo, FDis~BEF))
summary(mDis)
p<- coefficients(mDis)
plot(FDis~BEF, data=resultados_charcos_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDis (ponds)", main= "p-value: 0.5")
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
mDiv<-(lm(data= resultados_charcos_completo, FDiv~BEF))
summary(mDiv)
p<- coefficients(mDiv)
plot(FDiv~BEF, data=resultados_charcos_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDiv (ponds)", main= "p-value< 0.001")
curve((p[1]+p[2]*x), add=T, col="red",lwd=3)
mRao<-(lm(data= resultados_charcos_completo, RaoQ~BEF))
summary(mRao)
p<- coefficients(mRao)
plot(RaoQ~BEF, data=resultados_charcos_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "Rao (ponds)", main= "p-value: 0.08")
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)




# a nivel metacomm
par(mfrow=c(2,3), mar=c(5,5,1,1))
#layout(matrix(c(1, 2, 3, 4, 4, 5), 2, 3, byrow = TRUE))
mRic<-(lm(data= resultados_año_completo, FRic~R2.unexp))
summary(mRic)
p<- coefficients(mRic)
plot(FRic~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FRic ")
text(x=0.6,y=1.7e-06, "F(1,17): 9.3\nP: 0.013\nR2: 0.31", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mEve<-(lm(data= resultados_año_completo, FEve~R2.unexp))
summary(mEve)
p<- coefficients(mEve)
plot(FEve~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FEve")#, main= "p-value: 0.5")
text(x=0.68,y=0.5, "F(1,17): 0.07\nP: 0.5\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
mDis<-(lm(data= resultados_año_completo, FDis~R2.unexp))
summary(mDis)
p<- coefficients(mDis)
plot(FDis~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDis")#, main= "p-value< 0.001")
text(x=0.2,y=0.225, "F(1,17): 21\nP< 0.001\nR2: 0.53", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mDiv<-(lm(data= resultados_año_completo, FDiv~R2.unexp))
summary(mDiv)
p<- coefficients(mDiv)
plot(FDiv~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDiv")
text(x = 0.2, y = 0.84, labels = "F(1,17): 9.7\nP: 0.006\nR2: 0.33", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mRao<-(lm(data= resultados_año_completo, RaoQ~R2.unexp))
summary(mRao)
p<- coefficients(mRao)
plot(RaoQ~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "Rao")#, main= "p-value< 0.001")
text(x = 0.2, y = 0.055, labels = "F(1,17): 19\nP< 0.001\nR2: 0.5", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)


par(mfrow=c(2,3), mar=c(5,5,1,1))
mRic<-(lm(data= resultados_año_completo, FRic~BEF))
summary(mRic)
p<- coefficients(mRic)
plot(FRic~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FRic ")#, main= "p-value: 0.04")
text(x = 1.1, y = 0.5e-6, labels = "F(1,17): 4.9\nP: 0.04\nR2: 0.18", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="red",lwd=3)
mEve<-(lm(data= resultados_año_completo, FEve~BEF))
summary(mEve)
p<- coefficients(mEve)
plot(FEve~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FEve ")#, main= "p-value: 0.05")
text(x = 1, y = 0.625, labels = "F(1,17): 0.2\nP: 0.63\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
mDis<-(lm(data= resultados_año_completo, FDis~BEF))
summary(mDis)
p<- coefficients(mDis)
plot(FDis~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDis")#, main= "p-value: 0.8")
text(x = 1.1, y = 0.225, labels = "F(1,17): 0.02\nP: 0.9\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
mDiv<-(lm(data= resultados_año_completo, FDiv~BEF))
summary(mDiv)
p<- coefficients(mDiv)
plot(FDiv~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "FDiv")#, main= "p-value: 0.56")
text(x = 1.1, y = 0.84, labels = "F(1,17): 0.04\nP: 0.85\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
mRao<-(lm(data= resultados_año_completo, RaoQ~BEF))
summary(mRao)
p<- coefficients(mRao)
plot(RaoQ~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=19, col="lightblue3",
     cex.lab = 1.8, cex.axis = 1.3, ylab = "Rao")#, main= "p-value: 0.65")
text(x = 1.1, y = 0.054, labels = "F(1,17): 0.003\nP: 0.97\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)



#########################
# MODELO NULO
###

#sourse(Modelo_nulo)
read.csv("Modelo_nulo_3.csv")-> Modelo_nulo

library(ggplot2)
df <- data.frame(
  neu = R2.unexp,
  z_FEve = Modelo_nulo$z_FEve)
# Graficar con ggplot2
ggplot(df, aes(x = neu, y = z_FEve)) +
  geom_point(color = "blue", size = 2, alpha = 0.7) +  # Puntos más visibles
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "black") +  # Líneas de referencia
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Línea de regresión
  theme_pubr() +  # Tema más limpio
  labs(x = "Neutralidad", y = "Z FEve",
    title = "Relación entre Neutralidad y Z FEve") +
  theme(text = element_text(size = 14),  # Tamaño del texto
    plot.title = element_text(hjust = 0.5))  # Centrar el título


cbind(Modelo_nulo, R2.unexp, BEF=nn.1$BEF)-> Modelo_nulo

# a nivel metacomm
par(mfrow=c(2,3), mar=c(5,5,1,1))
mRic<-(lm(data= Modelo_nulo, z_FRic~R2.unexp))
summary(mRic)
p<- coefficients(mRic)
plot(cex=2, z_FRic~R2.unexp, data=Modelo_nulo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value FRic")#, main= "p-value: 0.19")
text(x=0.6,y=-1.5, "R2: 0.05", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")
mEve<-(lm(data= Modelo_nulo, z_FEve~R2.unexp))
summary(mEve)
p<- coefficients(mEve)
plot(cex=2, z_FEve~R2.unexp, data=Modelo_nulo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value FEve")#, main= "p-value: 0.029")
text(x=0.6,y=-2.1, "R2: 0.21", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")
mDis<-(lm(data= Modelo_nulo, z_FDis~R2.unexp))
summary(mDis)
p<- coefficients(mDis)
plot(cex=2, z_FDis~R2.unexp, data=Modelo_nulo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value FDis")#, main= "p-value: 0.04")
text(x=0.6,y=-1.6, "R2: 0.18", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")
mDiv<-(lm(data= Modelo_nulo, z_FDiv~R2.unexp))
summary(mDiv)
p<- coefficients(mDiv)
plot(cex=2, z_FDiv~R2.unexp, data=Modelo_nulo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value FDiv")#, main= "p-value: 0.15")
text(x=0.6,y=-2.3, "R2: 0.07", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")
mRao<-(lm(data= Modelo_nulo, RaoQ~R2.unexp))
summary(mRao)
p<- coefficients(mRao)
plot(cex=2, RaoQ~R2.unexp, data=Modelo_nulo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value Rao")#, main= "p-value: 0.017")
text(x=0.6,y=0.055, "R2: 0.26", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")


par(mfrow=c(2,3), mar=c(5,5,1,1))
mRic<-(lm(data= Modelo_nulo, z_FRic~BEF))
summary(mRic)
p<- coefficients(mRic)
plot(cex=2, z_FRic~BEF, data=Modelo_nulo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value FRic")#, main= "p-value: 0.15")
text(x=1,y=-1.6, "R2: 0.07", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")
mEve<-(lm(data= Modelo_nulo, z_FEve~BEF))
summary(mEve)
p<- coefficients(mEve)
plot(cex=2, z_FEve~BEF, data=Modelo_nulo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value FEve")#, main= "p-value: 0.96")
text(x=1,y=-2.2, "R2< 0.001", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")
mDis<-(lm(data= Modelo_nulo, z_FDis~BEF))
summary(mDis)
p<- coefficients(mDis)
plot(cex=2, z_FDis~BEF, data=Modelo_nulo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value FDis")#, main= "p-value: 0.51")
text(x=1,y=-1.5, "R2< 0.01", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")
mDiv<-(lm(data= Modelo_nulo, z_FDiv~BEF))
summary(mDiv)
p<- coefficients(mDiv)
plot(cex=2, z_FDiv~BEF, data=Modelo_nulo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2.2, cex.axis = 1.3, ylab = "Z-value FDiv")#, main= "p-value: 0.55")
text(x=1,y=-2.2, "R2< 0.001", cex=1.9, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")
mRao<-(lm(data= Modelo_nulo, RaoQ~BEF))
summary(mRao)
p<- coefficients(mRao)
plot(cex=2, RaoQ~BEF, data=Modelo_nulo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2.2, cex.axis = 1.8, ylab = "Z-value Rao")#, main= "p-value: 0.45")
text(x=1,y=0.055, "R2< 0.001", cex=1.7, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)
abline(h=-1.96,lty = 2, col = "black"); abline(h=1.96,lty = 2, col = "black")







# version final de lo q está un poco más arriba


# a nivel metacomm
#dev.new()
layout(matrix(1:15, ncol = 3, byrow = TRUE), widths = c(2, 1, 2), heights = rep(2, 5))
par(mar=c(5,5,2,1))


mRic<-(lm(data= resultados_año_completo, FRic~R2.unexp))
summary(mRic)
p<- coefficients(mRic)
plot(cex=2, FRic~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2, cex.axis = 1.8, ylab = "FRic ")
text(x=0.6,y=1.7e-06, "F(1,17): 9.3\nP: 0.013\nR2: 0.31", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)

plot.new()

mRic<-(lm(data= resultados_año_completo, FRic~BEF))
summary(mRic)
p<- coefficients(mRic)
plot(cex=2, FRic~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2, cex.axis = 1.8, ylab = "FRic ")#, main= "p-value: 0.04")
text(x = 1.1, y = 0.5e-6, labels = "F(1,17): 4.9\nP: 0.04\nR2: 0.18", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="red",lwd=3)



mEve<-(lm(data= resultados_año_completo, FEve~R2.unexp))
summary(mEve)
p<- coefficients(mEve)
plot(cex=2, FEve~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2, cex.axis = 1.8, ylab = "FEve")#, main= "p-value: 0.5")
text(x=0.68,y=0.5, "F(1,17): 0.07\nP: 0.5\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)

plot.new()

mEve<-(lm(data= resultados_año_completo, FEve~BEF))
summary(mEve)
p<- coefficients(mEve)
plot(cex=2, FEve~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2, cex.axis = 1.8, ylab = "FEve ")#, main= "p-value: 0.05")
text(x = 1, y = 0.625, labels = "F(1,17): 0.2\nP: 0.63\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)



mDis<-(lm(data= resultados_año_completo, FDis~R2.unexp))
summary(mDis)
p<- coefficients(mDis)
plot(cex=2, FDis~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2, cex.axis = 1.8, ylab = "FDis")#, main= "p-value< 0.001")
text(x=0.2,y=0.225, "F(1,17): 21\nP< 0.001\nR2: 0.53", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)

plot.new()

mDis<-(lm(data= resultados_año_completo, FDis~BEF))
summary(mDis)
p<- coefficients(mDis)
plot(cex=2, FDis~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2, cex.axis = 1.8, ylab = "FDis")#, main= "p-value: 0.8")
text(x = 1.1, y = 0.225, labels = "F(1,17): 0.02\nP: 0.9\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)



mDiv<-(lm(data= resultados_año_completo, FDiv~R2.unexp))
summary(mDiv)
p<- coefficients(mDiv)
plot(cex=2, FDiv~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2, cex.axis = 1.8, ylab = "FDiv")
text(x = 0.2, y = 0.84, labels = "F(1,17): 9.7\nP: 0.006\nR2: 0.33", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)

plot.new()

mDiv<-(lm(data= resultados_año_completo, FDiv~BEF))
summary(mDiv)
p<- coefficients(mDiv)
plot(cex=2, FDiv~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2, cex.axis = 1.8, ylab = "FDiv")#, main= "p-value: 0.56")
text(x = 1.1, y = 0.84, labels = "F(1,17): 0.04\nP: 0.85\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)



mRao<-(lm(data= resultados_año_completo, RaoQ~R2.unexp))
summary(mRao)
p<- coefficients(mRao)
plot(cex=2, RaoQ~R2.unexp, data=resultados_año_completo, xlab="Neutrality", bty="l", pch=16, col="red",
     cex.lab = 2, cex.axis = 1.8, ylab = "Rao")#, main= "p-value< 0.001")
text(x = 0.2, y = 0.055, labels = "F(1,17): 19\nP< 0.001\nR2: 0.5", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="darkgreen",lwd=3)

plot.new()

mRao<-(lm(data= resultados_año_completo, RaoQ~BEF))
summary(mRao)
p<- coefficients(mRao)
plot(cex=2, RaoQ~BEF, data=resultados_año_completo, xlab="BEF", bty="l", pch=18, col="navy",
     cex.lab = 2, cex.axis = 1.8, ylab = "Rao")#, main= "p-value: 0.65")
text(x = 1.1, y = 0.054, labels = "F(1,17): 0.003\nP: 0.97\nR2< 0.01", cex=1.4, pos = 1, offset = 0, xpd = TRUE)
curve((p[1]+p[2]*x), add=T, col="black",lwd=3)


