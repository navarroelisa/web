# Taller para el XVI Congreso de la SECEM (Granollers 6/12/2023)
# Introducción a la modelización jerárquica: detectabilidad imperfecta en modelos de ocupación y N-mixture
# Javier Fernández-López

########################
# MODELOS DE OCUPACIÓN #
########################

# Instalamos y cargamos las librerías que vamos a utilizar.
#install.packages("terra")
#install.packages("unmarked")
#install.packages("AICcmodavg")
library(terra)
library(unmarked)
library(AICcmodavg)

# Leemos nuestros datos en formato CVS
datos <- read.csv("https://raw.githubusercontent.com/jabiologo/web/master/tutorials/tallerSECEM_files/occuDatos.csv")
# Si hemos descargado previamente los datos, podemos indicar la ruta (la carpeta
# en nuestro ordenador donde se localizan los archivos) y cargarlos desde ahí.
# datos <- read.csv("mi/ruta/occuDatos.csv")

# Le pedimos a R que nos muestre las primeras filas de nuestros datos.
head(datos)

# Cargamos las capas raster correspondientes a nuestras covariables predictoras,
# en este caso elevación y porcentaje de vegetación arbustiva.
elev <- rast("https://github.com/jabiologo/web/raw/master/tutorials/tallerSECEM_files/elev.tif")
arbu <- rast("https://github.com/jabiologo/web/raw/master/tutorials/tallerSECEM_files/arbu.tif")
# Si hemos descargado previamente las capas, podemos indicar la ruta (la carpeta
# en nuestro ordenador donde se localizan los archivos) y cargarlos desde ahí.
# elev <- rast("mi/ruta/elev.tif")
# arbu <- rast("mi/ruta/arbu.tif")

names(elev) <- "elev"
names(arbu) <- "arbu"

# Vamos a graficar estas capas para hacernos una idea de cómo son. Además,
# colocaremos encima las localizaciones de nuestros dos modelos de cámaras.
par(mfrow = c(1,2))
plot(elev, main = "Elevación")
points(xyFromCell(elev, datos$id[datos$marca == "A"]), col = "darkred", pch = 19)
points(xyFromCell(elev, datos$id[datos$marca == "B"]), col = "darkblue", pch = 19)
plot(arbu, main = "Porcentaje arbusto")

# En este paso estandarizaremos nuestras variables para que su media sea 0 y su
# desviación estándar 1. Esta es una práctica habitual en modelización que ayuda
# a ajustar los modelos de forma más efectiva. La fórmula es:
# valor escalado = (valor de la cov - valor medio de la cov)/ SD de la cov
datos$elev <- scale(elev[datos$id])[1:100]
datos$arbu <- scale(arbu[datos$id])[1:100]

# Volvemos a inspeccionar las primeras filas de nuestro juego de datos.
head(datos)

# A continuación prepararemos los datos para que sean entendidos por el paquete
# unmarked. Para ello debemos colocar de una forma determinada nuestras
# observaciones (detecciones/no detecciones en este caso), nuestras covariables
# predictoras. Recordad que las covariables predictoras que pueden estar
# relacionadas con los sitios de muestreos (siteCovs) o con las ocasiones de
# muestreo (obsCovs).
obserCov <- list(dia = as.matrix(datos[,10:16]))
datosUM <- unmarkedFrameOccu(as.matrix(datos[,3:9]), siteCovs=datos[,c(2,17,18)], obsCovs=obserCov)

# A continuación ajustaremos el modelo de ocupación utilizando la marca de la 
# cámara trampa y los días desde que se colocó la cámara y el cebo como 
# covariables que afectan a la detectabilidad (proceso observacional), y 
# elevación y porcentaje de cobertura arbustiva como covariables que afectan a 
# la ocupación o probabilidad de presencia (proceso ecológico).
m1 <- occu(~ marca  + dia ~ elev + arbu, data = datosUM)

#Inspeccionamos los resultados del modelo.
summary(m1)

# Inspeccionamos los efectos de la covariable de la cámara  y del procentaje 
# de arbusto.
plotEffects(m1, type = "det", covariate="marca")
plotEffects(m1, type = "state", covariate="arbu")

# También podemos predecir de forma espacialmente explícita nuestro modelo 
# de ocupación para inspeccionar la probabilidad de presencia de nuestra especie
# sobre un mapa. Para ello primero estandarizaremos los rásters de la misma
# forma que lo hicimos con nuestras variables.
elevScaled <- (elev - 958.29) / 187.6693
arbuScaled <- (arbu - 49.78) / 7.041493
dataPred <- c(elevScaled, arbuScaled)
names(dataPred) <- c("elev", "arbu")

# Después, podemos utilizar la función predict con nuestro modelo y las nuevas
# capas raster estandarizadas y graficar los resultados.
pred <- predict(m1, dataPred, type = "state")
plot(pred)

# Por último podemos explorar la bondad de ajuste de nuestro modelo, esto es,
# cuanto de bien predice nuestras observaciones el modelo que hemos construido.
mb.gof.test(m1, nsim = 100)




#####################
# MODELOS N-MIXTURE #
#####################

# Si hemos cargado los paquetes necesarios anteriormente, no es necesario volver
# a cargarlos. Los paquetes que utilizaremos terra, unmarked y AICcmodavg
# Leemos nuestros datos en formato CVS.
datos <- read.csv("https://github.com/jabiologo/web/raw/master/tutorials/tallerSECEM_files/nmixDatos.csv")
# Si hemos descargado previamente los datos, podemos indicar la ruta (la carpeta
# en nuestro ordenador donde se localizan los archivos) y cargarlos desde ahí
# datos <- read.csv("mi/ruta/nmixDatos.csv").

# Visualizamos las primeras filas de nuestros datos.
head(datos)

# Cargamos las capas raster correspondientes a nuestras covariables predictoras,
# en este caso tipo de cobertura de vegetación y temperatura media de cada sitio.
cober <- rast("https://github.com/jabiologo/web/raw/master/tutorials/tallerSECEM_files/cober.tif")
tempe <- rast("https://github.com/jabiologo/web/raw/master/tutorials/tallerSECEM_files/tempe.tif")
# Si hemos descargado previamente las capas, podemos indicar la ruta (la carpeta
# en nuestro ordenador donde se localizan los archivos) y cargarlos desde ahí.
# cober <- rast("mi/ruta/cober.tif")
# tempe <- rast("mi/ruta/tempe.tif")
names(cober) <- "cober"
names(tempe) <- "tempe"

# Vamos a graficar estas capas para hacernos una idea de cómo son. Colocaremos 
# encima las localizaciones de nuestros conteos repetidos. Nótese que la 
# variable cobertura del suelo es categórica y cada númer corresponde a un tipo
# de cobertura:
# 1 = urbano
# 2 = agrícola
# 3 = transición
# 4 = bosque
par(mfrow = c(1,2))
plot(tempe, main = "Temperatura")
points(xyFromCell(tempe, datos$id), pch = 19, cex = 0.5)
plot(cober, main = "Cobertura del suelo")

# Estandarizamos la variable temperatura y convertimos a factor la variable de
# cobertura del suelo, ya que se trata de una variable categórica.
datos$tempe <- scale(tempe[datos$id])[1:100]
datos$cober <- as.factor(cober[datos$id][,1])

# Volvemos a visualizar las primeras filas de nuestro juego de datos.
head(datos)

# Como hicimos anteriormente, vamos a organizar nuestros datos de una manera
# que entienda el paquete unmarked.
obserCov <- list(tOc = datos[,5:7])
datosUM <- unmarkedFramePCount(y = as.matrix(datos[,2:4]), siteCovs=datos[,8:9], obsCovs=obserCov)

# A continuación ajustaremos un modelo N-mixture utilizando la temperatura de la
# ocasión de muestreo y la temperatura media del sitio como covariables que
# afectan a la detectabilidad y la temperatura media del sitio y el tipo de
# cobertura del suelo como variables que afectan a la abundancia de la especie.
m1 <- pcount(~ tOc + tempe ~ tempe + cober, data = datosUM, K = 200)

# Inspeccionamos los resultados del modelo
summary(m1)

# Inspeccionamos algunos efectos de algunas covariables
plotEffects(m1, type = "det", covariate="tempe")
plotEffects(m1, type = "state", covariate="tempe")

# Preparamos las variables para poder realizar las predicciones
tempeScaled <- (tempe[] - 26.22) / 5.02
dataPred <- data.frame(cbind(tempeScaled, cober[]))
names(dataPred) <- c("tempe", "cober")
dataPred$cober <- as.factor(dataPred$cober)

# Realizamos las predicciones del modelo y las graficamos
pred <- predict(m1, dataPred, type = "state")
predRas <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
predRas[] <- pred$Predicted
plot(predRas)

# Por último, podemos analizar la bondad de ajuste de nuestro modelo.
Nmix.gof.test(m1, nsim = 100)
