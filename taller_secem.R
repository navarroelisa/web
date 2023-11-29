
library(terra)
library(unmarked)
library(AICcmodavg)

################################################################################

rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))   stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

################################################################################
set.seed(1)
simgrid <- expand.grid(1:50, 1:50)
n <- nrow(simgrid)

distance <- as.matrix(dist(simgrid))
elev <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
elev[cellFromXY(elev,simgrid-0.5)] <- round(as.vector(rmvn(1, rep(1000, n), 50000 * exp(-(0.05) * distance))),0)
arbu <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
arbu[cellFromXY(arbu,simgrid-0.5)] <- round(as.vector(rmvn(1, rep(50, n), 50 * exp(-(0.5) * distance))),1)

psi <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
psi[] <- plogis( (-0.2 - 0.8*scale(elev) + 0.4*scale(arbu))[])
pres <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
pres[] <- rbinom(n, 1, psi[])
plot(pres)



nsitios <- 100
nnoches <- 7

df <- data.frame(id = sample(1:n, nsitios),
                 marca = sample(0:1, nsitios, replace = TRUE),
                 o1 = NA, o2 = NA, o3 = NA, o4 = NA, o5 = NA, o6 = NA, o7 = NA)
df <- cbind(pres[df$id],elev[df$id],arbu[df$id],df)
colnames(df)[1:3] <- c("pres","elev", "arbu")
plot(pres)
points(xyFromCell(pres, df$id))

for(i in 1:nsitios){
  cell <- df$id[i]
  for(j in 1:nnoches){
    marca <- df$marca[i] 
    p <- plogis(0.2+ 0.5*marca - 0.5*j - 0.6*as.numeric(scale(arbu)[cell]))
    df[i,5+j] <- rbinom(1, as.numeric(pres[cell]), p)
    #print(p)
  }
}

df$marca <- as.factor(df$marca)
levels(df$marca) <- c("A", "B")

#datos <- cbind(id=df$id, marca=df$marca, df[,6:12],
#               matrix(rep(1:7, 1, each =nsitios), nsitios,7))
#names(datos)[10:16] <- c("dia1","dia2","dia3","dia4","dia5","dia6","dia7")

datos <- cbind(id=df$id, marca=df$marca, elev=scale(df$elev), arbu=scale(df$arbu), df[,6:12],
               matrix(rep(1:7, 1, each =nsitios), nsitios,7))
names(datos)[12:18] <- c("dia1","dia2","dia3","dia4","dia5","dia6","dia7")

obserCov <- list(dia = as.matrix(datos[,12:18]))
datosUM <- unmarkedFrameOccu(as.matrix(datos[,5:11]), siteCovs=datos[,2:4], obsCovs=obserCov)
m1 <- occu(~ marca  + dia ~ elev + arbu, data = datosUM)
m2 <- occu(~ marca  + dia + arbu ~ elev + arbu, data = datosUM)
m3 <- occu(~ marca  + dia + arbu + elev ~ elev + arbu, data = datosUM)


plotEffects(m2, type = "det", covariate="arbu")

# scaled = (x-mean(x))/sd
elevScaled <- (elev - 958.29) / 187.6693
arbuScaled <- (arbu - 49.78) / 7.041493

dataPred <- c(elevScaled, arbuScaled)
names(dataPred) <- c("elev", "arbu")

pred <- predict(m2, dataPred, type = "state")
plot(pred)

mb.gof.test(m2, nsim = 50)

################################################################################

set.seed(2)
simgrid <- expand.grid(1:50, 1:50)
n <- nrow(simgrid)

distance <- as.matrix(dist(simgrid))
tempe <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
tempe[cellFromXY(tempe,simgrid-0.5)] <- round(as.vector(rmvn(1, rep(22, n), 22 * exp(-(0.08) * distance))),0)
cober <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
cober[cellFromXY(cober,simgrid-0.5)] <- round(as.vector(rmvn(1, rep(2, n), 0.3 * exp(-(0.05) * distance))),1)
cober[] <- round(cober[],0)
plot(cober)
cob1 <- cober == 1
cob2 <- cober == 2
cob3 <- cober == 3
cob4 <- cober == 4

lambda <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
lambda[] <- exp( ( -3*cob1 -2.2*cob2 + 0.4*cob3 + 2.8*cob4 + 0.8*scale(tempe))[])
abun <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
abun[] <- rpois(n, lambda[])
plot(abun)

nsitios <- 100
nnoches <- 3

df <- data.frame(id = sample(1:n, nsitios),
                 o1 = NA, o2 = NA, o3 = NA)
tiemp <- matrix(sample(1:3, nsitios*nnoches, replace = TRUE), 100, 3)
df <- cbind(abun[df$id],tempe[df$id],cober[df$id],df, tiemp)
colnames(df)[c(1:3,8:10)] <- c("abun","tempe","cober","t1","t2","t3")
points(xyFromCell(abun, df$id))

tt <- list(t1 = data.frame(t1 = df$t1==1,t2 = df$t1==2,t3 = df$t1==3),
           t2 = data.frame(t1 = df$t2==1,t2 = df$t2==2,t3 = df$t2==3),
           t3 = data.frame(t1 = df$t3==1,t2 = df$t3==2,t3 = df$t3==3))

for(i in 1:nsitios){
  cell <- df$id[i]
  for(j in 1:nnoches){
    p <- plogis(-1.7*tt[[j]][i,1] + 0.7*tt[[j]][i,2] + 1.2*tt[[j]][i,3]) 
    df[i,4+j] <- rbinom(1, as.numeric(abun[cell]), p)
    print(p)
  }
}


df$cober <- as.factor(df$cober)
levels(df$cober) <- c("urb", "agro", "tran", "bosq")
tOc <- df[8:10]
tOc[tOc[] == 1] <- "baja"
tOc[tOc[] == 2] <- "media"
tOc[tOc[] == 3] <- "alta"

datos <- cbind(id=df$id, cober=df$cober, tempe=scale(df$tempe), df[,5:7], tOc)

obserCov <- list(tOc = datos[,7:9])
datosUM <- unmarkedFramePCount(y = as.matrix(datos[,4:6]), siteCovs=datos[,2:3], obsCovs=obserCov)
m1 <- pcount(~ tOc ~ tempe + cober, data = datosUM)
m2 <- pcount(~ tempe ~ tempe + cober, data = datosUM)
m3 <- pcount(~ 1 ~ cober, data = datosUM)

plotEffects(m1, type = "det", covariate="tOc")
plotEffects(m1, type = "state", covariate="tempe")
plotEffects(m1, type = "state", covariate="cober")

# scaled = (x-mean(x))/sd
tempeScaled <- (tempe[] - 26.22) / 5.02

dataPred <- data.frame(cbind(tempeScaled, cober[]))
names(dataPred) <- c("tempe", "cober")
dataPred$cober[dataPred$cober == 1] <- "urb"
dataPred$cober[dataPred$cober == 2] <- "agro"
dataPred$cober[dataPred$cober == 3] <- "tran"
dataPred$cober[dataPred$cober == 4] <- "bosq"
dataPred$cober <- as.factor(dataPred$cober)


pred <- predict(m1, dataPred, type = "state")

predRas <- rast(nrows=50, ncols=50, xmin=0, xmax=50, ymin = 0, ymax = 50)
predRas[] <- pred$Predicted
plot(predRas)

Nmix.gof.test(m1, nsim = 50)






