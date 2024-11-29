library(ecospat)
library(magrittr)
library(raster)
library(ade4)
library(ape)
library(adehabitatHR)

#load climate data
chelsa_cur <- list.files("~/path-to-environmental-variables",
                         pattern = ".grd$", 
                         full.names = T)%>%stack

#read native and introduced populations' occurence data
ext_nat <- fread("native.csv")
ext_inv <- fread("introduced.csv")

#set study extent by using Minimum Convex Polygon
coordinates(ext_nat) <- c("x", "y")
proj4string(ext_nat) <- CRS( "+proj=longlat +datum=WGS84 +no_defs")
ext_nat <- mcp(ext_nat, percent = 100)
ext_nat_buf <- buffer(ext_nat,25000)

ext_inv <- ts_na
coordinates(ext_inv) <- c("x", "y")
proj4string(ext_inv) <- CRS( "+proj=longlat +datum=WGS84 +no_defs" )
ext_inv <- mcp(ext_inv, percent = 100)
ext_inv_buf <- raster::buffer(ext_inv,25000)

#crop the environmental data to the native and invasive geographical ranges
natEnvR <- crop(chelsa_cur, ext_nat_buf)%>%mask(ext_nat_buf)
invEnvR <- crop(chelsa_cur, ext_inv_buf)%>%mask(ext_inv_buf)

natEnvM <- getValues(natEnvR)
invEnvM <- getValues(invEnvR)

#remove missing values
natEnvM <- natEnvM [complete.cases(natEnvM ), ]
invEnvM <- invEnvM [complete.cases(invEnvM ), ]

#produce global environmental background data
globalEnvM <- rbind(natEnvM , invEnvM )

natEnv <- cbind(ext_nat,extract(natEnvR, ext_nat))
invEnv <- cbind(ext_inv,extract(invEnvR, ext_inv))

#Niche quantification
pca.clim <- dudi.pca(globalEnvM, center =T, scale =T, scannf = FALSE, nf = 2)
global.scores <- pca.clim$li

#map occur to 2D space with selected variables
nativeLS.scores <-
  suprow(pca.clim,
         data.frame(natEnv)[, colnames(globalEnvM)])$li   
nativeLS.scores <- na.omit(nativeLS.scores)%>%as.data.frame
invasiveLS.scores <-
  suprow(pca.clim,
         data.frame(invEnv)[, colnames(globalEnvM)])$li
invasiveLS.scores <- na.omit(invasiveLS.scores)%>%as.data.frame

nativeEnv.scores <- suprow(pca.clim, natEnvM)$li
invasiveEnv.scores <- suprow(pca.clim, invEnvM)$li

#calculate the occurrence density grid for both native and introduced population
nativeGrid <- ecospat.grid.clim.dyn(global.scores,
                                    nativeEnv.scores,
                                    nativeLS.scores)

invasiveGrid <- ecospat.grid.clim.dyn(global.scores,
                                      invasiveEnv.scores, 
                                      invasiveLS.scores)

#calculate niche overlap
ecospat.niche.overlap(nativeGrid, invasiveGrid, cor =T)

#visualizing niche categories, niche dynamics and climate analogy between ranges
ecospat.plot.niche.dyn(nativeGrid, invasiveGrid, quant = 0.05, interest = 2)
