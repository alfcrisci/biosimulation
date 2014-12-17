
########################################################################################################
# Install/ load packages

library(rgdal)
library(zoo)
library(xts)
library(lubridate)
library(deSolve)
library(plyr)
library(verification)

###########################################################################################################
# Setting della directory di lavoro


setwd("D:/lav_selmi_manfredi")

########################################################################################################
# load model's code

source("Rbiosim.006.r")


########################################################################################################
# Marco Selmi 's selection od sites
#########################################################################################################

sel_row=c(1,5,8,12,16,31,33,34,42,50)


#########################################################################################################

metadati_sim=read.csv("metadati_sim.csv")
matrix_2011=read.csv("selected_2011.csv",header=T)
matrix_2012=read.csv("selected_2012.csv",header=T)


##################################################################
# Instanziamento gli oggetti biologici in versione standard
#################################################################

alpha_l_seq=seq(0.5,1.8,0.1)
alpha_u_seq=c(0,seq(0.000001,0.000005,0.000001))
alpha_a_seq=c(0,seq(0.00001,0.00005,0.00001))


f_meteo_2011=list.files(path="dati_meteo_2011_location",full.names=T)
f_meteo_2012=list.files(path="dati_meteo_2012_location",full.names=T)

#############################################################################################################à
# Create results object
res_sim=list()


			
#############################################################################################################à
# Create results object

sim_index=1	
local_index=0		
for (loc_index in sel_row) {
local_index=local_index+1
########################################################################################################

metadati_meteo<-metameteo(nome=as.character(metadati_sim$Area[loc_index]),
			     rete="Modellistica",
			     tipo_dati="Simulazioni WRF",
			     tipo_stazione="Estrazione puntuale WRF",
			     fonte_dati="IBIMET",			     
			     lat=metadati_sim$Lat[loc_index],	
                 lon=metadati_sim$Lon[loc_index],
			     quota="local");


biocont<-biocontainer(nrecipients=50) #BS

biomet<-biometeo(f_meteo_2012[local_index],biocont,metadati_meteo);

biopop<-biopopulation(uova=100,
                      larve=0,
                      pupe=0,
                      adulti=0,
					  uova_diap=0)
					  
########################################################################################################
			  

for (l_index in alpha_l_seq) {

for (u_index in alpha_u_seq) {

for (a_index in alpha_a_seq) {

########################################################################################################
simul=list()					  
biopar<-bioparameters(alfa_l=l_index,
                      alfa_u=u_index,
					  alfa_a=a_index)
				  
########################################################################################################
simul$metadati=metadati_meteo
simul$biocont=biocont
#simul$biomet=biomet
simul$biopop=biopop		
simul$biopar=biopar	
simul$ts_sim=NA
simul$ts_obs=NA
simul$res=NA
simul$verify=NA
simul$success=0	

	
options(show.error.messages = FALSE)
simulazione<- try(rbiosim_run(biomet,biopop,biopar,biocont,volfix=1,soglia_pioggia=4,delta_ts=24,resamples=5))
options(show.error.messages = TRUE)
if (class(simulazione) == "try-error") {break}

						       ########################################################################################################
                               
                               zoo_popfin_xts=as.xts(zoo(simulazione$pop_fin_day$N_u/(biocont$nrecipients),as.Date(simulazione$biometeo$dates)))
                               zoo_eggs_xts=as.xts(zoo(matrix_2012[,local_index+1],as.Date(matrix_2012$Data)))
                               merged_ts=merge.xts( zoo_eggs_xts,zoo_popfin_xts,join = "inner")

                               ########################################################################################################
                               simul$ts_sim=zoo_popfin_xts
                               simul$ts_obs=zoo_eggs_xts
                               simul$res=merged_ts
                               simul$verify=verify(merged_ts$zoo_eggs_xts, merged_ts$zoo_popfin_xts, frcst.type = "cont", obs.type = "cont")
                               simul$success=1       
	

		rm(zoo_popfin_xts)
		rm(zoo_eggs_xts)
		rm(merged_ts)
		
        ###################################################################################################################
        res_sim=append(res_sim,simul)
		cat(paste(loc_index,sim_index,".....fatta!\n"))
        sim_index=sim_index+1
        ###################################################################################################################

}
}
}
}

saveRDS(res_sim,"simulation_selected.rds")
 
print("Fatto!")

 
########################################################################################################
# Reference
#  http://stackoverflow.com/questions/16442396/convert-daily-to-weekly-monthly-data-with-r
#  http://www.quantmod.com/examples/
#  http://www.r-bloggers.com/melbournes-weather-and-cross-correlations/


########################################################################################################
# Supplementary code


# A<- verify(obs, pred, frcst.type = "cont", obs.type = "cont")


# library(zoo)
# tt <- seq(Sys.Date(), by='day', length=365)
# vals <- data.frame(A=runif(365), B=rnorm(365), C=1:365)
#library(zoo)
# z <- read.zoo(df)

# classic graphics in separate and single plots
# plot(z)
# plot(z, screen = 1)

# lattice graphics in separate and single plots
# library(lattice)
# xyplot(z)
# xyplot(z, screen = 1)

# ggplot2 graphics in separate and single plots
# library(ggplot2)
# autoplot(z) + facet_free()
# autoplot(z, facet = NULL) z <- zoo(vals, tt)


# Now I define a function which extracts the year 
# and the number of the week (drop %Y if you don't need to distinguish between years):

# week <- function(x)format(x, '%Y.%W')

# You can use this function to aggregate the zoo object with mean (for example):

# aggregate(z, by=week, FUN=mean)

#ex1=biomodelday(biopop,biopar,biocont,tmed=24,twmed=21,deltatime=24)
#model_out=elab_day_sim(ex1,biopop,biocont, p100_su=1,n_resampling=5,times=12)
		
#merged_ts=merge.xts( zoo_eggs_xts,zoo_popfin_xts,join = "inner")						                      
#merged_ts=na.approx(merged_ts)
