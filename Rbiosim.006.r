#####################################################################################################################################################
# Costruttori di oggetto di una singola simulazione sito località e popolazione con dati relativi a stazione						   

#####################################################################################################################################################
# Caricamento funzioni

source("rbio_functions.r")						   
						   
############################################################################################
# Funzione metadati meteo default Pisa 
						   
metameteo<-function(nome="Pisa San Giusto",
			   rete="Aeronautica Militare",
			   tipo_dati="Osservazioni",
			   tipo_stazione="SYNOP",
			   fonte_dati="IBIMET CNR",			     
			   lat=43.0,	
               lon=11.0,
			   quota=40,
			   timeformat="giornaliero"
			    ){
  object <- list(nome=nome,
			     rete=rete,
			     tipo_dati=tipo_dati,
			     tipo_stazione=tipo_stazione,
				 fonte_dati=fonte_dati,
			     lat=lat,	
                 lon=lon,
			     quota=quota,
			     timeformat=timeformat
				 );
  class(object) <- "metameteo"
  attr(object,"nome") <- "Denominazione Stazione"
  attr(object,"rete") <- "Appartenenza di network osservativo"
  attr(object,"tipo_dati")<-"Tipologia di stazione"
  attr(object,"tipo_stazione")<-"Tipologia di stazione"
  attr(object,"fonte_dati")<-"Origine dati"
  attr(object,"lat")<-"Latitudine stazione Sessadecimali WGS 84"
  attr(object,"lon")<-"Longitudine stazione Sessadecimali WGS 84"	
  attr(object,"quota")<-"Quota stazione metri"	
  attr(object,"timeformat")<-"Nome del sito"	
  class(object) <- "metameteo"
  return(object)
}					   
##############################################################################################################################################
# Funzione oggetto  meteo

biometeo <- function(filemeteo,biocontainer,metameteo, w_model="watermodel.rda",errmodeltmax=0,errmodeltmin=0,errmodeltmed=0,S_dry=1,weigth_k=5)
                           {
						   
						    require(zoo);
						    
							if (class(biocontainer)!="biocontainer") { stop("Object metameteo argument must be of class biocontainer" ) } 
							if (class(metameteo)!="metameteo") { stop("Object metameteo argument must be of class metameteo" ) } 
							
							wmodel_tmed=readRDS("wmodel_tmed.rds")
							
							#################################################################################################
                            # Lettura ascii Text tab-formatted: anno	mese	giorno	tmed	tmax	tmin	urel	prec	
						    #################################################################################################
							
							if (is.data.frame(filemeteo)){ filedata=filemeteo}
							
                           	else {
                         	      filedata=read.table(filemeteo, header=TRUE, sep="\t",na.strings="NA", dec=".", strip.white=TRUE)
							     }
							
							filedata$tmax=filedata$tmax+errmodeltmax
							filedata$tmin=filedata$tmin+errmodeltmin
							filedata$tmed=filedata$tmed
							
							datadate=as.Date(paste(filedata$mese,filedata$giorno,filedata$anno,sep='/'), "%m/%d/%Y")
							
							jdays=jdday(datadate)
							
							#######################################################################################################################

							
							suncalc=suncalc(jdday(datadate),Lat=metameteo$lat,Lon=metameteo$lon)
							
							lm.filedata<- lm(tmed ~ tmax + tmin+as.factor(mese), data=filedata)
							
							tmed_estimate=lm.filedata$fitted.values+errmodeltmed
							
							#######################################################################################################################
                            # Water temperature estimation
							
							id <- which(!(filedata$mese %in% levels(wmodel_tmed$model[,5])))
                            filedata$mese[id]<-NA
							new_data=data.frame(mese=as.factor(filedata$mese),tmed=tmed_estimate,tmax=filedata$tmax,tmin=filedata$tmin)
							
							w_tmed=predict(wmodel_tmed,new_data)
							
							id.na<-which(is.na(w_tmed)) # check missing values
							w_tmed[id.na]<-tmed_estimate[id.na] # replace missing values with mean temperatures 
							
							#######################################################################################################################
                            # Water evaporation temperature estimation
							
							evap_obj<-evap_number(tmed_estimate,w_tmed,filedata$urel,L=biocontainer$D/100)
							
							tot_evap=evap_obj$mevday*biocontainer$nrecipients
							
							#######################################################################################################################
							# Pourcentage (%) of eggs diapausing
							
							p_diap=ore_schiusura(suncalc$daylength)
							
							#######################################################################################################################
							# Pourcentage (%) of eggs diapausing
							
							peso_pioggia=weigth_k
							
							#######################################################################################################################
							# Pourcentage (%) of eggs diapausing
							
							object <- list( tmed_est=tmed_estimate,
							                tmed=filedata$tmed, 
							                tmax=filedata$tmax,
											tmin=filedata$tmin,
											urel=filedata$urel,
											prec=filedata$prec,
											w_tmed=w_tmed,
											jdays=jdays,
											daylengths=suncalc$daylength,
											p_diap=p_diap,
											hsunrise=suncalc$sunrise,
											hsunset=suncalc$sunset,
											tot_evap=tot_evap,
											S_dry=S_dry,
											prevdrydays=drydaycons(filedata$prec,S=S_dry),						   
                                            weight_alpha=weigthrain(drydaycons(filedata$prec),S=peso_pioggia),
											weight_K=peso_pioggia,
											vol_iniziale=biocontainer$vol_iniziale,
											nrecipients=biocontainer$nrecipients,
											biocontainer=biocontainer,
										    stazione=metameteo$nome,
											lat=metameteo$lat,	
                                            lon=metameteo$lon,
			                                quota=metameteo$quota,
											dates=datadate,
											container_class=biocontainer$container_class,
											timeseries=as.zoo(cbind(w_tmed,new_data),datadate)
			                                )
                           class(object) <- "biometeo"
						   attr(object,"tmed_est") <- "Temperatura media aria stimata"
						   attr(object,"tmed") <- "Temperatura media aria"
                           attr(object,"tmax") <- "Temperatura massima aria"
                           attr(object,"tmin") <- "Temperatura minima aria"
                           attr(object,"urel") <- "Umidità relativa  aria"
                           attr(object,"prec") <- "Precipitazione"
                           attr(object,"w_tmed") <- "Temperatura acqua media stimata"
						   attr(object,"tot_evap") <- "Totale evaporazione globale in mg"
						   attr(object,"nrecipients") <- "Numero di trappole Breeding Site"
						   attr(object,"S_day") <- "Soglia in mm per giorno secco"
						   attr(object,"prevdrydays") <- "Serie giorni secchi consecutivi"
						   attr(object,"weight_K") <- "Fattore per alpha weigth"
						   attr(object,"weight_alpha") <- "Fattore per alpha carrying capacity"
						   attr(object,"vol_iniziale") <- "Volume iniziale"
						   attr(object,"biocontainer") <- "Oggetto biocontainer"
                           attr(object,"jdays") <- "Giorni giuliani"
						   attr(object,"daylengths") <- "Durate del giorno"
						   attr(object,"p_diap") <- "Probabilità schiusura diapausanti"
						   attr(object,"hsunrise") <- "Ora dell'alba in HH.DD"
						   attr(object,"hsunset") <- "Ora del tramonto in HH.DD"
						   attr(object,"stazione")<- "Stazione Nome"
						   attr(object,"tipo_contenitore")<- "Tipo_contenitore"
                           attr(object,"lat")<- "Latitudine stazione Sessadecimali WGS 84"
                           attr(object,"lon")<- "Longitudine stazione Sessadecimali WGS 84"	
                           attr(object,"quota")<- "Quota stazione meteo di riferimento metri"	
                           attr(object,"timeformat")<- "Formato dei dati meteo"
						   attr(object,"timeseries")<- "Oggetto timeseries dati"
						   return(object)						   					   
						   }
						   

############################################################################################
# Funzione oggetto popolazione

biopopulation<- function( uova=10000,
			   larve=0,	
               pupe=0,
			   adulti=0,
			   uova_diap=10,
			   egn=63,
			   ef=0.83,
			   lat=43.5,
			   lon=11.27,
			   quota=100,
			   ID=c("NA"),
			   nome_sito=c("NA")
			   ) {
                 object <- list(
			     uova=uova,
			     larve=larve,	
                 pupe=pupe,
			     adulti=adulti,
			     uova_diap=uova_diap,
                 egn=egn,
                 ef=ef,				 
			     lat=lat,
			     lon=lon,
			     quota=quota,
	             ID=ID,
				 nome_sito=nome_sito);
				 
  class(object) <- "biopopulation"
  attr(object,"uova") <- "Numero di uova attuale"
  attr(object,"larve") <- "Numero di larve attuale"
  attr(object,"pupe") <- "Numero di pupe attuale"
  attr(object,"adulti") <- "Numero di adulti depositanti uova attuale"
  attr(object,"uova_diap") <- "Uova Diapausanti"
  attr(object,"egn")<-"Numero uova medio di cova della popolazione"
  attr(object,"ef")<-"Fattore emergenza pupale della popolazione" 
  attr(object,"lat")<-"Latitudine"
  attr(object,"lon")<-"Longitudine"
  attr(object,"quota")<-"Quota"
  attr(object,"ID")<-"ID oggetto"
  attr(object,"nome_sito")<-"Nome del sito"			   			   
  class(object) <- "biopopulation"
  return(object)
}
############################################################################################
# Funzione oggetto contenitore
# vedi selmi trappola da 750 ml 0.75.lt
# va re-inizializzata ogni vola prima del lancio  di un ciclo di simulazione giornaliero
# trappola standard altezza 13 cm D 50

biocontainer<- function(
			   tipo="Trappola REDLAV ITA",
			   nrecipients=50,
			   container_class=1,
			   vol_iniziale_uni=1000,
			   perc_vol_iniz=0.75,
			   D=11,
			   sottovaso=0,
			   fraz_sottovaso=0,
			   sunexposure=0,
			   fraz_copertura=0,
			   lat=43.5,
			   lon=11.27,
			   quota=100,
			   ID=c("NA"),
			   nome_sito=c("NA")			   
			   ) {
			    
				 object <- list(
				 tipo=tipo,
			     container_class=container_class,
				 nrecipients=nrecipients,
			     D=D,
			     sottovaso=sottovaso,
			     fraz_sottovaso=fraz_sottovaso,
			     perc_vol_iniz=perc_vol_iniz,
				 vol_iniziale=vol_iniziale_uni*perc_vol_iniz*nrecipients,
				 sunexposure=sunexposure,
			     fraz_copertura=fraz_copertura,
				 lat=lat,
			     lon=lon,
			     quota=quota,
	             ID=ID,
				 nome_sito=nome_sito);
				 
  class(object) <- "biocontainer"
  attr(object,"Tipo") <- "Tipologia di container"
  attr(object,"container_class") <- "1 per parallelogramma, 2 per cilindro"
  attr(object,"D") <- "Diametro della superficie di base del cilindro in centimetri"
  attr(object,"nrecipients=") <- "Numero di trappole Breeding Sites"
  attr(object,"sottovaso") <- " Presenza Sottovaso 1 se SI, 0 se NO"
  attr(object,"fraz_sottovaso") <- "Frazione decimale di occupazione del vaso nel sottovaso"
  attr(object,"Perc_vol_iniz") <- "Frazione decimale di acqua iniziale sul volume totale(capacita)"
  attr(object,"sunexposure") <- "Frazione decimale di esposizione al sole (0 non esposto, 1 completamente esposto(<0.75)"
  attr(object,"fraz_copertura") <- "Frazione decimale di copertura del contenitore (0 se scoperto, 0.9 se completamente coperto)"
  attr(object,"riempi_svuota") <- "Frazione decimale di riempimento (da 0 a 1) o svuotamento (da -1 a 0) artificiale"
  attr(object,"capacita")<-"Capacità media del contenitore"
  attr(object,"vol_iniziale")<-"Volume acqua dell'insieme dei contenitori"
  attr(object,"tot_evap")<-"Volume acqua nel contenitore persa al giorno dai contenitori"
  attr(object,"vol_finale")<-"Volume acqua finale dei contenitore"
  attr(object,"lat")<-"latitudine"
  attr(object,"lon")<-"longitudine"
  attr(object,"quota")<-"Quota"
  attr(object,"ID")<-"ID oggetto"
  attr(object,"nome_sito")<-"Nome del sito"			   			   
  class(object) <- "biocontainer"
  return(object)
}
############################################################################################
# Funzione oggetto parametri della specie in esame


bioparameters <- function(alfa_l=1.5,
                          alfa_u=0,
                          alfa_a=0,
                         densita=70,
						 sexratio=0.5,
						 DR25_elr=0.24,
						 EA_elr=10798,
						 EI_elr=100000,
						 TI_elr=14184,
						 DR25_lpr=0.2088,
						 EA_lpr=26018,
						 EI_lpr=55990,
						 TI_lpr=304.6,
						 DR25_par=0.384,
						 EA_par=14931,
						 EI_par=-472379,
						 TI_par=148,
						 DR25_ovr1=0.216,
						 EA_ovr1=15725,
						 EI_ovr1=1756481,
						 TI_ovr1=447.2,
						 DR25_ovr2=0.372,
						 EA_ovr2=15725,
						 EI_ovr2=1756481,
						 TI_ovr2=447.2,
						 egn=63,
						 ef=0.83,
						 inib=0.63, 					 
						 nome_specie=c("Albopictus"),
			             genere_specie=c("Aedes"),
			             ordine_specie=c("Diptera"),
						 area_geo="Toscana",			    
			             nome_location=c("NA")
                        	) {
  object <- list(alfa_l=alfa_l,
                 alfa_u=alfa_u,
				 alfa_a=alfa_a,
                 densita=densita,
				 DR25_elr=DR25_elr,
				 EA_elr=EA_elr,
				 EI_elr=EI_elr,
				 TI_elr=TI_elr,
				 DR25_lpr=DR25_lpr,
				 EA_lpr=EA_lpr,
				 EI_lpr=EI_lpr,
				 TI_lpr=TI_lpr,
				 DR25_par=DR25_par,
				 EA_par=EA_par,
				 EI_par=EI_par,
				 TI_par=TI_par,
				 DR25_ovr1=DR25_ovr1,
				 EA_ovr1=EA_ovr1,
				 EI_ovr1=EI_ovr1,
				 TI_ovr1=TI_ovr1,
				 DR25_ovr2=DR25_ovr2,
				 EA_ovr2=EA_ovr2,
				 EI_ovr2=EI_ovr2,
				 TI_ovr2=TI_ovr2,
				 egn=egn,
				 ef=ef,
				 sexratio=sexratio,	 
				 inib=inib,					
				 nome_specie=nome_specie,
			     genere_specie=genere_specie,
			     ordine_specie=ordine_specie,
				 area_geo="Toscana",			    
				 nome_location=nome_location);
  class(object) <- "bioparameters"
  attr(object,"alfa_l") <- "Environmental Competition larva"
  attr(object,"alfa_u") <- "Environmental Competition eggs"
  attr(object,"alfa_a") <- "Environmental Competition adult"
   attr(object,"densita") <- "Densita critica larvale"
  attr(object,"sexratio") <- "Rapporto sessi"
  attr(object,"DR25_elr")<-"Tasso di sviluppo at 25°C Schoofield transferhl(elr) uovo->larva"
  attr(object,"EA_elr")<-"Entalpia di Attivazione Schoofield transferhl(elr) uovo->larva"
  attr(object,"EI_elr")<-"Entalpia di Inattivazione Schoofield transferhl(elr) uovo->larva"
  attr(object,"TI_elr")<-"Temperatura di Inattivazione Schoofield transferhl(elr) uovo->larva"
  attr(object,"DR25_lpr")<-"Tasso di sviluppo at 25°C Schoofield transferlp(lpr) larva->pupa"
  attr(object,"EA_lpr")<-"Entalpia di Attivazione transferlp(lpr) larva->pupa"
  attr(object,"EI_lpr")<-"Entalpia di Inattivazione transferlp(lpr)  larva->pupa"
  attr(object,"TI_lpr")<-"Temperatura di Inattivazione Schoofield transferlp(lpr) larva->pupa" 
  attr(object,"DR25_par")<-"Tasso di sviluppo at 25°C Schoofield transferpa(par)  pupa->dulto"
  attr(object,"EA_par")<-"Entalpia di Attivazione transferpa(par)  pupa->adulto"
  attr(object,"EI_par")<-"Entalpia di Inattivazione transferpa(par) pupa->adulto"
  attr(object,"TI_par")<-"Temperatura di Inattivazione transferpa(par)   pupa->adulto"
  attr(object,"DR25_ovr1")<-"Tasso di sviluppo at 25°C Schoofield ovideposizione(ovr1)  ciclo gonotropico per adulti depositanti uova"
  attr(object,"EA_ovr1")<-"Entalpia di Attivazione ovideposizione(ovr1)ciclo gonotropico per adulti depositanti uova"
  attr(object,"EI_ovr1")<-"Entalpia di Inattivazione ovideposizione(ovr1)  ciclo gonotropico per adulti depositanti uova"
  attr(object,"TI_ovr1")<-"Temperatura di Inattivazione ovideposizione(ovr1)  ciclo gonotropico per adulti depositanti uova"
  attr(object,"DR25_ovr2")<-"Tasso di sviluppo at 25°C Schoofield ovideposizione(ovr2)  ciclo gonotropico per adulti non depositanti uova"
  attr(object,"EA_ovr2")<-"Entalpia di Attivazione ovideposizione(ovr2)  ciclo gonotropico per adulti non depositanti uova"
  attr(object,"EI_ovr2")<-"Entalpia di Inattivazione ovideposizione(ovr2)  ciclo gonotropico per adulti non depositanti uova"
  attr(object,"TI_ovr2")<-"Temperatura di Inattivazione ovideposizione(ovr2) ciclo gonotropico per adulti non depositanti uova"
  attr(object,"egn")<-"Numero medio di uova depositate in un ciclo (83 per cesena e piacenza)"
  attr(object,"ef ")<-"Fattore di emergenza pupale"
  attr(object,"inib")<-"Fattore di inibizione"
  attr(object,"nome_specie")<-"Nome specie dell' organismo"
  attr(object,"genere_specie")<-"Genere della specie"
  attr(object,"ordine_specie")<-"Ordine della specie"
  attr(object,"area_geo")<-"Nome dell'area geografica"			    
  attr(object,"nome_location")<-"Nome del sito"
  class(object) <- "bioparameters"
  return(object)
}
############################################################################################
# Funzione oggetto popolazione

biopopulation<- function( uova=10000,
			   larve=0,	
               pupe=0,
			   adulti=0,
			   uova_diap=10,
			   egn=63,
			   ef=0.83,
			   lat=43.5,
			   lon=11.27,
			   quota=100,
			   ID=c("NA"),
			   nome_sito=c("NA")
			   ) {
                 object <- list(
			     uova=uova,
			     larve=larve,	
                 pupe=pupe,
			     adulti=adulti,
			     uova_diap=uova_diap,
                 egn=egn,
                 ef=ef,				 
			     lat=lat,
			     lon=lon,
			     quota=quota,
	             ID=ID,
				 nome_sito=nome_sito);
				 
  class(object) <- "biopopulation"
  attr(object,"uova") <- "Numero di uova attuale"
  attr(object,"larve") <- "Numero di larve attuale"
  attr(object,"pupe") <- "Numero di pupe attuale"
  attr(object,"adulti") <- "Numero di adulti depositanti uova attuale"
  attr(object,"uova_diap") <- "Uova Diapausanti"
  attr(object,"egn")<-"Numero uova medio di cova della popolazione"
  attr(object,"ef")<-"Fattore emergenza pupale della popolazione" 
  attr(object,"lat")<-"Latitudine"
  attr(object,"lon")<-"Longitudine"
  attr(object,"quota")<-"Quota"
  attr(object,"ID")<-"ID oggetto"
  attr(object,"nome_sito")<-"Nome del sito"			   			   
  class(object) <- "biopopulation"
  return(object)
}




biomodelday<-function(biopopulation,
               bioparameters, 			   
               biocontainer,
			   tmed=24,
               twmed=21,			   
			   State=c(L1=0,L3=0,L4=0,L5=0,L6=0,L7=0,L8=0,L10=0),
               deltatime=1440,
			   nome="Rbiosim"
			   ){
			   require(deSolve);
			   if (class(biopopulation)!="biopopulation") { stop("Object biopopulation argument must be of class biopopulation" )};
 			   if (class(bioparameters)!="bioparameters") { stop("Object bioparameters argument must be of class bioparameters" )};
  			   if (class(biocontainer)!="biocontainer") { stop("Object biocontainer argument must be of class biocontainer" )};
  			   
			   flag_inib=calc_flag_inib(biopopulation$larve,vol_acqua_finale=biocontainer$vol_iniziale,densita=bioparameters$densita);
			   inib=bioparameters$inib; # Condizione inibizione larve	
			   
			   #########################################################################################à
			   # Risoluzione temporale di integrazione
			   Time <- 0:deltatime;
			   
			   
			   #########################################################################################à
			   # calcolo dei tassi per riempire l'oggetto
			   
			   f_ovo_a=ovodepos_g(twmed); # Tacqua aria tasso deposizione uova adulti
			   f_trans_u2l=transferul(twmed,flag_inib); # Tacqua aria tasso uova2larve
               f_trans_l2p=transferlp(twmed); # Tacqua aria tasso larve2pupe
               f_trans_p2a=transferpa(twmed); # Tacqua aria tasso pupe2adulti
               ma=mort_adult(tmed);   # T° aria tasso mortalità adulti  
               mu=mort_uova(twmed); # T° aria tasso mortalità uova  
               mp=mort_pupe(twmed); # tasso mortalità pupale  
               ml=mort_lar(twmed); # tasso mortalità larvale 
			   alfa_lar=bioparameters$alfa_l/biocontainer$nrecipients; # carrying capacity delle larve
			   alfa_egg=bioparameters$alfa_u/biocontainer$nrecipients; # carrying capacity delle uova
			   alfa_ad=bioparameters$alfa_a; # carrying capacity delle larve
			   
			   #########################################################################################à
			                 
			    object <- list(nome=nome,
			                  states=State,
				              times =Time,
							  #########################################################################################################################################
				              # Parametri
							  Parameter <- c(
							  inib=bioparameters$inib, 			
							  f_ovo_a=f_ovo_a,
                              f_trans_u2l=f_trans_u2l,  
                              f_trans_l2p=f_trans_l2p, 
                              f_trans_p2a=f_trans_p2a,  
                              ma=ma,        
                              mu=mu,  
                              mp=mp,  
                              ml=ml,  
			                  alfa_l=alfa_lar, # parametro di competitivita ambientale carrying capacity delle larve
							  alfa_u=alfa_egg, # parametro di competitivita ambientale carrying capacity delle larve
			                  alfa_a=alfa_ad, # parametro di competitivita ambientale carrying capacity delle larve
			                  ef=bioparameters$ef, # fattore di emergenza pupale
						      egn=bioparameters$egn, # uova depositate ciclo
			                  N_ini_a=biopopulation$adulti, # adulti primipare  iniziali presenti
			                  N_ini_p=biopopulation$pupe, # pupe iniziali presenti
			                  N_ini_l=biopopulation$larve,  # larve iniziali presenti
			                  N_ini_u=biopopulation$uova,#uova iniziali presenti
			                  N_u_diap=biopopulation$uova_diap, #uova diapausanti presenti
			                  sexratio=bioparameters$sexratio,
							  deltatime=deltatime),# sex ratio simulazione
							  
							  #########################################################################################################################################
				              # Model Functionals
							  
							  
                aedesmodelday <- function(Time, State, Parameter){
                                       with(as.list(c(State, Parameter)), {
                                       L1=(f_ovo_a*( N_ini_a + (sexratio*L8) - L10 - L1 ) + (alfa_u*(N_ini_u)*(N_ini_u-1)))/ deltatime;                   # Uova deposte da adulti femmine
                                       L3=((mu*(N_ini_u+egn*L1-L3))/ deltatime);                                                                                  # Uova che muoiono
                                       L4=(f_trans_u2l*(N_ini_u+egn*L1-L3-L4))/deltatime;                                                                         # Uova che diventano larve 			
                                       L5=((f_trans_l2p*(N_ini_l+L4-L5-L6))/deltatime);                                                                       # Larve che diventano pupe						   
                                       L6=(ml*(N_ini_l+L4-L5-L6) + (alfa_l*(N_ini_l+L4-L5-L6)*(N_ini_l+L4-L5-L6)))/deltatime;                                 # Mortalità larve 
                                       L7=((mp+sexratio*f_trans_p2a)*(N_ini_p+L5-L7-L8)/deltatime);                                                           # Mortalità pupe
                                       L8=((f_trans_p2a *(sexratio*N_ini_p)+L5-L7-L8)+ (alfa_a*(N_ini_a)*(N_ini_a-1)))/deltatime;                             # Pupe che diventano adulti femmine
									   L10=(ma*(N_ini_a+(sexratio*L8)-L10-L1)/deltatime);                                                                     # Mortalità  di adulti femmine 
									   return(list(c(L1,L3,L4,L5,L6,L7,L8,L10)))})} ,
									   
									   
							        #########################################################################################################################################
				             	    aedesfunc_day <- ode(func = aedesmodelday , y = State, parms = Parameter,times = Time, method = rkMethod("rk45dp7"))
                            									  
				 );
  class(object) <- "biomodelday"
  attr(object,"nome") <- "Denominazione modello"
  attr(object,"Time") <- "Vettore iterazione giornaliero dipendete parametro deltat"
  attr(object,"parameter") <- "Parametri generali  modello"
  attr(object,"tmed") <- "Temperatura media aria"
  attr(object,"twmed") <- "Temperatura media acqua"
  attr(object,"inib") <- "Fattore di inibizione crescita larvale"
  attr(object,"flag_inib") <- "Presenza inibizione crescita larvale"
  attr(object,"aedesmodelday") <- "Instaziamento giornaliero modello"
  attr(object,"aedesfunc_day") <- "Funzione lancio modello"							  
  class(object) <- "biomodelday"
  return(object)
}
	
############################################################################################
# Funzione per la gestione degli output  dell' oggetto biomodelday. usata in biosim_run()
						   
elab_day_sim<- function(elab_obj,biopopulation,biocontainer, p100_su=1,n_resampling=25,times=1440)
                           {
						   
						    if (class(biopopulation)!="biopopulation") { stop("Object biopopulation argument must be of class biopopulation" )};
 						 	if (class(biocontainer)!="biocontainer") { stop("Object biocontainer argument must be of class biocontainer" )};
 						 	
						    pop_delay<-as.data.frame(cbind(
                            n1 = rpois(n_resampling,as.data.frame(elab_obj[[6]])$L1[times+1]),
                            n3 = rpois(n_resampling,as.data.frame(elab_obj[[6]])$L3[times+1]),
                            n4 = rpois(n_resampling,as.data.frame(elab_obj[[6]])$L4[times+1]),
                            n5 = rpois(n_resampling,as.data.frame(elab_obj[[6]])$L5[times+1]),
                            n6 = rpois(n_resampling,as.data.frame(elab_obj[[6]])$L6[times+1]),
                            n7 = rpois(n_resampling,as.data.frame(elab_obj[[6]])$L7[times+1]),
                            n8 = rpois(n_resampling,as.data.frame(elab_obj[[6]])$L8[times+1]),
                            n10 = rpois(n_resampling,as.data.frame(elab_obj[[6]])$L10[times+1])))
                           
                            pop_delay_means<-as.data.frame(t(apply(pop_delay, 2, mean)))
                            pop_delay_sd<-as.data.frame(t((apply(pop_delay, 2, sd))))

                            ######################################################################################à
                            # Inizializzazione delle popolazione
							egn=biopopulation$egn;
                            N_ini_a=biopopulation$adulti;
							N_ini_p=biopopulation$pupe;
							N_ini_l=biopopulation$larve;
							N_ini_u=biopopulation$uova;
							N_ini_u_diap=biopopulation$uova_diap;
							N_fin_a=biopopulation$adulti;
							N_fin_p=biopopulation$pupe;
							N_fin_l=biopopulation$larve;
							N_fin_u=biopopulation$uova;
							
							nrecipients=biocontainer$nrecipients;
							p100_su=p100_su/100;

							######################################
							# Adulti
							N_fin_a=0
                            N_new_a=pop_delay_means$n8-pop_delay_means$n10;              
							N_fin_a=N_ini_a+N_new_a;
                                                                                             # Adulti totali
                            if (N_fin_a < 1) { N_fin_a=1}
							
							######################################
							# Pupe
                            N_fin_p							
							N_new_p=pop_delay_means$n5-pop_delay_means$n7-pop_delay_means$n8;                     # Pupe  presenti
							N_fin_p=N_ini_p+N_new_p;                              
                            if (N_fin_p < 0) { N_fin_p=0;}
							
							######################################
							# Larve	
							N_fin_l=0
						    N_new_l=pop_delay_means$n4-pop_delay_means$n5-pop_delay_means$n6;                     # Larve  presenti
                            N_fin_l=N_ini_l+N_new_l;                                                                               
                            if (N_fin_l < 0) { N_fin_l=0;}
							
                            ######################################
							# Uova	
						    N_fin_u=0
							N_new_u=(egn*pop_delay_means$n1)-pop_delay_means$n3-pop_delay_means$n4+N_ini_u_diap*(p100_su);           # nuova uova presenti
                            N_fin_u=N_ini_u+N_new_u*(1-p100_su);
							
							if (N_fin_u < 0) { N_fin_u=0;}
							
                            ######################################
							# Uova normalizzate sul recipiente.	
						   
							N_ini_u_mean_rec=N_ini_u/nrecipients;                                                                   # uova medie attese
                            N_fin_u_mean_rec=N_fin_u/nrecipients;
                         
                            #############################################################################################################################################
                            # Per controllo nel ciclo per evitare l'estinzione.
                            # evito l'estinsione riportando ad 1 il valore della popolazione di uova qualore risulti <= 0
                            #############################################################################################################################################
                            
							######################################
							# Uova diapausanti presenti
                            N_fin_diap=0
							N_new_diap=N_new_u*(1-p100_su);
							N_fin_diap=N_new_diap-N_ini_u_diap*(p100_su);
							if (N_fin_diap <0) {N_fin_diap=0}
                            ######################################
						    # Variabilità delle popolazioni
                            
							sd_u=pop_delay_sd$n1;	
                            sd_l=pop_delay_sd$n4;	
                            sd_p=pop_delay_sd$n5;	
                            sd_a=pop_delay_sd$n8;
							
							######################################
						    # Variabilità delle popolazioni
                           
							day_new_final<-round(data.frame(N_new_u=N_new_u,
															N_new_l=N_new_l,
															N_new_p=N_new_p,
							                                N_new_a=N_new_a));
							
							day_state_final<-round(data.frame(N_fin_u=N_fin_u,
															  N_fin_l=N_fin_l,
															  N_fin_p=N_fin_p,
															  N_fin_a=N_fin_a));
							
							day_sd<-round(data.frame(sd_u=sd_u,
							                         sd_l=sd_l,
													 sd_p=sd_p,
													 sd_a=sd_a));
													 
						    res=list(day_new_final=day_new_final,
							         day_state_final=day_state_final,
									 day_sd=day_sd,
									 N_ini_u_mean_rec=as.numeric(N_ini_u_mean_rec),
									 N_fin_u_mean_rec=as.numeric(N_fin_u_mean_rec),
									 N_fin_diap=N_fin_diap)
                           
						   return( res);
					   
						   }	
						   
						   
#####################################################################################################################################################
						   						   
rbiosim_run <- function(biometeo,biopop,bioparameters,biocontainer,volfix=1,soglia_pioggia=4,delta_ts=1440,resamples=25)
                           {
						    if (class(biometeo)!="biometeo") { stop("Object biometeo argument must be of class biometeo" )};
 						    if (class(biopop)!="biopopulation") { stop("Object biopopulation argument must be of class biopopulation" )};
 						 	if (class(biocontainer)!="biocontainer") { stop("Object biocontainer argument must be of class biocontainer" )};
 						 	if (class(bioparameters)!="bioparameters") { stop("Object bioparameters argument must be of class bioparameters" )};
 						 	
							####################################################################
							# Preparazione degli oggetti utili per la simulazione
							res <- list();
							
							####################################################################
							# Carica gli oggetti nella simulazione 
							
							res$biometeo=biometeo;  # oggetto biometeo per le variabili di uso durante la simulazione
							
							res$bioparameters=bioparameters; # oggetto bioparametri della popolazione dei ditteri 
							
							res$biocontainer=biocontainer; # oggetto contenitori
							
							####################################################################
							# Stabilisci che i breeding sites sono a volume 
							
							res$volfix=volfix; # volume di acqua fisso 
							
							########################################################################################
							# Stabilisci la lunghezza del periodo di simulazione in giorni dall'oggetto biometeo
							# Riempi con lo stesso volume i giorni della simulazione ipotizzando un volume fisso
							
							res$days=length(biometeo$tmed_est); 
							
							res$vol_finale[1:res$days]=rep(biometeo$vol_iniziale,res$days); 
							
							########################################################################################
							# Crea una variabile volume corrente
							
							vol_current=biometeo$vol_iniziale;
							
							########################################################################################
							# Inizializza una lista per mettere i risultati
							
							temp=list(sim=list(pop_new_day=list(),
							              pop_fin_day=list(),
										  pop_sd_day=list()))
							
							
							#################################################################################################################################################à
							# Esegui la simulazione lungo i giorni dell'oggetto biometeo.
							
							for (k in 1:res$days) {
												   
												   #######################################################
												   
												   run_bioparameters=bioparameters;
												   
												   vol_current=vol_current-biometeo$tot_evap[k]
												   
												   if (biometeo$prec[k]>=soglia_pioggia) {vol_current=biometeo$vol_iniziale; 
												                                          run_bioparameters$alfa0=run_bioparameters$alfa0*biometeo$weight_alpha[k]}
												   
												   res$vol_finale[k]=vol_current
												   
												
												   if (volfix==1) {biocontainer$vol_iniziale=biometeo$vol_iniziale} else { biocontainer$vol_iniziale=vol_current;
												                                                                           run_bioparameters$alfa0=bioparameters$alfa0}
							                       
												   ##########################################################################################################
												   
												   biomodel=biomodelday(biopop,run_bioparameters,biocontainer,tmed=biometeo$tmed_est[k],twmed=biometeo$w_tmed[k],deltatime=delta_ts)
							                     
  												   restemp=elab_day_sim(biomodel,biopop,biocontainer,p100_su=biometeo$p_diap[k],n_resampling=resamples,times=delta_ts);
							                       
												   
												   ##########################################################################################################
												   
												   
												   temp$sim$pop_new_day[[k]]=restemp$day_new_final;
												   temp$sim$pop_fin_day[[k]]=restemp$day_state_final;
												   temp$sim$pop_sd_day[[k]]=restemp$day_sd;
												   temp$sim$N_fin_diap[[k]]=restemp$N_fin_diap;
												   
												   biopop$uova=as.numeric(restemp$day_state_final$N_fin_u)
												   biopop$larve=as.numeric(restemp$day_state_final$N_fin_l)
												   biopop$pupe=as.numeric(restemp$day_state_final$N_fin_p)
												   biopop$adulti=as.numeric(restemp$day_state_final$N_fin_a)
												   biopop$uova_diap=as.numeric(restemp$N_fin_diap)
																		
							}
 	                       
						   
						    #########################################################################################################
							#
							#
							
							 res$pop_new_day<-data.frame(matrix(unlist(temp$sim$pop_new_day), nrow=length(biometeo$tmed_est), byrow=T))
							 names(res$pop_new_day)<-c("N_u","N_l","N_p","N_a")
							 res$pop_fin_day<-data.frame(matrix(unlist(temp$sim$pop_fin_day), nrow=length(biometeo$tmed_est), byrow=T))
							 names(res$pop_fin_day)<-c("N_u","N_l","N_p","N_a")
							 res$pop_sd_day<-data.frame(matrix(unlist(temp$sim$pop_sd_day), nrow=length(biometeo$tmed_est), byrow=T))
							 names(res$pop_sd_day)=c("SD_u","SD_l","SD_p","SD_a")
							 res$N_fin_diap<-data.frame(N_fin_diap=unlist(temp$sim$N_fin_diap))
							 
							#########################################################################################################
							#
							
							return(res)
						   					   
						   }
						   




