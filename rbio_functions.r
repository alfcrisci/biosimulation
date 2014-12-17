#####################################################################################################################################################
# Funzioni per il lancio di una simulazione generata da file meteo.

drydaycons<-function(prec,S=1) {
                                   prec[prec<S]<-0
                                   prec[prec>=S]<-1
                                   prec.rle=rle(prec);
                                   res=NULL;
                                   for (k in 1:length(prec.rle$lengths)) { res=append(res,1:prec.rle$lengths[k])}
                                   res[which(prec>=1)]=0
                                   return(res)
						       }
						   
						   
weigthrain<-function(drylength,S=5) {
                                   res=NULL;
                                   res=round(1-exp(drylength/S-drylength),digits=2) 
								   return(res)
						            }	


############################################################################################
# Capacità dei contenitori # Le misure dimensionali sono in metri
 

capacita_cont<-function(SB,H,sottovaso,fraz_sottovaso)
{
 
 capacita=SB*H*0.001;  # in litri H altezza del contenitore

 if (sottovaso==1) {capacita=capacita*(1-fraz_sottovaso) }
	
 # capacità ridotta per la presenza del vaso

return(capacita);
}

############################################################################################
############################################################################################
# Funzioni psicrometriche utili per la biofisica dei contenitori e dei breeding site in generale

############################################################################################
# Pressione di vapore saturo

pvapsat<-function(ta=24)
{
   # Esce in hpa o millibar
   
   es=6.1078*10^((7.5*ta)/(237.7+ta));
  
  return(es);
}

############################################################################################
# Coefficente di diffusione massico di Fick  m2 al secondo vapore in funzione di T°

diffvap_t<-function(ta=24)
{  

   diff_t=(2.05*(10**-5)) * ((ta+273)/273)**2.072;
  
  return(diff_t);
}


############################################################################################
# Differenza di densita' di vapore su superficie a pelo libero

delta_d<-function(ta=24,twmed=21,rh=60)
{  
   R=461.48; #J K-l Kg-l.
	g_standard=9.80665;

   deltad=((100*pvapsat(twmed))/(R*(twmed+273)))-(100*(rh/100)*pvapsat(ta))/(R*(ta+273));

  return(deltad)
}

############################################################################################
# Differenza di pressione in pascal di vapore su superficie a pelo libero

delta_p<-function(ta=24,twmed=21,rh=60)
{  
    
   delta_p=((100*pvapsat(twmed))-(100*(rh/100)*pvapsat(ta)))

  return(delta_p)
}

############################################################################################
# Differenza di densita di vapore su superficie a pelo libero


delta_md<-function(ta=24,twmed=21,rh=60)
{   
    R=461.48; #J K-l Kg-l.
	g_standard=9.80665
    deltad=((100*pvapsat(twmed))/(R*(twmed+273)))+(100*(rh/100)*pvapsat(ta))/(R*(ta+273))/2;
    return(deltad)
}

############################################################################################
# Stima dei parametri utili per calcolare l'evaporazione di masse d'acqua in funzione della geometria dei contenitori

evap_number<-function(ta=24,twmed=21,vvent=0.5,rh=60,L)
{ 

A=3.1415*(L/2)**2;
mv=18.016; #kg/kmole
d_aria_20=0.8216; #kg/m3
cp_aria=1.005; #KJ/kg K° #calore specifico aria
visc_aria=1.583*(10**-5);# viscosità aria mi greco
g_standard=9.80665;
dilataz_aria=1/(ta+273.16);
Rzero=8314/mv; #J K-l kg-l.
entalpiavap=40.8; #kJ/mol
res=list()
res$reynolds=vvent*L/visc_aria;
res$raylegh=g_standard*((delta_d(ta,twmed,rh)*(L**3))/(delta_md(ta,twmed,rh)*visc_aria*diffvap_t(ta)));
res$schmidt=0.6;
res$sherwood=0.230*((res$schmidt)**(1/3))*res$raylegh**(0.321);
res$hm=((res$sherwood*diffvap_t(ta))/L);
res$mev=res$hm*A*(delta_p(ta,twmed,rh)/(461*(twmed+273))); # kg/sec
res$mevday=res$mev*3600*24*1000; #mg day
return(res)
}

 
############################################################################################
# Stima del coefficiente di trasporto massico  masse d'acqua in contenitori in funzione del T° aria, T° acqua, Umidità relativa e velocità del vento.


ktrasp_real<-function(ta=23,twmed=21,rh=50,vvent=0.5)
{  ktrasp_vap=3.4*(10**-8);
   k_trasp=ktrasp_vap*delta_p(ta,twmed,rh);
   return(k_trasp);
}


############################################################################################
# Stima del punto di rugiada

dewpoint<-function(ta=24,rh=50) {
	dpt=(17.271*(ta+273))/(237.7+(ta+273))-log(rh/100)
	return (dpt);
}

############################################################################################
# Stima del deficit di saturazione in hpa

degsat<-function(ta=24,rh=50) {

   # Calcola millibar che mancano alla saturazione data una t°.
   # Gli intervalli di secchezza sono fra 10 fino a 30. 
   # Sotto 10 è umido ma più cresce il deficit e più è secca l'aria.   
   
	mu=pvapsat(ta)*(1-rh/100);
	
	return (mu);
}



################################################################################################################
# T. MIN RAGGIUNTA DALL'ACQUA IN BASE ALLA TEMP. MAX E MIN DELL'ARIA E ALL'ESPOSIZIONE  # Focks CIMSIM 

 temperatura_min<-function( t_max_aria,  t_min_aria,  sunesp)
{
       # La temp. min raggiunta dall'acqua : t_min_acqua;
	   
       t_min_acqua=5.02-1.36*sunesp+0.81*t_min_aria+0.001*t_max_aria*t_max_aria; 
      
       return(t_min_acqua);
}

#####################################################################################################################################àààà
# Calcolo temperatura media della acqua fonte Merler

twmedacqua<-function(tmed=24,tmin=12,rel=60) {
t_med_acqua=NA;
if ((mese >= 3) && (mese<=10)) {
				t_med_acqua = (tmed * 0.642) + (tmin * 0.201) + (rel * 0.037) + 0.817;
			        }
			   else {
				t_med_acqua = tmed;
				 }
return(t_med_acqua); 
}

#es :twmedacqua(22,16,50,5)
#es :ef_iul = emer_factor_iul();
#es :ma_iul = mort_adult_iul()
	
################################################################################################################
# T. MAX RAGGIUNTA DALL'ACQUA IN BASE ALLA TEMP. MAX E MIN DELL'ARIA E ALL'ESPOSIZIONE SECONDO    #### Focks CIMSIM ####
# FUNZIONE CHE IN BASE ALL'UMIDTA' REL E ALL'ESPOSIZIONE AL SOLE MISURA IL VOLUME D'ACQUA PERSO PER EVAPORAZIONE

perdite_evap_focks<-function(SB,umi_rel,sunesp)

              {     
			  # Focks CIMSIM #    
              # cm d'acqua al giorno evaporati 
              # cm persi per l'evaporazione sono : cm_persi);
              # litri persi per evaporazione : litri_evap);

              cm_persi= 0.93+0.28*sunesp-0.01*umi_rel; 
              litri_evap= cm_persi*SB*0.001;

    return(litri_evap);
}

 temperatura_focks_max<-function(t_max_aria,t_min_aria,esp)
{
      # La temp. max raggiunta dall'acqua e' %f gradi Celsius\n" , t_max_acqua);
       t_max_acqua=15.03+0.27*t_min_aria+0.01*t_max_aria*t_max_aria+7.69*(esp^2); 
	   return(t_max_acqua);
}

###################################################################################################################
# Funzione di calcolo per il giorno successivo.

calcolo_volume_fin<-function(vol_acqua_iniz,litri_evap,fraz_copertura,riempi_svuota,capacita)
{
      # La frazione di copertura del contenitore  fraz_copertura;
      # La frazione di riempimento o svuotamento giornaliera :riempi_svuota;
      # Il volume finale di acqua e' litri vol_acqua_finale;

       vol_acqua_finale=vol_acqua_iniz-(litri_evap*(1-fraz_copertura))+vol_acqua_iniz*riempi_svuota;

      if (vol_acqua_finale>capacita) {vol_acqua_finale=capacita; } # controlla che l'acqua finale non superi la capacita. #

      if (vol_acqua_finale<0) {vol_acqua_finale=0; } # controlla che l'acqua finale non superi la capacita#
      return(vol_acqua_finale);
}

#####################################################################################################################
# Tasso di sviluppo generale Dati SOLARI OTERO 

 DR_calc<-function( temp, DR25,EA,EI,TI,t25=298)
			  {
	          temp = temp + 273.5;
	          DR = ((DR25 * (temp/t25) * exp((EA/R) * (1/t25 - 1/temp))) / (1 + exp((EI/R) * (1/TI - 1/temp))));
           	  return(DR);
             }
			 
				 
#####################################################################################################################################àààà
# Funzione per la modifica dei parametri dei tassi di transizioni fra stati di  popolazione Dati Merler

#######################################################################################################################
#  Mortality rates from Merler FBK



mort_uova<-function(t_med_acqua=24){
	
	me=1.220816
	if((t_med_acqua < 35.0)||(t_med_acqua > 15.0))  
      {
   	   me = 506.000 - (506.000*exp(-((t_med_acqua-25)/27.300)*((t_med_acqua-25)/27.300)*((t_med_acqua-25)/27.300)*((t_med_acqua-25)/27.300)*((t_med_acqua-25)/27.300)*((t_med_acqua-25)/27.300)));
      }
	  return(me)
	}

mort_lar<-function(t_med_acqua=24){
	
	mla = 0.029 + (858*exp(t_med_acqua-43.4));
	return(mla)
	}
	
	
mort_pupe<-function(t_med_acqua=24){
	# temperatura acqua
	
	mpu = 0.021 + (37.000*exp(t_med_acqua-36.800));	
	return(mpu)
	}



mort_adult<-function(tmed=24){
	# temperatura aria
	ma_mer = 0.031 + (95820*exp(tmed - 50.400));
	return(ma_mer)
	
}

mort_adult_deltasd<-function(sd=30){
	# temperatura aria
	delta_mer = 0.04*(sd-10)/(30-10)
	return(delta_mer)
	
}

ovodepos_g<-function(tmed=24)
                            { ovodepos_g1_day=(1 / ((0.046*tmed^2) - (2.770*tmed) + 45.300));
                              return(ovodepos_g1_day);
                            }
                            
						
							
							
transferlp <-function(tmed_acqua=24)
                            { transferlp_day=1 / ((0.120*tmed_acqua*tmed_acqua) - (6.600*tmed_acqua) + 98.000);
                             return(transferlp_day);
                             }
transferpa <-function(tmed_acqua=24)
                            { transferpa_day=1 / ((0.027*tmed_acqua*tmed_acqua) - (1.700*tmed_acqua) + 27.700);	
                             return(transferpa_day);
                             }
							 
transferul<-function(tmed_acqua=24,flag_I=0,inib=0.63)
                            { transferhl_day=1 / (6.900 - (4.000*exp(-(((tmed_acqua-20)/4.100)*((tmed_acqua-20)/4.100)))))
							  if(flag_I==1)
                                          { transferhl_day=transferhl_day*(1-inib);
                                          }
                             return(transferhl_day);
                             }
							 
							 
 

	  
#######################################################################################################################
# Definisco tassi di sviluppo attraverso Dati Juliano


emer_factor_iul<-function(temp=24){
	
	ma_oss=c(0.0087,0.0114,0.0124)
	t_oss2=c(22.0,24.0,26.0)		
	f <- splinefun(t_oss2, ma_oss)
	emer_iul=f(temp)
	return(emer_iul)
	
}

mort_adult_iul<-function(tmed=24){
	# temperatura aria
	ma_oss=c(0.0087,0.0114,0.0124)
	t_oss2=c(22.0,24.0,26.0)		
	f <- splinefun(t_oss2, ma_oss)
	ma_iul=f(tmed)
	return(ma_iul)
	
}

#######################################################################################################################
# Definisco un fattore indicante  inizione sviluppo  larvale
	
calc_flag_inib<-function(N_lar_ini,vol_acqua_finale=2,densita=70){
	
	flag_inib="0";
	if(N_lar_ini >= densita*vol_acqua_finale)
                { flag_inib="1";
                }  
	return(flag_inib)
	}	
	
#######################################################################################################################
# Definizione della carryng capacity larvale su set di contenitori:  carico massimo
	
calc_cmnl<-function(alfa0,nrecipients){
	# carryng capacity factor and BS number
	cmnl=(alfa0/nrecipients)
	return(cmnl)
	}	
	


###############################################################################################################################################
# Funzioni per uova diapausanti

mortality_diapause<-function(tmin=2)
{
	tmin=round(tmin,0);
	if(tmin>2){
		dead=0;
		return(dead);
	}
	else if(tmin> -10 && tmin <=2){
		dead=(-1.2346 * tmin) + 8.179;
		return(dead);
	}
	else if(tmin<=-10 && tmin > -14){
		dead= (-18.75 * tmin) - 166.67;
		return(dead);
	}
	else if(tmin< -14){
		dead=100;
		return(dead);
	}	
}
###########################################################################################################################################################
# Funzione  per calcolo della probabilità  schiusura di uova diapausanti in base alle ore di luce.

ore_schiusura<-function(x) {
                        1.000143/(1+exp((11.894826 -x)/0.580523 ))
                       }
###########################################################################################################################################################

###########################################################################################################################################################
# Funzioni per calcolo variabili di calendario

jdday_str<-function(date) {
                        strptime(date, "%m/%d/%Y")$yday+1
                       }
jdday<-function(dates) {
                        strptime(dates, "%Y-%m-%d")$yday+1
                       }
                       
jdday_diff<-function(data_i,data_f) {
                        strptime(data_f, "%m/%d/%Y")$yday+1-strptime(data_1, "%m/%d/%Y")$yday+1
                       }

rad<-function(x)pi*x/180					   
###########################################################################################################################################################
# Funzioni per calcolo effemeridi

suncalc<-function(d,Lat=43.77,Long=11.10){

  ## d is the day of year
  ## Lat is latitude in decimal degrees
  ## Long is longitude in decimal degrees (negative == West)
  
  ## This method is copied from:
  ## Teets, D.A. 2003. Predicting sunrise and sunset times.
  ## The College Mathematics Journal 34(4):317-321.
 
  ## At the default location the estimates of sunrise and sunset are within
  ## seven minutes of the correct times (http://aa.usno.navy.mil/data/docs/RS_OneYear.php)
  ## with a mean of 2.4 minutes error.

  ## Function to convert degrees to radians

  
  ##Radius of the earth (km)
  R=6378
  
  ##Radians between the xy-plane and the ecliptic plane
  epsilon=rad(23.45)

  ##Convert observer's latitude to radians
  L=rad(Lat)

  ## Calculate offset of sunrise based on longitude (min)
  ## If Long is negative, then the mod represents degrees West of
  ## a standard time meridian, so timing of sunrise and sunset should
  ## be made later.
  # es :suncalc(jdday("03/23/2012"),Lat=43.1,Long=10.75)

  timezone = -4*(abs(Long)%%15)*sign(Long)

  ## The earth's mean distance from the sun (km)
  r = 149598000

  theta = 2*pi/365.25*(d-80)

  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)

  t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
  
  ##a kludge adjustment for the radius of the sun
  that = t0+5 

  ## Adjust "noon" for the fact that the earth's orbit is not circular:
  n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)

  ## now sunrise and sunset are:
  sunrise = (n-that+timezone)/60
  sunset = (n+that+timezone)/60
  daylenght=sunset-sunrise
  
  return(list("sunrise" = sunrise,"sunset" = sunset,"daylength" = daylenght))
}

############################################################################################
# Funzione media mobile

running_mean<-function( X, fun=mean, width=min(length(X),20),
                     allow.fewer=FALSE,...)
{
  n<-length(X)

  from  <-  sapply( (1:n) - width + 1, function(x) max(x,1) )
  to    <-  1:n

  elements  <- apply(cbind(from,to), 1,function(x) seq(x[1], x[2]) )

  if(is.matrix(elements))
    elements  <- as.data.frame(elements)

  funct<-function(which,what,fun,...) fun(what[which],...)

  Xvar<-sapply(elements, funct, what=X, fun=fun, ...)
  names(Xvar) <- paste(from,to,sep=":")

  if(!allow.fewer)
    Xvar[1:(width-1)]  <- NA

  return(as.vector(Xvar))
}