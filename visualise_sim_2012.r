#########################################################
# Load libraries

library(rgl)
library(r2stl)
library(zoo)
library(xts)
library(fBasics)
library(plot3D)
library(geometry)
library(xtsExtra)

#########################################################
# Setup 

setwd("D:\\lav_selmi_manfredi")

source("utils.r")

sel_row=c(1,5,8,12,16,31,33,34,42,50)
metadati_sim=read.csv("metadati_sim.csv")

#########################################################
# Analisys of simulations

sim_results=readRDS("simulation_selected.rds")

l_sim_results=length(sim_results)/9

res_minimal=list()

# 1 Metadati
# 2 Trappola
# 3 popolazione
# 4 Parametri di simulazione
# 5 Serie annuale di simulazione
# 6 Serie Osservazioni
# 7 Confronto 
# 8 Verifica $MAE $ME sqrt($MSE)
# 9 successo

ind_meta=seq(1,length(sim_results),9)

res_sim=data.frame(nome_loc=rep(NA,length(ind_meta)),
            lat=rep(NA,length(ind_meta)),
		    lon=rep(NA,length(ind_meta)),
            alfa_l=rep(NA,length(ind_meta)),
		    alfa_u=rep(NA,length(ind_meta)),
		    alfa_a=rep(NA,length(ind_meta)),	
            MAE_sim=rep(NA,length(ind_meta)),
		    ME_sim=rep(NA,length(ind_meta)),
		    RMSE_sim=rep(NA,length(ind_meta)),
		    success_sim=rep(NA,length(ind_meta)),
		    ind_data_confronto=rep(NA,length(ind_meta))
		   )
		   


#########################################################


j=1
for (i in ind_meta) {
 res_sim$nome_loc[j]=as.character(sim_results[i]$metadati$nome)
 res_sim$lat=sim_results[i]$metadati$lat
 res_sim$lon=sim_results[i]$metadati$lon
 res_sim$alfa_l[j]=sim_results[i+3]$biopar$alfa_l
 res_sim$alfa_u[j]=sim_results[i+3]$biopar$alfa_u
 res_sim$alfa_a[j]=sim_results[i+3]$biopar$alfa_a
 res_sim$MAE_sim[j]=sim_results[i+7]$verify$MAE
 res_sim$ME_sim[j]=sim_results[i+7]$verify$ME
 res_sim$RMSE_sim[j]=sqrt(sim_results[i+7]$verify$MSE)
 res_sim$success_sim[j]=sim_results[i+8]
 res_sim$ind_data_confronto[j]=i+6
 j=j+1
 }

 
##########################################################################################################################
 
res_sim_location=split(res_sim,res_sim$nome_loc)

i=1
for ( i in seq(res_sim_location)) {

##########################################################################################################################
######################################################################################################################################
# Adulti

open3d()


plot_coords_l_a=xyz.coords(res_sim_location[[i]]$alfa_l, 
                           res_sim_location[[i]]$alfa_a,
						   res_sim_location[[i]]$ME_sim, 
						   xlab = "alpha_l", ylab = "alpha_a", 
						   zlab = paste("Mae Adulti",res_sim_location[[1]]$nome[1]))
						   
ME_la<- matrix(plot_coords_l_a$z, 
                nrow=length(unique(plot_coords_l_a$x)), 
				ncol=length(unique(plot_coords_l_a$y)))

zlim <- range(plot_coords_l_a$z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- heat.colors(zlen,alpha=0) # height color lookup table
color_ME <- colorlut[ plot_coords_l_a$z-zlim[1]+1 ] # assign colors to heights for each point

				
plot3d(plot_coords_l_a$x, 
        plot_coords_l_a$y,
	    plot_coords_l_a$z, 
	    aspect=c(1, 1, 0.5), 
	    col = color_ME, 
	     xlab = "alpha_l", 
	     ylab = "alpha_a", 
	    zlab = paste("Mae",res_sim_location[[1]]$nome[1]))

writeWebGL(dir = "webGL", filename = file.path("webGL", paste0(res_sim_location[[i]]$nome[1],"l_a_plot.html")))
	   
persp3d(x=unique(plot_coords_l_a$x),
                  y=unique(plot_coords_l_a$y),
				  ME_la, aspect=c(1, 1, 0.5), 
				  col = color_ME, 
				  xlab = "alpha_l", 
				  ylab = "alpha_a", 
				  zlab = paste("Mae",res_sim_location[[1]]$nome[1]))
				  
				  
writeWebGL(dir = "webGL", filename = file.path("webGL", paste0(res_sim_location[[i]]$nome[1],"l_a_persp.html")))


r2stl(x=unique(plot_coords_l_a$x), 
               y=unique(plot_coords_l_a$y),
			   ME_la, 
			   filename=paste0(res_sim_location[[i]]$nome[1],"l_a_surface.stl"), show.persp=FALSE)

######################################################################################################################################
# grid for plotting function & function for predicted surface from model &# predicted values in form that persp() can use

x <- seq(min(res_sim_location[[i]]$alfa_l), max(res_sim_location[[i]]$alfa_l), (diff(range(res_sim_location[[i]]$alfa_l))/50))
y <- seq(min(res_sim_location[[i]]$alfa_a), max(res_sim_location[[i]]$alfa_a), (diff(range(res_sim_location[[i]]$alfa_a))/50))
model = lm(formula = res_sim_location[[1]]$ME_sim ~ (res_sim_location[[i]]$alfa_l+res_sim_location[[i]]$alfa_a)^2 
                                                     + I(res_sim_location[[i]]$alfa_l^2) 
													 + I(res_sim_location[[i]]$alfa_a^2),  
													 data = res_sim_location[[i]])

f <- function(x, y) { cbind(1,x,y,x^2,y^2,x*y) %*% coef(model) }
z <- outer(x, y, f)

zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- heat.colors(zlen,alpha=0) # height color lookup table
color_ME <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point

persp3d(x, y, z, color=color_ME, alpha=0.75,aspect=c(1, 1, 0.5), back="lines")

writeWebGL(dir = "webGL", filename = file.path("webGL", paste0(res_sim_location[[i]]$nome[1],"l_a_qsurface.html")))

r2stl(x,y,z, filename=paste0(res_sim_location[[i]]$nome[1],"_l_a_qsurface.stl"),show.persp=FALSE)

#############################################################################################################################################
# uova


plot_coords_l_u=xyz.coords(res_sim_location[[i]]$alfa_l, 
                            res_sim_location[[i]]$alfa_u,
						    res_sim_location[[i]]$ME_sim, 
						    xlab = "alpha_l", ylab = "alpha_u", 
						    zlab = paste("Mae Uova",res_sim_location[[i]]$nome[1]))
						   
ME_lu<- matrix(plot_coords_l_u$z, 
                nrow=length(unique(plot_coords_l_a$x)), 
				ncol=length(unique(plot_coords_l_a$y)))

zlim <- range(plot_coords_l_u$z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- heat.colors(zlen,alpha=0) # height color lookup table
color_ME <- colorlut[ plot_coords_l_u$z-zlim[1]+1 ] # assign colors to heights for each point

				
plot3d(plot_coords_l_a$x, 
       plot_coords_l_a$y,
	   plot_coords_l_a$z, 
	   aspect=c(1, 1, 0.5), 
	   col = color_ME, 
	   xlab = "alpha_l", 
	   ylab = "alpha_u", 
	   zlab = paste("Mae",res_sim_location[[1]]$nome[1]))

writeWebGL(dir = "webGL", filename = file.path("webGL", paste0(res_sim_location[[i]]$nome[1],"l_u_plot.html")))
	   
persp3d(x=unique(plot_coords_l_a$x),
                  y=unique(plot_coords_l_a$y),
				  ME_la, aspect=c(1, 1, 0.5), 
				  col = color_ME, 
				  xlab = "alpha_l", 
				  ylab = "alpha_u", 
				  zlab = paste("Mae",res_sim_location[[1]]$nome[1]))

writeWebGL(dir = "webGL", filename = file.path("webGL", paste0(res_sim_location[[i]]$nome[1],"l_u_surface.html")))

r2stl(x=unique(plot_coords_l_a$x), 
               y=unique(plot_coords_l_a$y),
			   ME_lu, 
			   filename=paste0(res_sim_location[[i]]$nome[i],"l_u_surface.stl"), show.persp=FALSE)


######################################################################################################################################
######################################################################################################################################
# Final

temp_f=as.data.frame(res_sim_location[[i]])
write.csv(as.matrix(temp_f),file=gsub(" ","_",paste0(res_sim_location[[i]]$nome_loc[1],"_sim_err.csv")),row.names=FALSE)


id_min_RMSE=which.min(res_sim_location[[i]]$RMSE_sim)
temp_RMSE=res_sim_location[[i]][id_min_RMSE,]
temp_df=sim_results[[temp_RMSE$ind_data_confronto]]
names(temp_df)= c("Osservati","Stimati")

######################################################################################################################################

png(file=paste0(temp_RMSE$nome_loc,"_confronto_bestRMSE.png"),width = 1000, height = 700, units = "px", pointsize = 12)
xtsExtra::plot.xts(temp_df,
                   main = paste("Osservati e Stimati -",res_sim_location[[i]]$nome[i],"-","Parametri indici di competizione: Larve=",temp_RMSE$alfa_l,"Uova=",temp_RMSE$alfa_u,"Adulti=",temp_RMSE$alfa_a)
				   ,cex.axis = 1.2,cex.main = 2.5,
                   legend.loc = "bottomright", 
				   legend.pars = list(bty = "n",cex=2,horiz=TRUE),
				   legend.names = c("Osservati","Stimati")) 
				   Sys.sleep(1)
dev.off()


res_minimal=append(res_minimal,temp_RMSE)

#################################

id_min_ME=which.min(res_sim_location[[i]]$ME_sim)
temp_ME=res_sim_location[[i]][id_min_ME,]
temp_df=sim_results[[temp_ME$ind_data_confronto]]
names(temp_df)= c("Osservati","Stimati")
png(file=paste0(temp_RMSE$nome_loc,"_confronto_bestME.png"),width = 1000, height = 700, units = "px", pointsize = 12)

######################################################################################################################################

xtsExtra::plot.xts(temp_df,
                   main = paste("Osservati e Stimati -",res_sim_location[[i]]$nome[i],"-","Parametri indici di competizione: Larve=",temp_RMSE$alfa_l,"Uova=",temp_RMSE$alfa_u,"Adulti=",temp_RMSE$alfa_a)
				   ,cex.axis = 1.2,cex.main = 2.5,
                   legend.loc = "bottomright", 
				   legend.pars = list(bty = "n",cex=2,horiz=TRUE),
				   legend.names = c("Osservati","Stimati")) 
				   Sys.sleep(1)
dev.off()

res_minimal=append(res_minimal,temp_ME)

}






























#######################################################
# Reference

# http://www.r-bloggers.com/creating-3d-geographical-plots-in-r-using-rgl/
# http://trestletechnology.net/2013/10/interactive-3d-in-shiny-shinyrgl/
# http://life.bio.sunysb.edu/morph/morphmet/geomorph.pdf
# http://cran.r-project.org/web/packages/r2stl/r2stl.pdf
# http://www.r-bloggers.com/color-palettes-in-r/
# http://stackoverflow.com/questions/17258787/formating-of-persp3d-plot

#######################################################
# Supplementary code

# scene3d saves a copy of a scene to an R variable; writeWebGL, writePLY and writeSTL write the scene to a file in various other formats. 

# browseURL(paste(“file://”, writeWebGL(dir=file.path(tempdir(), “webGL”), width=700), sep=”"))
# writeWebGL(dir = "webGL", filename = file.path(dir, "index.html"), 
# template = system.file(file.path("WebGL", "template.html"), package = "rgl"),
# snapshot = TRUE, font = "Arial")
# ans=krigeInterp(res_sim_location[[1]]$alfa_l, res_sim_location[[1]]$alfa_a, res_sim_location[[1]]$MAE_sim, extrap = T)


# rgl.pop()
# nicer colored plot
# ylim <- range(calls)
# ylen <- ylim[2] - ylim[1] + 1
# col <- topo.colors(ylen)[ calls-ylim[1]+1 ]
# x = (1: nrow(calls))
# z = (1: ncol(calls))

# rgl.bg(sphere=FALSE, color=c("black"), lit=FALSE)
# rgl.viewpoint( theta = 300, phi = 30, fov = 170, zoom = 0.03)
# rgl.surface(x, z, calls, color = col, shininess = 10)
# rgl.bringtotop()