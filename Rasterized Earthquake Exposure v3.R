## Testing ability to create earthquakes as polygons

#Define an attenuation function with optimization on a single parameter:
# Given a Magnitude, Depth, and target Modified Mercalli Intensity, 
# what is the minimum radius such that all points within the radius experienced at least the given MMI?
## note that the attenuation function is monotonically decreasing.
## exponentiating the distance in the function reduces the range that optimize has to consider, and results in a higher true convergence rate. 

pacman::p_load(
  tidyverse,
  lubridate, #easy date handling
  tigris, #grab base county shapefiles
  rgeos, #base geospatial packages
  sp, #spatial polygons (and points, etc.)
  spdplyr, #dplyr functions for sp objects. Occasionally useful. 
  raster, #buffer function, and raster manipulations where useful
  animation, #print gifs!
  RColorBrewer #color palettes for plots
)

#specify the attenuation function and a penalty shape for optimization: 
#in this case, squared MMI error should be pretty precise for distance estimation. 
attenuate_optim<-function(param, M, D, MMI){
  x=exp(param[1])
  MMI_est = 11.72+2.36*(M-6)+0.1155*(M-6)^2-0.44*log10(sqrt((x/1000)^2+D^2+289))-.002044*sqrt((x/1000)^2+D^2+289)+2.31*ifelse(sqrt((x/1000)^2+D^2+289)<80, 0, log10(sqrt((x/1000)^2+D^2+289)/80)) -.479*M*log10(sqrt((x/1000)^2+D^2+289))
  fit = (MMI-MMI_est)^2
  return(fit)
}

##Select several radii: one for each MMI level in 2:5. 

Earthquake_raw<-read.csv("USGS Earthquake Data.csv")

Earthquake_raw<-Earthquake_raw %>%
  mutate(x=longitude, 
         y=latitude,
         yearmonth = parse_date_time2(gsub("-","",substr(as.character(Earthquake_raw$time),1,7)), 
                                      orders = c("%Y%m")))

##Grab a census shapefile to plot over
Southern_Counties<-counties(state = c('NE', 'TX', 'KS', 'OK', 'NM', 'AR', 'LA', 'CO', 'MO'))
Census_Projection<-proj4string(Southern_Counties)
#Obtain the minimizing parameter value for the objective function: this is our radius
#If the objective function is outside of our radius, this means that we didn't converge - set the radius to zero, as there is no radius for which the MMI was experienced.

for(i in 2:5){
  
  Earthquake_raw[,paste0("Optim_",i)] <-apply(Earthquake_raw, 1,
                                              function(x) optimize(f=attenuate_optim, M = as.numeric(x[["mag"]]),
                                                                   D=as.numeric(x[["depth"]]), MMI=i, interval = c(0,15))$objective
  )
  
  Earthquake_raw[,paste0("Radius_",i)]<-exp(apply(Earthquake_raw, 1, 
                                                  function(x) optimize(f=attenuate_optim, M = as.numeric(x[["mag"]]),
                                                                       D=as.numeric(x[["depth"]]), MMI=i, interval = c(0,15), tol = .01)$minimum
  ))
  
  Earthquake_raw[,paste0("Radius_",i)]<-ifelse(Earthquake_raw[,paste0("Optim_",i)]>.01,0,Earthquake_raw[,paste0("Radius_",i)])
}

#Define a spatial points dataframe for the earthquake data
Earthquake_SPDF<-
  SpatialPointsDataFrame(
    coords = cbind(Earthquake_raw$longitude, Earthquake_raw$latitude),
    data=Earthquake_raw,
    proj4string = CRS(Census_Projection)
  )

#Running all buffers at same time returns a single feature - we want one polygon per earthquake, thus the list.
#For each point, the ID slot of the polygon is initially 1: change this to the polygon ID itself.
buffer_store_2<-list()
buffer_store_3<-list()
buffer_store_4<-list()
buffer_store_5<-list()

for(i in 1:length(Earthquake_SPDF)){
  buffer_store_2[[i]]<-buffer(Earthquake_SPDF[i,], width = Earthquake_SPDF$Radius_2[i])
  slot(slot(buffer_store_2[[i]], "polygons")[[1]], "ID") = paste('"',i,'"', sep = "")
  
  buffer_store_3[[i]]<-buffer(Earthquake_SPDF[i,], width = Earthquake_SPDF$Radius_3[i]) 
  slot(slot(buffer_store_3[[i]], "polygons")[[1]], "ID") = paste('"',i,'"', sep = "")
  
  buffer_store_4[[i]]<-buffer(Earthquake_SPDF[i,], width = Earthquake_SPDF$Radius_4[i])
  slot(slot(buffer_store_4[[i]], "polygons")[[1]], "ID") = paste('"',i,'"', sep = "")
  
  buffer_store_5[[i]]<-buffer(Earthquake_SPDF[i,], width = Earthquake_SPDF$Radius_5[i])
  slot(slot(buffer_store_5[[i]], "polygons")[[1]], "ID") = paste('"',i,'"', sep = "")
}

###Prepping a grid for conversion of many polygons to raster
Earthquake_Extent<-bbox(Southern_Counties)

x<-seq(from = Earthquake_Extent[1,1], to = Earthquake_Extent[1,2], by = .05)
y<-seq(from = Earthquake_Extent[2,1], to = Earthquake_Extent[2,2], by = .05)
xy<-expand.grid(x=x, y=y)

#create a spatial grid defined by these evenly spaced points
grid_points<-SpatialPointsDataFrame(coords=xy, data=xy, proj4string=CRS(Census_Projection))
grid<-SpatialGrid(points2grid(grid_points))
proj4string(grid) <- CRS(Census_Projection)


#store the buffer results and select a MMI to move forward with. 

month_list<-seq(from=as.Date(min(Earthquake_raw$yearmonth))+months(1), 
                to=as.Date(max(Earthquake_raw$yearmonth)), 
                by="1 month")

##Create a function that takes the above estimates and transforms them into gif plots over time
#type is either incremental or cumulative
#buffer_set is one of the buffer stores from above
#output_name is the name of the output .gif file
#month_list is a list of months to evaluate at. A default of "the entire time period" is available above. 
buffer_create<-function(buffer_set, output_name, month_list, type){
  
  buffer_store<-buffer_set
  #consolidate all buffers into a single set.
  Earthquake_Buffers <- do.call("rbind",buffer_store) 
  raster_list<-list()
  
  ##For each month end where I have earthquake data, create a raster for a given MMI and store in a list with a date attribute
  #switching the <= for a == below gives the incremental change rather than the cumulative change
  for (i in seq_along(month_list)){
    if(type=="Cumulative"){
      Earthquake_Buffers_Date<-Earthquake_Buffers[which(Earthquake_raw$yearmonth <= as.Date(month_list[i]))]
    }else{
      Earthquake_Buffers_Date<-Earthquake_Buffers[which(Earthquake_raw$yearmonth == as.Date(month_list[i]))]
    }
    if(length(Earthquake_Buffers_Date)==0){
      next
    }
    #use over to list the polygons overlapping the grid. The length of the list of polygons is the number of earthquakes exposed!
    o <- over(grid, Earthquake_Buffers_Date, returnList=TRUE)
    ct <- sapply(o, length)
    raster_list[[which(month_list== month_list[i])]]<-SpatialGridDataFrame(grid, data=data.frame(ct=ct, date=month_list[i]))
    #optional progress tracking: this does take a while
    print(paste(round(i/max(seq_along(month_list)),2)*100, "Percent Done"))
  }
  
  # set breaks for the new plot: use 10 points, evenly spaced, and 0 as a special case (nulled out in color)
  my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
  
  ##obtain some max exposure for plotting.
  if(type=="Cumulative"){
    #obtain max for cumulative gif
    maxexp<-my.max(raster_list[[max(seq_along(month_list))]]@data$ct)
  }else{
    ##alternate method of determining max for incremental gif:
    for(i in seq_along(month_list)){
      if(i==1){maxexp<-0}
      maxexp<-max(my.max(raster_list[[i]]@data$ct), maxexp)
    }
  }
  #Set plot breakpoints. Using pseudo-log scale: discrete breaks for 0 and 1, and log over remaining interval. 
  #plot_breaks<-sort(unique(c(0,1,floor(maxexp / c(1, 1.5,1.75, 2, 2.75, 3:20)))))
  plot_breaks<-sort(unique(c(0,1,exp(seq(from=0, to=max(log(maxexp)%/%1,0)+1, by = 1)))))
  plot_breaks
  # set color scale such that exactly 0 exposure is transparent
  cols<-c(heat.colors(n=1, alpha=0), rev(heat.colors(n=length(plot_breaks)-2, alpha=1)))
  
  #storing the plot as a function to be called, so I can loop over multipart functions. 
  store_MMI<-function(grid){
    plot(grid, breaks = plot_breaks, col=cols, main = unique(grid@data$date))
    plot(Southern_Counties, add=TRUE)
  }
  
  plot_MMI<-function(x){store_MMI(raster_list[[x]])}
  
  saveGIF(
    for(i in seq_along(month_list)){
      plot_MMI(i)
    },
    movie.name = output_name, interval=.3
  )
}

buffer_create(buffer_store_2, "MMI2_Incremental.gif", month_list = month_list, type = "Incremental")
buffer_create(buffer_store_3, "MMI3_Incremental.gif", month_list = month_list, type = "Incremental")
buffer_create(buffer_store_4, "MMI4_Incremental.gif", month_list = month_list, type = "Incremental")
buffer_create(buffer_store_5, "MMI5_Incremental.gif", month_list = month_list, type = "Incremental")

buffer_create(buffer_store_2, "MMI2_Cumulative.gif", month_list = month_list, type = "Cumulative")
buffer_create(buffer_store_3, "MMI3_Cumulative.gif", month_list = month_list, type = "Cumulative")
buffer_create(buffer_store_4, "MMI4_Cumulative.gif", month_list = month_list, type = "Cumulative")
buffer_create(buffer_store_5, "MMI5_Cumulative.gif", month_list = month_list, type = "Cumulative")

