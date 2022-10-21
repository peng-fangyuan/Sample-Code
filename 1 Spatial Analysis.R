# Title: Electrification, Climate Adaptation, and Health in Africa
# Note: This is sample code about spatial sorting and matching for this project 

# Date: Qct 21, 2022
# Author: Fangyuan PENG

#************************************
## 1. Preliminary                 ----
#************************************

### 1.0 Load packages, these lines will install the packages if necessary 
rm(list=ls())
packagelist <- c("rgdal",  "sp", "sf","stats" ,"rlist", "texreg",
                 "lfe", "stargazer","plm","reshape","data.table","rgeos","ggplot2","dplyr","tidyr",
                 "foreign","geosphere", "ebal","WeightIt","cobalt","shiny")
newpackages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(newpackages)) install.packages(newpackages)
lapply(packagelist, function(i) require(i, character.only=TRUE))
lapply(packagelist, function(i) library(i, character.only=TRUE))
rm(packagelist, newpackages)
#install.packages("devtools")
dir <- list()
dir$root <- dirname(getwd())
dir$data <- paste0(dir$root,"/raw")


#*****************************************
# 2. Import and clean data ----
#*****************************************
## 2.1 Africa Map  ---
# africa map
map <-   readOGR(dsn = paste0(dir$data,"/afr_g2014_2013_2/afr_g2014_2013_2.shp"),layer = paste0("afr_g2014_2013_2"))
map_df<-as.data.frame(map)  #the prefecture-city level boundary

## 2.2 Village GIS ----

village_data <- read_dta(paste0(dir$data,"/village_GIS.dta"))
village_data_shp <-SpatialPointsDataFrame(cbind(village_data$longnum ,village_data$latnum),village_data,proj4string = crs(map))
# DHS Country code
DHS_ccode <- read_csv(paste0(dir$data,"/Country_code.csv"))


## 2.3 China aid----
China_aid <- read_csv(paste0(dir$data,"/AidDatasall.csv"))


# # exclude data without geolocation
China_aid <- China_aid %>% filter(!is.na( `geoJSON URL DL`) )

# # extract geolocation from website
 EUrl <-  China_aid$`geoJSON URL DL`
 Geolist <- lapply(1:length(EUrl),function(i){
   print(i)
   df <- geojson_sf(China_aid$`geoJSON URL DL`[i])
   return(df)
 })
 Geolist2 <- lapply(1:length(Geolist),function(i){
   print(i)
   df <- Geolist[[i]]%>%select(id,geometry)
   return(df)
 })
 geocol2<- bind_rows(Geolist2)
 China_aid$id <- China_aid$`AidData TUFF Project ID`
# #merge geolocation with raw table:
 China_aid2 <- geocol2 %>% merge(China_aid,by="id")
# Save data with geolocation
save(China_aid2,file = paste0(dir$root,"/chinaaid_geo.RData"))

#load(paste0(dir$root,"/chinaaid_geo.RData"))


## 2.4 Check extend ----

### 2.4.1 village Map----

map_sf <-  st_as_sf(map)
proj4string(village_data_shp) <- CRS("+init=epsg:4326") 
village_data_sf <-   st_as_sf(village_data_shp)

# Check whether location is correct through mapping
# Note: part of DHS Village Geo-location may have mistake 
mapview(village_data_shp)

# select africa country
sf::sf_use_s2(FALSE)
China_aid3 <-  
  st_filter( China_aid2,map_sf) 
m<-mapview(China_aid3)
mapshot(m, url = paste0(dir$root, "/dataset/map_africa.html"))



#*****************************************
# 3. Geo-Matching Method 1 : Buffer----
#*****************************************

# We plan to build up several buffers in different radius 
# which could help us identify the impact zone
# 10- 100
#  Simplest way
for (i in 1:10){
  aid_buffer<- st_buffer(China_aid2,i*10*10^3)
  aid_buffer <- aid_buffer %>%
    mutate(aidid=row_number())
  
  v_aid_match <- st_join(village_data_sf ,aid_buffer,left=TRUE)
  
  save(v_aid_match ,
       file =paste0(dir$root,"/dataset/v_aid_match_",as.character(i),"0km.RData"))
  
}

#*****************************************
# 4. Geo-Matching Method 2 : distance within country----
#*****************************************

# Consider different ring 
# village may fall in different projects sector
# 
# different buffer/different sectors

## 4.1 Get Center Point of each project(polygon)

nc_sp <- sf:::as_Spatial(China_aid2$geometry)
Aid_centroids <- gCentroid(nc_sp , byid = T)
Aid_centroids_p <- Aid_centroids@coords %>% as.data.frame()
aid_shp <-SpatialPointsDataFrame(cbind(Aid_centroids_p$x,Aid_centroids_p$y),Aid_centroids_p ,proj4string = crs(map))
aid_shp_sf <-  st_as_sf(aid_shp)
aid_shp_sf$id <- China_aid3$id
aid_place  <-  
  st_join(aid_shp_sf,map_sf,left=TRUE) 


## 4.2 Point to Point Matching
countrylist <- map_df %>% 
  distinct(ADM0_NAME) %>%
  arrange(ADM0_NAME) %>%
  as.list()
countrylist <- countrylist[[1]]
#length(countrylist)



dist_c<-lapply(1:length(countrylist),function(i){
  a<-subset.data.frame(aid_place ,aid_place$ADM0_NAME==countrylist[i])
  a2<-subset.data.frame(village_data_sf,village_data_sf$ADM0_NAME==countrylist[i])
  #i %in%no_m
  if (nrow(a)==0 |nrow(a2)==0){ 
    print(i)
  }else
  {
    b<-as.matrix(pointDistance(cbind(a$x,a$y),na.omit(cbind(a2$longnum,a2$latnum)),type='Euclidean',lonlat=TRUE) )
    if (ncol(b)==1){
      b<-t(b)
    }else{
      b <- b
    }
    b<-as.data.frame(b)
    colnames(b)<-na.omit(a2$dhsnum)
    rownames(b)<-a$id
    c<-melt(setDT(b,keep.rownames=TRUE),"rn")
    c$countryID=i
    #c$aidid <- c$rn
  }
  return(c)
})

dist_c1<-list.remove(dist_c ,no_in)
dist_c2<-do.call(rbind,dist_c1)
names(dist_c2)<-c("aidid","dhsnum","dist","CountryID")
dist_c2$aidid <-as.numeric( dist_c2$aidid ) 
dist_c2 <- dist_c2[order(dist_c2$aidid),]
rm(dist_c)
# merge with village_data, complete some village do not have bri projects
dist_c3<- dist_c2 %>%
  merge(village_data_2,by="dhsnum",all.y=T)
save(dist_c3,file= paste0(dir$root,"/dataset/dist_c3.RData"))

