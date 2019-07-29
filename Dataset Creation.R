ifelse(!require(pacman), install.packages("pacman"),"")
pacman::p_load(
  spdep, #spatial weighting matrices, etc. 
  tidycensus, #census data api calls require a key. Signup is really easy
  data.table, #Fantastic data reshaping capability
  RcppRoll, #used for rolling sums for cumulative lags. 
  quantmod, #useful for extracting economic data. Here, used for a few FRED datasets.
  geosphere, #good for distance calculations on latlong data. 
  dplyr, #duh.
  lubridate #easy date handling, though I go back and forth between formats
)

library(rgeos)#spatial manipulation. Seems to not like being called by p_load. 

adjacent_states<-c('TX', 'KS', 'OK','AR')


#census_api_key("getyourowncensuskeyitseasy", install=TRUE)

##Economic data: Case-Schiller house price index, 30 year treasury rate, and Housing starts. 
getSymbols(c("CSUSHPINSA", "GS30","HOUST"),src='FRED', auto.assign = T)
Economic_Data<-cbind(CSUSHPINSA, GS30,HOUST) %>%
  data.frame(.)%>%
  dplyr::mutate(YYYYMM = rownames(.),
                YYYYMM = parse_date_time2(YYYYMM, orders = "%Y-%m-%d"),
                YYYYMM.l = year(YYYYMM%m-%months(1))*100+month(YYYYMM%m-%months(1)))%>%
  dplyr::select(-YYYYMM)%>%
  filter(complete.cases(.))
###Pull census tract level data. For accurate measurement of earthquake exposure. 
###Plan is to calculate population weighted earthquake exposure across tract centroids, aggregated up to county level. 
#This should yield more accurate results than a standard "just take the county centroid" method". 
##Extract population data from tracts for the most recent 5-year ACS (with the given pop variable available). Remove tracts with no people in them.
#Returns a SpatialPolygonsDataFrame. 
tracts<-get_acs(geography = "tract", variables = "B01003_001", 
                state = adjacent_states, geometry = T)%>%
  filter(estimate>0)%>%
  as(., 'Spatial')

tract_centroids<-gCentroid(tracts, byid=TRUE, id=tracts$GEOID)
tract_coords<-data.frame(x=tract_centroids@coords[,1], y=tract_centroids@coords[,2])

counties<-get_acs(geography= 'county', variables = "B19013_001", state = adjacent_states, geometry = T) %>%
  as(., 'Spatial')

county_centroids<-gCentroid(counties, byid=TRUE, id=counties$GEOID)
county_coords<-data.frame(x=county_centroids@coords[,1], y=county_centroids@coords[,2], FIPS = rownames(county_centroids@coords))


Earthquake_raw<-read.csv("USGS Earthquake Data.csv") %>%
  mutate(yearmonth = parse_date_time2(gsub("-","",substr(as.character(time),1,7)), 
                                      orders = c("%Y%m")), 
         date = year(yearmonth)*100+month(yearmonth))

#County level seasonally adjusted data from https://www.zillow.com/research/data/
Home_Sales_SA_raw<-read.csv('Sale_Counts_Seas_Adj_County_Zillow.csv')
#remove periods from column names so dates can parse 
names(Home_Sales_SA_raw)<-gsub(pattern="[:.:]", "", names(Home_Sales_SA_raw))

##Zillow regions differ from census ones: crosswalk maps between them.
##downloaded from http://files.zillowstatic.com/research/public/CountyCrossWalk_Zillow.csv
Zillow_Crosswalk <- read.csv("CountyCrossWalk_Zillow.csv") %>% dplyr::select(FIPS, CountyRegionID_Zillow)

Home_Sales_SA<-Home_Sales_SA_raw %>%
  left_join(., Zillow_Crosswalk, by=c("RegionID"="CountyRegionID_Zillow")) %>%
  dplyr::select(-c(RegionName, StateName, SizeRank))%>%
  melt(., id=c("RegionID", "FIPS"))%>%
  mutate(YYYYMM = as.numeric(gsub("X", "",variable)),
         YYYY = as.numeric(substr(YYYYMM,1,4)),
         MM = as.numeric(substr(YYYYMM,5,6)),
         FIPS = sprintf("%05d", FIPS)) %>% 
  rename(Sales=value) %>%
  dplyr::select(-variable, -RegionID)

#County level seasonally adjusted for sale inventory, again from zillow
Home_Listings_SA_raw<-read.csv('MonthlyListings_SSA_AllHomes_County.csv')
names(Home_Listings_SA_raw)<-gsub(pattern="[:.:]", "", names(Home_Listings_SA_raw))

Home_Listings_SA<-Home_Listings_SA_raw %>%
  left_join(., Zillow_Crosswalk, by=c("RegionID"="CountyRegionID_Zillow")) %>%
  dplyr::select(-c(RegionName, StateName, SizeRank, RegionType))%>%
  melt(., id=c("RegionID", "FIPS"))%>%
  mutate(YYYYMM = as.numeric(gsub("X", "",variable)),
         YYYY = as.numeric(substr(YYYYMM,1,4)),
         MM = as.numeric(substr(YYYYMM,5,6)),
         FIPS = sprintf("%05d", FIPS)) %>% 
  rename(Listings=value) %>%
  dplyr::select(-variable, -RegionID)


#census estimates of total estimates. Valuation is midpoint of the year. Files are split for 2000-2010, 2011-2018 (with 2010 census in middle)
Housing_Units_2010<-read.csv("hu-est00int-tot.csv")%>% filter(SUMLEV==50) %>% mutate(GEO.id2 = paste(sprintf("%02d", STATE),sprintf("%03d", COUNTY), sep=""))
Housing_Units_2018<-read.csv("PEP_2018_PEPANNHU_with_ann.csv") %>% mutate(GEO.id2 = sprintf("%05d", GEO.id2))

#generic naming for Housing_Units: remove a 7 from the 2018 data, remove an _ from the 2010 data
names(Housing_Units_2018)<-tolower(gsub(pattern="t7", "t", names(Housing_Units_2018)))
names(Housing_Units_2010)<-tolower(gsub(pattern="_", "", names(Housing_Units_2010)))                              

#merge by geoid.2, which is state+county key. Take 2010 est from only the 2010 data. 
#gather to a long dataset, removing unnecessary variables. Note that there's some NA entries. 
Housing_Units<-Housing_Units_2018 %>% dplyr::select(-huest2010)%>%
  left_join(., Housing_Units_2010, by = "geo.id2") %>%
  dplyr::select(-c(geo.id, geo.display.label, hucen42010, hubase42010, x, sumlev, state, county, ctyname, huestbase2000, 
                   hucensus2010))%>%
  reshape(., direction='long', 
          varying =paste("huest", 2000:2018, sep=""),
          timevar = 'Year',
          idvar='geo.id2',
          sep='')

#pairwise distance matrix over properties and earthquakes. 
mat<-distm(tract_coords, 
           Earthquake_raw[,c('longitude', 'latitude')], 
           fun=distHaversine)

#creating vectors of values as prep for apply function
Dmat<-matrix(data=Earthquake_raw$depth, nrow=1, ncol=nrow(Earthquake_raw))
Mmat<-matrix(data=Earthquake_raw$mag, nrow=1, ncol=nrow(Earthquake_raw))

#For lagged models, need to convert to number of months rather than yyyymm
EQDatematL<-matrix(data=(as.numeric(substr(Earthquake_raw$date,1,4))*12+as.numeric(substr(Earthquake_raw$date,5,6))), nrow=1, ncol=nrow(Earthquake_raw))

#for attenuation, apply the attenuation function over the matrix. use by column, as earthquake attributes change by quake. 
start<-Sys.time()
attenuate<-function(x){
  11.72+2.36*(M-6)+0.1155*(M-6)^2-0.44*log10(sqrt((x/1000)^2+D^2+289))-.002044*sqrt((x/1000)^2+D^2+289)+2.31*ifelse(sqrt((x/1000)^2+D^2+289)<80, 0, log10(sqrt((x/1000)^2+D^2+289)/80)) -.479*M*log10(sqrt((x/1000)^2+D^2+289))
}
for (j in 1:ncol(mat)){
  D<- Dmat[,j]
  M<- Mmat[,j]
  mat[,j]<-sapply(mat[,j],attenuate)
}
time.taken<-Sys.time()-start
time.taken


##For each centroid for each date that I have data for, calculate the number of earthquakes above a threshold as of X months. 
##Can iterate over this for each of the 5 MMI levels considered, and then convert to a long dataset by tract/month. 
##From tract/month dataset, merge back on the initial population dataset. 
##For each earthquake metric, county, and month, take a weighted average of the exposure. Note that county is first 5 digits of geoid (2 state, 3 county)
##With properly applied lagged sums, this becomes our modeling dataset. 
##For simplicity of sums, iterating over a filtered version of the attenuation matrix, and only tracking incremental gains.

#initialize final storage dataframe
County_Exposure<-data.frame(county=factor(), Month_EQ=numeric(), MMI=numeric(), EQ=numeric())

for (Month in (2004*12+11):(2019*12+5)){
  Month_EQ<-rep(Month, length(tract_centroids))
  mat_month<-mat[,which(EQDatematL==Month)]
  for (MMI in 2:6){
    assign(paste0("tract_EQ",MMI, sep=""),
           apply(mat_month, 1, function(x) sum(x>=MMI)))
  }
  County_Temp<-data.frame(county=substr(tracts@data$GEOID,1,5),
                          population = tracts@data$estimate) %>%
    cbind(tract_EQ2, tract_EQ3, tract_EQ4, tract_EQ5, tract_EQ6, Month_EQ) %>%
    group_by(county, Month_EQ) %>%
    summarise(EQ2 = sum(tract_EQ2*population)/sum(population),
              EQ3 = sum(tract_EQ3*population)/sum(population),
              EQ4 = sum(tract_EQ4*population)/sum(population),
              EQ5 = sum(tract_EQ5*population)/sum(population),
              EQ6 = sum(tract_EQ6*population)/sum(population)
    ) %>%
    data.frame() %>% 
    reshape(., direction='long', 
            varying =c("EQ2", "EQ3","EQ4", "EQ5", "EQ6" ),
            timevar = 'MMI',
            idvar='county',
            sep='')
  
  County_Exposure<-rbind(County_Exposure,County_Temp)
}
remove(County_Temp)
remove(Dmat)
remove(mat)
##With this lagged dataset, can calculate lagged responses: 
##End goal is a dataset of lagged exposure variables by county and month: 
##For the final dcast, data is long on County+Month_EQ, wide on lags and MMI
## at the end, convert my stupid date convention back to something reasonable.
##Note that each lag is an accumulation of past n months, excluding the current month (so roll_sum of 3 - current = sum of 2 lags)

#Temporary storage for testing. 
CE_Temp<-County_Exposure
County_Exposure<-CE_Temp

County_Exposure<-County_Exposure %>%
  group_by(county, MMI) %>%
  arrange(county, MMI, desc(Month_EQ))%>%
  mutate(EQ_L1 = roll_sum(EQ, 2, fill=NA)-EQ,
         EQ_L2 = roll_sum(EQ, 3, fill=NA)-EQ,
         EQ_L3 = roll_sum(EQ, 4, fill=NA)-EQ, 
         EQ_L6 = roll_sum(EQ, 7, fill=NA)-EQ,
         EQ_L9 = roll_sum(EQ, 10, fill=NA)-EQ,
         EQ_L12 = roll_sum(EQ, 13, fill=NA)-EQ,
         EQ_L15 = roll_sum(EQ, 16, fill=NA)-EQ,
         EQ_L18 = roll_sum(EQ, 19, fill=NA)-EQ,
         EQ_L21 = roll_sum(EQ, 22, fill=NA)-EQ,
         EQ_L24 = roll_sum(EQ, 25, fill=NA)-EQ,
         EQ_L36 = roll_sum(EQ, 37, fill=NA)-EQ)%>%
  dplyr::select(-EQ)%>%
  arrange(county, MMI, Month_EQ) %>%
  melt(., id.vars = c("county", "Month_EQ", "MMI")) %>%
  dcast(., county+Month_EQ~variable+MMI) %>%
  mutate(temp = ifelse(Month_EQ%%12==0,12,Month_EQ%%12),
         Month_EQ = as.numeric(paste(floor(Month_EQ/12), sprintf("%02d", temp), sep=""))
  ) %>%
  dplyr::select(-temp)

##Crafting the "final" dataset: 
##Sales from Zillow are taken at year and month combinations. Using this as a base table. 
##Listings are also from Zillow, same granularity. Not the base, as sparse timespan covered.
##Housing units from census are taken at year points
##Earthquake Exposures are taken at year and month combinations.
##Economic data are taken at year and month combinations, but merged on as lags. 
##FIPS codes are used to remove states not in the study area.
data(fips_codes)
#Time defined as months from first month of sales data. 
start<-min(Home_Sales_SA$YYYYMM)

Main_Data <- Home_Sales_SA %>%
  left_join(., Home_Listings_SA, by = c("YYYY", "YYYYMM", "MM", "FIPS")) %>% 
  full_join(., Housing_Units, by = c("YYYY"= "Year", "FIPS" = "geo.id2")) %>%
  left_join(., County_Exposure, by = c("YYYYMM" = "Month_EQ", "FIPS" = "county")) %>%
  inner_join(.,fips_codes %>% 
               mutate(FIPS = paste(state_code, county_code, sep=""))%>%
               filter(state %in% adjacent_states), by = c("FIPS" = "FIPS")) %>% 
  left_join(., Economic_Data, by = c("YYYYMM"="YYYYMM.l"))%>%
  inner_join(., county_coords, by="FIPS")%>%
  mutate(Time = (floor(YYYYMM/100)*12+(YYYYMM-floor(YYYYMM/100)*100))-
           (floor(start/100)*12+(start-floor(start/100)*100)),
         obstv = as.yearmon(paste(YYYY, MM, sep="-"))
  ) #%>%
#dplyr::filter(complete.cases(.)) %>%

#create a spatialpolygons dataframe, because it could be useful, who knows. 
Main_Data_SPDF<- SpatialPointsDataFrame(cbind(Main_Data$x, Main_Data$y), data= Main_Data)
