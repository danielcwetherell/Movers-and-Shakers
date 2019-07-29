#memory.limit(size=4000000)
ifelse(!require(pacman), install.packages("pacman"),"")
pacman::p_load(
  glmnet, #elastic net regression
  rethinking, # A wrapper for STAN (hamiltonian MC engine)
  animation, #Writing fun counterfactual plot GIFs
  gam #for GAMs. If that wasn't obvious
)

Main_Data.pois<-Main_Data %>%
  mutate(Time2 = Time^2,
         Time3 = Time^3,
         Time4 = Time^4,
         Time5 = Time^5,
         HOUST.log = log(HOUST),
         GS30.log = log(GS30),
         CSUSHPINSA.log = log(CSUSHPINSA))%>%
  filter(complete.cases(Listings, huest, EQ_L36_6))


Var_list<-Main_Data.pois %>% 
  dplyr::select(starts_with("EQ"), HOUST.log, GS30.log, CSUSHPINSA.log, FIPS, starts_with("Time")) %>%
  colnames()
formula.pois<-as.formula(paste("Listings ~ ", paste(Var_list,collapse="+"), " + offset(log(huest))"))

model.pois<-glm(formula.pois, family = "poisson", data=Main_Data.pois)
summary(model.pois)
#There's no guarantee that earthquakes matter, so fit a model without the EQ variables and see if it worsens fit
Var_list.noeq<-Main_Data.pois %>% 
  dplyr::select(HOUST.log, GS30.log, CSUSHPINSA.log, FIPS, starts_with("Time")) %>%
  colnames()
formula.pois.noeq<-as.formula(paste("Listings ~ ", paste(Var_list.noeq,collapse="+"), " + offset(log(huest))"))
model.pois.noeq<-glm(formula.pois.noeq, family = "poisson", data=Main_Data.pois)

#Run a chi-squared test to check whether information from earthquakes improves model by more than noise/chance
anova(model.pois, model.pois.noeq, test = "Chisq")
#Neat, looks okay. 

#as a counterfactual exercise: remove seismicity from the region and estimate how the model responses change
#estimate with EQ - estimate without EQ = net change in sales rate from seismicity
##Yes, there's still some natural seismicity in the region, so this isn't a perfect exercise. It should be adequate, though. 
Main_Data.pois$Listings.Predicted<-model.pois$fitted.values
Main_Data.pois$List.Rate <- Main_Data.pois$Listings.Predicted/Main_Data.pois$huest

Main_Data.pois_EQ<- Main_Data.pois %>% 
  data.frame() %>% 
  mutate_at(vars(contains("EQ")), list(~ ifelse(is.na(.), ., 0)))

Main_Data.pois$withoutEQ<-predict(object=model.pois, newdata= Main_Data.pois_EQ, type="response")

Main_Data.pois$List.Rate.Change <- (Main_Data.pois$Listings.Predicted - Main_Data.pois$withoutEQ)/Main_Data.pois$huest

#normalize change in the sales rate to the response standard deviation within the county over the modeled time periods. 
Main_Data.pois<- Main_Data.pois %>% 
  group_by(FIPS) %>% 
  mutate(List.Rate.SD = sd(Listings/huest),
         List.Rate.Mean = mean(Listings/huest),
         List.Rate.Change.S=(Listings.Predicted - withoutEQ)/huest/List.Rate.SD,
         List.Rate.Change.M= (Listings.Predicted - withoutEQ)/huest/List.Rate.Mean)


##create a gif showing progression over time, in terms of SDs of the population sales rate. 
##Runtime isn't super because of the internal merge of county polygons, with the requirement
saveGIF(
  for(Year in unique(Main_Data.pois$YYYYMM)){
    temp<-merge(counties, Main_Data.pois %>% filter(YYYYMM==Year) %>% dplyr::select(FIPS, YYYYMM, List.Rate.Change.S), by.x = "GEOID", by.y = "FIPS")
    print(
      spplot(temp, zcol = "List.Rate.Change.S", col.regions = c("red", "orange", "white", "white", "blue", "purple"),
             main=list(label=paste("List Rate #SD Change from Earthquakes", Year),cex=1.5), at = c(-2,-1,-.25,0,.25,1,2))
    )
  },
  movie.name = "poisSD_Predictions.gif", interval=.3
)

saveGIF(
  for(Year in unique(Main_Data.pois$YYYYMM)){
    temp<-merge(counties, Main_Data.pois %>% filter(YYYYMM==Year) %>% dplyr::select(FIPS, YYYYMM, List.Rate.Change.M), by.x = "GEOID", by.y = "FIPS")
    print(
      spplot(temp, zcol = "List.Rate.Change.M", col.regions = c("red", "orange", "white", "white", "blue", "purple"),
             main=list(label=paste("List Rate % of FIPS Mean Change from Earthquakes", Year),cex=1), at = c(-.25,-.1,-.05,0,.05,.1,.25))
    )
  },
  movie.name = "poisMU_Predictions.gif", interval=.3
)


##Example regression: Listings as a percentage of total housing units
##Fitting polynomial time trend to capture general market, and measuring perturbations from earthquakes and some economic indicators.
##capturing logged economic data as well (for % change treatment)
##note that the complete cases for the EQ lag determines how many older years are kept, with a select_if removing variables that have NA at later lags
Main_Data.glmnet<-Main_Data %>%
  mutate(Time2 = Time^2,
         Time3 = Time^3,
         Time4 = Time^4,
         Time5 = Time^5,
         HOUST.log = log(HOUST),
         GS30.log = log(GS30),
         CSUSHPINSA.log = log(CSUSHPINSA))%>%
  filter(complete.cases(Listings, huest, EQ_L21_6))%>%
  select_if(~ !any(is.na(.)))
#extracting all earthquake and economic variables, time, and new powers of time
#using model.matrix to one hot encode FIPS codes. 
FIPS_Matrix<- model.matrix( ~ Main_Data.glmnet$FIPS-1, data=factor(Main_Data$FIPS))

Main_Data.glmnet_x<-as.matrix(cbind(Main_Data.glmnet 
                                    %>% dplyr::select(starts_with("EQ"),starts_with("Time"), HOUST.log, GS30.log, CSUSHPINSA.log)
                                    ,FIPS_Matrix)
                              )
##Running an elastic net regression to pare down on earthquake exposure metrics. CV method uses (default) 10 fold cross-validation to select lambda.
##Note that there's some hyperparameter optimization to be done around alpha (and glmnet doesn't by default search over alphas) 
##One thing on my mind is that I can't think of a strong argument to penalize the FIPS parameters (which I'm just thinking of as fixed effects)
###A better approach might have been to train the model on residuals, or to write a custom function that only penalizes certain betas. 
Model.glmnet<-cv.glmnet(y=as.vector(Main_Data.glmnet$Listings),x=Main_Data.glmnet_x, offset = as.vector(log(Main_Data.glmnet$huest)), 
       family = "poisson", alpha = .9)

##Extract predicted new listings on the input data, using the lambda that minimized the objective function of cv.glmnet 
Main_Data.glmnet$Listings.Predicted<-predict(object=Model.glmnet, newx = Main_Data.glmnet_x, 
                            newoffset = as.vector(log(Main_Data.glmnet$huest)),
                            type = "response", "lambda.min")[,1]

#as a counterfactual exercise: remove seismicity from the region and estimate how the model responses change
#estimate with EQ - estimate without EQ = net change in sales rate from seismicity
##Yes, there's still some natural seismicity in the region, so this isn't a perfect exercise. It should be adequate, though. 

Main_Data.glmnet_EQ<- Main_Data.glmnet_x %>% 
  data.frame() %>% 
  mutate_at(vars(contains("EQ")), list(~ ifelse(is.na(.), ., 0))) %>%
  as.matrix(.)

Main_Data.glmnet$List.Rate <- Main_Data.glmnet$Listings.Predicted/Main_Data.glmnet$huest

Main_Data.glmnet$withoutEQ<-predict(object=Model.glmnet, newx = Main_Data.glmnet_EQ, 
              newoffset = as.vector(log(Main_Data.glmnet$huest)),
              type = "response", "lambda.min")[,1]

Main_Data.glmnet$withEQ<-predict(object=Model.glmnet, newx = Main_Data.glmnet_x, 
              newoffset = as.vector(log(Main_Data.glmnet$huest)),
              type = "response", "lambda.min")[,1]

Main_Data.glmnet$List.Rate.Change <- (withEQ - withoutEQ)/Main_Data.glmnet$huest
#normalize change in the sales rate to the response standard deviation within the county over the modeled time periods. 
Main_Data.glmnet<- Main_Data.glmnet %>% 
  group_by(FIPS) %>% 
  mutate(List.Rate.SD = sd(Listings/huest),
         List.Rate.Mean = mean(Listings/huest),
         List.Rate.Change.S=(Listings.Predicted - withoutEQ)/huest/List.Rate.SD,
         List.Rate.Change.M= (Listings.Predicted - withoutEQ)/huest/List.Rate.Mean)

##create a gif showing progression over time, in terms of SDs of the population sales rate. 
##Runtime isn't super because of the internal merge of county polygons, with the requirement
saveGIF(
  for(Year in unique(Main_Data.glmnet$YYYYMM)){
    temp<-merge(counties, Main_Data.glmnet %>% filter(YYYYMM==Year) %>% dplyr::select(FIPS, YYYYMM, List.Rate.Change.S), by.x = "GEOID", by.y = "FIPS")
    print(
        spplot(temp, zcol = "List.Rate.Change.S", col.regions = c("red", "orange", "white", "white", "blue", "purple"),
           main=list(label=paste("List Rate #SD Change from Earthquakes", Year),cex=1.5), at = c(-2,-1,-.5,0,.5,1,2))
    )
  },
  movie.name = "GlmnetSD_Predictions.gif", interval=.3
)

saveGIF(
  for(Year in unique(Main_Data.glmnet$YYYYMM)){
    temp<-merge(counties, Main_Data.glmnet %>% filter(YYYYMM==Year) %>% dplyr::select(FIPS, YYYYMM, List.Rate.Change.M), by.x = "GEOID", by.y = "FIPS")
    print(
      spplot(temp, zcol = "List.Rate.Change.M", col.regions = c("red", "orange", "white", "white", "blue", "purple"),
             main=list(label=paste("List Rate % of FIPS Mean Change from Earthquakes", Year),cex=1), at = c(-.25,-.1,-.05,0,.05,.1,.25))
    )
  },
  movie.name = "GlmnetMU_Predictions.gif", interval=.3
)

##Point of interest: Can we detect spatial autocorrelation? 
##Counterfactual above tends to show evidence that listings increase in regions adjacent to declines. 
##I don't have a strong prior on this.
weight<-distm(cbind(Main_Data.glmnet$x, Main_Data.glmnet$y),
              cbind(Main_Data.glmnet$x, Main_Data.glmnet$y),
              fun=distHaversine)

##Calculate Moran's I, a measure of spatial autocorrelation, based on the residual list rates
##A spatiotemporal distance matrix would be more appropriate, but that relies on some selection of a relative time vs. distance weighting.
ape::Moran.I(as.vector((Main_Data.glmnet$Listings.Predicted-Main_Data.glmnet$Listings)/Main_Data.glmnet$huest), weight=weight, scaled = TRUE, na.rm = FALSE,
        alternative = "two.sided")
##Here, none observed. 
##Other approaches to Moran's I use adjacency based matrices. Again, though, you need some selection on time vs. distance weighting.
##One approach would be to select a weight set that minimizes some objective/penalty function. Might be worth looking into

##At worst with the elastic net approach, I have a sense of which earthquake variables to keep in other modeling.
##Extract the earthquake variables that weren't removed by the lasso
#Side note: the predict.cv.glmnet function can return this list as well. Realized that after I wrote the below section.

##Going to try a non-spatial panel regression method here.
useful_vars<-rownames(coef(Model.glmnet, "lambda.min"))[which(coef(Model.glmnet, "lambda.min")!=0 & substr(rownames(coef(Model.glmnet, "lambda.min")),1,2)=="EQ")][-1]
additional_vars<-paste("HOUST.log", "GS30.log", "CSUSHPINSA.log", "log(huest)", sep = "+")
formula<-as.formula(paste("Listings ~", paste(useful_vars, collapse = "+"),"+",additional_vars, sep = " "))

##Unfortunately, having tried a few times to get the offset to work, it doesn't seem to work with this package. 
##The parameter estimate for log(huest) is close to 1, though. 
model.pglm<-pglm(formula=formula, data=Main_Data.glmnet, family = "poisson", index=c("FIPS","Time"), model = "pooling", effect="twoways")
summary(model.pglm)


##Trying out a GAM approach, breaking from the GLM approaches above. 
##One major worry out of the GLM approach is a nonlinearity in how we think about earthquake exposure
##Smoothed approaches might yield more accurate impressions of behavior

Main_Data.gam<-Main_Data %>%
  mutate(HOUST.log = log(HOUST),
         GS30.log = log(GS30),
         CSUSHPINSA.log = log(CSUSHPINSA))%>%
  filter(complete.cases(Listings, huest, EQ_L36_6))

#Extract variable list for splines
Var_list<-Main_Data.pois %>% 
  ungroup() %>% 
  dplyr::select(starts_with("EQ"), HOUST.log, GS30.log, CSUSHPINSA.log) %>%
  colnames() 
##Wrap each of the above contiunuous parameters in s() to estimate splines. 
Var_list<-paste("s(", Var_list,")", sep = "")

formula.gam<-as.formula(paste("Listings ~ ", paste(Var_list,collapse="+"), "+ s(Time,df=6) + FIPS + offset(log(huest))"))

model.gam<-gam(formula.gam, data=Main_Data.gam, family = poisson(link=log), sp = )

summary(model.gam)

#as a counterfactual exercise: remove seismicity from the region and estimate how the model responses change
#estimate with EQ - estimate without EQ = net change in sales rate from seismicity
##Yes, there's still some natural seismicity in the region, so this isn't a perfect exercise. It should be adequate, though. 
Main_Data.gam$Listings.Predicted<-model.gam$fitted.values
Main_Data.gam$List.Rate <- Main_Data.gam$Listings.Predicted/Main_Data.gam$huest

Main_Data.gam_EQ<- Main_Data.gam %>% 
  data.frame() %>% 
  mutate_at(vars(contains("EQ")), list(~ ifelse(is.na(.), ., 0)))

Main_Data.gam$withoutEQ<-predict(object=model.gam, newdata= Main_Data.gam_EQ, type="response")

Main_Data.gam$List.Rate.Change <- (Main_Data.gam$Listings.Predicted - Main_Data.gam$withoutEQ)/Main_Data.gam$huest

#normalize change in the sales rate to the response standard deviation within the county over the modeled time periods. 
Main_Data.gam<- Main_Data.gam %>% 
  group_by(FIPS) %>% 
  mutate(List.Rate.SD = sd(Listings/huest),
         List.Rate.Mean = mean(Listings/huest),
         List.Rate.Change.S=(Listings.Predicted - withoutEQ)/huest/List.Rate.SD,
         List.Rate.Change.M= (Listings.Predicted - withoutEQ)/huest/List.Rate.Mean)


##create a gif showing progression over time, in terms of SDs of the population sales rate. 
##Runtime isn't super because of the internal merge of county polygons, with the requirement
saveGIF(
  for(Year in unique(Main_Data.gam$YYYYMM)){
    temp<-merge(counties, Main_Data.gam %>% filter(YYYYMM==Year) %>% dplyr::select(FIPS, YYYYMM, List.Rate.Change.S), by.x = "GEOID", by.y = "FIPS")
    print(
      spplot(temp, zcol = "List.Rate.Change.S", col.regions = c("red", "orange", "white", "white", "blue", "purple"),
             main=list(label=paste("List Rate #SD Change from Earthquakes", Year),cex=1.5), at = c(-2,-1,-.25,0,.25,1,2))
    )
  },
  movie.name = "GamSD_Predictions.gif", interval=.3
)

saveGIF(
  for(Year in unique(Main_Data.gam$YYYYMM)){
    temp<-merge(counties, Main_Data.gam %>% filter(YYYYMM==Year) %>% dplyr::select(FIPS, YYYYMM, List.Rate.Change.M), by.x = "GEOID", by.y = "FIPS")
    print(
      spplot(temp, zcol = "List.Rate.Change.M", col.regions = c("red", "orange", "white", "white", "blue", "purple"),
             main=list(label=paste("List Rate % of FIPS Mean Change from Earthquakes", Year),cex=1), at = c(-.25,-.1,-.05,0,.05,.1,.25))
    )
  },
  movie.name = "GamMU_Predictions.gif", interval=.3
)
