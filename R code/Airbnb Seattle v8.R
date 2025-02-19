rm(list = ls())
###############################
#Script Structure:
#1. Packages
# Data import
# Data preparation
# Preliminary data visual analysis
# Baseline Regressions
# Regressions with Smoothing
# Cross-validation
# Diagnostic plot to check for quantile effects
# Price surface prediction

##############################
######### Packages ###########
##############################

library(car)
library(dplyr)
library(raster)
library(sp)
library(terra)
library(mgcv)
library(smoothLUR)
library(corrplot)
library(ggplot2)
library(quantreg)
library(lmtest)
library(sandwich)
library(tripack)
library(akima)
library(arulesViz)
library(gridExtra)
library(piecewiseSEM)
library(caret)
library(yarrr)

###########################################
########        Data import          ######
###########################################

setwd('/Users/andreinevskii/Documents/Data analysis/R/Airbnb Seattle')
dat_raw <- read.csv('seattle_01.csv')
dat <- as.data.frame(dat_raw, row.names = dat_raw$room_id, col.names = colnames(dat_raw))[1:15]

col_names <- colnames(dat)
nums <- unlist(lapply(dat, is.numeric), use.names = FALSE)  
col_names_num <- col_names[nums][c(4:11)]
col_names_categ <- c('bedrooms', 'bathrooms', 'room_type', 'no_reviws')

##########################################
#######     Data preparation    ##########
##########################################

#1. NAs in Bathrooms (only 2 vslues missing) - replace with median
dat$bathrooms[is.na(dat$bathrooms)] <- median(dat$bathrooms, na.rm = T)

# 2. NAs in bathrooms and overall_satisfaction - drop observations with NAs and #reviews > 0
#replace NAs in satisfaction in other cases with 0, create a dummy 'no reviews'
ind_rating <- is.na(dat$overall_satisfaction)
dat_mod <- dat
dat_mod$no_reviws <- as.numeric(is.na(dat_mod$overall_satisfaction) | dat_mod$reviews == 0)
dat_mod$overall_satisfaction[ind_rating & dat_mod$reviews == 0] <- 0
dat_clean <- na.exclude(dat_mod)

# 3.Multicollinearity of bedrooms and accommodates:omit accommodates
# (in linear and quantile regression its coefficient is not significant)
# cor_matr <- cor(dat_clean[col_names_num])
# #Correlations
# corrplot(cor_matr, type = "upper", order = "hclust", 
#          tl.col = "black", tl.srt = 45)
# ggplot(dat, aes(x=bedrooms, y=accommodates)) + geom_point() + geom_smooth(method=loess)

# 4. Splite 'room_type' into binary variables "Private room", "Entire home/apt","Shared room" 
# dat_clean$private_room <- as.numeric(dat_clean$room_type == "Private room")
# dat_clean$entire_home <- as.numeric(dat_clean$room_type == "Entire home/apt")
# dat_clean$shared_room <- as.numeric(dat_clean$room_type == "Shared room")

rm(dat_mod, dat, dat_raw)

##########################################
##### Preliminary visual analysis ########
##########################################

summary(dat_clean)

par(mfrow=c(2,4))
for (i in col_names_num){
  hist(dat_clean[[i]], main = i, breaks = 25, col = 'peachpuff', border="black",
       prob = TRUE)
  # lines(density(dat_clean[[i]]), # density plot
  #       lwd = 2, # thickness of line
  #       col = "chocolate3")
}
par(mfrow=c(1,1))

par(mfrow=c(1,4))
for (i in col_names_categ){
  boxplot(dat_clean$price ~ dat_clean[[i]], main = i)
}
par(mfrow=c(1,1))

par(mfrow=c(3,3))
for (i in col_names_num){
  plot(dat_clean$price ~ dat_clean[[i]], main = i)
}
par(mfrow=c(1,1))

##########################################
#####    Baseline Regressions     #######
##########################################

#Linear regression
lin_reg_short <- lm(price ~  reviews + overall_satisfaction 
              + bedrooms + bathrooms +  longitude + latitude , data = dat_clean)
summary(lin_reg_short) 
crPlots(lin_reg, terms = ~ bedrooms + bathrooms + reviews + overall_satisfaction, order = 6,
        smooth = TRUE)
coeftest(lin_reg, vcov = vcovHC(lin_reg, type = "HC0"))

#Quantile regression
q_reg_med <- rq(price ~ room_type + reviews + overall_satisfaction 
                + bedrooms + bathrooms +  longitude + latitude + no_reviws, 
                data = dat_clean, tau = 0.5)
summary(q_reg_med)

##########################################
###### Regressions with Smoothing ########
##########################################

#1. Reference linear model with smoothing only lat and long
lin_smooth <-  gam(price ~ s(longitude,latitude) + reviews
                   + bedrooms + bathrooms + overall_satisfaction, data = dat_clean,
                   method="ML")
summary(lin_smooth)
plot(lin_smooth, page = 1, residuals = TRUE)

#2. Conditional mean model allowing for smoothing of other parameters

smoothB <- smoothLUR(data = dat_clean,
                     x = c("longitude", "latitude", "reviews", 
                           'bedrooms', 'bathrooms', 'overall_satisfaction')
                     ,spVar1 = "longitude"
                     ,spVar2 = "latitude"
                     ,y = "price"
                     ,thresh = 0.95)
summary(smoothB)

#3. Additive Quantile Regression Smoothing
smoothQ_m <- rqss(price ~ qss(cbind(longitude,latitude), lambda = .8) + reviews
                        + bedrooms + bathrooms + overall_satisfaction
                        , tau = 0.5, data = dat_clean)
# (Intercept)          34.410966   3.968111   8.672  < 2e-16 ***
#   reviews              -0.055639   0.007634  -7.288 3.49e-13 ***
#   bedrooms             28.112094   1.037775  27.089  < 2e-16 ***
#   bathrooms            21.493507   1.808981  11.882  < 2e-16 ***
#   overall_satisfaction  0.882691   0.347642   2.539   0.0111 *  
# EDF Lambda Penalty F value Pr(>F)    
# cbind(longitude, latitude) 11131    0.8   11138   23.63 <2e-16 ***

smoothQ_0.1 <- rqss(price ~ qss(cbind(longitude,latitude), lambda = .8)  + reviews
                    + bedrooms + bathrooms + overall_satisfaction
                    , tau = 0.1, data = dat_clean)

smoothQ_0.9 <- rqss(price ~ qss(cbind(longitude,latitude), lambda = .8) + reviews
                    + bedrooms + bathrooms + overall_satisfaction
                    , tau = 0.9, data = dat_clean)

smoothQ_0.33 <- rqss(price ~ qss(cbind(longitude,latitude), lambda = .8)  + reviews
                     + bedrooms + bathrooms + overall_satisfaction
                     , tau = 0.33, data = dat_clean)

smoothQ_0.66 <- rqss(price ~ qss(cbind(longitude,latitude), lambda = .8)  + reviews
                     + bedrooms + bathrooms + overall_satisfaction
                     , tau = 0.66, data = dat_clean)

####################################################
##########       Cross-validation     ##############
####################################################

#0. Some in-sample metrics for comparison
# summary(smoothB) #R-sq.(adj) =  0.51 
# summary(lin_smooth) #R-sq.(adj) =  0.377 
# summary(smoothQ_m) #
# fit_q_m <- predict(smoothQ_m, newdata = dat_clean)
# R_2_q <- cor(fit_q_m, dat_clean$price) ^2 #0.3951699
# R_2_q_adj <- 1 - ( ((1 - R_2_q) * (nrow(dat_clean) - 1))/(nrow(dat_clean) - 1 - 5)) #0.3947315

# AIC(smoothB) # 77465.67
# AIC(lin_smooth) # 79109.07
# AIC(smoothQ_m) # 73610.93

#1.Prepare objects for results

n	<- nrow(dat_clean)
R	<- 100	# number of replications

PR2.output	<- array(
  data = NA, 
  dim = c(R, 4), 
  dimnames = list(NULL, c("MLR", "GAM1", "GAM2", "QR"))
)
MSE.output	<- array(
  data = NA, 
  dim = c(R, 4), 
  dimnames = list(NULL, c("MLR", "GAM1", "GAM2", "QR"))
)
MAE.output	<- array(
  data = NA, 
  dim = c(R, 4), 
  dimnames = list(NULL, c("MLR", "GAM1", "GAM2", "QR"))
)
AIC.output <- array(
  data = NA, 
  dim = c(R, 4), 
  dimnames = list(NULL, c("MLR", "GAM1", "GAM2", "QR"))
)

perc	<- 80
n.est	<- round(n*perc/100)

id.est	<- array(NA, c(R, n.est))
id.val	<- array(NA, c(R, n - n.est))

set.seed(42)
for(r in 1:R){
  sam <- sample(1:n, replace = FALSE)
  id.est[r, ]	<- sam[1:n.est]
  id.val[r, ]	<- sam[(n.est+1):n]
}

#2.Metrics computation for R iterations
for(r in 1:R){
  dat.est			<- dat_clean[id.est[r, ], ]
  dat.val			<- dat_clean[id.val[r, ], ]
  
  mlr.est <- lm(price ~  reviews + overall_satisfaction 
                + bedrooms + bathrooms +  longitude + latitude , data = dat.est)
  mlr.fit.est <- fitted(mlr.est)
  AIC.output[r, "MLR"]	<- AIC(mlr.est)
  PR2.output[r, "MLR"]	<- (cor(dat.est$price, mlr.fit.est))^2
  mlr.fit.val		<- predict(mlr.est, newdata = dat.val)
  MSE.output[r, "MLR"]	<- mean((mlr.fit.val - dat.val$price)^2)
  MAE.output[r, "MLR"]	<- mean(abs(mlr.fit.val - dat.val$price))
  
  gam1.est <- gam(price ~ s(longitude,latitude) + reviews
                  + bedrooms + bathrooms + overall_satisfaction, data = dat.est,
                  method="ML")
  gam1.fit.est <- fitted(gam1.est)
  AIC.output[r, "GAM1"]	<- AIC(gam1.est)
  PR2.output[r, "GAM1"]	<- (cor(dat.est$price, gam1.fit.est))^2
  gam1.fit.val		<- predict(gam1.est, newdata = dat.val)
  MSE.output[r, "GAM1"]	<- mean((gam1.fit.val - dat.val$price)^2)
  MAE.output[r, "GAM1"]	<- mean(abs(gam1.fit.val - dat.val$price))
    
  gam2.est <- smoothLUR(data = dat.est,
                        x = c("longitude", "latitude", "reviews", 
                              'bedrooms', 'bathrooms', 'overall_satisfaction')
                        ,spVar1 = "longitude", spVar2 = "latitude"
                        ,y = "price", thresh = 0.95)
  gam2.fit.est <- fitted(gam2.est)
  AIC.output[r, "GAM2"]	<- AIC(gam2.est)
  PR2.output[r, "GAM2"]	<- (cor(dat.est$price, gam2.fit.est))^2
  gam2.fit.val		<- predict(gam2.est, newdata = dat.val)
  MSE.output[r, "GAM2"]	<- mean((gam2.fit.val - dat.val$price)^2)
  MAE.output[r, "GAM2"]	<- mean(abs(gam2.fit.val - dat.val$price))
  
  qr.est <- rqss(price ~ qss(cbind(longitude,latitude), lambda = .8) + reviews
                 + bedrooms + bathrooms + overall_satisfaction
                 , tau = 0.5, data = dat.est)
  qr.fit.est <- fitted(qr.est)
  AIC.output[r, "QR"]	<- AIC(qr.est)
  PR2.output[r, "QR"]	<- (cor(dat.est$price, qr.fit.est))^2
  xy.unique <- unique(dat.est[,c("longitude", "latitude")])
  T <- with(data=xy.unique, expr=tri.mesh(longitude,latitude))
  tri.m <- tri.mesh(x=T$x, y=T$y, duplicate="error")
  ndum <- nrow(dat.val)
  xd <- dat.val$longitude
  yd <- dat.val$latitude
  con.idx <- in.convex.hull(tri.obj=T, x=xd, y=yd)
  qr.fit.val		<- predict(qr.est, newdata = dat.val[con.idx,])
  MSE.output[r, "QR"]	<- mean((qr.fit.val - dat.val[con.idx,'price'])^2)
  MAE.output[r, "QR"]	<- mean(abs(qr.fit.val - dat.val[con.idx,'price']))
}

# CV_output <- as.data.frame(cbind(AIC.output, PR2.output, MSE.output, MAE.output))
# names(CV_output) <- paste( c(rep('AIC',4), rep('PR2',4), rep('MSE',4), rep('MAE',4)),
#                            c("MLR", "GAM1", "GAM2", "QR"))
# write.csv(CV_output, file = 'CV_output.csv')

#3. Output
#3.1. Mean values
round(apply(PR2.output, 2, mean),3)
round(apply(MSE.output, 2, mean),0)
round(apply(MAE.output, 2, mean),0)
round(apply(AIC.output, 2, mean),0)

#3.2.Visualization of empirical distribution functions
f.edf <- function(
    f.measure
){
  f.temp	<- get(paste(f.measure, ".output", sep = ""))
  
  plot(		# for "MLR"
    sort(f.temp[, 1]), (1:R)/R, 
    type = "l", main = paste("Empirical distribution function of", f.measure), 
    xlab = f.measure, ylab = "Empirical distribution function",
    xlim = c(min(f.temp), max(f.temp)), ylim = c(0, 1), col = 1, lwd = 3
  )
  
  lines(		# for "GAM1"
    sort(f.temp[, 2]), (1:R)/R, col = 2, lwd = 3
  )
  
  lines(		# for "GAM2"
    sort(f.temp[, 3]), (1:R)/R, col = 3, lwd = 3
  )
  
  lines(		# for "QR"
    sort(f.temp[, 4]), (1:R)/R, col = 4, lwd = 3
  )
  
  abline(v = apply(f.temp, 2, min), col = 1:4, lwd = 1) # Minima of the EDFs
  abline(v = apply(f.temp, 2, max), col = 1:4, lwd = 1) # Maxima of the EDFs
  abline(h = 0.5, col = "darkgrey", lwd = 2, lty = 2)
  
  legend("bottomright", c("MLR", "GAM1", "GAM2", "QR"), col = 1:4, lwd = 4, bg = "white")
}

par(mfrow = c(2, 2), cex = 0.7)
f.edf("PR2")
f.edf("MSE")
f.edf("MAE")
f.edf("AIC")
par(mfrow = c(1, 1), cex = 1)

####################################################
# Diagnostic plot to check for quantile effects  ###
####################################################
#From linear reference model
coef_mean <- summary(lin_smooth)$p.coeff
se_mean <- summary(lin_smooth)$se
#Coefficients from quantile regressions
coef_0.1 <- summary(smoothQ_0.1)$coef
coef_0.33 <- summary(smoothQ_0.33)$coef
coef_median <- summary(smoothQ_m)$coef
coef_0.66 <- summary(smoothQ_0.66)$coef
coef_0.9 <- summary(smoothQ_0.9)$coef
#replace with sd

quant_coef <- rbind(coef_0.1[,1], coef_0.33[,1], coef_median[,1], coef_0.66[,1], coef_0.9[,1])
quantile <- c(.1,.33, .5, .66, .9)
quant_coef <- as.data.frame(cbind(quantile, quant_coef))
names(quant_coef)
# write.csv(quant_coef, file='quant_coef.csv')
quant_coef <- read.csv('quant_coef.csv')

quant_se <- rbind(coef_0.1[,2], coef_0.33[,2], coef_median[,2], coef_0.66[,2], coef_0.9[,2])
quantile <- c(.1,.33, .5, .66, .9)
quant_se <- as.data.frame(cbind(quantile, quant_se))
write.csv(quant_se, file='quant_se.csv')
# quant_se <- read.csv('quant_se.csv')

upper_bond <- quant_coef[,c(3:6)] + 1.96 * quant_se[,c(3:6)]
lower_bond <- quant_coef[,c(3:6)] - 1.96 * quant_se[,c(3:6)]

#Visualization
reviews_pl <- ggplot(data = quant_coef, aes(x = quantile, y = reviews)) +
  geom_point(col = 'red') +
  geom_line(col = 'red') +
  geom_hline(yintercept = coef_mean[2], col = 'blue') +
  geom_hline(yintercept = c(coef_mean[2] - (1.96 * se_mean[2]), coef_mean[2] + (1.96 * se_mean[2])),
             linetype='dotted', col = 'blue') +
  geom_ribbon(aes(ymin=upper_bond[,1],ymax=upper_bond[,1]),alpha=0.5)

bedrooms_pl <- ggplot(data = quant_coef, aes(x = quantile, y = bedrooms)) +
  geom_point(col = 'red') +
  geom_line(col = 'red') +
  geom_hline(yintercept = coef_mean[3], col = 'blue') +
  geom_hline(yintercept = c(coef_mean[3] - (1.96 * se_mean[3]), coef_mean[3] + (1.96 * se_mean[3])),
             linetype='dotted', col = 'blue')
geom_ribbon(aes(ymin=upper_bond[,1],ymax=upper_bond[,1]),alpha=0.5)

bathrooms_pl <- ggplot(data = quant_coef, aes(x = quantile, y = bathrooms)) +
  geom_point(col = 'red') +
  geom_line(col = 'red') +
  geom_hline(yintercept = coef_mean[4], col = 'blue') +
  geom_hline(yintercept = c(coef_mean[4] - (1.96 * se_mean[4]), coef_mean[4] + (1.96 * se_mean[4])),
             linetype='dotted', col = 'blue')
geom_ribbon(aes(ymin=upper_bond[,1],ymax=upper_bond[,1]),alpha=0.5)

overall_satisfaction_pl <- ggplot(data = quant_coef, aes(x = quantile, y = overall_satisfaction)) +
  geom_point(col = 'red') +
  geom_line(col = 'red') +
  geom_hline(yintercept = coef_mean[5], col = 'blue') +
  geom_hline(yintercept = c(coef_mean[5] - (1.96 * se_mean[5]), coef_mean[5] + (1.96 * se_mean[5])),
             linetype='dotted', col = 'blue')
geom_ribbon(aes(ymin=upper_bond[,1],ymax=upper_bond[,1]),alpha=0.5)

grid.arrange(bathrooms_pl, bedrooms_pl, overall_satisfaction_pl, reviews_pl, ncol = 2)
rm(overall_satisfaction_pl,bathrooms_pl, bedrooms_pl, reviews_pl)

##########################################
########## Seattle map grid ##############
##########################################

#data set area
bottom_left_dat <- c(min(dat_clean$latitude), min(dat_clean$longitude))
top_right_dat <- c(max(dat_clean$latitude), max(dat_clean$longitude))
dat_area <- rast(xmin=bottom_left_dat[2], xmax=top_right_dat[2], 
                 ymin=bottom_left_dat[1], ymax=top_right_dat[1],
                 res = c(0.0001086422, 0.0001086422))
# lat_step * 111 #1 cell is (12 x 12) ~140 m2
r <- rast("pm_heat_index_f_ranger-2.tif") * 0 + 1 
r_mod <- project(r,  'EPSG:4326') #change coordinate references
r_mod_cropped <- crop(r_mod, dat_area) #crop the initial map to match the data set
plot(r_mod_cropped, col = 'orange', main = 'Properties from the sample', cex=1.5)
points(dat_clean$longitude, dat_clean$latitude, 
       pch = 16,
       col = transparent("darkorange4", trans.val = .8))

####################################################
##########  Price surface prediction  ##############
####################################################

#Conditional mean prediction

#1. Create a raster object with layers of predictors
#1.1 longitude and latitude
lon_seq <- seq(from = bottom_left_dat[2], to = top_right_dat[2], length.out = ncol(r_mod_cropped))
lat_seq <- seq(from = bottom_left_dat[1], to = top_right_dat[1], length.out = nrow(r_mod_cropped))

lat_vec <- c()
for (i in 1:length(lat_seq)) {
  lat_vec <- append(lat_vec, rep(lat_seq[i], ncol(r_mod_cropped)), after = length(lat_vec))
}
lon_vec <- rep(lon_seq, nrow(r_mod_cropped)) #data vectors to fill the layers of lon and lat

longitude <- rast(xmin=bottom_left_dat[2], xmax=top_right_dat[2], 
                  ymin=bottom_left_dat[1], ymax=top_right_dat[1],
                  res = c(0.0001086422, 0.0001086422))
longitude <- resample(longitude, r_mod_cropped) #align extents
names(longitude) <- 'longitude'
longitude <- setValues(longitude, lon_vec)

latitude <- rast(xmin=bottom_left_dat[2], xmax=top_right_dat[2], 
                 ymin=bottom_left_dat[1], ymax=top_right_dat[1],
                 res = c(0.0001086422, 0.0001086422))
latitude <- resample(latitude, r_mod_cropped) #align extents
names(latitude) <- 'latitude'
latitude <- setValues(latitude, lat_vec)

#1.2 Other variables as all-zero layers: reviews, bedrooms, bathrooms, accommodates
#alternatively set them eqal to median?
vec_rev <- rep(median(dat_clean$reviews), ncell(r_mod_cropped))
reviews <- setValues(longitude, vec_rev)
names(reviews) <- 'reviews'

vec_bed <- rep(median(dat_clean$bedrooms), ncell(r_mod_cropped))
bedrooms <- setValues(longitude, vec_bed)
names(bedrooms) <- 'bedrooms'

vec_bath <- rep(median(dat_clean$bathrooms), ncell(r_mod_cropped))
bathrooms <- setValues(longitude, vec_bath)
names(bathrooms) <- 'bathrooms'

vec_sat <- rep(median(dat_clean$overall_satisfaction), ncell(r_mod_cropped))
overall_satisfaction <- setValues(longitude, vec_sat)
names(overall_satisfaction) <- 'overall_satisfaction'

vec_no_rew <- rep(median(dat_clean$no_reviws), ncell(r_mod_cropped))
no_reviws <- setValues(longitude, vec_no_rew)
names(no_reviws) <- 'no_reviws'

vec_room_t <- rep('Entire home/apt', ncell(r_mod_cropped))
room_type <- setValues(longitude, vec_room_t)
names(room_type) <- 'room_type'

vec_ent_home <- rep(median(dat_clean$entire_home), ncell(r_mod_cropped))
entire_home <- setValues(longitude, vec_ent_home)
names(entire_home) <- 'entire_home'

vec_pr_room <- rep(median(dat_clean$private_room), ncell(r_mod_cropped))
private_room <- setValues(longitude, vec_pr_room)
names(private_room) <- 'private_room'

vec_sh_room <- rep(median(dat_clean$shared_room), ncell(r_mod_cropped))
shared_room <- setValues(longitude, vec_sh_room)
names(shared_room) <- 'shared_room'

lon_lat_surf <- stack(x = c(latitude, longitude,reviews, bedrooms, 
                            bathrooms, overall_satisfaction,no_reviws, 
                            shared_room, private_room, entire_home, room_type)) #add layers to a single object
plot(lon_lat_surf)

#2. Predict a price surface for the locations grid and all other variables equal
#to the median values
#2.1. for SmoothLUR model
price_predict <- predict(object = lon_lat_surf, model = smoothA)
names(price_predict) <- 'price'
price_predict_rast <- terra::rast(price_predict) #convert to SpatRaster
# write.csv(as.data.frame(price_predict), file='price_meanmeth_medianflat.csv')
# price_predict <- read.csv('price_meanmeth_medianflat.csv')

#2.2. for GAM model
price_predict_gam <- predict(object = lon_lat_surf, model = lin_smooth)
names(price_predict_gam) <- 'price'
price_predict_gam_rast <- terra::rast(price_predict_gam)

#3. Combine the map and the price surface
#3.1. for SmoothLUR model
price_surf <- cover(r_mod_cropped, price_predict_rast,  values=1)
plot(price_surf, main = 'Price surface of Seattle Airbnb rent')

#3.2. for GAM model
price_surf_gam <- cover(r_mod_cropped, price_predict_gam_rast,  values=1)
plot(price_surf_gam, main = 'Price surface of Seattle Airbnb rent')

#Tau-quantile prediction
#1. Logical vector to check if an observation is within the convex hull
set.seed(123)
new_dat <- as.data.frame(lon_lat_surf) #convert a multilayer raster object into a data frame
xy.unique <- unique(dat_clean[,c("longitude", "latitude")]) #unique spatial points
T <- with(data=xy.unique, expr=tri.mesh(longitude,latitude))
tri.m <- tri.mesh(x=T$x, y=T$y, duplicate="error")
ndum <- nrow(new_dat)
xd <- new_dat$longitude
yd <- new_dat$latitude
con.idx <- in.convex.hull(tri.obj=T, x=xd, y=yd)
table(con.idx) * 100 / nrow(new_dat) #~6.5% loss

#2. Prediction
#2.1 Median
nd <- new_dat[con.idx,]
pred_q_m_med <- predict.rqss(object = smoothQ_m, newdata = nd)

predictions <- as.data.frame(pred_q_m_med)
names(predictions) <- c('median_medflad')
write.csv(predictions, file='price_prediction_medflat_medmod.csv')
# predictions <- read.csv(file = 'price_prediction_medflat_medmod.csv')

#2.2 For tau = 0.1
pred_q_0.1 <- predict.rqss(object = smoothQ_0.1, newdata = nd)

predictions_0.1 <- as.data.frame(pred_q_0.1)
names(predictions_0.1) <- c('0.1_medflad')
write.csv(predictions_0.1, file='price_prediction_medflat_medq01.csv')

#2.3 For tau = 0.9
pred_q_0.9 <- predict.rqss(object = smoothQ_0.9, newdata = nd)

predictions_0.9 <- as.data.frame(pred_q_0.9)
names(predictions_0.9) <- c('0.9_medflad')
write.csv(predictions_0.9, file='price_prediction_medflat_medq09.csv')

#3.Back to raster objects
#3.1 Median
median_rast <- new_dat[,c('longitude','latitude')]
median_rast$price <- mean(pred_q_m_med)
median_rast[con.idx,'price'] <- pred_q_m_med
price_predict_q_m_med_rast <- setValues(longitude, median_rast$price) #convert to SpatRaster
names(price_predict_q_m_med_rast) <- 'price'
price_predict_q_m_med_rast <- cover(r_mod_cropped, price_predict_q_m_med_rast,  values=1)

#3.2 For tau = 0.1
tau0.1_rast <- new_dat[,c('longitude','latitude')]
tau0.1_rast$price <- mean(pred_q_0.1)
tau0.1_rast[con.idx,'price'] <- pred_q_0.1
price_predict_q_0.1_rast <- setValues(longitude, tau0.1_rast$price) #convert to SpatRaster
names(price_predict_q_0.1_rast) <- 'price'
price_surf_q_0.1 <- cover(r_mod_cropped, price_predict_q_0.1_rast,  values=1)

#3.3 For tau = 0.9
tau0.9_rast <- new_dat[,c('longitude','latitude')]
tau0.9_rast$price <- mean(pred_q_0.9)
tau0.9_rast[con.idx,'price'] <- pred_q_0.9
price_predict_q_0.9_rast <- setValues(longitude, tau0.9_rast$price) #convert to SpatRaster
names(price_predict_q_0.9_rast) <- 'price'
price_surf_q_0.9 <- cover(r_mod_cropped, price_predict_q_0.9_rast,  values=1)

#4. Vizualization
#4.1 Price surfaces
par(mfrow=c(1,3)) 
plot(price_surf_gam, main = 'Price surface of Airbnb in Seattle - GAM I', cex=1.5 )
plot(price_surf, main = 'Price surface of Airbnb in Seattle - GAM II',cex=1.5 )
plot(price_predict_q_m_med_rast, main = 'Price surface of Airbnb in Seattle - QR', cex=1.5 )
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot(price_surf_q_0.1, main = 'Tau = 0.1 Price surface'
    , cex=1.5)
plot(price_predict_q_m_med_rast, main = 'Median Price surface'
    , cex=1.5)
plot(price_surf_q_0.9, main = 'Tau = 0.9 Price surface'
    , cex=1.5)
par(mfrow=c(1,1))



