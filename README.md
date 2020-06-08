# SSNutils: Miscellaneous functions for stream networks
## Generating artificial stream network data

![Alt text](https://github.com/EdgarSantos-Fernandez/SSNutils/blob/master/temp_ggplot.png?raw=true "Title")
Fig: Observation and prediction sites in a stream network.

Generating artificial stream network data

```{r, warning=F, message=F}
library('Rcpp')
library('ggplot2')
library('dplyr')
library('reshape2')
library('spdep')
library('RColorBrewer')
library('ggrepel')
library('viridis')
library("SSN")
library('SSNdesign')
```


In this tutorial we generate artificial stream networks/data, fit spatial models, and make predictions using the R packages ```SSN``` and ```SSNdesign``` [@SSN; @SSNdesign].

# Spatial data
## Simulated stream network and data

### Creating spatial stream network objects 

The function ```SSN::createSSN()``` generates artificial branching networks and (observation/predition) points.

*Arguments:*

- ```n```: is the vector. length = number of networks. values = number of line segments

- ```obsDesign```: sampling design for the observations. Arguments: binomialDesign, systematicDesign, poissonDesign, etc.

- ```predDesign```: sampling design for the predictions. Arguments: binomialDesign, systematicDesign, poissonDesign, etc.

- ```path```: path where the .ssn object will be saved.

See the other arguments in ```?SSN::createSSN()```.

Let us simulate two networks with 100 segments each. 
We generate 60 and 40 observation points on the networks and 40 and 20 prediction sites respectively. 

```{r, warning=F, message=F, cache = TRUE}
set.seed(1)

path <- "./ex1.ssn" # where to save your SSN object.

rewrite <- T
if(rewrite == T){ # if the layer already exist and we want to rewrite it

  if(file.exists(path) == T){
    unlink(path, recursive=TRUE)
    file.remove(path)}
  
ssn <- createSSN(
  n = c(100, 100), # lenght = number of networks / # segments to grow
  obsDesign = binomialDesign(c(60, 40)),
  predDesign = binomialDesign(c(40, 20)),
  importToR = TRUE,
  path = path,
  treeFunction = iterativeTreeLayout
)
}
```

Plotting can be done directly using the included ```plot``` (```SSN::plot.SpatialStreamNetwork()```) function. 

```{r, warning=F, message=F, cache = TRUE}
col <- brewer.pal(9,'YlGnBu')[6]

plot(
  ssn,
  lwdLineCol = "addfunccol",
  lwdLineEx = 8,
  lineCol = col,
  col = 1,
  pch = 16,
  xlab = "x-coordinate",
  ylab = "y-coordinate"
)

# adding the predictions sites
plot(ssn, PredPointsID = "preds", add = T, pch = 16, col = "#E41A1C")
```



### Simulating data on the sites

We first calculate the hydrologic distance matrix and simulate some continuous and categorical covariates:
```{r, warning=F, message=F, cache = TRUE}
createDistMat(ssn, "preds", o.write=TRUE, amongpred = TRUE)

rawDFobs <- getSSNdata.frame(ssn, "Obs")
rawDFpred <- getSSNdata.frame(ssn, "preds")

# generating continous covariates
rawDFobs[,"X1"] <- rnorm(length(rawDFobs[,1]))
rawDFpred[,"X1"] <- rnorm(length(rawDFpred[,1]))
rawDFobs[,"X2"] <- rnorm(length(rawDFobs[,1]))
rawDFpred[,"X2"] <- rnorm(length(rawDFpred[,1]))

# generating categorical covariates
rawDFobs[,"F1"] <- as.factor(sample.int(4,length(rawDFobs[,1]),
                                        replace = TRUE))
rawDFpred[,"F1"] <- as.factor(sample.int(4,length(rawDFpred[,1]),
                                         replace = TRUE))

# random effects
rawDFobs[,"RE1"] <- as.factor(sample(1:3,length(rawDFobs[,1]),
                                     replace = TRUE))
rawDFobs[,"RE2"] <- as.factor(sample(1:4,length(rawDFobs[,1]),
                                     replace = TRUE))
rawDFpred[,"RE1"] <- as.factor(sample(1:3,length(rawDFpred[,1]),
                                      replace = TRUE))
rawDFpred[,"RE2"] <- as.factor(sample(1:4,length(rawDFpred[,1]),
                                      replace = TRUE))
names(rawDFobs)
names(rawDFpred)

# updating the observation and prediction data on the SSN object
ssn <- putSSNdata.frame(rawDFobs,ssn, Name = 'Obs')
ssn <- putSSNdata.frame(rawDFpred,ssn, Name = 'preds')
```

Now using the function ```SSNdesign::SimulateOnSSN()``` to generate the response variable that will be stored in the column ```Sim_Values```.
```{r, warning=F, message=F, cache = TRUE}
set.seed(2020)
sim.out <- SimulateOnSSN(
  ssn,
  ObsSimDF = rawDFobs,
  PredSimDF = rawDFpred,
  PredID = "preds",
  formula = ~ X1 + X2 + F1,
  coefficients = c(10, 1, 0, -2, 0, 2),
  CorModels = c(
    "LinearSill.tailup", 
    "RE1", 
    "RE2"),
  use.nugget = TRUE,
  CorParms = c(3, 10,  1, .5, .1),#partial sill, range, random effects variance components ordered alphanumerically and nugget 
  addfunccol = "addfunccol"
)
```

The figure below shows the values of the response variable at the observation sites.

```{r, warning=F, message=F, cache = TRUE}
sim.ssn <- sim.out$ssn.object

plot(
  sim.ssn,"Sim_Values",
  lwdLineCol = "addfunccol",
  lwdLineEx = 8,
  lineCol = col,
  pch = 16,
  xlab = "x-coordinate",
  ylab = "y-coordinate"
)
```


We fit the spatial mixed effect stream network model using the function ```glmssn()```, which requires a formula object as in ```lm()```, an SSN object, 
the spatial autocorrelation structure for stream networks e.g. tail-up with linear sill and the variable used to compute the spatial weights (```addfunccol```). 
We then, use the prediction sites to estimate the response variable.  
See @SSN for more details.

```{r, warning=F, message=F, cache = TRUE}
simDFobs <- getSSNdata.frame(sim.ssn, "Obs")
simDFpred <- getSSNdata.frame(sim.ssn, "preds")

simpreds <- simDFpred[,"Sim_Values"]
simDFpred[,"Sim_Values"] <- NA
sim.ssn <- putSSNdata.frame(simDFpred, sim.ssn, "preds")

# obtaining the predictions 
glmssn.out <- glmssn(Sim_Values ~ X1 + X2 + F1, 
                     sim.ssn,
                     CorModels = c("LinearSill.tailup", "RE1", "RE2"), 
                     addfunccol = "addfunccol")
summary(glmssn.out)

# predictions
glmssn.pred <- predict(glmssn.out, "preds")
predDF <- getSSNdata.frame(glmssn.pred, "preds")
```

How well the prediction sites are recovered using the tail up model?
```{r, warning=F, message=F, cache = TRUE, fig.height = 4, fig.width = 4}
ggplot() +
  geom_point(aes(x = simpreds,
                 y = predDF[, "Sim_Values"])) +
                 xlab("True") +
                 ylab("Predicted") +
                 coord_fixed() +  
                 theme_bw()
```


### Plotting it using ggplot
The SSN object can be also visualized using ggplot2.
But, we need to extract to a data frame the segments from the SSN object. 


```{r, warning=F, message=F, cache = TRUE}
collapse <- function(t){
  df_all <- NULL
  for (i in 1:length(t@lines)){
    df <- data.frame(t@lines[[i]]@Lines[[1]]@coords)
    df$slot <- t@lines[[i]]@ID
    df$computed.afv <-  t@data$computed.afv[i]
    df_all<- rbind(df, df_all)
    
    df_all$slot <- as.numeric(as.character(df_all$slot))
  }
  df_all <-  dplyr::arrange(df_all, slot)
  df_all
}

sim.ssn.df <- collapse(sim.ssn)

obs_data <- getSSNdata.frame(sim.ssn, "Obs")
pred_data <- getSSNdata.frame(sim.ssn, "preds")
pred_data$Sim_Values <- simpreds

obs_data_coord <- data.frame(ssn@obspoints@SSNPoints[[1]]@point.coords)
obs_data_coord$locID <- factor(1:nrow(obs_data_coord))
pred_data_coord <- data.frame(ssn@predpoints@SSNPoints[[1]]@point.coords)
pred_data_coord$locID <- factor((nrow(obs_data_coord)+1):(nrow(obs_data_coord)+nrow(pred_data_coord)))

obs_data <- obs_data %>% left_join(obs_data_coord, by = c('locID'))
obs_data$point <- 'Obs'
pred_data <- pred_data %>% left_join(pred_data_coord, by = c('locID'))
pred_data$point <- 'Pred'
obs_pred_points <- rbind(obs_data, pred_data)

sim.ssn.df$addfunccol <- rep(sim.ssn@data$addfunccol, each = 2)
sim.ssn.df$addfunccol_cat <- cut(sim.ssn.df$addfunccol, 
                                 breaks = seq(min(sim.ssn.df$addfunccol),
                                            max(sim.ssn.df$addfunccol),
                                            length.out=6),
                                 labels = 1:5,
                                 include.lowest = T)

ggplot(sim.ssn.df) + 
  geom_path(aes(X1,X2, group = slot, size = addfunccol_cat), lineend = 'round', linejoin = 'round', col = col)+
  geom_point(data = obs_pred_points, aes(x = coords.x1, y = coords.x2, col = Sim_Values, shape = point))+
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  scale_color_viridis(option = 'C')+
  scale_shape_manual(values = c(16,15))+
  xlab("x-coordinate") +
  ylab("y-coordinate")+
  coord_fixed()+
  theme_bw()
```


![Alt text](https://github.com/EdgarSantos-Fernandez/SSNutils/blob/master/sim_ggplot.png?raw=true "Title")
Fig: Artificial stream network.



## Using an existing stream network topology to simulate data.
We use a stream network dataset with temperature measurements from the
Middle Fork in the Salmon River, Idaho, USA.
The full dataset can be found at https://www.fs.fed.us/rm/boise/AWAE/projects/SSN_STARS/software_data.html


```{r, warning=F, message=F, cache = TRUE}
paths <- "MedSpacetime.ssn"
t  <- importSSN(paths,
                o.write = TRUE)
```

The dataset is composed of 400 observations from 80 unique locations.

```{r, warning=F, message=F, cache = TRUE, fig.height = 5, fig.width = 5}
plot(
  t,
  lwdLineCol = "afvArea",
  lwdLineEx = 10,
  lineCol = col,
  pch = 16,
  xlab = "x-coordinate (m)",
  ylab = "y-coordinate (m)",
  asp = 1
)
```

We will generate new observation and prediction locations over the existing river network. We overwrite the existing data points using a binomial design with 100 observations and 50 prediction points.

This takes a few mins:
```{r, warning=F, message=F, cache = TRUE, fig.height = 5, fig.width = 5}
t2 <- generateSites(t, 
                    obsDesign = binomialDesign(100),
                    predDesign = binomialDesign(50),
                    o.write = TRUE) 
```

```{r, warning=F, message=F, cache = TRUE, fig.height = 5, fig.width = 5}
plot(
  t2,
  lwdLineCol = "afvArea",
  lwdLineEx = 8,
  lineCol = col,
  col = 1,
  pch = 16,
  xlab = "x-coordinate",
  ylab = "y-coordinate"
)

plot(t2, PredPointsID = "preds", add = T, pch = 15, col = "#E41A1C")
```

Simulating data on the SSN object:

```{r, warning=F, message=F, cache = TRUE}
createDistMat(t2, "preds", o.write=TRUE, amongpred = TRUE)

t3 <- additive.function(t2, "h2oAreaKm2", "computed.afv")

raw_t_obs <- getSSNdata.frame(t3, "Obs")
raw_t_pred <- getSSNdata.frame(t3, "preds")

# generating continous covariates
raw_t_obs[,"X1"] <- rnorm(length(raw_t_obs[,1]))
raw_t_pred[,"X1"] <- rnorm(length(raw_t_pred[,1]))
raw_t_obs[,"X2"] <- rnorm(length(raw_t_obs[,1]))
raw_t_pred[,"X2"] <- rnorm(length(raw_t_pred[,1]))

# generating categorical covariates
raw_t_obs[,"F1"] <- as.factor(sample.int(4,length(raw_t_obs[,1]),replace = TRUE))
raw_t_pred[,"F1"] <- as.factor(sample.int(4,length(raw_t_pred[,1]),replace = TRUE))

# random effects
raw_t_obs[,"RE1"] <- as.factor(sample(1:3,length(raw_t_obs[,1]),replace = TRUE))
raw_t_obs[,"RE2"] <- as.factor(sample(1:4,length(raw_t_obs[,1]),replace = TRUE))
raw_t_pred[,"RE1"] <- as.factor(sample(1:3,length(raw_t_pred[,1]),replace = TRUE))
raw_t_pred[,"RE2"] <- as.factor(sample(1:4,length(raw_t_pred[,1]),replace = TRUE))

```

```{r, warning=F, message=F, cache = TRUE}
t3 <- putSSNdata.frame(raw_t_obs,t3, Name = 'Obs')
t3 <- putSSNdata.frame(raw_t_pred,t3, Name = 'preds')
head(t3@data)
head(getSSNdata.frame(t3)[, c("computed.afv")])
```


Simulated values at the locations:
```{r, warning=F, message=F, cache = TRUE, fig.height = 5.5, fig.width = 5.5}
set.seed(2020)
# simulated values will be added in the column "Sim_Values"
sim_t <- SimulateOnSSN(
  t3,
  ObsSimDF = raw_t_obs,
  PredSimDF = raw_t_pred,
  PredID = "preds",
  formula = ~ X1 + X2 + F1,
  coefficients = c(10, 1, 0, -2, 0, 2),
  CorModels = c(
    "LinearSill.tailup",
    "RE1", 
    "RE2"),
  use.nugget = TRUE,
  CorParms = c(3, 10,  1, .5, .1),#partial sill, range, random effects variance components ordered alphanumerically and nugget #only tail up  2, 10, 1, 5,
  addfunccol = "computed.afv"
)

sim.ssn <- sim_t$ssn.object


plot(
  sim.ssn,"Sim_Values",
  lwdLineCol = "afvArea",
  lwdLineEx = 8,
  lineCol = col,
  pch = 16,
  xlab = "x-coordinate",
  ylab = "y-coordinate",
  asp = 1
)
```


Using the network model to obtain predictions: 
```{r, warning=F, message=F, cache = TRUE}
simDFobs <- getSSNdata.frame(sim.ssn, "Obs")
simDFpred <- getSSNdata.frame(sim.ssn, "preds")

# saving the simulated values in the prediction points
simpreds <- simDFpred[,"Sim_Values"]
# setting the prediction points to NA
simDFpred[,"Sim_Values"] <- NA 
sim.ssn <- putSSNdata.frame(simDFpred, sim.ssn, "preds")
```


```{r, warning=F, message=F, cache = TRUE, fig.height = 5, fig.width = 5}
glmssn.out <- glmssn(Sim_Values ~ X1 + X2 + F1, 
                     sim.ssn,
                     CorModels = c("LinearSill.tailup", "RE1", "RE2"), 
                     addfunccol = "computed.afv")

summary(glmssn.out)

glmssn.pred <- predict(glmssn.out, "preds")
predDF <- getSSNdata.frame(glmssn.pred, "preds")

```

How well the prediction sites are recovered using the model?
```{r, warning=F, message=F, cache = TRUE, fig.height = 4, fig.width = 4}
ggplot() +
  geom_point(aes(x = simpreds,
                 y = predDF[, "Sim_Values"])) +
                 xlab("True") +
                 ylab("Predicted") +
                 coord_fixed() +  
                 theme_bw()
```

Plotting the observation and prediction points using ggplot:

```{r, warning=F, message=F, cache = TRUE, fig.height = 7, fig.width = 7}
sim.ssn.df <- collapse(sim.ssn)

obs_data <- getSSNdata.frame(sim.ssn, "Obs")
pred_data <- getSSNdata.frame(sim.ssn, "preds")
pred_data$Sim_Values <- simpreds

obs_data_coord <- data.frame(sim.ssn@obspoints@SSNPoints[[1]]@point.coords)
obs_data_coord$locID <- factor(1:nrow(obs_data_coord))
pred_data_coord <- data.frame(sim.ssn@predpoints@SSNPoints[[1]]@point.coords)
pred_data_coord$locID <- factor((nrow(obs_data_coord)+1):(nrow(obs_data_coord)+nrow(pred_data_coord)))

obs_data <- obs_data %>% left_join(obs_data_coord, by = c('locID'))
obs_data$point <- 'Obs'
pred_data <- pred_data %>% left_join(pred_data_coord, by = c('locID'))
pred_data$point <- 'Pred'
obs_pred_points <- rbind(obs_data, pred_data)

sim.ssn.df$addfunccol_cat <- cut(sim.ssn.df$computed.afv, 
                                 breaks = seq(min(sim.ssn.df$computed.afv),
                                              max(sim.ssn.df$computed.afv),
                                              length.out=6),
                                 labels = 1:5,
                                 include.lowest = T)

ggplot(sim.ssn.df) + 
  geom_path(aes(X1,X2, group = slot, size = addfunccol_cat), lineend = 'round', linejoin = 'round', col = col)+ 
  geom_point(data = obs_pred_points, aes(x = coords.x1, y = coords.x2, col = Sim_Values, shape = point))+
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  scale_color_viridis(option = 'C')+
  scale_shape_manual(values = c(16,15))+
  xlab("x-coordinate") +
  ylab("y-coordinate")+
  coord_fixed()+
  theme_bw()

```


# Spatio-temporal stream network data [It will be added soon]

```{r, warning=F, message=F, cache = TRUE}
df_t <- collapse(t)
t_data <- t@obspoints@SSNPoints[[1]]@point.data
t_data$date <- as.Date(t_data$date_)

ggplot(df_t) + geom_path(aes(X1,X2, group = slot))+
  geom_point(data = t_data, aes(x = NEAR_X, y = NEAR_Y, col = Daily_mn))+
  facet_wrap(~date)+
  scale_color_viridis()+
  coord_fixed()+
  theme_bw()
```






# Bibliography
Jay M. Ver Hoef, Erin E. Peterson, David Clifford, Rohan Shah (2014). SSN: An R Package for Spatial Statistical Modeling on
  Stream Networks. Journal of Statistical Software, 56(3), 1-43. URL http://www.jstatsoft.org/v56/i03/.
Alan Pearse,James McGree, Nicholas Som, Jay Ver Hoef, and Erin Peterson. 2020. SSNdesign: An R Package for Optimal and Adaptive Sampling Designs on Stream Networks.
