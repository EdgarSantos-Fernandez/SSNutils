# SSNutils: Miscellaneous functions for stream networks

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
```


# Spatial stream network data 

The function ```SSN::createSSN()``` generates artificial branching networks and (observation/predition) points.

*Arguments:*

- ```n```: is the vector. length = number of networks. values = number of line segments

- ```obsDesign```: sampling design for the observations. Arguments: binomialDesign, systematicDesign, poissonDesign, etc.

- ```predDesign```: sampling design for the predictions. Arguments: binomialDesign, systematicDesign, poissonDesign, etc.

- ```path```: path where the .ssn object will be saved.

See other arguments in ?```SSN::createSSN()```.

```{r, warning=F, message=F}
set.seed(1)

path <- "./ex1.ssn" # where to save your SSN object.

rewrite <- T
if(rewrite == T){ # if the layer already exist and we want to rewrite it

  if(file.exists(path) == T){
    unlink(path, recursive=TRUE)
    file.remove(path)}
  
ssn1 <- createSSN(
  n = c(60, 30), # lenght = number of networks / # segments to grow
  obsDesign = binomialDesign(c(20, 10)),
  predDesign = binomialDesign(c(5, 3)),
  importToR = T,
  path = path,
  treeFunction = iterativeTreeLayout
)
}
```


```{r, warning=F, message=F}
col <- brewer.pal(9,'YlGnBu')[6]
col_points <- brewer.pal(8,'Paired')[c(2,4)]

plot(
  ssn1,
  lwdLineCol = "addfunccol",
  lwdLineEx = 8,
  lineCol = col,
  cex = 1,
  col = "#33A02C",
  xlab = "x-coordinate",
  ylab = "y-coordinate",
  pch = 1
)
# adding the predictions sites
plot(ssn1, PredPointsID = "preds", add = T, col = "#E41A1C")
```


### Using ggplot:

```{r, warning=F, message=F}
collapse <- function(t){ # will collapse the .ssn line segments into a data.frame.
  df_all <- NULL
  for (i in 1:length(t@lines)){
    df <- data.frame(t@lines[[i]]@Lines[[1]]@coords)
    df$slot <- t@lines[[i]]@ID
    df_all<- rbind(df, df_all)
    
    df_all$slot <- as.numeric(as.character(df_all$slot))
  }
  df_all
}

ssn1_df <- collapse(ssn1)

# getting from the .ssn object the obs and pred data
obs_data <- data.frame(ssn1@obspoints@SSNPoints[[1]]@point.coords)
pred_data <- data.frame(ssn1@predpoints@SSNPoints[[1]]@point.coords)


ggplot(ssn1_df) + 
  geom_path(aes(X1,X2, group = slot), col = col)+
  geom_point(data = obs_data, aes(x = coords.x1, y = coords.x2),col = "#33A02C")+
  geom_point(data = pred_data, aes(x = coords.x1, y = coords.x2), col = "#E41A1C")+
  xlab("x-coordinate") +
  ylab("y-coordinate")+
  coord_fixed()+
  theme_bw()
```
