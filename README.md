#agRank: the Package for Rank Aggregation in Phenotypic Selection

------

##Introduction
In phenotypic selection, we usually want to measure some traits which are difficult for researchers to measure quantitatively, like the taste and the ease of cooking. In this situation, we may adopt ``crowdsourcing'' methodology to incorporate farmers' wise. The general setup is that each farmer receives three varieties and ranks the varieties according to some trait from the best to the worst. After receiving all the partial rankings from farmers, we may aggregate them into an overall consensus ranking so that the overall ranking can reflect the level of the traits we ask farmers to measure. Meanhile, what we have before the actual field trials is the genotype data. Since the additive relationship matrix can be derived from the genotype data and this matrix measures the simularity among all the varieties, we can incorporate the information in this matrix into the rank aggregation process. We use a simple linear model, the Thurstone model, Bradley-Terry model and Plackett-Luce model to do the rank aggregation. MAP estimators are derived to incorporate genotypic information. Various experimental designs are also explored to improve the performance of those models. For details, please check the vignettes as well as the paper [RankAg_Paper.pdf](https://github.com/shuxiaoc/agRank/blob/master/RankAg_Paper.pdf) attached.

##Main Components
There are 4 main components in this package. Firstly, the simulation results are included in this package. Secondly, there are four sampling functions which can sample *complete* rankings from Bradley-Terry model, Plackett-Luce model, Thurstone model and Mallows model. Thirdly, the function of generating experimental designs (allocation of varieties to farmers) is included. And lastly, rank aggregation functions using linear model, Bradley-Terry model, Plackett-Luce model, Thurstone model and Multinomial Preference Model are included. There is also one wrapper function for those rank aggregation functions. 

Note: The names of the experimental designs in this package are different with the function names in the paper. In the paper, only KM4, KM6, GD3, GD6 and GD9 are included, and they are namely as KM1, KM2, GD1, GD2 and GD3 respectively. Meanwhile, the Multinomial Preference model is not included in the paper since it's relatively slow and cannot defeat the other models in terms of accuracy.

##Installation
Package *devtools* is needed for installation. See the example below:
```{r}
install.packages('devtools')
library(devtools)
#if choose to build vignettes
install_github('shuxiaoc/agRank', build_vignettes = TRUE)
```


