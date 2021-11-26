# VWMmodels
2021-11-26 by Ru-Yuan Zhang

A set of computational models (mostly probablistic models) of visual working memory



## Installnation

To run this, you need to add below files to your matlab path

* Ru-Yuan Zhang's matlab utility functions, https://github.com/ruyuanzhang/RZutil
* bads optimization toolbox, https://github.com/lacerbi/bads



## Model specifications

Currently we have 7 variants of models

* IL: item-limit model
* MIX: mixture model
* SA: slot-plus-averaging model
* EP: equal presicion model
* VP: variable precision model
* VPcap: variable-precision-plus-capacity model
* cosSA: slots-plus-averaging-with-cosine-precision model

 The total number of model is thus 7 x 2 = 14.



All model related files are located in 'modelfiles' directory. We label the model as 'VWM\_<model><nt>_config' format.  <model> is the model name, <nt> indicates whether to includes non-target swap error. For example,  VWM\_SA\_* indicates the slot-plus-averaging model without non-target swaping error, and VWM\_SAnt\_* indicates vice versa.



Each model has three files:

* VWM\_<model><nt>_config. Configurations for model fiting
* VWM\_<model><nt>_fit. The main function to preprocess data and conduct fit multiple times
* VWM\_<model><nt>_nll.  Negative log likelihood function (this is the key function of a model)



## Data preparation

Suppose you have run an experiment including <nTrials> trials. 

Basically you need to construct a data structure
