
### Gradient Boosting Trees with Squared Error Epsilon Insensitive Loss

The write-up is located [here](http://kenlau177.github.io/cpsc546/writeup.pdf)

#### File Descriptions:

*main.R:* Data processing, model fitting, and verification.

*functions:*
- brt_fit.R : functions to fit the brt model to data
- loss_fn.R : functions corresponding to loss functions, gradients, hessians
- predict2.rpart.R : functions to modify terminal node 
- estimates from tree and predict new data
- process_data_fn.R : functions to process and clean up data
- subsample.R : sub sampling function
- verif_fn.R : verification function to calculate MSPE for certain quantiles of the data

*visualizations:*
- visualizations.R : script to generate figures in paper

