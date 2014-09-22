
write up:
writeup.pdf

Organization of scripts:

main.R : 
Runs processing of data, model fitting, and validation and plots.

functions:
brt_fit.R - functions to fit the brt model to data
loss_fn.R - functions corresponding to loss functions, 
gradients, hessians
predict2.rpart.R - functions to modify terminal node 
estimates from tree and predict new data
process_data_fn.R - functions to process and clean up data
subsample.R - sub sampling function
verif_fn.R - verification function to calculate MSPE for 
certain quantiles of the data

visualizations:
visualizations.R - script to generate figures in paper

