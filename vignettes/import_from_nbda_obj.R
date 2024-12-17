library(STbayes)
library(NBDA)
#load example NBDAdata object from Tutorial 4.1 from Hasenjager et al. 2021
nbdaData_cTADA <- STbayes::nbdaData_cTADA
#Fit cTADA model using tadaFit:
model_constant<-NBDA::tadaFit(nbdaData_cTADA)
data.frame(Variable=model_constant@varNames,MLE=model_constant@outputPar,SE=model_constant@se)

#import into STbayes
data_list = import_NBDA_STb(nbdaData_cTADA)

#generate STAN model from input data
model_obj = generate_STb_model(data_list)

#fit model
fit = fit_STb(data_list, model_obj, chains = 2, cores = 2, iter=2000, control = list(adapt_delta=0.99) )

#check summary ~ the same as NBDA estimates. the priors could be adjusted to be less skeptical of the large s value
STb_summary(fit, depth=2)
