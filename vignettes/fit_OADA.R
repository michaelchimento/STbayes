library(NBDA)
library(STbayes)
library(ggplot2)

#### OADA fit ####
# Let's pretend you only have the order of acquisition rather than times.
#load example NBDAdata object from Tutorial 1.1 from Hasenjager et al. 2021
nbda_object <- STbayes::tutorial1_1
model_social <-oadaFit(nbda_object)
data.frame(Variable=model_social@varNames,MLE=model_social@outputPar,SE=model_social@se)

# import into Stbayes
data_list_user = import_NBDA_STb(nbda_object)
# generate the STAN model
model_social = generate_STb_model_OADA(data_list_user, gq=T)
#write(model_obj, file = "../data/STAN_example_OADA.stan")

# fit the model
fit_social = fit_STb(data_list_user, model_social, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99))
# peep the summary, it's similar estimate to NBDA package
STb_summary(fit_social)

#compare with asocial model (s=0)
# generate the STAN model
model_asocial = generate_STb_model_OADA_asocial(data_list_user, gq=T)
#write(model_obj, file = "../data/STAN_example_OADA_asocial.stan")

# fit the model
fit_asocial = fit_STb(data_list_user, model_asocial, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99))

#there are no parameters but we can compare WAIC
extract_WAIC(fit_social) #143
extract_WAIC(fit_asocial) #149
#delta_WAIC ~ 5


