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
# generate the STAN model, make sure to specify "order" as datatype rather than "time" (default behavior)
STb_data=data_list_user

#default uniform prior results in divergences, so set a slightly more informative prior
model_social = generate_STb_model(data_list_user, data_type="order", gq=T, model_type="full", priors=list(log_s="normal(1,3)"))
write(model_social, file = "../data/STAN_example_OADA.stan")
model_social
# fit the model
fit_social = fit_STb(data_list_user, model_social, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99))
# peep the summary, it's similar estimate to NBDA package
STb_summary(fit_social)

#compare with asocial model (s=0)
# generate the STAN model
model_asocial = generate_STb_model(data_list_user, data_type="order", gq=T, model_type="asocial")
write(model_asocial, file = "../data/STAN_example_OADA_asocial.stan")

# fit the model
fit_asocial = fit_STb(data_list_user, model_asocial, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99), algorithm = "Fixed_param")

#there are no parameters in asocial model but we can compare elpd metrics
loo_output = STb_compare(fit_social, fit_asocial, method="waic")
loo_output$comparison

