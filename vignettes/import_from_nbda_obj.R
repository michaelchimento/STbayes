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
model_obj = generate_STb_model(data_list, gq=T, est_acqTime = T)

#fit model
fit_social = fit_STb(data_list, model_obj, chains = 2, cores = 2, iter=2000, control = list(adapt_delta=0.99) )

#check summary ~ the same as NBDA estimates. the priors could be adjusted to be less skeptical of the large s value
STb_summary(fit_social, depth=1)

df_acq_times = extract_acqTime(fit_social, data_list)

library(ggplot2)
ggplot(df_acq_times, aes(x = observed_time, y = mean_time)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower_hpd, ymax = upper_hpd), width = 0.2) +
    facet_wrap(~trial, scales = "free_x") +
    labs(
        title = "Estimated Acquisition Times with 95% HPD Intervals",
        x = "Observed time",
        y = "Estimated time"
    ) +
    theme_minimal()

# Compare with asocial model (no s param)
model_obj = generate_STb_asocial_model(data_list)
fit_asocial = fit_STb(data_list, model_obj, chains = 2, cores = 2, iter=2000, control = list(adapt_delta=0.99) )
STb_summary(fit_asocial, depth=1)

AIC(fit_social)



