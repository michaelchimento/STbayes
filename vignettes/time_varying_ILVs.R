library(STbayes)

#### load tutorial object: cTADA with timevarying effects ####
nbda_object <- STbayes::tutorial4_2

# creates two TV ILVs, treatment and maleTV.
# despite sex being a constant, once imported into NBDA alongside of a TV variable, it is treated as TV as well
# thus, when importing from an NBDA object, it will get treated as TV in STbayes
# However, if you import your own data using import_user_STb, this is not necessary
data_list <- import_NBDA_STb(nbda_object = nbda_object)

model_tada <- generate_STb_model(STb_data = data_list)
model_tada
STb_fit_tada <- fit_STb(data_list,
  model_tada,
  chains = 5,
  cores = 5,
  iter = 2000,
  control = list(adapt_delta = 0.99)
)

STb_summary(STb_fit_tada)

# compare summary to nbda
# Fit the model
nbda_fit <- tadaFit(nbda_object)
# Display the output
data.frame(Variable = nbda_fit@varNames, MLE = nbda_fit@outputPar, SE = nbda_fit@se)

#### load tutorial object: OADA with timevarying effects ####
nbda_object <- STbayes::tutorial2_2

# creates two TV ILVs, treatment and maleTV.
# despite sex being a constant, once imported into NBDA alongside of a TV variable, it is treated as TV as well
# thus, when importing from an NBDA object, it will get treated as TV in STbayes
# However, if you import your own data using import_user_STb, this is not necessary
data_list <- import_NBDA_STb(nbda_object = nbda_object)

model_oada <- generate_STb_model_OADA(STb_data = data_list)
model_oada
STb_fit_oada <- fit_STb(data_list,
  model_oada,
  chains = 5,
  cores = 5,
  iter = 2000,
  control = list(adapt_delta = 0.99)
)

STb_summary(STb_fit_oada)

# compare summary to nbda
# Fit the model
nbda_fit <- oadaFit(nbda_object)
# Display the output
data.frame(Variable = nbda_fit@varNames, MLE = nbda_fit@outputPar, SE = nbda_fit@se)
