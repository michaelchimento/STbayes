library(STbayes)

#### import from bisonR ####
# load example bisonr fit object
bisonr_fit = STbayes::bisonr_fit

#network has 10 individuals, create mock diffusion data
event_data <- data.frame(
    trial = 1,
    id = c(1:10),
    time = sample(1:101, 10, replace = FALSE),
    t_end = 100
)

#create data_list as usual
data_list = import_user_STb(event_data, networks=bisonr_fit)
model = generate_STb_model(data_list)
model_oada = generate_STb_model(data_list, data_type = "order") #distributions can also be used in oada
write(model, file="../inst/extdata/STAN_example_edge_uncertainty.stan")

#the fit will be garbage because it's made up, but works
fit = fit_STb(data_list, model)
STb_summary(fit)

#### Import from STRAND package ####
# load example STRAND results object
# this is not the STRAND fit itself, but created using STRAND::summarize_strand_results(fit)
strand_results = STbayes::strand_results_obj
data_list = import_user_STb(event_data, networks=strand_results)
model = generate_STb_model(data_list)
fit = fit_STb(data_list, model)
STb_summary(fit)

