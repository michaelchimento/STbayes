library(STbayes)

#### import from bisonR ####
# load example bisonr fit object
bisonr_fit = STbayes::bisonr_fit

# convert to STb networks format (long)
# draws defines how many samples to take from posterior.
networks = extract_bisonr_edgeweights(bisonr_fit, draws=100)

# this particular fit was modeling edgeweights centered at 0, but
# edgeweights must be positive for STbayes, so rescale between 0 and 1 here
networks$value = scales::rescale(networks$value) #networks can now be used in import_user_STb following normal workflow

#network has 10 individuals, create mock diffusion data
event_data <- data.frame(
    trial = 1,
    id = c(1:10),
    time = sample(1:101, 10, replace = FALSE),
    t_end = 100
)

#create data_list as usual
data_list = import_user_STb(event_data, networks)

# STb detects that you've entered posterior distributions as edgeweights automatically
# it will generate a model wherein each iteration, model will marginalize LL over S=100 draws
model = generate_STb_model(data_list)
model_oada = generate_STb_model(data_list, data_type = "order") #distributions can also be used in oada
write(model, file="../data/STAN_example_edge_uncertainty.stan")

#the fit will be garbage because it's made up, but works
fit = fit_STb(data_list, model_oada)
STb_summary(fit)

#### Import from STRAND package ####
# load example STRAND results object
# this is not the STRAND fit itself, but created using STRAND::summarize_strand_results(fit)
strand_results = STbayes::strand_results_obj

#the rest is the same as above
networks = extract_strand_edgeweights(strand_results, draws=100)
event_data <- data.frame(
    trial = 1,
    id = c(1:10),
    time = sample(1:101, 10, replace = FALSE),
    t_end = 100
)
data_list = import_user_STb(event_data, networks)
model = generate_STb_model(data_list)
fit = fit_STb(data_list, model)
STb_summary(fit)

