library(STbayes)

# load example bisonr fit object
bisonr_fit = STbayes::bisonr_fit

# convert to STb networks format (long)
# draws defines how many samples to take from posterior.
networks = extract_bisonr_edgeweights(bisonr_fit, draws=100)

# this particular fit was modeling edgeweights centered at 0, but
# edgeweights must be positive for STbayes.
networks$value = networks$value - min(networks$value) #networks can now be used in import_user_STb following normal workflow

#network has 10 individuals, create mock diffusion data
diffusion_data <- data.frame(
    trial = 1,
    id = c(1:10),
    time = sample(1:100, 10, replace = FALSE),
    max_time = 101
)

#create data_list as usual
data_list = import_user_STb(diffusion_data, networks)

# STb detects that you've entered posterior distributions as edgeweights automatically
# it will generate a model wherein each iteration, model will marginalize LL over S=100 draws
model = generate_STb_model(data_list)
#write(model, file="../data/STAN_example_edge_uncertainty.stan")

#the fit will be garbage because it's made up, but works
fit = fit_STb(data_list, model)
STb_summary(fit)
