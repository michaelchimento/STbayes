library(STbayes)
#very mock data
event_data <- data.frame(
    trial = rep(1:2, each = 3),
    id = LETTERS[1:6],
    time = c(0, 1, 2, 0, 1, 5),
    t_end = c(3, 3, 3, 4, 4, 4)
)
networks <- data.frame(
    trial = rep(1:2, each = 3),
    from = c("A", "A", "B", "D", "D", "E"),
    to = c("B", "C", "C", "E", "F", "F"),
    kin = c(1, 0, 1, 0, 1, 1),
    inverse_distance = c(0, 1, .5, .25, .1, 0)
)
ILV_c <- data.frame(
    id = LETTERS[1:6],
    age = c(-1, -2, 0, 1, 2, 3), # continuous variables should be normalized
    sex = c(0, 1, 1, 0, 1, 0), # Factor ILVs must be input as numeric
    weight = c(0.5, .25, .3, 0, -.2, -.4)
)
ILV_tv <- data.frame(
    trial = c(rep(1, each = 9),rep(2, each = 9)),
    id = c(rep(LETTERS[1:3], each=3), rep(LETTERS[4:6], each=3)),
    # these times correspond to the inter-acquisition periods
    #e.g. 1 is from [t_0 to t_1), 2 is [t_1 to t_2), 3 = [t_2 to t_3 or t_end] if censored inds. present)
    time = c(rep(1:3, times = 3), rep(1:3, times=3)),
    #ensure the variable is summarizing these inter-acquisition time periods
    dist_from_resource = rnorm(18)
)

#### explicitly setting which variables are additive and multiplicative ####
STb_data <- import_user_STb(
    event_data = event_data,
    networks = networks,
    ILV_c = ILV_c,
    ILV_tv = ILV_tv,
    ILVi = c("age", "dist_from_resource"), # estimate effects of constant ILV 'age' and time-varying ILV 'dist_from_resource' on asocial learning rates
    ILVs = c("sex"), # Use only 'sex' for social learning
    ILVm = c("weight") # Use weight for multiplicative effect on asocial and social learning
)

model_obj = generate_STb_model(STb_data)
write(model_obj, file = "../data/STAN_example_ILV.stan")

fit = fit_STb(STb_data, model_obj)

#### if not explicitly set, these variables will be additive unconstrained ####
STb_data <- import_user_STb(
    event_data = event_data,
    networks = networks,
    ILV_c = ILV_c,
    ILV_tv = ILV_tv
)
model_obj = generate_STb_model(STb_data)
fit = fit_STb(STb_data, model_obj)

#### set varying effects ####
model_obj = generate_STb_model(STb_data, veff_ID = c("weight"))
fit = fit_STb(STb_data, model_obj)

