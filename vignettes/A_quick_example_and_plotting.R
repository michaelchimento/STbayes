library(STbayes)

#import some simulated example data
event_data = STbayes::event_data
edge_list = STbayes::edge_list

#format data
data_list_user = import_user_STb(event_data, edge_list)

#generate STAN model from input data
model_obj = generate_STb_model(data_list_user, gq=T)
cat(model_obj)
# fit model
full_fit = fit_STb(data_list_user,
                   model_obj,
                   parallel_chains=5,
                   chains = 5,
                   cores = 5,
                   iter= 5000)

# check estimates
STb_summary(full_fit, digits=3)
#get estimated times
acqdata = extract_acqTime(full_fit, data_list_user)

#plot estimated times versus observed times w/ residuals
ggplot(acqdata, aes(x = observed_time, y = median_time)) +
    geom_segment(
        aes(x = observed_time, xend = observed_time, y = median_time, yend = observed_time), # connect predicted to slope line
        color = "red",
        alpha = 0.2
    ) +
    geom_point(alpha = 0.6, size=2) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    facet_wrap(~trial, scales = "free_x") +
    labs(
        title = "Estimated acquisition time with residuals",
        x = "Observed time",
        y = "Estimated time"
    ) +
    theme_minimal()


#plot estimated times versus observed times w/ HPD
ggplot(acqdata, aes(x = observed_time, y = mean_time)) +
    geom_abline(slope = 1, intercept = 0,             # Line y = x
                color = "red", linetype = "dotted", size = 1) +
    geom_pointrange(aes(ymin = lower_hpd, ymax = upper_hpd), size=.8) +
    facet_wrap(~trial, scales = "free_x") +
    labs(
        title = "Estimated acquisition time with 95% HPD Intervals",
        x = "Observed time",
        y = "Estimated time"
    ) +
    theme_minimal()

ggsave("../docs/estimates_residuals.png", width=8, height=8, units="cm", scale=1.5)


