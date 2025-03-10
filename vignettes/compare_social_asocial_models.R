library(STbayes)
library(ggplot2)

diffusion_data = STbayes::diffusion_data
edge_list = STbayes::edge_list

# format data
data_list_user = import_user_STb(diffusion_data, edge_list)

# reusable function to generate and fit a model
generate_and_fit_model <- function(data, model_type, chains = 5, cores = 5, iter = 2000, control = list(adapt_delta = 0.99)) {
    model = generate_STb_model(data, est_acqTime = TRUE, model_type=model_type)
    fit = fit_STb(data, model, chains, cores, iter, control)
    return(fit)
}

# generate, fit, and summarize models
full_fit = generate_and_fit_model(data_list_user, "full")
asocial_fit = generate_and_fit_model(data_list_user, "asocial")

#check estimates
STb_summary(full_fit, digits = 4)
STb_summary(asocial_fit, digits = 4)

# extract WAIC and labels
WAIC_full = extract_WAIC(full_fit)
WAIC_null = extract_WAIC(asocial_fit)
label_full = paste0("WAIC=", round(WAIC_full[[5]]), "+/-", round(WAIC_full[[6]]))
label_null = paste0("WAIC=", round(WAIC_null[[5]]), "+/-", round(WAIC_null[[6]]))


# reusable function for plotting
plot_acq_time <- function(fit, data, title, label) {
    acqdata = extract_acqTime(fit, data)
    p = ggplot(acqdata, aes(x = observed_time, y = mean_time)) +
        annotate("text", x = 100, y = 350, label = label) +
        geom_segment(
            aes(x = observed_time, xend = observed_time, y = mean_time, yend = observed_time),
            color = "red",
            alpha = 0.2
        ) +
        geom_point(alpha = 0.6, size = 2) +
        geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
        labs(
            title = title,
            x = "Observed time",
            y = "Estimated time"
        ) +
        theme_minimal()
    return(p)
}


# plot estimated times
p1 = plot_acq_time(asocial_fit, data_list_user, "Asocial (null) model estimates", label_null)
p2 = plot_acq_time(full_fit, data_list_user, "Full model estimates", label_full)

library(ggpubr)
ggarrange(p1, p2)
ggsave("../data/compare_social_asocial.png", width=10, height=5, units="cm", scale=1.5)
