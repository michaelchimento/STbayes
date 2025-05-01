library(STbayes)
library(ggplot2)

event_data = STbayes::event_data
edge_list = STbayes::edge_list

# format data
data_list_user = import_user_STb(event_data, edge_list)

# generate, fit, and summarize models
model = generate_STb_model(data_list_user, est_acqTime = TRUE, model_type="full")
cat(model)
full_fit = fit_STb(data_list_user,
                   model,
                   chains = 5,
                   cores = 5,
                   parallel_chains=5,
                   iter = 2000,
                   control = list(adapt_delta = 0.99))
STb_summary(full_fit, depth=2)
model = generate_STb_model(data_list_user, est_acqTime = TRUE, model_type="asocial")
asocial_fit = fit_STb(data_list_user,
                      model,
                      chains = 5,
                      cores = 5,
                      parallel_chains=5,
                      iter = 2000,
                      control = list(adapt_delta = 0.99))

STb_summary(full_fit, digits = 4)
STb_summary(asocial_fit, digits = 4)

# convenient workflow for model comparison.
# You can either use WAIC
waic_output = STb_compare(full_fit, asocial_fit, method="waic")
waic_output$loo_objects
waic_output$comparison

# Or LOO-PSIS
loo_output = STb_compare(full_fit, asocial_fit, method="loo-psis")
loo_output$comparison

# LOO-PSIS gives pareto-k diagnostics
pareto_df = as.data.frame(loo_output$pareto_diagnostics)
ggplot(pareto_df, aes(x=observation, y=pareto_k, color=model))+
    geom_point() +
    scale_color_viridis_d(begin=0.2, end=0.7)+
    geom_hline(yintercept = 0.7, linetype="dashed", color="orange")+
    geom_hline(yintercept = 1, linetype="dashed", color="red")+
    labs(x="observation", y="pareto k value", title="Pareto-k diagnostics")+
    theme_minimal()

# if values fall above 0.7, consider using k-fold validation. values here seem fine.
# a note that kfold validation means manually making k folds of data, fitting model to each fold and
# computing kfold elpd.

# Extract comparison results for plotting
comparison_df <- as.data.frame(loo_output$comparison)
comparison_df$model <- rownames(comparison_df)

# Plot comparison
# Best predictive model will be at top. If SE of difference crosses zero,
# you should be less certain that the best fit model actually provides a better fit
ggplot(comparison_df, aes(x = reorder(model, elpd_diff), y = elpd_diff)) +
    geom_point(size = 3) +  #mean elpd difference
    geom_errorbar(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff), width = 0.2) + #SE of elpd diff
    coord_flip() +
    theme_minimal() +
    labs(x = "Model", y = "ELPD Difference", title = "Model Comparison") +
    theme(axis.text = element_text(size = 12))

# extract Loo-ic and labels
looic_full <- loo_output$loo_objects$full_fit$estimates["looic",]
looic_asocial <- loo_output$loo_objects$asocial_fit$estimates["looic",]
label_full = paste0("LOO-PSIS=", round(looic_full["Estimate"]), "+/-", round(looic_full["SE"]))
label_null = paste0("LOO-PSIS=", round(looic_asocial["Estimate"]), "+/-", round(looic_asocial["SE"]))

# reusable function for plotting
plot_acq_time <- function(fit, data, title, label) {
    acqdata = extract_acqTime(fit, data)
    p = ggplot(acqdata, aes(x = observed_time, y = median_time)) +
        annotate("text", x = 100, y = 350, label = label) +
        geom_segment(
            aes(x = observed_time, xend = observed_time, y = median_time, yend = observed_time),
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
p1 = plot_acq_time(asocial_fit, data_list_user, "Asocial model estimates", label_null)
p2 = plot_acq_time(full_fit, data_list_user, "Full model estimates", label_full)

library(ggpubr)
ggarrange(p1, p2)
ggsave("../docs/compare_social_asocial.png", width=10, height=5, units="cm", scale=2.5)
