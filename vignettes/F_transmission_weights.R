library(dplyr)
library(STbayes)

event_data = STbayes::event_data
networks = STbayes::edge_list

n_distinct(event_data$time) #45 unique events

t_weights = data.frame(
    trial = 1,
    time = rep(1:45, each=50),
    id = rep(1:50, times=45),
    t_weight = rep(exp(rnorm(50)), times=45)
)

#this generates a single transmission weight per individual
n_distinct(t_weights$t_weight)
#but there's nothing stopping you from inputting dynamic transmission weights.

# format data
data_list_user = import_user_STb(event_data = event_data,
                                 networks=edge_list,
                                 t_weights = t_weights)

#proceed as usual
