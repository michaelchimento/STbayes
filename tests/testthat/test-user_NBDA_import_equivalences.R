test_that("import_user and import_nbda result in element-wise equivalent data lists", {
    tie_vec = STbayes::event_data %>%
        dplyr::arrange(time) %>%
        dplyr::group_by(time, .drop = T) %>%
        dplyr::mutate(tie=ifelse(dplyr::n()>1,1,0)) %>%
        dplyr::pull(tie)
    seed_vec = STbayes::event_data %>%
        dplyr::arrange(time) %>%
        dplyr::group_by(time, .drop = T) %>%
        dplyr::mutate(tie=ifelse(dplyr::n()>1,1,0),
               seed=ifelse(time==0,1,0)) %>%
        dplyr::pull(seed)
    g = igraph::graph_from_edgelist(as.matrix(STbayes::edge_list[1:2]), directed = FALSE)
    adj_matrix = igraph::as_adjacency_matrix(g,attr=NULL, sparse = FALSE)
    dim(adj_matrix) = c(50,50,1)
    nbdaData_object = NBDA::nbdaData(label="sim_data",
                 assMatrix = adj_matrix,
                 orderAcq = STbayes::event_data$id,
                 timeAcq = STbayes::event_data$time,
                 endTime = 402,
                 ties = tie_vec,
                 demons = seed_vec)
    data_list_nbda = import_NBDA_STb(nbdaData_object, network_names = c("assoc"))
    data_list_user = import_user_STb(STbayes::event_data, STbayes::edge_list)
    mismatches <- which(sapply(seq_along(data_list_user), function(i) !identical(data_list_nbda[[i]], data_list_user[[i]])))
    expect(length(mismatches) == 0, paste("Mismatches found at indices:", paste(mismatches, collapse = ", ")))
})
