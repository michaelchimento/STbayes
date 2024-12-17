data_list_nbda = import_NBDA_STb(d, network_names = c("assoc"))
data_list_user = import_user_STb(diffusion_data, edge_list)

for (name in names(data_list_user)){
    print(name)
    print(sum(data_list_user[[name]] != data_list_nbda[[name]]))
}
