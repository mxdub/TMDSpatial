#' Reformating simulate_MC() output to abundances matrix
#'
#' Reformats simulate_MC() output to 3D matrix (species x patch x time)
#'
#' @param sim_output Output from simulate_MC()
#'
#' @return 3D matrix (species x patch x time)
#'
#' @examples
#' # output = simulate_MC(5,5)
#' # abundances = sim_to_matrix(output)
#'
#' @export
#'
sim_to_matrix = function(sim_output){
  x = sim_output$dynamics.df %>% select(-env_niche_breadth, -max_r, -optima, -env)
  x_ = x %>%spread(time, N)
  tmp = array(0, dim = c(length(unique(x$species)),
                         length(unique(x$patch)),
                         length(unique(x$time))))
  for(s in 1:length(unique(x$species))){
    tmp[s,,] = x_ %>% filter(species == s)%>%arrange(patch)%>%select(-patch, -species)%>%as.matrix()
  }
  tmp
}
