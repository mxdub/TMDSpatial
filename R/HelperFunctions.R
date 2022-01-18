#' Reformating simulate_MC() output to abundances matrix
#'
#' Reformats simulate_MC() output to 3D matrix (species x patch x time)
#'
#' @param sim_output Output from simulate_MC()
#'
#' @return 3D abundances matrix (species x patch x time)
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

#' From abundances matrix to occupancies matrix
#'
#' Turns abundances matrix to occupancies matrix  (species x patch x time)
#'
#' @param abundances Abundances matrix
#'
#' @return 3D occupancies matrix (species x patch x time)
#'
#' @examples
#' # occupancies = abund_to_occ(abundances)
#'
#' @export
#'
abund_to_occ = function(abundances){
  abundances[abundances>0]=1
  abundances
}

#' Plots occupancies along time
#'
#' Simply plots proportion of occupied patch by species as function of time
#'
#' @param occupancies Occupancies matrix
#' @param
#'
#' @return none
#'
#' @examples
#' # plots_occupancies(occupancies)
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
plots_occupancies = function(occupancies){
  if( !all(unique(as.vector(occupancies)) <= 1) ){
    stop("Matrix is not a occurences matrix (values > 1)")
  }
  occupancies = apply(occupancies, c(1,3), sum) / dim(occupancies)[2]
  occupancies = occupancies %>% as_tibble() %>% rowid_to_column("species") %>%
    pivot_longer(-species) %>%
    mutate(name = str_replace(name, "V", "")) %>%
    mutate_at(.vars = c("name"), as.numeric)

  ggplot(occupancies%>%filter(name>0), aes(x=name, y=value, color = as.factor(species)))+
    geom_line()+
    scale_y_continuous(limits=c(0,1))+
    labs(x = "Time", y = "Occupancy", color = "Species")+
    theme_bw()
}


