#' Plant Community Assembly Simulation (PCAS)
#'
#' Simulate a plant community assembly
#'
#' @param L The grid size
#' @param jc_strength The strength of the Janzen-Connell (negative density dependence) effect. When 0, it represents a neutral model
#' @param sig_disp The dispersal width
#' @param tend The simulating time
#' @param output The directory and file name of the output
#' @param num_cores The number of the cores for parallel computing
#' @param torus With bundary or not
#'
#' @return None
#'
#' @export
PlantCommunitySim <-
  function(L,
           jc_strength,
           sig_disp,
           tend,
           output,
           speciation_rate,
           num_cores = 4,
           torus = TRUE) {
    numCores <- min(num_cores, parallel::detectCores())

    community <-
      matrix(0, nrow = L, ncol = L)  # initialize the grid
    ini.row <- sample(seq(0, L), size = 1)
    ini.col <- sample(seq(0, L), size = 1)
    community[ini.row, ini.col] <- 1
    abundance_vector <- c(L ^ 2 - 1, 1)
    coordinate_bank <-
      cbind(rep(c(1:L), L), rep(c(1:L), each = L))

    for (t in c(1:tend)) {
      # sample the dead individual
      dead.row <- sample(seq(0, L), size = 1)
      dead.col <- sample(seq(0, L), size = 1)
      dead.species <- community[dead.row, dead.col]
      # reduce the corresponding abundance
      abundance_vector[dead.species + 1] <-
        abundance_vector[dead.species + 1] - 1
      # sample the birth or speciation event
      sample_birth_speciation_event <-
        sample(c(1, 2),
               size = 1,
               prob = c(1 - speciation_rate, speciation_rate))
      # the birth event
      if (sample_birth_speciation_event == 1) {
        # calculate the colonization probability of individuals
        colonization_probability_cluster <-
          parallel::mclapply(
            1:L ^ 2,
            PlantCommunity::colonization_probability,
            coordinate_bank = coordinate_bank,
            coordinate2 = c(dead.row, dead.col),
            jc_strength = jc_strength,
            abundance_vector = abundance_vector,
            L = L,
            sig_disp = sig_disp,
            community = community
          )
        colonization_probability_vector <-
          unlist(colonization_probability_cluster)
        # sample which individual gives the birth to the vacant site
        sample.parent.to.give.birth <-
          sample(community[1:L ^ 2], size = 1, prob = colonization_probability_vector)
        community[dead.row, dead.col] <- sample.parent.to.give.birth
        # update the corresponding abundance
        abundance_vector[sample.parent.to.give.birth + 1] <-
          abundance_vector[sample.parent.to.give.birth + 1] + 1

      } else{
        community[dead.row, dead.col] <- length(abundance_vector)
        abundance_vector <- c(abundance_vector, 1)
      }



    }
  }
