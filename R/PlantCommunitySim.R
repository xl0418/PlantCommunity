#' @export
PlantCommunitySim <-
  function(L,
           jc_strength,
           sig_disp,
           tend,
           output,
           num_cores = 4,
           torus = TRUE) {
    numCores <- min(num_cores, parallel::detectCores())

    community <-
      matrix(0, nrow = L, ncol = L)  # initialize the grid
    ini.row <- sample(seq(0, L), size = 1)
    ini.col <- sample(seq(0, L), size = 1)
    community[ini.row, ini.col] <- 1
    abundance_vector <- c(L ^ 2 - 1, 1)

    for (t in c(1:tend)) {
      dead.row <- sample(seq(0, L), size = 1)
      dead.col <- sample(seq(0, L), size = 1)
      dead.species <- community[dead.row, dead.col]
      abundance_vector[dead.species + 1] <-
        abundance_vector[dead.species + 1] - 1
      coordinate_bank <-
        cbind(rep(c(1:L), each = L), rep(c(1:L), L))

      colonization_probability_cluster <-
        parallel::mclapply(
          1:L ^ 2,
          colonization_probability,
          coordinate_bank = coordinate_bank,
          coordinate2 = c(dead.row, dead.col),
          jc_strength = jc_strength,
          abundance_vector = abundance_vector,
          L = L,
          sig_disp = sig_disp,
          community = community
        )
    }
  }
