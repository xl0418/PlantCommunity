#' @export
colonization_probability <-
  function(bank_index,
           coordinate_bank,
           coordinate2,
           jc_strength,
           abundance_vector,
           L,
           sig_disp,
           community,
           torus = TRUE) {
    coordinate1 <- coordinate_bank[bank_index,]
    if (all(coordinate1 == coordinate2)) return(0)
    else {
      if (torus) {
        x.dist <- min((coordinate1[1] - coordinate2[1]) ^ 2,
                      (
                        L - max(coordinate1[1], coordinate2[1]) + min(coordinate1[1], coordinate2[1])
                      ) ^ 2)
        y.dist <- min((coordinate1[2] - coordinate2[2]) ^ 2,
                      (
                        L - max(coordinate1[2], coordinate2[2]) + min(coordinate1[2], coordinate2[2])
                      ) ^ 2)
        distance <- sqrt(x.dist + y.dist)
      } else {
        distance <- sqrt(sum(coordinate1 - coordinate2) ^ 2)
      }
      abundance <- abundance_vector[community[coordinate1[1], coordinate1[2]] + 1]
      probability <-
        (1 - jc_strength * abundance / L ^ 2) * (1 / sqrt(2 * pi) / sig_disp * exp(-distance ^
                                                                                     2 / 2 / sig_disp ^ 2))
      return(probability)
    }

  }
