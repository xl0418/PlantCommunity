test_that("self_excluded", {
  expect_identical(colonization_probability(bank_index = 1,
                                            coordinate_bank = t(as.matrix(c(1,1))),
                                            coordinate2 = c(1, 1),
                                            jc_strength = 0.5,
                                            abundance_vector = c(100),
                                            L = 10,
                                            sig_disp = 1,
                                            community = matrix(0, 10, 10),
                                            torus = TRUE), 0)
})


