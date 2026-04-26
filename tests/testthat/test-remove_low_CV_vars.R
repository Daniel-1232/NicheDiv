test_that("remove.low.CV.vars returns expected list elements", {
  sp1_occ <- data.frame(id = 1:5, var1 = c(1, 2, 3, 4, 5), var2 = c(1, 1, 1, 1, 1))
  sp2_occ <- data.frame(id = 1:5, var1 = c(2, 3, 4, 5, 6), var2 = c(1, 1, 1, 1, 1))
  sp1_bg  <- sp1_occ
  sp2_bg  <- sp2_occ
  out <- remove.low.CV.vars(
    Sp1.occurrence.data = sp1_occ,
    Sp2.occurrence.data = sp2_occ,
    Sp1.background.data = sp1_bg,
    Sp2.background.data = sp2_bg,
    exclude.cols = "id",
    CV.threshold = 0.01,
    verbose = FALSE
  )
  expect_true(is.list(out))
  expect_true(all(c(
    "occurrence_Sp1", "occurrence_Sp2", "background.Sp1", "background.Sp2",
    "dropped_non_numeric", "dropped_NA_only", "dropped_lowCV", "kept_variables"
  ) %in% names(out)))
})
