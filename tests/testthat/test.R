context("Functions")

test_that("merge_ibaq throws error without valid input", {
  expect_error(merge_ibaq("test_data", test_pep))
  expect_error(merge_ibaq(test_data, "test_pep"))
  expect_error(merge_ibaq(test_data[,-(30)], test_pep))
  expect_error(merge_ibaq(test_data[,-(31)], test_pep))
  expect_error(merge_ibaq(test_data[,-(15:20)], test_pep))
  expect_error(merge_ibaq(test_data, test_pep[,-(14)]))
  expect_error(merge_ibaq(test_data, test_pep[,-(6)]))
})

test_that("merge_ibaq returns a data.frame", {
  expect_is(merge_ibaq(test_data, test_pep), "data.frame")
  expect_is(merge_ibaq(tibble::as_tibble(test_data), test_pep), "data.frame")
  expect_is(merge_ibaq(test_data, tibble::as_tibble(test_pep)), "data.frame")
})

test_that("merge_ibaq returns an object with the rigth dimensions and columns", {
  result <- merge_ibaq(test_data, test_pep)
  expect_equal(grep("iBAQ", colnames(result)), 4:9)
  expect_equal(dim(result), c(359,10))
})

test_that("get_stoichiometry throws error without valid input", {
  expect_error(get_stoichiometry("test_dep", test_ibaq, "GFP_vs_WT", "Rbbp4", 1))
  expect_error(get_stoichiometry(test_dep, "test_ibaq", "GFP_vs_WT", "Rbbp4", 1))
  expect_error(get_stoichiometry(test_dep, test_ibaq, GFP_vs_WT, "Rbbp4", 1))
  expect_error(get_stoichiometry(test_dep, test_ibaq, "GFP_vs_WT", Rbbp4, 1))
  expect_error(get_stoichiometry(test_dep, test_ibaq, "GFP_vs_WT", "Rbbp4", "1"))
  expect_error(get_stoichiometry(test_dep, test_ibaq[,-(2:3)], "GFP_vs_WT", "Rbbp4", 1))
  expect_error(get_stoichiometry(test_dep, test_ibaq[,-(4:9)], "GFP_vs_WT", "Rbbp4", 1))

  test_ibaq_sign_error1 <- test_dep
  SummarizedExperiment::rowData(test_ibaq_sign_error1) <- SummarizedExperiment::rowData(test_ibaq_sign_error1)[,-(1)]
  expect_error(get_stoichiometry(test_ibaq_sign_error1, test_ibaq, "GFP_vs_WT", "Rbbp4", 1))

  test_ibaq_sign_error2 <- test_dep
  SummarizedExperiment::rowData(test_ibaq_sign_error2) <- SummarizedExperiment::rowData(test_ibaq_sign_error2)[,-c(31,32)]
  expect_error(get_stoichiometry(test_ibaq_sign_error2, test_ibaq, "GFP_vs_WT", "Rbbp4", 1))

  test_ibaq_sign_error3 <- test_dep
  SummarizedExperiment::rowData(test_ibaq_sign_error3) <- SummarizedExperiment::rowData(test_ibaq_sign_error3)[,-(35)]
  expect_error(get_stoichiometry(test_ibaq_sign_error3, test_ibaq, "GFP_vs_WT", "Rbbp4", 1))

})

test_that("get_stoichiometry returns a data.frame", {
  expect_is(get_stoichiometry(test_dep, test_ibaq, "GFP_vs_WT", "Rbbp4", 1), "data.frame")
})

test_that("get_stoichiometry returns an object with the rigth dimensions and columns", {
  result <- get_stoichiometry(test_dep, test_ibaq, "GFP_vs_WT", "Rbbp4", 1)
  expect_equal(result$stoichiometry[result$name == "Rbbp4"], 1)
  expect_equal(dim(result), c(5,4))
})

test_that("plot_stoichiometry throws error without valid input", {
  expect_error(plot_stoichiometry("test_stoi", 0.001, NULL))
  expect_error(plot_stoichiometry(test_stoi, "0.001", NULL))
  expect_error(plot_stoichiometry(test_stoi, 0.001, "0.05"))
  expect_error(plot_stoichiometry(test_stoi[,-(1)], 0.001, NULL))
  expect_error(plot_stoichiometry(test_stoi[,-(2)], 0.001, NULL))
  expect_error(plot_stoichiometry(test_stoi[,-(3)], 0.001, NULL))
  expect_error(plot_stoichiometry(test_stoi[,-(4)], 0.001, NULL))
})

test_that("plot_stoichiometry returns a ggplot object", {
  expect_is(plot_stoichiometry(test_stoi, 0.001), "ggplot")
  expect_is(plot_stoichiometry(test_stoi, 0.001, 0.5), "ggplot")
})

test_that("plot_ibaq throws error without valid input", {
  expect_error(plot_ibaq("test_dep", "GFP_vs_WT", 3))
  expect_error(plot_ibaq(test_dep, GFP_vs_WT, 3))
  expect_error(plot_ibaq(test_dep, "GFP_vs_WT", "3"))
  expect_error(plot_ibaq(test_dep, "test", 3))

  test_ibaq_sign_error1 <- test_dep
  SummarizedExperiment::rowData(test_ibaq_sign_error1) <- SummarizedExperiment::rowData(test_ibaq_sign_error1)[,-(1)]
  expect_error(plot_ibaq(test_ibaq_sign_error1, "GFP_vs_WT", 3))

  test_ibaq_sign_error2 <- test_dep
  SummarizedExperiment::rowData(test_ibaq_sign_error2) <- SummarizedExperiment::rowData(test_ibaq_sign_error2)[,-c(31,32)]
  expect_error(plot_ibaq(test_ibaq_sign_error2, "GFP_vs_WT", 3))

  test_ibaq_sign_error3 <- test_dep
  SummarizedExperiment::rowData(test_ibaq_sign_error3) <- SummarizedExperiment::rowData(test_ibaq_sign_error3)[,-(16:21)]
  expect_error(plot_ibaq(test_ibaq_sign_error3, "GFP_vs_WT", 3))

  test_ibaq_sign_error4 <- test_dep
  SummarizedExperiment::rowData(test_ibaq_sign_error4) <- SummarizedExperiment::rowData(test_ibaq_sign_error4)[,-(35)]
  expect_error(plot_ibaq(test_ibaq_sign_error4, "GFP_vs_WT", 3))

})

test_that("plot_ibaq returns a ggplot object", {
  expect_is(plot_ibaq(test_dep, "GFP_vs_WT", 3), "ggplot")
})

test_that("iBAQ throws error without valid input", {
  expect_error(iBAQ("test_result", test_pep, "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(test_result, "test_pep", "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(test_result, test_pep, GFP_vs_WT, "Rbbp4"))
  expect_error(iBAQ(test_result, test_pep, "GFP_vs_WT", Rbbp4))
  expect_error(iBAQ(test_result[-(1)], test_pep, "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(test_result, test_pep[,-(6)], "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(test_result, test_pep[,-(14)], "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(test_result, test_pep, "test", "Rbbp4"))
  expect_error(iBAQ(test_result, test_pep, "GFP_vs_WT", "test"))

  result_error <- test_result
  result_error$data <- result_error$data[,-(16:21)]
  expect_error(iBAQ(result_error, test_pep, "GFP_vs_WT", "Rbbp4"))

  result_error2 <- test_result
  SummarizedExperiment::rowData(result_error2$dep) <- SummarizedExperiment::rowData(result_error2$dep)[,-(1)]
  expect_error(iBAQ(result_error2, test_pep, "GFP_vs_WT", "Rbbp4"))

  result_error3 <- test_result
  SummarizedExperiment::rowData(result_error3$dep) <- SummarizedExperiment::rowData(result_error3$dep)[,-c(31,35)]
  expect_error(iBAQ(result_error3, test_pep, "GFP_vs_WT", "Rbbp4"))
})

test_that("iBAQ returns a data.frame", {
  expect_is(iBAQ(test_result, test_pep, "GFP_vs_WT", "Rbbp4", level = 1), "data.frame")
  expect_is(iBAQ(test_result, test_pep, "GFP_vs_WT", "Rbbp4", level = 1L), "data.frame")
  expect_is(iBAQ(test_result, tibble::as_tibble(test_pep), "GFP_vs_WT", "Rbbp4", level = 1), "data.frame")
})

test_that("run_app throws error without valid input", {
	expect_error(run_app("test"))
})


