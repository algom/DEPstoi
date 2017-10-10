.onLoad <- function(libname = find.package("DEPstoi"), pkgname = "DEPstoi"){

  # CRAN Note avoidance
  if(getRversion() >= "2.15.1")
    utils::globalVariables(
      c( # functions.R globalVariables
        "ID", ".", "Peptides", "Protein.group.IDs",
        "Unique..Groups.", "condition", "iBAQ_value",
        "ibaq", "label", "lfc", "mean_ctrl",
        "name", "stoichiometry", "ymax", "ymin"
        )
    )
  invisible()
}
