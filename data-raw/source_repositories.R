library(ctv)
library(purrr)
library(dplyr)

ctv_pkgs <- read.ctv("Epidemiology.md") |>
  pluck("packagelist", "name")

source_from_desc <- tools::CRAN_package_db() |>
  filter(Package %in% ctv_pkgs) |> 
  transmute(
    package = Package,
    github_repo = dplyr::coalesce(
      stringr::str_match(URL, "https://github.com/([^/]+/[^/,>]+)")[, 2],
      stringr::str_match(BugReports, "https://github.com/([^/]+/[^/,]+)/issues")[, 2]
    )
  )

source_from_manual <- tribble(
  ~ package, ~ github_repo,
  "coarseDataTools", "nickreich/coarseDataTools",
  "colmozzie", "thiyangt/colmozzie",
  "covid19france", "Covid19R/covid19france",
  "covid19us", "aedobbyn/covid19us",
  "epiR", "mstevenson888/epiR",
  "EpiReport", "EU-ECDC/EpiReport",
  "episensr", "dhaine/episensr",
  "etm", "mclements/etm",
  "mediation", "kosukeimai/mediation",
  "outbreaker2", "reconhub/outbreaker2",
  "riskCommunicator", "jgrembi/riskCommunicator",
  "shinySIR", "SineadMorris/shinySIR",
  "TransPhylo", "xavierdidelot/TransPhylo",
  "tsiR", "adbecker/tsiR",
  "excessmort", "rafalab/excessmort",
  "EpiForsk", "Laksafoss/EpiForsk",
  "tidyrates", "rfsaldanha/tidyrates",
  "genSEIR", "selcukorkmaz/genSEIR",
  # "dataquieR", "https://gitlab.com/libreumg/dataquier",
  "linelistBayes", "cmilando/linelistBayes",
  "rplanes", "signaturescience/rplanes",
  "EVI", "ku-awdc/EVI"
)

anti_join(source_from_desc, source_from_manual, by = "package") |>
  full_join(source_from_manual) |>
  dplyr::arrange(package) |>
  dplyr::filter(
    package != "epimdr" # declared repo does not contain source code
  ) |>
  write.csv("data/source_repositories.csv", row.names = FALSE)
