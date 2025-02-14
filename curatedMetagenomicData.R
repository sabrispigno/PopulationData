library(curatedMetagenomicData)
curatedMetagenomicData("AsnicarF_20.+")
curatedMetagenomicData("AsnicarF_2017.relative_abundance", dryrun = FALSE, rownames = "short")

#'
#' @param dataType (default: "relative_abundance")
#' Data type, passed on to \code{\link[curatedMetagenomicData]{returnSamples}}.
#' @param counts (Default: FALSE)
#' Whether to convert to count-like data by multiplying through by read depth.
#' Passed on to \code{\link[curatedMetagenomicData]{returnSamples}}.
#'
#' @return
#' The SummarizedExperiment or TreeSummarizedExperiment containing all of cMD for the data type requested.
#' Calling this function has the side-effect of also writing two csv files, one for the assay data of this
#' object and one for the colData.
#'
#' @details This function also removes control samples from the NielsenHB_2014 study,
#' which are duplicated in the LeChatelierE_2013 study. The cMD version number is put in the filename.
#'
#' @export
#'
#' @importFrom curatedMetagenomicData returnSamples
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#'    tse <- dataDump()
#'    tse
#'    dir(pattern = "\\.csv")
#' }
dataDump <- function(dataType = "relative_abundance", counts = FALSE) {
  dups <-
    curatedMetagenomicData::sampleMetadata$study_name == "NielsenHB_2014" &
    grepl(pattern = "^MH",
          x = curatedMetagenomicData::sampleMetadata$sample_id)
  
  x <-
    returnSamples(curatedMetagenomicData::sampleMetadata[!dups, ],
                  dataType = dataType,
                  counts = counts)
  
  vn <- as.character(packageVersion("curatedMetagenomicData"))  # Versione del pacchetto
  metadatafile <- paste0("metadata-", vn, ".csv")
  datafile <- paste0(dataType, "-", vn, ".csv")
  
  utils::write.csv(assay(x), file = datafile)
  utils::write.csv(colData(x), file = metadatafile)
  
  return(x)
}

results <- dataDump()
list.files(pattern = "\\.csv")

## ----message = FALSE----------------------------------------------------------
library(dplyr)
library(DT)

## ----collapse = TRUE----------------------------------------------------------
sampleMetadata |>
  filter(study_name == "AsnicarF_2017") |>
  select(where(~ !any(is.na(.x)))) |>
  slice(1:10) |>
  select(1:10) |>
  datatable(options = list(dom = "t"), extensions = "Responsive")