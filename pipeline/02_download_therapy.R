# Download per-patient drug and radiation tables (BCR Biotab) from GDC for all
# 33 TCGA projects, and stack them into two raw tables.
# Usage: Rscript pipeline/02_download_therapy.R   (from repo root)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
})

out_drug <- "data/raw/therapy_drug.tsv"
out_rad  <- "data/raw/therapy_radiation.tsv"
if (file.exists(out_drug) && file.exists(out_rad)) { message("skip therapy"); quit() }

gdc_dir <- "data/raw/GDCdata"
projects <- grep("^TCGA-", getGDCprojects()$project_id, value = TRUE)

na_tokens <- c("", "NA", "[Not Available]", "[Not Applicable]", "[Unknown]",
               "[Not Evaluated]", "[Discrepancy]", "[Completed]", "[Not Reported]")

drug <- list(); rad <- list()
for (p in sort(projects)) {
  message("== ", p)
  q <- tryCatch(GDCquery(project = p, data.category = "Clinical",
                         data.type = "Clinical Supplement", data.format = "BCR Biotab"),
                error = function(e) NULL)
  if (is.null(q)) next
  files <- getResults(q)$file_name
  want <- grep("clinical_(drug|radiation)_", files, value = TRUE)
  if (!length(want)) next
  q$results[[1]] <- q$results[[1]][q$results[[1]]$file_name %in% want, ]
  GDCdownload(q, directory = gdc_dir, method = "api", files.per.chunk = 10)
  tabs <- GDCprepare(q, directory = gdc_dir)
  for (nm in names(tabs)) {
    # Biotab files have 2 extra header rows (CDE names / ids) at the top
    d <- as.data.table(tabs[[nm]])[-(1:2)]
    for (col in names(d)) set(d, which(d[[col]] %in% na_tokens), col, NA)
    d[, `:=`(project = sub("^TCGA-", "", p), source_table = nm)]
    if (grepl("_drug_", nm)) drug[[nm]] <- d else rad[[nm]] <- d
  }
}

fwrite(rbindlist(drug, fill = TRUE), out_drug, sep = "\t")
fwrite(rbindlist(rad,  fill = TRUE), out_rad,  sep = "\t")
message("drug rows: ", sum(sapply(drug, nrow)), "  radiation rows: ", sum(sapply(rad, nrow)))
