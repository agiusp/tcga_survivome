# Download raw source files listed in sources.csv, plus cBioPortal PanCancer Atlas
# clinical attributes (MSI, FGA, hypoxia, ancestry). Re-runs skip files already present.
# Usage: Rscript pipeline/01_download.R   (from repo root)

suppressPackageStartupMessages({
  library(data.table)
  library(httr2)
})

raw_dir <- "data/raw"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# ---- GDC / Xena flat files -------------------------------------------------
src <- fread("pipeline/sources.csv")
for (i in seq_len(nrow(src))) {
  dest <- file.path(raw_dir, src$filename[i])
  if (file.exists(dest) && file.size(dest) > 0) {
    message("skip  ", src$key[i]); next
  }
  message("fetch ", src$key[i])
  request(src$url[i]) |>
    req_retry(max_tries = 3) |>
    req_timeout(600) |>
    req_perform(path = dest)
}

# ---- cBioPortal PanCancer Atlas clinical data ------------------------------
api <- "https://www.cbioportal.org/api"
get_clin <- function(study, type) {
  request(sprintf("%s/studies/%s/clinical-data", api, study)) |>
    req_url_query(clinicalDataType = type, projection = "SUMMARY",
                  pageSize = 10000000) |>
    req_retry(max_tries = 3) |>
    req_perform() |> resp_body_json(simplifyVector = TRUE) |>
    as.data.table()
}
pancan_studies <- function() {
  studies <- request(paste0(api, "/studies")) |>
    req_url_query(projection = "SUMMARY") |>
    req_perform() |> resp_body_json(simplifyVector = TRUE)
  grep("_tcga_pan_can_atlas_2018$", studies$studyId, value = TRUE)
}

for (type in c("SAMPLE", "PATIENT")) {
  dest <- file.path(raw_dir, sprintf("cbioportal_pancan_%s.tsv", tolower(type)))
  if (file.exists(dest)) { message("skip  cbioportal ", type); next }
  ids <- pancan_studies()
  message("cBioPortal ", type, ": ", length(ids), " PanCanAtlas studies")
  rows <- rbindlist(lapply(ids, function(s) get_clin(s, type)), fill = TRUE)
  keys <- if (type == "SAMPLE") "studyId + patientId + sampleId" else "studyId + patientId"
  wide <- dcast(rows, as.formula(paste(keys, "~ clinicalAttributeId")),
                value.var = "value")
  fwrite(wide, dest, sep = "\t")
}

message("done")
