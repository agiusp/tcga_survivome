# Harmonize all raw sources into one row per patient.
# Outputs:
#   data/processed/survivome.parquet     - analysis table (site loads this)
#   data/processed/survivome.csv.gz      - same, for R/Python users
#   data/processed/data_dictionary.json  - column metadata for the UI
#   data/processed/qc_report.md          - coverage summary
# Usage: Rscript pipeline/03_build.R   (from repo root)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(arrow)
  library(jsonlite)
})

raw <- function(f) file.path("data/raw", f)
out_dir <- "data/processed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

NA_TOKENS <- c("", "NA", "[Not Available]", "[Not Applicable]", "[Unknown]",
               "[Not Evaluated]", "[Discrepancy]", "[Completed]", "[Not Reported]",
               "#N/A")

clean_chr <- function(x) { x <- trimws(as.character(x)); x[x %in% NA_TOKENS] <- NA; x }
patient_of <- function(b) substr(b, 1, 12)
sample_of  <- function(b) substr(b, 1, 15)
sample_code <- function(b) as.integer(substr(b, 14, 15))

# Tumor sample preference when a patient has several: primary solid tumor (01),
# primary blood cancer (03, LAML), metastatic (06), recurrent (02/05), other tumor.
SAMPLE_PRIORITY <- c(`1` = 1, `3` = 2, `9` = 3, `6` = 4, `2` = 5, `5` = 6, `4` = 7, `7` = 8, `8` = 9)

# Keep one tumor sample per patient from a sample-keyed table.
# `id` is a column holding a sample- or aliquot-level barcode.
one_per_patient <- function(dt, id) {
  dt <- copy(dt)
  dt[, sample_id := sample_of(get(id))]
  dt[, code := sample_code(sample_id)]
  dt <- dt[code %in% as.integer(names(SAMPLE_PRIORITY))]
  dt[, prio := SAMPLE_PRIORITY[as.character(code)]]
  setorder(dt, prio)
  dt[, patient_id := patient_of(sample_id)]
  dt <- dt[!duplicated(patient_id)]
  dt[, c("code", "prio") := NULL]
  if (id != "sample_id") dt[, (id) := NULL]
  dt
}

num <- function(x) suppressWarnings(as.numeric(x))

# ---------------------------------------------------------------- clinical (CDR)
cdr <- as.data.table(read_excel(raw("TCGA-CDR-SupplementalTableS1.xlsx"),
                                sheet = "TCGA-CDR", na = NA_TOKENS, guess_max = 20000))
cdr[, `...1` := NULL]
chr_cols <- names(cdr)[sapply(cdr, is.character)]
cdr[, (chr_cols) := lapply(.SD, clean_chr), .SDcols = chr_cols]

roman_stage <- function(s) {
  s <- toupper(s)
  out <- fcase(grepl("^STAGE IV", s), "IV",
               grepl("^STAGE III", s), "III",
               grepl("^STAGE II", s), "II",
               grepl("^STAGE I", s), "I",
               default = NA_character_)
  out
}

clin <- cdr[, .(
  patient_id   = bcr_patient_barcode,
  cancer_type  = type,
  age          = num(age_at_initial_pathologic_diagnosis),
  sex          = tolower(gender),
  race         = tolower(race),
  stage_detail = fcoalesce(ajcc_pathologic_tumor_stage, clinical_stage),
  stage_source = fifelse(!is.na(roman_stage(ajcc_pathologic_tumor_stage)), "pathologic",
                 fifelse(!is.na(roman_stage(clinical_stage)), "clinical", NA_character_)),
  stage        = fcoalesce(roman_stage(ajcc_pathologic_tumor_stage), roman_stage(clinical_stage)),
  histology    = histological_type,
  grade        = histological_grade,
  dx_year      = as.integer(num(initial_pathologic_dx_year)),
  menopause    = menopause_status,
  tumor_status = tolower(tumor_status),
  residual_tumor,
  OS  = as.integer(num(OS)),  OS_time  = num(OS.time),
  DSS = as.integer(num(DSS)), DSS_time = num(DSS.time),
  DFI = as.integer(num(DFI)), DFI_time = num(DFI.time),
  PFI = as.integer(num(PFI)), PFI_time = num(PFI.time),
  redacted = !is.na(Redaction)
)]
clin[, age_group := cut(age, c(-Inf, 40, 50, 60, 70, Inf),
                        labels = c("<40", "40-49", "50-59", "60-69", "70+"), right = FALSE)]
clin[, age_group := as.character(age_group)]

# ------------------------------------------------------ cBioPortal patient-level
cbp <- fread(raw("cbioportal_pancan_patient.tsv"), na.strings = NA_TOKENS)
cbp <- cbp[, .(patient_id = patientId,
               genetic_ancestry = GENETIC_ANCESTRY_LABEL,
               ethnicity = tolower(ETHNICITY),
               neoadjuvant_tx = tolower(HISTORY_NEOADJUVANT_TRTYN),
               hypoxia_buffa  = num(BUFFA_HYPOXIA_SCORE),
               hypoxia_ragnum = num(RAGNUM_HYPOXIA_SCORE),
               hypoxia_winter = num(WINTER_HYPOXIA_SCORE))]
cbp <- cbp[!duplicated(patient_id)]

# ------------------------------------------------------- cBioPortal sample-level
cbs <- fread(raw("cbioportal_pancan_sample.tsv"), na.strings = NA_TOKENS)
cbs <- one_per_patient(cbs[, .(sampleId,
                               msi_mantis   = num(MSI_SCORE_MANTIS),
                               msi_sensor   = num(MSI_SENSOR_SCORE),
                               fga          = num(FRACTION_GENOME_ALTERED),
                               tmb          = num(TMB_NONSYNONYMOUS),
                               tumor_break_load = num(TBL_SCORE),
                               cbio_sample_type = SAMPLE_TYPE)], "sampleId")
# This sample choice defines the canonical tumor sample shown in the UI.
setnames(cbs, "sample_id", "sample_id")

# ------------------------------------------------------------- sample QC flags
qc <- fread(raw("merged_sample_quality_annotations.tsv"), na.strings = NA_TOKENS)
qc[, sample_id := sample_of(aliquot_barcode)]
qc_flag <- qc[, .(qc_do_not_use = any(Do_not_use %in% c("True", TRUE)),
                  qc_pathology_excluded = any(num(AWG_excluded_because_of_pathology) %in% 1)),
              by = sample_id]

# ---------------------------------------------------------------- mutation load
ml <- fread(raw("mutation-load_updated.txt"))
ml <- one_per_patient(ml[, .(Tumor_Sample_ID,
                             tmb_mc3 = num(`Non-silent per Mb`),
                             tmb_mc3_silent = num(`Silent per Mb`))], "Tumor_Sample_ID")

# ------------------------------------------------------------------------- HRD
hrd <- fread(raw("TCGA.HRD_withSampleID.txt"))
hrd <- one_per_patient(hrd[, .(sampleID, hrd_score = num(HRD), hrd_loh = num(`hrd-loh`),
                               hrd_tai = num(ai1), hrd_lst = num(lst1))], "sampleID")

# ------------------------------------------------------------ aneuploidy + arms
# file uses bare CR line endings
an_txt <- gsub("\r\n?", "\n", readChar(raw("PANCAN_ArmCallsAndAneuploidyScore_092817.txt"),
                                       file.size(raw("PANCAN_ArmCallsAndAneuploidyScore_092817.txt")),
                                       useBytes = TRUE))
an <- fread(text = an_txt, na.strings = NA_TOKENS)
arm_cols <- setdiff(names(an), c("Sample", "Type", "Aneuploidy Score"))
arm_names <- paste0("arm_", sub(" .*", "", arm_cols))
setnames(an, arm_cols, arm_names)
setnames(an, "Aneuploidy Score", "aneuploidy_score")
an[, Type := NULL]
an <- one_per_patient(an, "Sample")

# ------------------------------------------------------------ segment burden
seg <- fread(raw("seg_based_scores.tsv"))
seg <- one_per_patient(seg[, .(Sample, cna_n_segments = num(n_segs),
                               cna_frac_altered = num(frac_altered))], "Sample")

# ------------------------------------------------------------------ ABSOLUTE
ab <- fread(raw("TCGA_mastercalls.abs_tables_JSedit.fixed.txt"), na.strings = NA_TOKENS)
ab <- ab[`call status` == "called"]
ab <- one_per_patient(ab[, .(array, purity = num(purity), ploidy = num(ploidy),
                             genome_doublings = num(`Genome doublings`),
                             subclonal_genome_frac = num(`Subclonal genome fraction`))], "array")

# --------------------------------------------------------- immune features
leuk <- fread(raw("TCGA_all_leuk_estimate.masked.20170107.tsv"), header = FALSE,
              col.names = c("type", "aliquot", "leukocyte_fraction"))
leuk <- one_per_patient(leuk[, .(aliquot, leukocyte_fraction = num(leukocyte_fraction))], "aliquot")

imm <- fread(raw("Subtype_Immune_Model_Based.txt.gz"))
imm[, immune_subtype := sub(".*\\((Immune )?(C[1-6])\\).*", "\\2", Subtype_Immune_Model_Based)]
imm[, immune_subtype_label := sub(" \\(Immune C[1-6]\\)", "", Subtype_Immune_Model_Based)]
imm <- one_per_patient(imm[, .(sample, immune_subtype, immune_subtype_label)], "sample")

tcr <- fread(raw("mitcr_sampleStatistics_20160714.tsv"), na.strings = NA_TOKENS)
bcr <- fread(raw("tcga.pancan.igh.div.txt"))
bcr <- merge(bcr, tcr[, .(CGHub_analysis_id, AliquotBarcode)],
             by.x = "sample", by.y = "CGHub_analysis_id")
bcr <- one_per_patient(bcr[, .(AliquotBarcode, bcr_shannon = num(shannon_entropy_by_clus),
                               bcr_clonotypes = num(num_unique_clus))], "AliquotBarcode")
tcr <- one_per_patient(tcr[, .(AliquotBarcode, tcr_shannon = num(shannon),
                               tcr_clonotypes = num(numClones))], "AliquotBarcode")

vir <- fread(raw("viral.tsv"))
vir <- one_per_patient(vir[, .(AliquotBarcode, viral_hpv = num(HPV), viral_hbv = num(HBV),
                               viral_hcv = num(HCV), viral_ebv = num(EBV),
                               viral_mcv = num(MCV))], "AliquotBarcode")

# ------------------------------------------------------- molecular subtypes
sub <- fread(raw("TCGASubtype.20170308.tsv"), na.strings = NA_TOKENS)
sub <- sub[, .(patient_id = patient_of(pan.samplesID),
               subtype_selected = Subtype_Selected, subtype_mrna = Subtype_mRNA,
               subtype_dna_meth = Subtype_DNAmeth, subtype_cna = Subtype_CNA,
               subtype_integrative = Subtype_Integrative)]
sub <- sub[!duplicated(patient_id)]

# ------------------------------------------------------------------- therapy
therapy <- NULL
if (file.exists(raw("therapy_drug.tsv"))) {
  drug <- fread(raw("therapy_drug.tsv"), na.strings = NA_TOKENS)
  rad  <- fread(raw("therapy_radiation.tsv"), na.strings = NA_TOKENS)
  drug[, cls := tolower(pharmaceutical_therapy_type)]
  has <- function(pat) drug[grepl(pat, cls), unique(bcr_patient_barcode)]
  drug_projects <- unique(drug$project)
  therapy <- data.table(patient_id = clin$patient_id, cancer_type = clin$cancer_type)
  flag <- function(ids) fifelse(therapy$patient_id %in% ids, "yes",
                         fifelse(therapy$cancer_type %in% drug_projects, "no record", NA_character_))
  therapy[, `:=`(
    tx_any_drug     = flag(unique(drug$bcr_patient_barcode)),
    tx_chemo        = flag(has("chemo")),
    tx_hormone      = flag(has("hormone")),
    tx_targeted     = flag(has("targeted")),
    tx_immunotherapy= flag(has("immuno|vaccine")),
    tx_radiation    = fifelse(patient_id %in% rad$bcr_patient_barcode, "yes",
                      fifelse(cancer_type %in% unique(rad$project), "no record", NA_character_))
  )]
  therapy[, cancer_type := NULL]
}

# ------------------------------------------------------------------- merge
dt <- merge(clin, cbs, by = "patient_id", all.x = TRUE)
dt <- merge(dt, qc_flag, by = "sample_id", all.x = TRUE)
for (tab in list(cbp, sub, therapy)) if (!is.null(tab)) dt <- merge(dt, tab, by = "patient_id", all.x = TRUE)
for (tab in list(ml, hrd, an, seg, ab, leuk, imm, tcr, bcr, vir)) {
  dt <- merge(dt, tab[, !"sample_id"], by = "patient_id", all.x = TRUE)
}
dt[is.na(qc_do_not_use), qc_do_not_use := FALSE]
dt[is.na(qc_pathology_excluded), qc_pathology_excluded := FALSE]

# --------------------------------------------------------- derived classes
dt[, msi_status := fifelse(is.na(msi_mantis), NA_character_,
                    fifelse(msi_mantis >= 0.4, "MSI-H", "MSS"))]
dt[, tmb_class := fifelse(is.na(tmb), NA_character_, fifelse(tmb >= 10, "TMB-H (>=10)", "TMB-L (<10)"))]
dt[, hrd_class := fifelse(is.na(hrd_score), NA_character_, fifelse(hrd_score >= 42, "HRD-H (>=42)", "HRD-L (<42)"))]
dt[, wgd := fifelse(is.na(genome_doublings), NA_character_,
                    fifelse(genome_doublings >= 1, "WGD", "no WGD"))]

setcolorder(dt, c("patient_id", "sample_id", "cancer_type"))
setorder(dt, cancer_type, patient_id)

# ------------------------------------------------------------- dictionary
# group, label, type (numeric|categorical|boolean|survival_event|survival_time), unit, source
D <- function(col, group, label, type, source, unit = NULL, desc = NULL)
  list(column = col, group = group, label = label, type = type, unit = unit,
       source = source, description = desc)
S_CDR <- "TCGA-CDR (Liu et al. 2018, Cell)"
S_CBIO <- "cBioPortal PanCancer Atlas 2018"
S_IMM <- "Thorsson et al. 2018, Immunity"
dict <- list(
  D("patient_id", "id", "Patient", "id", S_CDR),
  D("sample_id", "id", "Tumor sample", "id", S_CBIO, desc = "Canonical tumor sample used for sample-level features (primary > blood > metastatic > recurrent)."),
  D("cancer_type", "clinical", "Cancer type", "categorical", S_CDR),
  D("age", "clinical", "Age at diagnosis", "numeric", S_CDR, "years"),
  D("age_group", "clinical", "Age group", "categorical", "derived"),
  D("sex", "clinical", "Sex", "categorical", S_CDR),
  D("race", "clinical", "Race (self-reported)", "categorical", S_CDR),
  D("ethnicity", "clinical", "Ethnicity", "categorical", S_CBIO),
  D("genetic_ancestry", "clinical", "Genetic ancestry", "categorical", S_CBIO),
  D("stage", "clinical", "Stage (I-IV)", "categorical", S_CDR, desc = "AJCC pathologic stage collapsed to I-IV; clinical stage used when pathologic is missing."),
  D("stage_detail", "clinical", "Stage (detailed)", "categorical", S_CDR),
  D("stage_source", "clinical", "Stage source", "categorical", "derived"),
  D("histology", "clinical", "Histological type", "categorical", S_CDR),
  D("grade", "clinical", "Histological grade", "categorical", S_CDR),
  D("dx_year", "clinical", "Year of diagnosis", "numeric", S_CDR),
  D("menopause", "clinical", "Menopause status", "categorical", S_CDR),
  D("tumor_status", "clinical", "Tumor status at last contact", "categorical", S_CDR),
  D("residual_tumor", "clinical", "Residual tumor", "categorical", S_CDR),
  D("neoadjuvant_tx", "therapy", "Neoadjuvant treatment history", "categorical", S_CBIO),
  D("tx_any_drug", "therapy", "Any systemic therapy recorded", "categorical", "GDC BCR Biotab drug tables", desc = "'no record' means nothing was recorded, not confirmed untreated."),
  D("tx_chemo", "therapy", "Chemotherapy", "categorical", "GDC BCR Biotab drug tables"),
  D("tx_hormone", "therapy", "Hormone therapy", "categorical", "GDC BCR Biotab drug tables"),
  D("tx_targeted", "therapy", "Targeted therapy", "categorical", "GDC BCR Biotab drug tables"),
  D("tx_immunotherapy", "therapy", "Immunotherapy / vaccine", "categorical", "GDC BCR Biotab drug tables"),
  D("tx_radiation", "therapy", "Radiation therapy", "categorical", "GDC BCR Biotab radiation tables"),
  D("OS", "survival", "Overall survival event", "survival_event", S_CDR),
  D("OS_time", "survival", "Overall survival time", "survival_time", S_CDR, "days"),
  D("DSS", "survival", "Disease-specific survival event", "survival_event", S_CDR),
  D("DSS_time", "survival", "Disease-specific survival time", "survival_time", S_CDR, "days"),
  D("DFI", "survival", "Disease-free interval event", "survival_event", S_CDR),
  D("DFI_time", "survival", "Disease-free interval time", "survival_time", S_CDR, "days"),
  D("PFI", "survival", "Progression-free interval event", "survival_event", S_CDR),
  D("PFI_time", "survival", "Progression-free interval time", "survival_time", S_CDR, "days"),
  D("tmb", "mutation", "TMB (nonsynonymous)", "numeric", S_CBIO, "mut/Mb", desc = "cBioPortal TMB from MC3 calls. Primary TMB column; about 1.3x the PanCanAtlas mutation-load values (smaller exome denominator)."),
  D("tmb_mc3", "mutation", "TMB (PanCanAtlas mutation load)", "numeric", "MC3 (Ellrott et al. 2018)", "mut/Mb", desc = "Non-silent mutations per Mb from the PanCanAtlas mutation-load file. Lower coverage (e.g. CESC 53%)."),
  D("tmb_mc3_silent", "mutation", "Silent mutation rate", "numeric", "MC3 (Ellrott et al. 2018)", "mut/Mb"),
  D("tmb_class", "mutation", "TMB class", "categorical", "derived", desc = "TMB-H at >=10 mut/Mb (FDA pembrolizumab threshold), applied to the cBioPortal TMB column."),
  D("msi_mantis", "mutation", "MSI score (MANTIS)", "numeric", "Bonneville et al. 2017 via cBioPortal"),
  D("msi_sensor", "mutation", "MSI score (MSIsensor)", "numeric", "Bonneville et al. 2017 via cBioPortal"),
  D("msi_status", "mutation", "MSI status", "categorical", "derived", desc = "MSI-H when MANTIS score >= 0.4. Not available for metastatic samples (most SKCM)."),
  D("hrd_score", "dna_repair", "HRD score", "numeric", "Knijnenburg et al. 2018, Cell Rep", desc = "Sum of HRD-LOH, TAI and LST."),
  D("hrd_loh", "dna_repair", "HRD-LOH", "numeric", "Knijnenburg et al. 2018, Cell Rep"),
  D("hrd_tai", "dna_repair", "Telomeric allelic imbalance", "numeric", "Knijnenburg et al. 2018, Cell Rep"),
  D("hrd_lst", "dna_repair", "Large-scale state transitions", "numeric", "Knijnenburg et al. 2018, Cell Rep"),
  D("hrd_class", "dna_repair", "HRD class", "categorical", "derived", desc = "HRD-H at score >= 42 (Telli et al. 2016)."),
  D("aneuploidy_score", "chromosomal_instability", "Aneuploidy score", "numeric", "Taylor et al. 2018, Cancer Cell", "arms"),
  D("fga", "chromosomal_instability", "Fraction genome altered", "numeric", S_CBIO),
  D("cna_frac_altered", "chromosomal_instability", "CNA fraction altered (segments)", "numeric", S_IMM),
  D("cna_n_segments", "chromosomal_instability", "CNA segment count", "numeric", S_IMM),
  D("tumor_break_load", "chromosomal_instability", "Tumor break load", "numeric", S_CBIO),
  D("ploidy", "chromosomal_instability", "Ploidy (ABSOLUTE)", "numeric", "ABSOLUTE (PanCanAtlas)"),
  D("genome_doublings", "chromosomal_instability", "Genome doublings", "numeric", "ABSOLUTE (PanCanAtlas)"),
  D("wgd", "chromosomal_instability", "Whole-genome doubling", "categorical", "derived"),
  D("subclonal_genome_frac", "chromosomal_instability", "Subclonal genome fraction", "numeric", "ABSOLUTE (PanCanAtlas)"),
  D("purity", "tumor", "Tumor purity (ABSOLUTE)", "numeric", "ABSOLUTE (PanCanAtlas)"),
  D("hypoxia_buffa", "tumor", "Hypoxia score (Buffa)", "numeric", S_CBIO),
  D("hypoxia_ragnum", "tumor", "Hypoxia score (Ragnum)", "numeric", S_CBIO),
  D("hypoxia_winter", "tumor", "Hypoxia score (Winter)", "numeric", S_CBIO),
  D("subtype_selected", "subtype", "Molecular subtype", "categorical", "PanCanAtlas subtypes"),
  D("subtype_mrna", "subtype", "mRNA subtype", "categorical", "PanCanAtlas subtypes"),
  D("subtype_dna_meth", "subtype", "DNA methylation subtype", "categorical", "PanCanAtlas subtypes"),
  D("subtype_cna", "subtype", "CNA subtype", "categorical", "PanCanAtlas subtypes"),
  D("subtype_integrative", "subtype", "Integrative subtype", "categorical", "PanCanAtlas subtypes"),
  D("immune_subtype", "immune", "Immune subtype", "categorical", S_IMM, desc = "C1 wound healing, C2 IFN-gamma dominant, C3 inflammatory, C4 lymphocyte depleted, C5 immunologically quiet, C6 TGF-beta dominant. Primary tumors only; DLBC and THYM not assigned."),
  D("immune_subtype_label", "immune", "Immune subtype (name)", "categorical", S_IMM),
  D("leukocyte_fraction", "immune", "Leukocyte fraction", "numeric", S_IMM),
  D("tcr_shannon", "immune", "TCR diversity (Shannon)", "numeric", S_IMM),
  D("tcr_clonotypes", "immune", "TCR clonotypes", "numeric", S_IMM),
  D("bcr_shannon", "immune", "BCR diversity (Shannon)", "numeric", S_IMM),
  D("bcr_clonotypes", "immune", "BCR clonotypes", "numeric", S_IMM),
  D("viral_hpv", "viral", "HPV reads", "numeric", S_IMM, "RPM"),
  D("viral_hbv", "viral", "HBV reads", "numeric", S_IMM, "RPM"),
  D("viral_hcv", "viral", "HCV reads", "numeric", S_IMM, "RPM"),
  D("viral_ebv", "viral", "EBV reads", "numeric", S_IMM, "RPM"),
  D("viral_mcv", "viral", "MCPyV reads", "numeric", S_IMM, "RPM"),
  D("qc_do_not_use", "qc", "AWG 'do not use' flag", "boolean", "PanCanAtlas sample QC"),
  D("qc_pathology_excluded", "qc", "Excluded on pathology review", "boolean", "PanCanAtlas sample QC"),
  D("redacted", "qc", "Redacted case", "boolean", S_CDR),
  D("cbio_sample_type", "qc", "Canonical sample type", "categorical", S_CBIO)
)
for (a in arm_names) dict[[length(dict) + 1]] <-
  D(a, "arm_level_cna", paste("Arm", sub("arm_", "", a)), "categorical", "Taylor et al. 2018, Cancer Cell",
    desc = "-1 loss, 0 neutral, 1 gain.")
for (a in arm_names) dt[, (a) := as.character(get(a))]

missing_dict <- setdiff(names(dt), sapply(dict, `[[`, "column"))
if (length(missing_dict)) warning("columns without dictionary entry: ", paste(missing_dict, collapse = ", "))
dict <- Filter(function(d) d$column %in% names(dt), dict)
for (i in seq_along(dict)) dict[[i]]$coverage <- round(mean(!is.na(dt[[dict[[i]]$column]])), 3)

# ------------------------------------------------------------------ write
write_parquet(dt, file.path(out_dir, "survivome.parquet"), compression = "zstd")
fwrite(dt, file.path(out_dir, "survivome.csv.gz"))
write_json(list(built = format(Sys.time(), "%Y-%m-%d"), n_patients = nrow(dt), columns = dict),
           file.path(out_dir, "data_dictionary.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")

# ---------------------------------------------------------------- QC report
cov <- dt[, lapply(.SD, function(x) round(100 * mean(!is.na(x)))),
          by = cancer_type,
          .SDcols = c("stage", "OS", "PFI", "tmb", "msi_mantis", "hrd_score",
                      "aneuploidy_score", "purity", "immune_subtype",
                      intersect("tx_any_drug", names(dt)))]
ev <- dt[, .(n = .N, OS_events = sum(OS == 1, na.rm = TRUE), PFI_events = sum(PFI == 1, na.rm = TRUE),
             median_followup_mo = round(median(OS_time, na.rm = TRUE) / 30.44, 1),
             tx_drug_recorded = if ("tx_any_drug" %in% names(dt)) sum(tx_any_drug == "yes", na.rm = TRUE) else NA_integer_),
         by = cancer_type]
knitr_table <- function(x) {
  h <- paste("|", paste(names(x), collapse = " | "), "|")
  s <- paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|")
  r <- apply(x, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  c(h, s, r)
}
md <- c("# Survivome build QC", "",
        sprintf("Built %s. %d patients, %d columns.", Sys.Date(), nrow(dt), ncol(dt)), "",
        "## Events and follow-up by cancer type", "", knitr_table(ev), "",
        "## Feature coverage (% of patients)", "", knitr_table(cov))
writeLines(md, file.path(out_dir, "qc_report.md"))

message(sprintf("wrote %d patients x %d columns", nrow(dt), ncol(dt)))
