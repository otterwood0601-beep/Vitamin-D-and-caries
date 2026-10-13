default_config <- function() {
  list(
    output_dir = NULL, cache_dir = file.path(getwd(), "nhanes_cache"),
    split_seed = 123L, cv_seed = 456L, shap_seed = 20261013L,
    bootstrap_reps = 500L, shap_nsim = 512L,
    expected_main_n = 3089L, expected_main_events = 592L,
    expected_day1_n = 2969L, expected_day1_events = 578L,
    check_sample_counts = TRUE, allow_existing_output = FALSE,
    download_retries = 3L, timeout_seconds = 300L,
    stages = c("download", "data", "regression", "rcs", "day1_regression",
               "lasso", "models", "calibration", "shap", "subgroups", "tables", "figures")
  )
}

load_required_packages <- function(stages) {
  pkgs <- c("nhanesA", "dplyr", "tibble", "readr", "survey", "Hmisc",
            "jsonlite", "ggplot2", "patchwork", "writexl")
  if (any(c("lasso", "models") %in% stages)) {
    pkgs <- c(pkgs, "glmnet", "RColorBrewer")
  }
  if ("models" %in% stages) {
    pkgs <- c(pkgs, "caret", "gbm", "lightgbm", "xgboost", "catboost",
              "nnet", "ranger", "pROC", "kernlab", "rpart", "naivebayes")
  }
  if ("shap" %in% stages) pkgs <- c(pkgs, "fastshap", "shapviz")
  pkgs <- unique(pkgs)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install required R packages before running: ",
         paste(missing, collapse = ", "),
         ". CatBoost installation: https://catboost.ai/docs/en/installation/r-installation-binary-installation")
  }
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(pkgs)
}

make_context <- function(config) {
  if (is.null(config$output_dir)) {
    config$output_dir <- file.path(getwd(), "code_results",
                                 format(Sys.time(), "run_%Y%m%d_%H%M%S"))
  }
  if (dir.exists(config$output_dir) &&
      length(list.files(config$output_dir, all.files = TRUE, no.. = TRUE)) &&
      !isTRUE(config$allow_existing_output)) {
    stop("Output directory is not empty. Choose a new directory, or explicitly set allow_existing_output=TRUE.")
  }
  ctx <- list(config = config, root = config$output_dir, cache = config$cache_dir)
  for (nm in c("data", "tables", "rcs", "lasso", "models", "shap", "subgroups", "figures")) {
    ctx[[nm]] <- file.path(ctx$root, nm)
  }
  for (p in unlist(ctx[setdiff(names(ctx), "config")], use.names = FALSE)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  ctx$memory <- new.env(parent = emptyenv())
  ctx
}

nhanes_table_name <- function(component, period) {
  if (!period %in% c("I", "P")) stop("Only cycles I and P are supported.")
  if (period == "I") paste0(component, "_I") else paste0("P_", component)
}

fetch_nhanes <- function(table, ctx) {
  if (!grepl("^[A-Za-z0-9_]+$", table)) stop("Invalid NHANES table name.")
  if (exists(table, envir = ctx$memory, inherits = FALSE)) {
    return(get(table, envir = ctx$memory, inherits = FALSE))
  }
  cache_file <- file.path(ctx$cache, paste0(table, ".rds"))
  if (file.exists(cache_file)) {
    obj <- readRDS(cache_file)
    if (!is.list(obj) || !identical(obj$table, table) ||
        !identical(obj$provider, "nhanesA") || !is.data.frame(obj$data)) {
      stop("Unrecognized cache file: ", cache_file)
    }
    value <- obj$data
  } else {
    last_error <- "No data returned"
    value <- NULL
    for (attempt in seq_len(ctx$config$download_retries)) {
      message("Downloading ", table, " using nhanesA (attempt ", attempt, ")")
      value <- tryCatch(
        nhanesA::nhanes(table, translated = FALSE, includelabels = FALSE,
                       cleanse_numeric = FALSE),
        error = function(e) { last_error <<- conditionMessage(e); NULL }
      )
      if (is.data.frame(value) && nrow(value) > 0L) break
      if (attempt < ctx$config$download_retries) Sys.sleep(attempt)
    }
    if (!is.data.frame(value) || !nrow(value)) {
      stop("NHANES download failed for ", table, ": ", last_error)
    }
    value <- as.data.frame(value)
    obj <- list(table = table, provider = "nhanesA", translated = FALSE,
                retrieved_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                package_version = as.character(utils::packageVersion("nhanesA")),
                data = value)
    saveRDS(obj, cache_file)
  }
  assign(table, value, envir = ctx$memory)
  value
}

download_inputs <- function(ctx) {
  components <- c("DEMO", "BMX", "SMQ", "ALQ", "PAQ", "OHQ", "DIQ",
                  "BPQ", "MCQ", "KIQ_U", "DR1TOT", "DSQTOT", "OHXDEN",
                  "DS1TOT", "DS1IDS")
  tables <- c(unlist(lapply(c("I", "P"), function(p) {
    vapply(components, nhanes_table_name, character(1), period = p)
  }), use.names = FALSE), "DSPI", "DSII", "DSBI")
  rows <- lapply(tables, function(table) {
    z <- fetch_nhanes(table, ctx)
    f <- file.path(ctx$cache, paste0(table, ".rds"))
    data.frame(table = table, rows = nrow(z), columns = ncol(z),
               provider = "nhanesA", translated = FALSE,
               cache_file = basename(f), md5 = unname(tools::md5sum(f)))
  })
  write.csv(do.call(rbind, rows), file.path(ctx$root, "download_manifest.csv"), row.names = FALSE)
  invisible(tables)
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Run prerequisite stages first. Missing: ", paste(missing, collapse = "; "))
}

finite_complete <- function(df, vars) {
  ok <- complete.cases(df[, vars, drop = FALSE])
  nums <- vars[vapply(df[, vars, drop = FALSE], is.numeric, logical(1))]
  if (length(nums)) ok <- ok & Reduce("&", lapply(df[nums], is.finite))
  ok
}

save_plot_pair <- function(plot, stem, width, height) {
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height,
                  units = "in", dpi = 300, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
                  units = "in", device = grDevices::cairo_pdf, bg = "white")
}

prepare_data <- function(ctx) {
  work_dir <- ctx$data
  read_component <- function(component, period, variables) {
    x <- fetch_nhanes(nhanes_table_name(component, period), ctx)
    if (!"SEQN" %in% names(x)) stop(component, " lacks SEQN")
    keep <- intersect(unique(c("SEQN", variables)), names(x))
    out <- x[, keep, drop = FALSE]
    if (anyDuplicated(out$SEQN)) stop(component, " contains duplicate SEQN")
    out
  }
  tooth_ids <- sprintf("%02d", setdiff(1:32, c(1, 16, 17, 32)))

  derive_oral <- function(period) {
    x <- fetch_nhanes(nhanes_table_name("OHXDEN", period), ctx)
    tc_vars <- paste0("OHX", tooth_ids, "TC")
    csc_vars <- paste0("OHX", tooth_ids, "CSC")
    stopifnot(all(tc_vars %in% names(x)), all(csc_vars %in% names(x)))

    tc <- as.data.frame(lapply(x[tc_vars], as.numeric))
    csc <- as.data.frame(lapply(x[csc_vars], function(v) {
      z <- as.character(v)
      z[is.na(z)] <- ""
      z
    }))

    present <- rowSums(do.call(cbind, lapply(tc, function(v) v %in% c(1, 2, 5))), na.rm = TRUE)
    observed <- rowSums(!is.na(tc))

    surface_positive <- do.call(cbind, lapply(csc, grepl, pattern = "-[4-8]", perl = TRUE))
    event <- rowSums(surface_positive, na.rm = TRUE) > 0

    tibble(
      SEQN = as.numeric(x$SEQN),
      oral_exam_status = as.numeric(x$OHDDESTS),
      oral_evaluable = as.numeric(x$OHDDESTS) %in% c(1, 2),
      tooth_status_observed_n_28 = observed,
      present_teeth_28 = present,
      missing_teeth_28 = 28 - present,
      untreated_surface_caries_28 = as.integer(event),
      untreated_carious_teeth_28 = rowSums(surface_positive, na.rm = TRUE)
    )
  }

  domain_minutes <- function(yes_no, days, minutes) {
    yes_no <- as.numeric(yes_no)
    days <- as.numeric(days)
    minutes <- as.numeric(minutes)
    out <- rep(NA_real_, length(yes_no))
    out[!is.na(yes_no) & yes_no == 2] <- 0
    valid <- !is.na(yes_no) & yes_no == 1 & !is.na(days) & days %in% 1:7 &
      !is.na(minutes) & minutes >= 1 & minutes < 7777
    out[valid] <- days[valid] * minutes[valid]
    out
  }

  derive_period <- function(period) {
    demo <- read_component("DEMO", period, c(
      "RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "DMDMARTL",
      "DMDMARTZ", "INDFMPIR", "SDMVSTRA", "SDMVPSU"
    ))
    bmx <- read_component("BMX", period, "BMXBMI")
    smq <- read_component("SMQ", period, c("SMQ020", "SMQ040"))
    alq <- read_component("ALQ", period, c(
      "ALQ101", "ALQ110", "ALQ111", "ALQ120Q", "ALQ120U", "ALQ121", "ALQ130"
    ))
    paq <- read_component("PAQ", period, c(
      "PAQ605", "PAQ610", "PAD615", "PAQ620", "PAQ625", "PAD630",
      "PAQ635", "PAQ640", "PAD645", "PAQ650", "PAQ655", "PAD660",
      "PAQ665", "PAQ670", "PAD675"
    ))
    ohq <- read_component("OHQ", period, c("OHQ030", "OHQ870"))
    diq <- read_component("DIQ", period, "DIQ010")
    bpq <- read_component("BPQ", period, "BPQ020")
    mcq <- read_component("MCQ", period, c(
      "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F", "MCQ160L"
    ))
    kiq <- read_component("KIQ_U", period, "KIQ022")
    diet <- read_component("DR1TOT", period, c(
      "DR1DRSTZ", "DR1TKCAL", "DR1TSUGR", "DR1TVD", "DR1TCHOL",
      "DR1TTFAT", "DR1TFIBE", "DR1TPROT", "DR1TCALC", "DR1TPHOS",
      "WTDRD1", "WTDRD1PP"
    ))
    dsq <- read_component("DSQTOT", period, c("DSD010", "DSQTVD"))
    oral <- derive_oral(period)

    d <- Reduce(function(x, y) left_join(x, y, by = "SEQN"),
                list(demo, bmx, smq, alq, paq, ohq, diq, bpq, mcq, kiq, diet, dsq, oral))
    for (nm in c("DMDMARTL", "DMDMARTZ", "ALQ120Q", "ALQ120U", "ALQ121",
                 "ALQ130", "WTDRD1", "WTDRD1PP")) {
      if (!nm %in% names(d)) d[[nm]] <- NA_real_
    }

    freq_map <- c(`0` = 0, `1` = 365, `2` = 300, `3` = 182, `4` = 104,
                  `5` = 52, `6` = 30, `7` = 12, `8` = 9, `9` = 4.5, `10` = 1.5)
    freq_days <- if (period == "I") {
      ifelse(d$ALQ120Q == 0, 0, as.numeric(d$ALQ120Q) *
               ifelse(d$ALQ120U == 1, 52,
                      ifelse(d$ALQ120U == 2, 12, ifelse(d$ALQ120U == 3, 1, NA_real_))))
    } else {
      unname(freq_map[as.character(d$ALQ121)])
    }
    avg_drinks_all_days <- freq_days * as.numeric(d$ALQ130) / 365
    activity_equivalent_minutes <-
      2 * domain_minutes(d$PAQ605, d$PAQ610, d$PAD615) +
      domain_minutes(d$PAQ620, d$PAQ625, d$PAD630) +
      domain_minutes(d$PAQ635, d$PAQ640, d$PAD645) +
      2 * domain_minutes(d$PAQ650, d$PAQ655, d$PAD660) +
      domain_minutes(d$PAQ665, d$PAQ670, d$PAD675)

    marital_raw <- if (period == "P") d$DMDMARTZ else d$DMDMARTL
    cvd_vars <- c("MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F")
    cvd_mat <- as.data.frame(lapply(d[cvd_vars], as.numeric))
    cvd_yes <- rowSums(cvd_mat == 1, na.rm = TRUE) > 0
    cvd_all_no <- rowSums(cvd_mat == 2, na.rm = TRUE) == length(cvd_vars)

    d %>%
      mutate(
        source_period = period,
        age = as.numeric(RIDAGEYR),
        sex = case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female", TRUE ~ NA_character_),
        race = case_when(
          RIDRETH3 == 1 ~ "Mexican American", RIDRETH3 == 2 ~ "Other Hispanic",
          RIDRETH3 == 3 ~ "Non-Hispanic White", RIDRETH3 == 4 ~ "Non-Hispanic Black",
          RIDRETH3 == 6 ~ "Non-Hispanic Asian", RIDRETH3 == 7 ~ "Other/multiracial",
          TRUE ~ NA_character_
        ),
        education = case_when(
          DMDEDUC2 == 1 ~ "Less than 9th grade",
          DMDEDUC2 %in% 2:5 ~ "9th grade or above", TRUE ~ NA_character_
        ),
        marital = if (period == "P") case_when(
          marital_raw == 1 ~ "Married/living with partner",
          marital_raw == 2 ~ "Widowed/divorced/separated",
          marital_raw == 3 ~ "Never married", TRUE ~ NA_character_
        ) else case_when(
          marital_raw %in% c(1, 6) ~ "Married/living with partner",
          marital_raw %in% c(2, 3, 4) ~ "Widowed/divorced/separated",
          marital_raw == 5 ~ "Never married", TRUE ~ NA_character_
        ),
        pir = as.numeric(INDFMPIR),
        bmi = as.numeric(BMXBMI),
        smoking = case_when(
          SMQ020 == 2 ~ "Never smoker",
          SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
          SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker", TRUE ~ NA_character_
        ),
        alcohol = case_when(
          ALQ130 %in% c(777, 999) ~ NA_character_,
          is.na(ALQ130) ~ "None",
          sex == "Male" & avg_drinks_all_days > 2 ~ "Heavy",
          sex == "Female" & avg_drinks_all_days > 1 ~ "Heavy",
          TRUE ~ "Light to moderate"
        ),
        physical_activity = case_when(
          is.na(activity_equivalent_minutes) ~ NA_character_,
          activity_equivalent_minutes >= 150 ~ ">=150 min/week", TRUE ~ "<150 min/week"
        ),
        floss_days = ifelse(OHQ870 %in% 0:7, as.numeric(OHQ870), NA_real_),
        dental_visit = case_when(
          OHQ030 %in% c(1, 2) ~ "Within the last year",
          OHQ030 %in% c(3, 4, 5) ~ ">1 year but <=5 years",
          OHQ030 %in% c(6, 7) ~ ">5 years/never", TRUE ~ NA_character_
        ),
        diabetes = case_when(DIQ010 == 1 ~ "Yes", DIQ010 %in% c(2, 3) ~ "No", TRUE ~ NA_character_),
        hypertension = case_when(BPQ020 == 1 ~ "Yes", BPQ020 == 2 ~ "No", TRUE ~ NA_character_),
        cvd = ifelse(cvd_yes, "Yes", ifelse(cvd_all_no, "No", NA_character_)),
        ckd = case_when(KIQ022 == 1 ~ "Yes", KIQ022 == 2 ~ "No", TRUE ~ NA_character_),
        liver_disease = case_when(MCQ160L == 1 ~ "Yes", MCQ160L == 2 ~ "No", TRUE ~ NA_character_),
        diet_recall_reliable = DR1DRSTZ == 1,
        dietary_vitd_ug = as.numeric(DR1TVD),
        supplement_question_valid = DSD010 %in% c(1, 2),
        supplement_vitd_ug = ifelse(supplement_question_valid & is.na(DSQTVD), 0, as.numeric(DSQTVD)),
        energy_kcal = as.numeric(DR1TKCAL), sugar_g = as.numeric(DR1TSUGR),
        cholesterol_mg = as.numeric(DR1TCHOL), fat_g = as.numeric(DR1TTFAT),
        fiber_g = as.numeric(DR1TFIBE), protein_g = as.numeric(DR1TPROT),
        calcium_mg = as.numeric(DR1TCALC), phosphorus_mg = as.numeric(DR1TPHOS),
        diet_weight = as.numeric(if (period == "P") WTDRD1PP else WTDRD1),
        SDMVSTRA = as.numeric(SDMVSTRA), SDMVPSU = as.numeric(SDMVPSU)
      ) %>%
      select(
        SEQN, source_period, age, sex, race, education, marital, pir, bmi,
        smoking, alcohol, physical_activity, floss_days, dental_visit,
        diabetes, hypertension, cvd, ckd, liver_disease, diet_recall_reliable,
        dietary_vitd_ug, supplement_question_valid, supplement_vitd_ug,
        energy_kcal, sugar_g, cholesterol_mg, fat_g, fiber_g, protein_g,
        calcium_mg, phosphorus_mg, oral_evaluable, present_teeth_28,
        missing_teeth_28, untreated_surface_caries_28, untreated_carious_teeth_28,
        tooth_status_observed_n_28, SDMVSTRA, SDMVPSU, diet_weight
      )
  }

  period_data <- bind_rows(derive_period("I"), derive_period("P"))

  factor_levels <- list(
    source_period = c("I", "P"),
    sex = c("Female", "Male"),
    race = c("Mexican American", "Other Hispanic", "Non-Hispanic White",
             "Non-Hispanic Black", "Non-Hispanic Asian", "Other/multiracial"),
    education = c("Less than 9th grade", "9th grade or above"),
    marital = c("Married/living with partner", "Never married", "Widowed/divorced/separated"),
    smoking = c("Never smoker", "Former smoker", "Current smoker"),
    alcohol = c("None", "Light to moderate", "Heavy"),
    physical_activity = c("<150 min/week", ">=150 min/week"),
    dental_visit = c("Within the last year", ">1 year but <=5 years", ">5 years/never"),
    diabetes = c("No", "Yes"), hypertension = c("No", "Yes"),
    cvd = c("No", "Yes"), ckd = c("No", "Yes"), liver_disease = c("No", "Yes")
  )
  for (v in names(factor_levels)) period_data[[v]] <- factor(period_data[[v]], levels = factor_levels[[v]])

  model_covariates <- c(
    "source_period", "age", "sex", "race", "education", "marital", "pir", "bmi",
    "smoking", "alcohol", "physical_activity", "floss_days", "dental_visit",
    "diabetes", "hypertension", "cvd", "ckd", "liver_disease", "present_teeth_28",
    "energy_kcal", "sugar_g", "calcium_mg", "cholesterol_mg", "fat_g", "fiber_g",
    "protein_g", "phosphorus_mg"
  )
  required <- c("untreated_surface_caries_28", "dietary_vitd_ug", "supplement_vitd_ug",
                model_covariates, "diet_weight", "SDMVSTRA", "SDMVPSU")

  finite_complete <- function(df, vars) {
    ok <- complete.cases(df[, vars, drop = FALSE])
    nums <- vars[vapply(df[, vars, drop = FALSE], is.numeric, logical(1))]
    if (length(nums)) ok <- ok & Reduce(`&`, lapply(df[nums], is.finite))
    ok
  }

  flow <- tibble(
    step = c(
      "All I and P records", "Public age code >=60 years", "Oral examination evaluable",
      "At least one natural tooth (28-tooth definition)", "Reliable day-1 dietary recall",
      "Food and supplement vitamin D available", "Survey design variables available",
      "All model covariates complete"
    )
  )
  idx1 <- rep(TRUE, nrow(period_data))
  idx2 <- idx1 & is.finite(period_data$age) & period_data$age >= 60
  idx3 <- idx2 & period_data$oral_evaluable & period_data$untreated_surface_caries_28 %in% 0:1
  idx4 <- idx3 & period_data$present_teeth_28 > 0
  idx5 <- idx4 & period_data$diet_recall_reliable
  idx6 <- idx5 & period_data$supplement_question_valid & is.finite(period_data$dietary_vitd_ug) &
    is.finite(period_data$supplement_vitd_ug)
  idx7 <- idx6 & is.finite(period_data$diet_weight) & period_data$diet_weight > 0 &
    !is.na(period_data$SDMVSTRA) & !is.na(period_data$SDMVPSU)
  idx8 <- idx7 & finite_complete(period_data, model_covariates)
  indices <- list(idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8)
  flow$remaining_n <- vapply(indices, sum, integer(1))
  flow$excluded_at_step_n <- c(0L, -diff(flow$remaining_n))

  analysis <- period_data[idx8, , drop = FALSE] %>%
    mutate(
      analysis_weight = ifelse(source_period == "I", diet_weight * (2 / 5.2), diet_weight * (3.2 / 5.2)),
      survey_strata = interaction(source_period, SDMVSTRA, drop = TRUE),
      survey_psu = interaction(source_period, SDMVPSU, drop = TRUE),
      outcome = factor(untreated_surface_caries_28, levels = c(0, 1), labels = c("No", "Yes"))
    )

  design <- svydesign(ids = ~survey_psu, strata = ~survey_strata, weights = ~analysis_weight,
                      nest = TRUE, data = analysis)

  fmt_p <- function(p) {
    if (!is.finite(p)) return("")
    if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  }

  weighted_mean_se <- function(des, var, digits) {
    st <- svymean(as.formula(paste0("~", var)), des, na.rm = TRUE)
    sprintf(paste0("%.", digits, "f +/- %.", digits, "f"), unname(coef(st)[1]), unname(SE(st)[1]))
  }

  weighted_median_iqr <- function(des, var, digits) {
    st <- svyquantile(as.formula(paste0("~", var)), des, quantiles = c(.25, .50, .75),
                      ci = FALSE, na.rm = TRUE)
    q <- unname(as.numeric(coef(st)))
    sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), q[2], q[1], q[3])
  }

  weighted_n_percent <- function(des, var, level) {
    z <- as.numeric(as.character(des$variables[[var]]) == level)
    des2 <- update(des, .indicator = z)
    st <- svymean(~.indicator, des2, na.rm = TRUE)
    n <- sum(as.character(des$variables[[var]]) == level, na.rm = TRUE)
    sprintf("%d (%.2f)", n, 100 * unname(coef(st)[1]))
  }

  p_cont <- function(var, type) {
    tryCatch({
      if (type == "median") {
        unname(svyranktest(as.formula(paste0(var, "~outcome")), design, test = "wilcoxon")$p.value)
      } else {
        fit <- svyglm(as.formula(paste0(var, "~outcome")), design)
        unname(coef(summary(fit))[2, 4])
      }
    }, error = function(e) NA_real_)
  }

  p_cat <- function(var) {
    tryCatch(unname(svychisq(as.formula(paste0("~", var, "+outcome")), design, statistic = "F")$p.value),
             error = function(e) NA_real_)
  }

  continuous_specs <- tribble(
    ~variable, ~label, ~type, ~digits,
    "dietary_vitd_ug", "Dietary vitamin D, Mean +/- S.E., ug/day", "mean", 2,
    "supplement_vitd_ug", "Supplement vitamin D, Mean +/- S.E., ug/day", "mean", 2,
    "age", "Age, Mean +/- S.E., years", "mean", 2,
    "bmi", "BMI, Mean +/- S.E., kg/m2", "mean", 2,
    "pir", "PIR, Mean +/- S.E.", "mean", 2,
    "floss_days", "Flosser/interdental device use, Mean +/- S.E., days/week", "mean", 2,
    "present_teeth_28", "Present natural teeth, Mean +/- S.E., n", "mean", 2,
    "missing_teeth_28", "Missing teeth, Mean +/- S.E., n", "mean", 2,
    "energy_kcal", "Energy intake, Mean +/- S.E., kcal/day", "mean", 2,
    "sugar_g", "Total sugar, Mean +/- S.E., g/day", "mean", 2,
    "cholesterol_mg", "Cholesterol, Mean +/- S.E., mg/day", "mean", 2,
    "fat_g", "Total fat, Mean +/- S.E., g/day", "mean", 2,
    "fiber_g", "Dietary fiber, Mean +/- S.E., g/day", "mean", 2,
    "protein_g", "Protein, Mean +/- S.E., g/day", "mean", 2,
    "calcium_mg", "Calcium, Mean +/- S.E., mg/day", "mean", 2,
    "phosphorus_mg", "Phosphorus, Mean +/- S.E., mg/day", "mean", 2
  )

  categorical_specs <- list(
    list(variable = "sex", label = "Gender (%)", levels = setNames(factor_levels$sex, factor_levels$sex)),
    list(variable = "race", label = "Race/ethnicity (%)", levels = setNames(factor_levels$race, factor_levels$race)),
    list(variable = "education", label = "Education (%)", levels = setNames(factor_levels$education, factor_levels$education)),
    list(variable = "marital", label = "Marital status (%)", levels = setNames(factor_levels$marital, factor_levels$marital)),
    list(variable = "smoking", label = "Smoking status (%)", levels = setNames(factor_levels$smoking, factor_levels$smoking)),
    list(variable = "alcohol", label = "Alcohol use (%)", levels = setNames(factor_levels$alcohol, factor_levels$alcohol)),
    list(variable = "physical_activity", label = "Physical activity (%)", levels = setNames(factor_levels$physical_activity, factor_levels$physical_activity)),
    list(variable = "dental_visit", label = "Last dentist visit (%)", levels = setNames(factor_levels$dental_visit, factor_levels$dental_visit)),
    list(variable = "diabetes", label = "Diabetes (%)", levels = setNames(factor_levels$diabetes, factor_levels$diabetes)),
    list(variable = "hypertension", label = "Hypertension (%)", levels = setNames(factor_levels$hypertension, factor_levels$hypertension)),
    list(variable = "cvd", label = "Broad CVD (%)", levels = setNames(factor_levels$cvd, factor_levels$cvd)),
    list(variable = "ckd", label = "CKD (%)", levels = setNames(factor_levels$ckd, factor_levels$ckd)),
    list(variable = "liver_disease", label = "Liver disease (%)", levels = setNames(factor_levels$liver_disease, factor_levels$liver_disease))
  )

  design_no <- subset(design, outcome == "No")
  design_yes <- subset(design, outcome == "Yes")

  rows <- list()
  for (i in seq_len(nrow(continuous_specs))) {
    s <- continuous_specs[i, ]
    formatter <- if (s$type == "median") weighted_median_iqr else weighted_mean_se
    rows[[length(rows) + 1]] <- tibble(
      Variables = s$label,
      Total = formatter(design, s$variable, s$digits),
      No = formatter(design_no, s$variable, s$digits),
      Yes = formatter(design_yes, s$variable, s$digits),
      `p-value` = fmt_p(p_cont(s$variable, s$type)),
      row_type = "continuous"
    )
  }

  for (s in categorical_specs) {
    rows[[length(rows) + 1]] <- tibble(
      Variables = s$label, Total = "", No = "", Yes = "",
      `p-value` = fmt_p(p_cat(s$variable)), row_type = "group"
    )
    for (raw_level in names(s$levels)) {
      display_level <- unname(s$levels[[raw_level]])
      rows[[length(rows) + 1]] <- tibble(
        Variables = display_level,
        Total = weighted_n_percent(design, s$variable, raw_level),
        No = weighted_n_percent(design_no, s$variable, raw_level),
        Yes = weighted_n_percent(design_yes, s$variable, raw_level),
        `p-value` = "", row_type = "level"
      )
    }
  }

  table1 <- bind_rows(rows)
  n_total <- nrow(analysis)
  n_no <- sum(analysis$untreated_surface_caries_28 == 0)
  n_yes <- sum(analysis$untreated_surface_caries_28 == 1)

  metadata <- tibble(
    key = c("n_total", "n_no", "n_yes", "design_df", "I_n", "P_n",
            "weighted_event_prevalence", "surface_positive_teeth_max",
            "outcome_definition", "supplement_definition", "weight_definition"),
    value = c(
      n_total, n_no, n_yes, degf(design),
      sum(analysis$source_period == "I"), sum(analysis$source_period == "P"),
      sprintf("%.6f", unname(coef(svymean(~untreated_surface_caries_28, design))[1])),
      max(analysis$untreated_carious_teeth_28),
      "Any untreated carious surface among 28 non-third-molar teeth",
      "DSQTVD, average daily vitamin D from supplements in the prior 30 days; valid nonusers coded 0",
      "I WTDRD1*(2/5.2); P WTDRD1PP*(3.2/5.2)"
    )
  )

  if (isTRUE(ctx$config$check_sample_counts)) {
    stopifnot(nrow(analysis) == ctx$config$expected_main_n,
              sum(analysis$untreated_surface_caries_28) == ctx$config$expected_main_events)
  }
  write_excel_csv(table1, file.path(work_dir, "table1_values.csv"))
  write_excel_csv(flow, file.path(work_dir, "sample_flow.csv"))
  write_excel_csv(metadata, file.path(work_dir, "metadata.csv"))
  write_excel_csv(analysis, file.path(work_dir, "main_analysis.csv"))
  saveRDS(analysis, file.path(work_dir, "analysis_complete_case.rds"))
  saveRDS(list(period_data = period_data, factor_levels = factor_levels,
               model_covariates = model_covariates, continuous_specs = continuous_specs,
               categorical_specs = categorical_specs),
          file.path(work_dir, "data_context.rds"))
  writexl::write_xlsx(list(Table_1 = table1, Sample_flow = flow, Metadata = metadata),
                      file.path(ctx$tables, "Table_1.xlsx"))
  message("Main cohort: ", nrow(analysis), "; events: ", sum(analysis$untreated_surface_caries_28))
  invisible(analysis)
}

fit_main_regression <- function(ctx) {
  require_files(file.path(ctx$data, "analysis_complete_case.rds"))
  work_dir <- ctx$tables
  analysis <- readRDS(file.path(ctx$data, "analysis_complete_case.rds"))
  design <- svydesign(
    ids = ~survey_psu,
    strata = ~survey_strata,
    weights = ~analysis_weight,
    nest = TRUE,
    data = analysis
  )

  demographic_covariates <- c(
    "age", "sex", "race", "education", "marital", "pir"
  )
  all_covariates <- c(
    demographic_covariates,
    "bmi", "smoking", "alcohol", "physical_activity", "floss_days",
    "dental_visit", "diabetes", "hypertension", "cvd", "ckd",
    "liver_disease", "present_teeth_28", "energy_kcal", "sugar_g",
    "cholesterol_mg", "fat_g", "fiber_g", "protein_g", "calcium_mg",
    "phosphorus_mg"
  )

  exposures <- tribble(
    ~variable, ~label,
    "dietary_vitd_ug", "Dietary vitamin D (per 1 ug/day)",
    "supplement_vitd_ug", "Supplement vitamin D (per 1 ug/day)"
  )

  fit_model <- function(exposure, covariates = character()) {
    formula <- reformulate(c(exposure, covariates), response = "untreated_surface_caries_28")
    svyglm(formula, design = design, family = quasibinomial("logit"))
  }

  extract_effect <- function(fit, exposure) {
    beta <- unname(coef(fit)[exposure])
    standard_error <- sqrt(unname(vcov(fit)[exposure, exposure]))
    degrees_freedom <- degf(design)
    critical_value <- qt(0.975, degrees_freedom)
    tibble(
      estimate = exp(beta),
      conf_low = exp(beta - critical_value * standard_error),
      conf_high = exp(beta + critical_value * standard_error),
      p_value = 2 * pt(abs(beta / standard_error), df = degrees_freedom, lower.tail = FALSE)
    )
  }

  format_p <- function(value) {
    if (!is.finite(value)) return("")
    if (value < 0.001) "<0.001" else sprintf("%.3f", value)
  }

  format_or_ci <- function(estimate, low, high) {
    sprintf("%.3f (%.3f-%.3f)", estimate, low, high)
  }

  rows <- lapply(seq_len(nrow(exposures)), function(index) {
    exposure <- exposures$variable[index]
    crude <- extract_effect(fit_model(exposure), exposure)
    model1 <- extract_effect(fit_model(exposure, demographic_covariates), exposure)
    model2 <- extract_effect(fit_model(exposure, all_covariates), exposure)
    tibble(
      Variables = exposures$label[index],
      Crude_p = format_p(crude$p_value),
      Crude_OR_CI = format_or_ci(crude$estimate, crude$conf_low, crude$conf_high),
      Model1_p = format_p(model1$p_value),
      Model1_OR_CI = format_or_ci(model1$estimate, model1$conf_low, model1$conf_high),
      Model2_p = format_p(model2$p_value),
      Model2_OR_CI = format_or_ci(model2$estimate, model2$conf_low, model2$conf_high)
    )
  })

  results <- bind_rows(rows)
  details <- bind_rows(lapply(seq_len(nrow(exposures)), function(index) {
    exposure <- exposures$variable[index]
    bind_rows(
      mutate(extract_effect(fit_model(exposure), exposure), model = "Crude"),
      mutate(extract_effect(fit_model(exposure, demographic_covariates), exposure), model = "Model 1"),
      mutate(extract_effect(fit_model(exposure, all_covariates), exposure), model = "Model 2")
    ) |>
      mutate(
        variable = exposure,
        label = exposures$label[index],
        n = nrow(analysis),
        events = sum(analysis$untreated_surface_caries_28 == 1),
        design_df = degf(design),
        .before = 1
      )
  }))

  write_excel_csv(results, file.path(work_dir, "logistic_results_wide.csv"))
  write_excel_csv(details, file.path(work_dir, "logistic_results_details.csv"))

  invisible(details)
}

prepare_day1_exposure <- function(ctx) {
  require_files(file.path(ctx$data, c("data_context.rds", "analysis_complete_case.rds")))
  context <- readRDS(file.path(ctx$data, "data_context.rds"))
  base <- context$period_data
  model_covariates <- context$model_covariates
  out <- ctx$rcs
  pi <- fetch_nhanes("DSPI", ctx)
  ii <- fetch_nhanes("DSII", ctx)
  bi <- fetch_nhanes("DSBI", ctx)
  bi_d<-bi[grepl('\\bVITAMIN D\\b|\\bCHOLECALCIFEROL\\b|\\bERGOCALCIFEROL\\b|\\bCALCIFEDIOL\\b',bi$DSDBCNAM,ignore.case=TRUE),]
  vitd_ids<-unique(c(ii$DSDPID[toupper(trimws(ii$DSDINGR))=='VITAMIN D'],ii$DSDPID[ii$DSDIID%in%bi_d$DSDIID]))
  totals<-list();products<-list()
  for(p in c('I','P')){
   t<-fetch_nhanes(nhanes_table_name('DS1TOT',p),ctx)
   a<-fetch_nhanes(nhanes_table_name('DS1IDS',p),ctx)
   if(!'DSDPID'%in%names(a))a$DSDPID<-pi$DSDPID[match(a$DSDSUPID,pi$DSDSUPID)]
   cat(p,'product names:',paste(names(a),collapse=', '),'\n')
   a$source_period<-p
   a$label_known<-a$DSDPID%in%ii$DSDPID
   a$label_d<-a$DSDPID%in%vitd_ids
   a$unresolved<-!a$label_known|(a$label_d & !is.finite(a$DS1IVD))
   a$vitd_known<-ifelse(is.finite(a$DS1IVD),as.numeric(a$DS1IVD),ifelse(a$label_known&!a$label_d,0,NA_real_))
   products[[p]]<-a[,c('SEQN','source_period','DSDPID','DS1IVD','label_known','label_d','unresolved','vitd_known')]
   agg<-a%>%group_by(SEQN)%>%summarise(products=n(),unresolved=any(unresolved),reconstructed=sum(vitd_known),.groups='drop')
   t<-left_join(t,agg,by='SEQN')
   t$source_period<-p
   none<-t$DS1DS%in%2 & t$DS1AN%in%2
   taken<-t$DS1DS%in%1 | t$DS1AN%in%1
   flags_valid<-t$DS1DS%in%c(1,2)&t$DS1AN%in%c(1,2)
   t$status<-case_when(!flags_valid~'Missing supplement/antacid recall',none~'Confirmed nonuser on day 1',taken&is.na(t$products)~'User without product record',taken&t$unresolved~'Product or vitamin D dose unresolved',taken&!t$unresolved~'All products resolved',TRUE~'Unresolved')
   t$day1_vitd_ug<-ifelse(none,0,ifelse(flags_valid&t$status=='All products resolved',t$reconstructed,NA_real_))
   t$raw_total<-as.numeric(t$DS1TVD)
   good<-is.finite(t$day1_vitd_ug)&is.finite(t$raw_total)
   stopifnot(all(abs(t$day1_vitd_ug[good]-t$raw_total[good])<.011))
   stopifnot(!anyDuplicated(t$SEQN))
   totals[[p]]<-t[,c('SEQN','source_period','DS1DS','DS1AN','DS1DSCNT','DS1ANCNT','products','status','raw_total','day1_vitd_ug')]
  }
  t<-bind_rows(totals);a<-bind_rows(products)
  d<-left_join(base,t,by=c('SEQN','source_period'))
  eligible<-with(d,is.finite(age)&age>=60&oral_evaluable%in%TRUE&present_teeth_28>0&diet_recall_reliable%in%TRUE&is.finite(dietary_vitd_ug)&is.finite(diet_weight)&diet_weight>0&!is.na(SDMVSTRA)&!is.na(SDMVPSU))&finite_complete(d,model_covariates)
  eligible[is.na(eligible)]<-FALSE
  d$eligible_before_day1<-eligible
  d$included<-eligible&is.finite(d$day1_vitd_ug)
  d$analysis_weight<-ifelse(d$source_period=='I',d$diet_weight*(2/5.2),d$diet_weight*(3.2/5.2))
  d$survey_strata<-interaction(d$source_period,d$SDMVSTRA,drop=TRUE)
  d$survey_psu<-interaction(d$source_period,d$SDMVPSU,drop=TRUE)
  d$x<-d$day1_vitd_ug;d$y<-d$untreated_surface_caries_28
  write.csv(d[eligible,c('SEQN','source_period','status','raw_total','day1_vitd_ug','included')],file.path(out,'day1_exposure_audit.csv'),row.names=FALSE)
  write.csv(a,file.path(out,'day1_product_audit.csv'),row.names=FALSE)
  saveRDS(d,file.path(out,'full_design_data.rds'))
  z<-d[d$included,];stopifnot(!anyDuplicated(z$SEQN),all(z$x>=0))
  saveRDS(z,file.path(out,'day1_analysis.rds'))
  old<-readRDS(file.path(ctx$data,'analysis_complete_case.rds'))
  print(table(d$status[eligible],useNA='always'));print(quantile(z$x,c(0,.1,.5,.9,.95,.99,1)));print(table(z$source_period));print(table(z$x==0,z$y))
  s<-list(eligible_before_exposure=sum(eligible),n=nrow(z),events=sum(z$y),zero_n=sum(z$x==0),positive_n=sum(z$x>0),lost_vs_old=sum(!old$SEQN%in%z$SEQN),added_vs_old=sum(!z$SEQN%in%old$SEQN),quantiles=quantile(z$x,c(0,.1,.5,.9,.95,.99,1)),status=as.list(table(d$status[eligible])))
  write_json(s,file.path(out,'cleaning_summary.json'),pretty=TRUE,auto_unbox=TRUE);print(s)
  if (isTRUE(ctx$config$check_sample_counts)) {
    stopifnot(nrow(z) == ctx$config$expected_day1_n,
              sum(z$y) == ctx$config$expected_day1_events)
  }
  write.csv(z, file.path(out, "day1_analysis.csv"), row.names = FALSE)
  invisible(z)
}

fit_dose_response <- function(ctx) {
  out <- ctx$rcs
  require_files(file.path(out,"full_design_data.rds"))
  set.seed(20260929L)
  all<-readRDS(file.path(out,'full_design_data.rds'))
  all<-all[is.finite(all$analysis_weight)&all$analysis_weight>0&!is.na(all$survey_strata)&!is.na(all$survey_psu),]
  d<-all[all$included,]

  knots<-rcspline.eval(d$x,nk=3,knots.only=TRUE)
  bn<-paste0('r',seq_len(length(knots)-1))
  B<-rcspline.eval(all$x,knots=knots,inclx=TRUE)
  all[bn]<-B;all$x10<-all$x/10
  full<-svydesign(ids=~survey_psu,strata=~survey_strata,weights=~analysis_weight,nest=TRUE,data=all)
  ds<-subset(full,included);df<-degf(ds);crit<-qt(.975,df)
  rcs<-svyglm(reformulate(bn,'y'),design=ds,family=quasibinomial())
  po<-as.numeric(regTermTest(rcs,reformulate(bn),df=df)$p)
  pn<-as.numeric(regTermTest(rcs,reformulate(bn[-1]),df=df)$p)
  lin<-svyglm(y~x10,design=ds,family=quasibinomial())
  effect<-function(fit,L,label,unit){b<-sum(L*coef(fit));se<-sqrt(as.numeric(L%*%vcov(fit)%*%L));data.frame(model=label,unit=unit,OR=exp(b),lower=exp(b-crit*se),upper=exp(b+crit*se),p=2*pt(abs(b/se),df,lower.tail=FALSE))}
  linear<-rbind(effect(lin,c(0,.1),'Linear','per 1 ug/day'),effect(lin,c(0,1),'Linear','per 10 ug/day'))
  write.csv(linear,file.path(out,'linear_effects.csv'),row.names=FALSE)
  xmax<-as.numeric(quantile(d$x,.99));grid<-sort(unique(c(seq(0,xmax,length.out=601),knots)))
  bg<-rcspline.eval(grid,knots=knots,inclx=TRUE);br<-as.numeric(rcspline.eval(0,knots=knots,inclx=TRUE))
  delta<-sweep(bg,2,br,'-');b<-coef(rcs)[bn];V<-vcov(rcs)[bn,bn]
  eta<-as.numeric(delta%*%b);se<-sqrt(rowSums((delta%*%V)*delta))
  write.csv(data.frame(dose=grid,OR=exp(eta),low=exp(eta-crit*se),high=exp(eta+crit*se)),file.path(out,'rcs_curve.csv'),row.names=FALSE)

  lo<-sort(d$x[d$x>0])[50];hi<-as.numeric(quantile(d$x,.95))
  eligible<-function(k)sum(d$x<=k)>=100&&sum(d$x>k)>=100&&sum(d$y[d$x<=k])>=10&&sum(d$y[d$x>k])>=10&&sum(d$x>0&d$x<=k)>=50
  scan<-sort(unique(c(lo,hi,seq(lo,hi,length.out=90),d$x[d$x>=lo&d$x<=hi])))
  scan<-scan[vapply(scan,eligible,logical(1))];lo<-min(scan);hi<-max(scan)
  ux<-sort(unique(d$x));ix<-match(d$x,ux)
  aggregate_w<-function(w){sw<-as.numeric(rowsum(w,ix,reorder=FALSE)[as.character(seq_along(ux)),1]);sy<-as.numeric(rowsum(w*d$y,ix,reorder=FALSE)[as.character(seq_along(ux)),1]);list(w=sw,y=ifelse(sw>0,sy/sw,0),sy=sy)}
  softplus<-function(z)pmax(z,0)+log1p(exp(-abs(z)))
  fitk<-function(k,ag,full=FALSE){X<-cbind(1,ux/10,pmax(ux-k,0)/10);f<-suppressWarnings(glm.fit(X,ag$y,weights=ag$w,family=quasibinomial(),control=glm.control(maxit=100,epsilon=1e-9)));if(!f$converged||any(!is.finite(coef(f))))return(if(full)NULL else Inf);val<-2*sum(ag$w*softplus(f$linear.predictors)-ag$sy*f$linear.predictors);if(full)list(k=k,objective=val,coef=unname(coef(f)))else val}
  searchk<-function(w,save_scan=FALSE){ag<-aggregate_w(w/sum(w)*nrow(d));v<-vapply(scan,fitk,numeric(1),ag=ag);i<-which.min(v);candidates<-list(fitk(scan[i],ag,TRUE));loc<-which(v<=c(Inf,head(v,-1))&v<=c(tail(v,-1),Inf));loc<-head(loc[order(v[loc])],3);for(i in loc){ab<-scan[c(max(1,i-1),min(length(scan),i+1))];if(diff(ab)>1e-6){o<-optimize(fitk,ab,ag=ag,tol=.001);candidates[[length(candidates)+1]]<-fitk(o$minimum,ag,TRUE)}};ans<-candidates[[which.min(vapply(candidates,function(a)a$objective,numeric(1)))]];if(save_scan)write.csv(data.frame(k=scan,objective=v),file.path(out,'breakpoint_profile.csv'),row.names=FALSE);ans}
  best<-searchk(d$analysis_weight,TRUE);k<-best$k
  ds<-update(ds,hinge10=pmax(x-k,0)/10)
  seg<-svyglm(y~x10+hinge10,design=ds,family=quasibinomial())
  stopifnot(max(abs(unname(coef(seg))-best$coef))<1e-5)
  segfx<-rbind(effect(seg,c(0,.1,0),'Below breakpoint','per 1 ug/day'),effect(seg,c(0,.1,.1),'Above breakpoint','per 1 ug/day'),effect(seg,c(0,1,0),'Below breakpoint','per 10 ug/day'),effect(seg,c(0,1,1),'Above breakpoint','per 10 ug/day'))
  write.csv(segfx,file.path(out,'segmented_effects_conditional.csv'),row.names=FALSE)
  G<-cbind(0,grid/10,pmax(grid-k,0)/10);e<-as.numeric(G%*%coef(seg));s<-sqrt(rowSums((G%*%vcov(seg))*G))
  write.csv(data.frame(dose=grid,OR=exp(e),low=exp(e-crit*s),high=exp(e+crit*s)),file.path(out,'segmented_curve_conditional.csv'),row.names=FALSE)

  previous<-c(0,5,62.5);bp<-rcspline.eval(d$x,knots=previous,inclx=TRUE)
  dp<-update(ds,old1=bp[,1],old2=bp[,2]);fp<-svyglm(y~old1+old2,design=dp,family=quasibinomial())

  pos<-subset(ds,x>0);pk<-as.numeric(quantile(pos$variables$x,c(.1,.5,.9)))
  pb<-rcspline.eval(pos$variables$x,knots=pk,inclx=TRUE);pos<-update(pos,pos1=pb[,1],pos2=pb[,2]);pf<-svyglm(y~pos1+pos2,design=pos,family=quasibinomial())
  summary<-list(exposure='DS1TVD: day-1 24-hour supplement vitamin D',n=nrow(d),events=sum(d$y),zero_n=sum(d$x==0),positive_n=sum(d$x>0),years='2015-March 2020',adjustment='No covariates; day-1 dietary weights, strata and PSU; covariate-complete sample maintained',design_df=df,knots=knots,knots_method='Hmisc nk=3 with automatic repeated-quantile fallback; actual knots recorded',p_overall=po,p_nonlinear=pn,linear_effects=linear,breakpoint=k,search_range=c(lo,hi),conditional_segment_effects=segfx,previous_knots=previous,previous_knots_p_overall=as.numeric(regTermTest(fp,~old1+old2,df=df)$p),previous_knots_p_nonlinear=as.numeric(regTermTest(fp,~old2,df=df)$p),positive_only_knots=pk,positive_only_p_overall=as.numeric(regTermTest(pf,~pos1+pos2,df=degf(pos))$p),positive_only_p_nonlinear=as.numeric(regTermTest(pf,~pos2,df=degf(pos))$p),plot_max=xmax,fit_range=range(d$x),reference=0)
  write_json(summary,file.path(out,'summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
  print(summary);flush.console()
  boot_n<-as.integer(ctx$config$bootstrap_reps)
  if (boot_n > 0L) {

  reps<-as.svrepdesign(full,type='subbootstrap',replicates=boot_n,mse=TRUE)
  rw<-weights(reps,type='analysis')[all$included,,drop=FALSE]
  boot<-matrix(NA_real_,boot_n,3,dimnames=list(NULL,c('k','below_beta10','above_beta10')))
  for(j in seq_len(boot_n)){a<-tryCatch(searchk(rw[,j]),error=function(e)NULL);if(!is.null(a))boot[j,]<-c(a$k,a$coef[2],sum(a$coef[2:3]));if(j%%25==0){cat('Bootstrap',j,'/',boot_n,'\n');flush.console();write.csv(boot,file.path(out,'bootstrap_progress.csv'),row.names=FALSE)}}
  ok<-complete.cases(boot);stopifnot(mean(ok)>=.95)
  summary$breakpoint_bootstrap_ci<-as.numeric(quantile(boot[ok,'k'],c(.025,.975)))
  summary$bootstrap_success<-sum(ok);summary$bootstrap_n<-boot_n
  summary$bootstrap_at_lower<-sum(boot[ok,'k']<lo+.01);summary$bootstrap_at_upper<-sum(boot[ok,'k']>hi-.01)
  summary$bootstrap_slope_intervals_per1<-exp(apply(boot[ok,2:3,drop=FALSE]/10,2,quantile,c(.025,.975)))
  write.csv(boot,file.path(out,'bootstrap_estimates.csv'),row.names=FALSE)
  write_json(summary,file.path(out,'summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
  }
  saveRDS(list(rcs=rcs,linear=lin,segmented=seg,summary=summary,positive_rcs=pf,previous_knots_rcs=fp),file.path(out,'models.rds'))
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));print(summary)
  lrt <- survey::regTermTest(seg, ~hinge10, method = "LRT", df = df,
                            lrt.approximation = "satterthwaite")
  write_json(list(breakpoint = k, p_value = as.numeric(lrt$p),
                  design_df = df, conditional_on_selected_breakpoint = TRUE,
                  accounts_for_breakpoint_search = FALSE),
             file.path(out, "conditional_lrt.json"), auto_unbox = TRUE, pretty = TRUE)
  table3 <- rbind(linear[linear$unit == "per 1 ug/day", ],
                   segfx[segfx$unit == "per 1 ug/day", ])
  table3$exposure <- "Day-1 supplement vitamin D"
  table3$n <- nrow(d)
  table3$breakpoint <- k
  table3$conditional_LRT_p <- as.numeric(lrt$p)
  write.csv(table3, file.path(ctx$tables, "Table_3.csv"), row.names = FALSE)
  invisible(summary)
}

fit_correlation_lasso <- function(ctx) {
  work_dir <- ctx$data
  output_dir <- ctx$lasso
  require_files(file.path(work_dir,"analysis_complete_case.rds"))
  analysis <- readRDS(file.path(work_dir, "analysis_complete_case.rds"))

  stopifnot(
    nrow(analysis) == 3089L,
    sum(analysis$untreated_surface_caries_28 == 1L) == 592L,
    !anyNA(analysis[, c("untreated_surface_caries_28", "dietary_vitd_ug", "supplement_vitd_ug")])
  )

  predictor_names <- c(
    "dietary_vitd_ug", "supplement_vitd_ug",
    "age", "sex", "race", "education", "marital", "pir",
    "bmi", "smoking", "alcohol", "physical_activity", "floss_days",
    "dental_visit", "diabetes", "hypertension", "cvd", "ckd",
    "liver_disease", "present_teeth_28", "energy_kcal", "sugar_g",
    "cholesterol_mg", "fat_g", "fiber_g", "protein_g", "calcium_mg",
    "phosphorus_mg"
  )

  display_labels <- c(
    dietary_vitd_ug = "Dietary vitamin D",
    supplement_vitd_ug = "Supplement vitamin D",
    age = "Age",
    sex = "Gender",
    race = "Race/ethnicity",
    education = "Education",
    marital = "Marital status",
    pir = "PIR",
    bmi = "BMI",
    smoking = "Smoking",
    alcohol = "Alcohol use",
    physical_activity = "Physical activity",
    floss_days = "Flosser/interdental device use",
    dental_visit = "Last dentist visit",
    diabetes = "Diabetes",
    hypertension = "Hypertension",
    cvd = "Broad CVD",
    ckd = "CKD",
    liver_disease = "Liver disease",
    present_teeth_28 = "Present natural teeth",
    energy_kcal = "Energy intake",
    sugar_g = "Total sugar",
    cholesterol_mg = "Cholesterol",
    fat_g = "Total fat",
    fiber_g = "Dietary fiber",
    protein_g = "Protein",
    calcium_mg = "Calcium",
    phosphorus_mg = "Phosphorus"
  )

  as_code <- function(x) {
    if (is.factor(x)) return(as.integer(x))
    if (is.logical(x)) return(as.integer(x))
    as.numeric(x)
  }

  heatmap_order <- c(
    "dietary_vitd_ug", "supplement_vitd_ug",
    "energy_kcal", "sugar_g", "cholesterol_mg", "fat_g", "fiber_g",
    "protein_g", "calcium_mg", "phosphorus_mg",
    "sex", "age", "race", "education", "marital", "pir", "bmi",
    "smoking", "alcohol", "physical_activity", "floss_days", "dental_visit",
    "diabetes", "hypertension", "cvd", "ckd", "liver_disease",
    "present_teeth_28"
  )

  heatmap_data <- as.data.frame(lapply(analysis[heatmap_order], as_code))
  names(heatmap_data) <- unname(display_labels[heatmap_order])
  cor_matrix <- cor(heatmap_data, method = "spearman", use = "complete.obs")

  matrix_long <- expand.grid(
    row = rownames(cor_matrix),
    column = colnames(cor_matrix),
    stringsAsFactors = FALSE
  )
  matrix_long$row_index <- match(matrix_long$row, rownames(cor_matrix))
  matrix_long$column_index <- match(matrix_long$column, colnames(cor_matrix))
  matrix_long$correlation <- mapply(
    function(r, c) cor_matrix[r, c],
    matrix_long$row_index,
    matrix_long$column_index
  )
  matrix_long$triangle <- ifelse(
    matrix_long$row_index < matrix_long$column_index,
    "upper",
    ifelse(matrix_long$row_index > matrix_long$column_index, "lower", "diagonal")
  )
  matrix_long$fill_value <- ifelse(
    matrix_long$triangle %in% c("upper", "diagonal"),
    matrix_long$correlation,
    NA_real_
  )
  matrix_long$number_label <- ifelse(
    matrix_long$triangle == "lower",
    sprintf("%.2f", matrix_long$correlation),
    ""
  )

  axis_levels <- colnames(cor_matrix)
  matrix_long$column <- factor(matrix_long$column, levels = axis_levels)
  matrix_long$row <- factor(matrix_long$row, levels = rev(axis_levels))

  heatmap_plot <- ggplot(matrix_long, aes(x = column, y = row)) +
    geom_tile(aes(fill = fill_value), color = "#C8C8C8", linewidth = 0.28) +
    geom_text(aes(label = number_label), family = "Arial", size = 2.95, color = "black") +
    scale_fill_gradient2(
      low = "#B2182B", mid = "#FFFFFF", high = "#053061",
      midpoint = 0, limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.2),
      na.value = "#FFFFFF", name = NULL
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(clip = "off") +
    theme_void(base_family = "Arial") +
    theme(
      axis.text.x = element_text(
        angle = 90, hjust = 0, vjust = 0.5, size = 9.7,
        color = "black", margin = margin(b = 5)
      ),
      axis.text.y = element_text(
        size = 9.7, color = "black", hjust = 1,
        margin = margin(r = 5)
      ),
      legend.position = "right",
      legend.key.height = unit(4.3, "cm"),
      legend.key.width = unit(0.42, "cm"),
      legend.text = element_text(size = 8.5, color = "black"),
      plot.margin = margin(t = 18, r = 18, b = 12, l = 12)
    )

  heatmap_file <- file.path(output_dir, "Figure_3_Spearman.png")
  ggsave(
    filename = heatmap_file,
    plot = heatmap_plot,
    width = 19.0,
    height = 16.8,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  write_excel_csv(
    data.frame(Variable = rownames(cor_matrix), cor_matrix, check.names = FALSE),
    file.path(output_dir, "Spearman_correlations.csv")
  )

  predictor_data <- analysis[predictor_names]
  lasso_formula <- reformulate(predictor_names)
  terms_object <- terms(lasso_formula, data = predictor_data)
  model_matrix_full <- model.matrix(terms_object, data = predictor_data)
  assignment_full <- attr(model_matrix_full, "assign")
  term_labels <- attr(terms_object, "term.labels")

  intercept_column <- which(colnames(model_matrix_full) == "(Intercept)")
  X <- model_matrix_full[, -intercept_column, drop = FALSE]
  column_assignments <- assignment_full[-intercept_column]
  column_to_term <- setNames(term_labels[column_assignments], colnames(X))
  y <- as.integer(analysis$untreated_surface_caries_28)

  X <- scale(X, center = TRUE, scale = TRUE)

  set.seed(123)
  fold_id <- integer(length(y))
  for (class_value in sort(unique(y))) {
    indices <- which(y == class_value)
    fold_id[indices] <- sample(rep(seq_len(10L), length.out = length(indices)))
  }

  full_model <- glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = 1,
    standardize = FALSE,
    intercept = TRUE
  )

  cv_model <- cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = 1,
    lambda = full_model$lambda,
    foldid = fold_id,
    type.measure = "deviance",
    standardize = FALSE,
    intercept = TRUE,
    keep = TRUE
  )

  if (.Platform$OS.type == "windows") {
    windowsFonts(Arial = windowsFont("Arial"))
  }

  cv_file <- file.path(output_dir, "Figure_4A_LASSO_CV.png")
  png(cv_file, width = 4500, height = 3000, res = 300, bg = "white")
  par(
    family = "Arial", mar = c(5.5, 5.5, 3.4, 1.2),
    mgp = c(3.2, 1.0, 0), cex.axis = 1.18, cex.lab = 1.32,
    font.lab = 1, las = 1
  )
  plot(
    cv_model,
    main = "",
    xlab = "Log(Lambda)",
    ylab = "Deviance"
  )
  cv_usr <- par("usr")
  text(
    x = cv_usr[1] - 0.065 * diff(cv_usr[1:2]),
    y = cv_usr[4] + 0.075 * diff(cv_usr[3:4]),
    labels = "A", xpd = NA, cex = 2.15, font = 2, family = "Arial"
  )
  dev.off()

  path_colors <- colorRampPalette(brewer.pal(8, "Set2"))(ncol(X))
  path_file <- file.path(output_dir, "Figure_4B_LASSO_path.png")
  png(path_file, width = 4500, height = 3000, res = 300, bg = "white")
  par(
    family = "Arial", mar = c(5.5, 5.5, 3.4, 1.2),
    mgp = c(3.2, 1.0, 0), cex.axis = 1.18, cex.lab = 1.32,
    font.lab = 1, las = 1
  )
  plot(
    full_model,
    xvar = "lambda",
    label = FALSE,
    main = "",
    xlab = "Log(Lambda)",
    ylab = "Coefficients",
    lwd = 2,
    col = path_colors
  )
  path_usr <- par("usr")
  text(
    x = path_usr[1] - 0.065 * diff(path_usr[1:2]),
    y = path_usr[4] + 0.075 * diff(path_usr[3:4]),
    labels = "B", xpd = NA, cex = 2.15, font = 2, family = "Arial"
  )
  dev.off()

  combined_file <- file.path(output_dir, "Figure_4_LASSO.png")
  png(combined_file, width = 4500, height = 6100, res = 300, bg = "black")
  layout(matrix(c(1, 2, 3), ncol = 1), heights = c(1, 0.035, 1))
  par(
    family = "Arial", mar = c(5.5, 5.5, 3.4, 1.2),
    mgp = c(3.2, 1.0, 0), cex.axis = 1.18, cex.lab = 1.32,
    font.lab = 1, las = 1, bg = "white"
  )
  plot(cv_model, main = "", xlab = "Log(Lambda)", ylab = "Deviance")
  mtext("A", side = 3, line = 1.15, adj = 0, cex = 2.15, font = 2, family = "Arial")
  par(mar = c(0, 0, 0, 0), bg = "black")
  plot.new()
  par(
    family = "Arial", mar = c(5.5, 5.5, 3.4, 1.2),
    mgp = c(3.2, 1.0, 0), cex.axis = 1.18, cex.lab = 1.32,
    font.lab = 1, las = 1, bg = "white"
  )
  plot(
    full_model, xvar = "lambda", label = FALSE, main = "",
    xlab = "Log(Lambda)", ylab = "Coefficients",
    lwd = 2, col = path_colors
  )
  mtext("B", side = 3, line = 1.15, adj = 0, cex = 2.15, font = 2, family = "Arial")
  dev.off()

  extract_selection <- function(cv_fit, s_name) {
    coefficient_matrix <- as.matrix(coef(cv_fit, s = s_name))
    coefficient_names <- rownames(coefficient_matrix)[
      coefficient_matrix[, 1] != 0 & rownames(coefficient_matrix) != "(Intercept)"
    ]
    selected_terms <- unique(unname(column_to_term[coefficient_names]))
    selected_terms <- predictor_names[predictor_names %in% selected_terms]
    data.frame(
      rule = s_name,
      lambda = if (s_name == "lambda.min") cv_fit$lambda.min else cv_fit$lambda.1se,
      conceptual_variable = selected_terms,
      display_name = unname(display_labels[selected_terms]),
      stringsAsFactors = FALSE
    )
  }

  selection_min <- extract_selection(cv_model, "lambda.min")
  selection_1se <- extract_selection(cv_model, "lambda.1se")
  selection_all <- rbind(selection_min, selection_1se)

  write_excel_csv(selection_all, file.path(output_dir, "LASSO_selected_variables.csv"))

  coefficient_table <- do.call(
    rbind,
    lapply(c("lambda.min", "lambda.1se"), function(s_name) {
      coefficient_matrix <- as.matrix(coef(cv_model, s = s_name))
      keep <- coefficient_matrix[, 1] != 0 & rownames(coefficient_matrix) != "(Intercept)"
      coefficient_names <- rownames(coefficient_matrix)[keep]
      data.frame(
        rule = s_name,
        lambda = if (s_name == "lambda.min") cv_model$lambda.min else cv_model$lambda.1se,
        coefficient = coefficient_names,
        conceptual_variable = unname(column_to_term[coefficient_names]),
        display_name = unname(display_labels[column_to_term[coefficient_names]]),
        estimate = as.numeric(coefficient_matrix[keep, 1]),
        stringsAsFactors = FALSE
      )
    })
  )
  write_excel_csv(coefficient_table, file.path(output_dir, "LASSO_nonzero_coefficients.csv"))

  summary_lines <- c(
    sprintf("Analysis sample: n=%d; events=%d", length(y), sum(y == 1L)),
    sprintf("Candidate conceptual variables: %d", length(predictor_names)),
    sprintf("Design-matrix coefficients: %d", ncol(X)),
    "Cross-validation: stratified 10-fold; seed=123; type.measure=deviance",
    "Predictors: z-standardized before fitting; intercept retained.",
    "LASSO is unweighted, consistent with the manuscript's machine-learning analyses.",
    "",
    sprintf("lambda.min = %.9f", cv_model$lambda.min),
    sprintf("lambda.min selected conceptual variables: %d", nrow(selection_min)),
    paste(unname(selection_min$display_name), collapse = "; "),
    sprintf("lambda.min non-zero design-matrix coefficients: %d", sum(coefficient_table$rule == "lambda.min")),
    "",
    sprintf("lambda.1se = %.9f", cv_model$lambda.1se),
    sprintf("lambda.1se selected conceptual variables: %d", nrow(selection_1se)),
    paste(unname(selection_1se$display_name), collapse = "; "),
    sprintf("lambda.1se non-zero design-matrix coefficients: %d", sum(coefficient_table$rule == "lambda.1se"))
  )
  writeLines(summary_lines, file.path(output_dir, "LASSO_summary.txt"), useBytes = TRUE)

  saveRDS(
    list(
      cv_model = cv_model,
      full_model = full_model,
      X_columns = colnames(X),
      column_to_term = column_to_term,
      selection_min = selection_min,
      selection_1se = selection_1se,
      cor_matrix = cor_matrix,
      fold_id = fold_id
    ),
    file.path(output_dir, "LASSO_analysis.rds")
  )

  cat(paste(summary_lines, collapse = "\n"), "\n")
  cat("\nOutput directory:\n", output_dir, "\n", sep = "")
  save_plot_pair(heatmap_plot, file.path(output_dir, "Figure_3_Spearman"), 19, 16.8)
  grDevices::cairo_pdf(file.path(output_dir, "Figure_4_LASSO.pdf"), width = 10, height = 12)
  old_par <- par(mfrow = c(2, 1), mar = c(5, 5, 3, 1), family = "sans")
  plot(cv_model, xlab = "Log(Lambda)", ylab = "Deviance")
  mtext("A", side = 3, adj = 0, font = 2)
  plot(full_model, xvar = "lambda", xlab = "Log(Lambda)", ylab = "Coefficients",
       col = path_colors, lwd = 1.5)
  mtext("B", side = 3, adj = 0, font = 2)
  par(old_par)
  dev.off()
  invisible(selection_all)
}

fit_prediction_models <- function(ctx) {
  data_file <- file.path(ctx$data, "analysis_complete_case.rds")
  output_dir <- ctx$models
  require_files(c(data_file, file.path(ctx$lasso, "LASSO_analysis.rds")))
  split_seed <- ctx$config$split_seed
  cv_seed <- ctx$config$cv_seed
  split_ratio <- 0.80
  cv_folds <- 10L

  analysis <- readRDS(data_file)
  outcome_name <- "untreated_surface_caries_28"
  predictor_names <- c(
    "supplement_vitd_ug", "age", "sex", "race", "pir", "smoking",
    "alcohol", "floss_days", "dental_visit", "diabetes", "sugar_g",
    "fiber_g"
  )

  stopifnot(
    nrow(analysis) == 3089L,
    sum(analysis[[outcome_name]] == 1L) == 592L,
    !anyNA(analysis[c(outcome_name, predictor_names)])
  )

  selected <- readRDS(file.path(ctx$lasso, "LASSO_analysis.rds"))$selection_1se$conceptual_variable
  stopifnot(setequal(predictor_names, selected))
  y <- as.integer(analysis[[outcome_name]])
  analysis_id <- if ("SEQN" %in% names(analysis)) analysis$SEQN else seq_len(nrow(analysis))

  stratified_split <- function(y, train_fraction, seed) {
    set.seed(seed)
    train_index <- unlist(lapply(sort(unique(y)), function(class_value) {
      class_index <- which(y == class_value)
      sample(
        class_index,
        size = floor(length(class_index) * train_fraction),
        replace = FALSE
      )
    }), use.names = FALSE)
    sort(train_index)
  }

  train_index <- stratified_split(y, split_ratio, split_seed)
  test_index <- setdiff(seq_len(nrow(analysis)), train_index)

  design_formula <- reformulate(predictor_names)
  design_matrix <- model.matrix(design_formula, data = analysis)
  design_matrix <- design_matrix[, colnames(design_matrix) != "(Intercept)", drop = FALSE]
  storage.mode(design_matrix) <- "double"

  x_train_raw <- design_matrix[train_index, , drop = FALSE]
  x_test_raw <- design_matrix[test_index, , drop = FALSE]
  y_train <- y[train_index]
  y_test <- y[test_index]

  train_center <- colMeans(x_train_raw)
  train_scale <- apply(x_train_raw, 2, sd)
  train_scale[!is.finite(train_scale) | train_scale == 0] <- 1
  x_train_scaled <- sweep(sweep(x_train_raw, 2, train_center, "-"), 2, train_scale, "/")
  x_test_scaled <- sweep(sweep(x_test_raw, 2, train_center, "-"), 2, train_scale, "/")

  y_train_factor <- factor(
    ifelse(y_train == 1L, "Caries", "NoCaries"),
    levels = c("Caries", "NoCaries")
  )

  set.seed(cv_seed)
  cv_fold_id <- as.integer(createFolds(y_train_factor, k = cv_folds, list = FALSE))
  cv_validation_folds <- lapply(seq_len(cv_folds), function(k) which(cv_fold_id == k))
  cv_training_folds <- lapply(
    cv_validation_folds,
    function(validation_index) setdiff(seq_along(y_train), validation_index)
  )
  fold_names <- paste0("Fold", seq_len(cv_folds))
  names(cv_validation_folds) <- fold_names
  names(cv_training_folds) <- fold_names

  cv_control <- trainControl(
    method = "cv",
    number = cv_folds,
    index = cv_training_folds,
    indexOut = cv_validation_folds,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    allowParallel = FALSE
  )

  model_names <- c(
    "Logistic Regression", "Ridge Logistic Regression", "Elastic Net",
    "GBM", "LightGBM", "XGBoost", "CatBoost", "KNN", "SVM",
    "Random Forest", "Decision Tree", "Naive Bayes", "AdaBoost", "MLP",
    "Extra Trees"
  )

  prob_train <- setNames(vector("list", length(model_names)), model_names)
  prob_test <- setNames(vector("list", length(model_names)), model_names)
  model_objects <- setNames(vector("list", length(model_names)), model_names)
  best_parameter_rows <- list()
  tuning_details <- list()

  progress <- function(text) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), text))
    flush.console()
  }

  paste_parameters <- function(x) {
    if (is.null(x) || ncol(x) == 0L) return("None")
    paste(paste(names(x), as.character(x[1, , drop = TRUE]), sep = "="), collapse = "; ")
  }

  record_parameters <- function(model, parameters, metric, score) {
    best_parameter_rows[[length(best_parameter_rows) + 1L]] <<- data.frame(
      Model = model,
      Best_parameters = parameters,
      Selection_metric = metric,
      Selection_score = as.numeric(score),
      stringsAsFactors = FALSE
    )
  }

  record_caret_fit <- function(model, fit) {
    record_parameters(
      model,
      paste_parameters(fit$bestTune),
      "10-fold CV AUC",
      max(fit$results$ROC, na.rm = TRUE)
    )
    details <- fit$results
    details$Model <- model
    tuning_details[[model]] <<- details
  }

  progress("1/15 Logistic Regression")
  logistic_train <- data.frame(outcome = y_train, x_train_raw, check.names = TRUE)
  logistic_test <- data.frame(x_test_raw, check.names = TRUE)
  logistic_fit <- glm(outcome ~ ., data = logistic_train, family = binomial())
  prob_train[["Logistic Regression"]] <- as.numeric(predict(logistic_fit, logistic_train, type = "response"))
  prob_test[["Logistic Regression"]] <- as.numeric(predict(logistic_fit, logistic_test, type = "response"))
  model_objects[["Logistic Regression"]] <- logistic_fit
  record_parameters("Logistic Regression", "None", "Not tuned", NA_real_)

  progress("2/15 Ridge Logistic Regression")
  set.seed(cv_seed)
  ridge_fit <- cv.glmnet(
    x_train_scaled, y_train,
    family = "binomial", alpha = 0, foldid = cv_fold_id,
    type.measure = "auc", standardize = FALSE
  )
  prob_train[["Ridge Logistic Regression"]] <- as.numeric(predict(
    ridge_fit, x_train_scaled, s = "lambda.min", type = "response"
  ))
  prob_test[["Ridge Logistic Regression"]] <- as.numeric(predict(
    ridge_fit, x_test_scaled, s = "lambda.min", type = "response"
  ))
  model_objects[["Ridge Logistic Regression"]] <- ridge_fit
  record_parameters(
    "Ridge Logistic Regression",
    sprintf("alpha=0; lambda=%.8g", ridge_fit$lambda.min),
    "10-fold CV AUC",
    max(ridge_fit$cvm, na.rm = TRUE)
  )
  tuning_details[["Ridge Logistic Regression"]] <- data.frame(
    Model = "Ridge Logistic Regression",
    alpha = 0,
    lambda = ridge_fit$lambda,
    CV_AUC = ridge_fit$cvm,
    CV_AUC_SE = ridge_fit$cvsd
  )

  progress("3/15 Elastic Net")
  elastic_alpha_grid <- c(0.25, 0.50, 0.75)
  elastic_fits <- vector("list", length(elastic_alpha_grid))
  elastic_summary <- vector("list", length(elastic_alpha_grid))
  for (i in seq_along(elastic_alpha_grid)) {
    set.seed(cv_seed + i)
    fit <- cv.glmnet(
      x_train_scaled, y_train,
      family = "binomial", alpha = elastic_alpha_grid[i], foldid = cv_fold_id,
      type.measure = "auc", standardize = FALSE
    )
    elastic_fits[[i]] <- fit
    elastic_summary[[i]] <- data.frame(
      Model = "Elastic Net",
      alpha = elastic_alpha_grid[i],
      lambda = fit$lambda.min,
      CV_AUC = max(fit$cvm, na.rm = TRUE)
    )
  }
  elastic_summary <- do.call(rbind, elastic_summary)
  best_elastic_index <- which.max(elastic_summary$CV_AUC)
  elastic_fit <- elastic_fits[[best_elastic_index]]
  elastic_alpha <- elastic_summary$alpha[best_elastic_index]
  prob_train[["Elastic Net"]] <- as.numeric(predict(
    elastic_fit, x_train_scaled, s = "lambda.min", type = "response"
  ))
  prob_test[["Elastic Net"]] <- as.numeric(predict(
    elastic_fit, x_test_scaled, s = "lambda.min", type = "response"
  ))
  model_objects[["Elastic Net"]] <- elastic_fit
  record_parameters(
    "Elastic Net",
    sprintf("alpha=%.2f; lambda=%.8g", elastic_alpha, elastic_fit$lambda.min),
    "10-fold CV AUC",
    elastic_summary$CV_AUC[best_elastic_index]
  )
  tuning_details[["Elastic Net"]] <- elastic_summary

  progress("4/15 GBM")
  gbm_grid <- expand.grid(
    n.trees = c(100L, 300L, 500L),
    interaction.depth = c(1L, 2L, 3L),
    shrinkage = c(0.01, 0.05),
    n.minobsinnode = 10L
  )
  set.seed(cv_seed)
  gbm_fit <- train(
    x = x_train_raw,
    y = y_train_factor,
    method = "gbm",
    metric = "ROC",
    trControl = cv_control,
    tuneGrid = gbm_grid,
    distribution = "bernoulli",
    verbose = FALSE
  )
  prob_train[["GBM"]] <- as.numeric(predict(gbm_fit, x_train_raw, type = "prob")[, "Caries"])
  prob_test[["GBM"]] <- as.numeric(predict(gbm_fit, x_test_raw, type = "prob")[, "Caries"])
  model_objects[["GBM"]] <- gbm_fit
  record_caret_fit("GBM", gbm_fit)

  progress("5/15 LightGBM")
  lgb_train <- lgb.Dataset(x_train_raw, label = y_train, free_raw_data = FALSE)
  lgb_grid <- expand.grid(
    num_leaves = c(7L, 15L),
    learning_rate = c(0.03, 0.08),
    stringsAsFactors = FALSE
  )
  lgb_tuning <- vector("list", nrow(lgb_grid))
  lgb_params_list <- vector("list", nrow(lgb_grid))
  for (i in seq_len(nrow(lgb_grid))) {
    params <- list(
      objective = "binary",
      metric = "auc",
      num_leaves = lgb_grid$num_leaves[i],
      learning_rate = lgb_grid$learning_rate[i],
      min_data_in_leaf = 20L,
      feature_fraction = 0.8,
      bagging_fraction = 0.8,
      bagging_freq = 1L,
      lambda_l2 = 1,
      num_threads = 1L,
      verbosity = -1L,
      seed = cv_seed,
      deterministic = TRUE,
      force_col_wise = TRUE
    )
    set.seed(cv_seed + i)
    cv_fit <- lgb.cv(
      params = params,
      data = lgb_train,
      nrounds = 500L,
      nfold = cv_folds,
      folds = cv_validation_folds,
      verbose = -1L,
      early_stopping_rounds = 30L,
      eval_train_metric = FALSE
    )
    auc_values <- as.numeric(unlist(cv_fit$record_evals$valid$auc$eval))
    lgb_tuning[[i]] <- data.frame(
      Model = "LightGBM",
      num_leaves = lgb_grid$num_leaves[i],
      learning_rate = lgb_grid$learning_rate[i],
      best_iteration = cv_fit$best_iter,
      CV_AUC = max(auc_values, na.rm = TRUE)
    )
    lgb_params_list[[i]] <- params
  }
  lgb_tuning <- do.call(rbind, lgb_tuning)
  best_lgb_index <- which.max(lgb_tuning$CV_AUC)
  lgb_fit <- lgb.train(
    params = lgb_params_list[[best_lgb_index]],
    data = lgb_train,
    nrounds = lgb_tuning$best_iteration[best_lgb_index],
    verbose = -1L
  )
  prob_train[["LightGBM"]] <- as.numeric(predict(lgb_fit, x_train_raw))
  prob_test[["LightGBM"]] <- as.numeric(predict(lgb_fit, x_test_raw))
  model_objects[["LightGBM"]] <- lgb_fit
  record_parameters(
    "LightGBM",
    sprintf(
      "num_leaves=%d; learning_rate=%.2f; iterations=%d",
      lgb_tuning$num_leaves[best_lgb_index],
      lgb_tuning$learning_rate[best_lgb_index],
      lgb_tuning$best_iteration[best_lgb_index]
    ),
    "10-fold CV AUC",
    lgb_tuning$CV_AUC[best_lgb_index]
  )
  tuning_details[["LightGBM"]] <- lgb_tuning

  progress("6/15 XGBoost")
  xgb_train <- xgb.DMatrix(x_train_raw, label = y_train)
  xgb_test <- xgb.DMatrix(x_test_raw, label = y_test)
  xgb_grid <- expand.grid(
    max_depth = c(2L, 4L),
    eta = c(0.03, 0.10),
    stringsAsFactors = FALSE
  )
  xgb_tuning <- vector("list", nrow(xgb_grid))
  xgb_params_list <- vector("list", nrow(xgb_grid))
  for (i in seq_len(nrow(xgb_grid))) {
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = xgb_grid$max_depth[i],
      eta = xgb_grid$eta[i],
      min_child_weight = 1,
      subsample = 0.8,
      colsample_bytree = 0.8,
      lambda = 1,
      nthread = 1L,
      seed = cv_seed
    )
    set.seed(cv_seed + i)
    cv_fit <- xgb.cv(
      params = params,
      data = xgb_train,
      nrounds = 500L,
      nfold = cv_folds,
      folds = cv_validation_folds,
      verbose = 0,
      maximize = TRUE,
      early_stopping_rounds = 30L
    )
    best_iteration <- which.max(cv_fit$evaluation_log$test_auc_mean)
    xgb_tuning[[i]] <- data.frame(
      Model = "XGBoost",
      max_depth = xgb_grid$max_depth[i],
      eta = xgb_grid$eta[i],
      best_iteration = best_iteration,
      CV_AUC = cv_fit$evaluation_log$test_auc_mean[best_iteration]
    )
    xgb_params_list[[i]] <- params
  }
  xgb_tuning <- do.call(rbind, xgb_tuning)
  best_xgb_index <- which.max(xgb_tuning$CV_AUC)
  xgb_fit <- xgb.train(
    params = xgb_params_list[[best_xgb_index]],
    data = xgb_train,
    nrounds = xgb_tuning$best_iteration[best_xgb_index],
    verbose = 0
  )
  prob_train[["XGBoost"]] <- as.numeric(predict(xgb_fit, xgb_train))
  prob_test[["XGBoost"]] <- as.numeric(predict(xgb_fit, xgb_test))
  model_objects[["XGBoost"]] <- xgb_fit
  record_parameters(
    "XGBoost",
    sprintf(
      "max_depth=%d; eta=%.2f; iterations=%d",
      xgb_tuning$max_depth[best_xgb_index],
      xgb_tuning$eta[best_xgb_index],
      xgb_tuning$best_iteration[best_xgb_index]
    ),
    "10-fold CV AUC",
    xgb_tuning$CV_AUC[best_xgb_index]
  )
  tuning_details[["XGBoost"]] <- xgb_tuning

  progress("7/15 CatBoost")
  cat_train_pool <- catboost.load_pool(x_train_raw, label = y_train)
  cat_test_pool <- catboost.load_pool(x_test_raw)
  cat_grid <- expand.grid(
    depth = c(4L, 6L),
    learning_rate = c(0.03, 0.08),
    stringsAsFactors = FALSE
  )
  cat_tuning <- vector("list", nrow(cat_grid))
  cat_params_list <- vector("list", nrow(cat_grid))
  for (i in seq_len(nrow(cat_grid))) {
    params <- list(
      loss_function = "Logloss",
      eval_metric = "AUC",
      iterations = 500L,
      depth = cat_grid$depth[i],
      learning_rate = cat_grid$learning_rate[i],
      l2_leaf_reg = 3,
      random_seed = cv_seed,
      thread_count = 1L,
      logging_level = "Silent",
      allow_writing_files = FALSE
    )
    set.seed(cv_seed + i)
    cv_result <- catboost.cv(
      pool = cat_train_pool,
      params = params,
      fold_count = cv_folds,
      partition_random_seed = cv_seed,
      shuffle = TRUE,
      stratified = TRUE,
      early_stopping_rounds = 30L
    )
    best_iteration <- which.max(cv_result$test.AUC.mean)
    cat_tuning[[i]] <- data.frame(
      Model = "CatBoost",
      depth = cat_grid$depth[i],
      learning_rate = cat_grid$learning_rate[i],
      best_iteration = best_iteration,
      CV_AUC = cv_result$test.AUC.mean[best_iteration]
    )
    params$iterations <- best_iteration
    cat_params_list[[i]] <- params
  }
  cat_tuning <- do.call(rbind, cat_tuning)
  best_cat_index <- which.max(cat_tuning$CV_AUC)
  cat_fit <- catboost.train(cat_train_pool, params = cat_params_list[[best_cat_index]])
  prob_train[["CatBoost"]] <- as.numeric(catboost.predict(
    cat_fit, cat_train_pool, prediction_type = "Probability"
  ))
  prob_test[["CatBoost"]] <- as.numeric(catboost.predict(
    cat_fit, cat_test_pool, prediction_type = "Probability"
  ))
  model_objects[["CatBoost"]] <- cat_fit
  record_parameters(
    "CatBoost",
    sprintf(
      "depth=%d; learning_rate=%.2f; iterations=%d",
      cat_tuning$depth[best_cat_index],
      cat_tuning$learning_rate[best_cat_index],
      cat_tuning$best_iteration[best_cat_index]
    ),
    "10-fold CV AUC",
    cat_tuning$CV_AUC[best_cat_index]
  )
  tuning_details[["CatBoost"]] <- cat_tuning

  progress("8/15 KNN")
  set.seed(cv_seed)
  knn_fit <- train(
    x = x_train_scaled,
    y = y_train_factor,
    method = "knn",
    metric = "ROC",
    trControl = cv_control,
    tuneGrid = data.frame(k = c(5L, 9L, 15L, 25L, 35L))
  )
  prob_train[["KNN"]] <- as.numeric(predict(knn_fit, x_train_scaled, type = "prob")[, "Caries"])
  prob_test[["KNN"]] <- as.numeric(predict(knn_fit, x_test_scaled, type = "prob")[, "Caries"])
  model_objects[["KNN"]] <- knn_fit
  record_caret_fit("KNN", knn_fit)

  progress("9/15 SVM")
  set.seed(cv_seed)
  svm_fit <- train(
    x = x_train_scaled,
    y = y_train_factor,
    method = "svmRadial",
    metric = "ROC",
    trControl = cv_control,
    tuneLength = 5L
  )
  prob_train[["SVM"]] <- as.numeric(predict(svm_fit, x_train_scaled, type = "prob")[, "Caries"])
  prob_test[["SVM"]] <- as.numeric(predict(svm_fit, x_test_scaled, type = "prob")[, "Caries"])
  model_objects[["SVM"]] <- svm_fit
  record_caret_fit("SVM", svm_fit)

  progress("10/15 Random Forest")
  candidate_mtry <- unique(pmax(1L, pmin(
    ncol(x_train_raw),
    c(floor(sqrt(ncol(x_train_raw))), floor(ncol(x_train_raw) / 3), floor(ncol(x_train_raw) / 2))
  )))
  rf_grid <- expand.grid(
    mtry = candidate_mtry,
    splitrule = "gini",
    min.node.size = c(1L, 5L, 10L),
    stringsAsFactors = FALSE
  )
  set.seed(cv_seed)
  rf_fit <- train(
    x = x_train_raw,
    y = y_train_factor,
    method = "ranger",
    metric = "ROC",
    trControl = cv_control,
    tuneGrid = rf_grid,
    num.trees = 1000L,
    importance = "none",
    num.threads = 1L
  )
  prob_train[["Random Forest"]] <- as.numeric(predict(rf_fit, x_train_raw, type = "prob")[, "Caries"])
  prob_test[["Random Forest"]] <- as.numeric(predict(rf_fit, x_test_raw, type = "prob")[, "Caries"])
  model_objects[["Random Forest"]] <- rf_fit
  record_caret_fit("Random Forest", rf_fit)

  progress("11/15 Decision Tree")
  set.seed(cv_seed)
  dt_fit <- train(
    x = x_train_raw,
    y = y_train_factor,
    method = "rpart",
    metric = "ROC",
    trControl = cv_control,
    tuneGrid = data.frame(cp = c(0.001, 0.005, 0.01, 0.02, 0.05))
  )
  prob_train[["Decision Tree"]] <- as.numeric(predict(dt_fit, x_train_raw, type = "prob")[, "Caries"])
  prob_test[["Decision Tree"]] <- as.numeric(predict(dt_fit, x_test_raw, type = "prob")[, "Caries"])
  model_objects[["Decision Tree"]] <- dt_fit
  record_caret_fit("Decision Tree", dt_fit)

  progress("12/15 Naive Bayes")
  nb_grid <- expand.grid(
    laplace = c(0, 1),
    usekernel = c(TRUE, FALSE),
    adjust = c(0.5, 1.0),
    stringsAsFactors = FALSE
  )
  set.seed(cv_seed)
  nb_fit <- train(
    x = x_train_scaled,
    y = y_train_factor,
    method = "naive_bayes",
    metric = "ROC",
    trControl = cv_control,
    tuneGrid = nb_grid
  )
  prob_train[["Naive Bayes"]] <- as.numeric(predict(nb_fit, x_train_scaled, type = "prob")[, "Caries"])
  prob_test[["Naive Bayes"]] <- as.numeric(predict(nb_fit, x_test_scaled, type = "prob")[, "Caries"])
  model_objects[["Naive Bayes"]] <- nb_fit
  record_caret_fit("Naive Bayes", nb_fit)

  progress("13/15 AdaBoost")
  adaboost_grid <- expand.grid(
    interaction.depth = c(1L, 2L),
    shrinkage = c(0.01, 0.05),
    n.minobsinnode = 10L,
    stringsAsFactors = FALSE
  )
  adaboost_train_data <- data.frame(y = y_train, x_train_raw, check.names = TRUE)
  adaboost_test_data <- data.frame(x_test_raw, check.names = TRUE)
  adaboost_models <- vector("list", nrow(adaboost_grid))
  adaboost_tuning <- vector("list", nrow(adaboost_grid))
  for (i in seq_len(nrow(adaboost_grid))) {
    set.seed(cv_seed + i)
    fit <- gbm(
      y ~ ., data = adaboost_train_data,
      distribution = "adaboost",
      n.trees = 1500L,
      interaction.depth = adaboost_grid$interaction.depth[i],
      shrinkage = adaboost_grid$shrinkage[i],
      n.minobsinnode = adaboost_grid$n.minobsinnode[i],
      bag.fraction = 1,
      train.fraction = 1,
      cv.folds = cv_folds,
      keep.data = FALSE,
      verbose = FALSE,
      n.cores = 1L
    )
    best_iteration <- gbm.perf(fit, method = "cv", plot.it = FALSE)
    adaboost_models[[i]] <- fit
    adaboost_tuning[[i]] <- data.frame(
      Model = "AdaBoost",
      interaction.depth = adaboost_grid$interaction.depth[i],
      shrinkage = adaboost_grid$shrinkage[i],
      n.minobsinnode = adaboost_grid$n.minobsinnode[i],
      best_iteration = best_iteration,
      minimum_CV_error = min(fit$cv.error, na.rm = TRUE)
    )
  }
  adaboost_tuning <- do.call(rbind, adaboost_tuning)
  best_adaboost_index <- which.min(adaboost_tuning$minimum_CV_error)
  adaboost_fit <- adaboost_models[[best_adaboost_index]]
  adaboost_n_trees <- adaboost_tuning$best_iteration[best_adaboost_index]
  prob_train[["AdaBoost"]] <- as.numeric(predict(
    adaboost_fit,
    newdata = adaboost_train_data[, -1, drop = FALSE],
    n.trees = adaboost_n_trees,
    type = "response"
  ))
  prob_test[["AdaBoost"]] <- as.numeric(predict(
    adaboost_fit,
    newdata = adaboost_test_data,
    n.trees = adaboost_n_trees,
    type = "response"
  ))
  model_objects[["AdaBoost"]] <- adaboost_fit
  record_parameters(
    "AdaBoost",
    sprintf(
      "interaction.depth=%d; shrinkage=%.2f; iterations=%d",
      adaboost_tuning$interaction.depth[best_adaboost_index],
      adaboost_tuning$shrinkage[best_adaboost_index],
      adaboost_n_trees
    ),
    "10-fold CV error",
    adaboost_tuning$minimum_CV_error[best_adaboost_index]
  )
  tuning_details[["AdaBoost"]] <- adaboost_tuning

  progress("14/15 MLP")
  mlp_grid <- expand.grid(
    size = c(3L, 5L, 8L, 12L),
    decay = c(0, 0.001, 0.01, 0.1)
  )
  set.seed(cv_seed)
  mlp_fit <- train(
    x = x_train_scaled,
    y = y_train_factor,
    method = "nnet",
    metric = "ROC",
    trControl = cv_control,
    tuneGrid = mlp_grid,
    trace = FALSE,
    maxit = 1000L,
    MaxNWts = 10000L
  )
  prob_train[["MLP"]] <- as.numeric(predict(mlp_fit, x_train_scaled, type = "prob")[, "Caries"])
  prob_test[["MLP"]] <- as.numeric(predict(mlp_fit, x_test_scaled, type = "prob")[, "Caries"])
  model_objects[["MLP"]] <- mlp_fit
  record_caret_fit("MLP", mlp_fit)

  progress("15/15 Extra Trees")
  extra_grid <- expand.grid(
    mtry = candidate_mtry,
    splitrule = "extratrees",
    min.node.size = c(1L, 5L, 10L),
    stringsAsFactors = FALSE
  )
  set.seed(cv_seed)
  extra_fit <- train(
    x = x_train_raw,
    y = y_train_factor,
    method = "ranger",
    metric = "ROC",
    trControl = cv_control,
    tuneGrid = extra_grid,
    num.trees = 1000L,
    num.random.splits = 10L,
    importance = "none",
    num.threads = 1L
  )
  prob_train[["Extra Trees"]] <- as.numeric(predict(extra_fit, x_train_raw, type = "prob")[, "Caries"])
  prob_test[["Extra Trees"]] <- as.numeric(predict(extra_fit, x_test_raw, type = "prob")[, "Caries"])
  model_objects[["Extra Trees"]] <- extra_fit
  record_caret_fit("Extra Trees", extra_fit)

  clip_probability <- function(x) pmin(pmax(as.numeric(x), 1e-6), 1 - 1e-6)
  for (model in model_names) {
    stopifnot(
      length(prob_train[[model]]) == length(y_train),
      length(prob_test[[model]]) == length(y_test),
      all(is.finite(prob_train[[model]])),
      all(is.finite(prob_test[[model]]))
    )
    prob_train[[model]] <- clip_probability(prob_train[[model]])
    prob_test[[model]] <- clip_probability(prob_test[[model]])
  }

  predictions_train <- data.frame(
    id = analysis_id[train_index],
    row_index = train_index,
    Dataset = "Training",
    Outcome = y_train,
    check.names = FALSE
  )
  predictions_test <- data.frame(
    id = analysis_id[test_index],
    row_index = test_index,
    Dataset = "Testing",
    Outcome = y_test,
    check.names = FALSE
  )
  for (model in model_names) {
    predictions_train[[model]] <- prob_train[[model]]
    predictions_test[[model]] <- prob_test[[model]]
  }
  predictions <- rbind(predictions_train, predictions_test)

  get_roc <- function(y_true, probability) {
    roc(
      response = y_true,
      predictor = probability,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )
  }

  calculate_confusion <- function(y_true, probability, threshold, model, dataset) {
    predicted <- as.integer(probability >= threshold)
    tn <- sum(y_true == 0L & predicted == 0L)
    fp <- sum(y_true == 0L & predicted == 1L)
    fn <- sum(y_true == 1L & predicted == 0L)
    tp <- sum(y_true == 1L & predicted == 1L)
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
    accuracy <- (tp + tn) / length(y_true)
    balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
    f1 <- if (is.finite(precision) && is.finite(sensitivity) && (precision + sensitivity) > 0) {
      2 * precision * sensitivity / (precision + sensitivity)
    } else {
      NA_real_
    }

    metrics <- data.frame(
      Model = model,
      Dataset = dataset,
      N = length(y_true),
      Events = sum(y_true == 1L),
      Threshold = threshold,
      TN = tn,
      FP = fp,
      FN = fn,
      TP = tp,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Precision = precision,
      NPV = npv,
      Accuracy = accuracy,
      Balanced_accuracy = balanced_accuracy,
      F1 = f1,
      stringsAsFactors = FALSE
    )

    counts <- data.frame(
      Model = model,
      Dataset = dataset,
      Actual = c("No caries", "No caries", "Caries", "Caries"),
      Predicted = c("No caries", "Caries", "No caries", "Caries"),
      Count = c(tn, fp, fn, tp),
      stringsAsFactors = FALSE
    )
    counts$Row_total <- ave(counts$Count, counts$Actual, FUN = sum)
    counts$Row_percent <- counts$Count / counts$Row_total

    list(metrics = metrics, counts = counts)
  }

  training_thresholds <- setNames(numeric(length(model_names)), model_names)
  performance_rows <- list()
  confusion_rows <- list()
  roc_curve_rows <- list()

  for (model in model_names) {
    train_roc <- get_roc(y_train, prob_train[[model]])
    test_roc <- get_roc(y_test, prob_test[[model]])

    threshold_result <- coords(
      train_roc,
      x = "best",
      best.method = "youden",
      ret = "threshold",
      transpose = FALSE
    )
    threshold <- as.numeric(threshold_result)[1]
    training_thresholds[[model]] <- threshold

    for (dataset_name in c("Training", "Testing")) {
      y_true <- if (dataset_name == "Training") y_train else y_test
      probability <- if (dataset_name == "Training") prob_train[[model]] else prob_test[[model]]
      roc_object <- if (dataset_name == "Training") train_roc else test_roc
      auc_ci <- as.numeric(ci.auc(roc_object, method = "delong"))

      confusion_result <- calculate_confusion(
        y_true,
        probability,
        threshold,
        model,
        dataset_name
      )
      metric_row <- confusion_result$metrics
      metric_row$AUC <- as.numeric(auc(roc_object))
      metric_row$AUC_CI_low <- auc_ci[1]
      metric_row$AUC_CI_high <- auc_ci[3]
      performance_rows[[length(performance_rows) + 1L]] <- metric_row
      confusion_rows[[length(confusion_rows) + 1L]] <- confusion_result$counts

      roc_coordinates <- coords(
        roc_object,
        x = "all",
        ret = c("specificity", "sensitivity"),
        transpose = FALSE
      )
      roc_curve_rows[[length(roc_curve_rows) + 1L]] <- data.frame(
        Model = model,
        Dataset = dataset_name,
        One_minus_specificity = 1 - roc_coordinates$specificity,
        Sensitivity = roc_coordinates$sensitivity,
        stringsAsFactors = FALSE
      )
    }
  }

  performance <- do.call(rbind, performance_rows)
  confusion_long <- do.call(rbind, confusion_rows)
  roc_curves <- do.call(rbind, roc_curve_rows)
  performance$Model <- factor(performance$Model, levels = model_names)
  performance$Dataset <- factor(performance$Dataset, levels = c("Training", "Testing"))
  performance <- performance[order(performance$Dataset, performance$Model), ]
  performance$Model <- as.character(performance$Model)
  performance$Dataset <- as.character(performance$Dataset)

  decision_curve_data <- function(y_true, probability, model, dataset, thresholds) {
    n <- length(y_true)
    rows <- lapply(thresholds, function(threshold) {
      predicted <- probability >= threshold
      tp <- sum(predicted & y_true == 1L)
      fp <- sum(predicted & y_true == 0L)
      data.frame(
        Model = model,
        Dataset = dataset,
        Threshold = threshold,
        NetBenefit = tp / n - fp / n * threshold / (1 - threshold),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }

  dca_thresholds <- seq(0.01, 0.80, by = 0.01)
  dca_rows <- list()
  for (dataset_name in c("Training", "Testing")) {
    y_true <- if (dataset_name == "Training") y_train else y_test
    prevalence <- mean(y_true == 1L)
    for (model in model_names) {
      probability <- if (dataset_name == "Training") prob_train[[model]] else prob_test[[model]]
      dca_rows[[length(dca_rows) + 1L]] <- decision_curve_data(
        y_true,
        probability,
        model,
        dataset_name,
        dca_thresholds
      )
    }
    dca_rows[[length(dca_rows) + 1L]] <- data.frame(
      Model = "Treat all",
      Dataset = dataset_name,
      Threshold = dca_thresholds,
      NetBenefit = prevalence - (1 - prevalence) * dca_thresholds / (1 - dca_thresholds)
    )
    dca_rows[[length(dca_rows) + 1L]] <- data.frame(
      Model = "Treat none",
      Dataset = dataset_name,
      Threshold = dca_thresholds,
      NetBenefit = 0
    )
  }
  dca_curves <- do.call(rbind, dca_rows)

  split_summary <- do.call(rbind, lapply(
    list(Full = seq_len(nrow(analysis)), Training = train_index, Testing = test_index),
    function(index) {
      data.frame(
        N = length(index),
        Events = sum(y[index] == 1L),
        Event_percent = mean(y[index] == 1L),
        Non_events = sum(y[index] == 0L),
        stringsAsFactors = FALSE
      )
    }
  ))
  split_summary$Dataset <- rownames(split_summary)
  rownames(split_summary) <- NULL
  split_summary <- split_summary[, c("Dataset", "N", "Events", "Event_percent", "Non_events")]

  best_parameters <- do.call(rbind, best_parameter_rows)
  best_parameters$Model <- factor(best_parameters$Model, levels = model_names)
  best_parameters <- best_parameters[order(best_parameters$Model), ]
  best_parameters$Model <- as.character(best_parameters$Model)

  test_rank <- performance[performance$Dataset == "Testing", ]
  test_rank <- test_rank[order(-test_rank$AUC), ]
  test_rank$AUC_rank <- seq_len(nrow(test_rank))
  test_rank <- test_rank[, c(
    "AUC_rank", "Model", "AUC", "AUC_CI_low", "AUC_CI_high",
    "Sensitivity", "Specificity", "Precision", "NPV", "Accuracy",
    "Balanced_accuracy", "F1", "Threshold", "TN", "FP", "FN", "TP"
  )]

  auc_gap <- merge(
    performance[performance$Dataset == "Training", c("Model", "AUC")],
    performance[performance$Dataset == "Testing", c("Model", "AUC")],
    by = "Model",
    suffixes = c("_Training", "_Testing")
  )
  auc_gap$AUC_gap <- auc_gap$AUC_Training - auc_gap$AUC_Testing
  auc_gap <- auc_gap[order(-auc_gap$AUC_gap), ]

  model_colors <- c(
    "Logistic Regression" = "#E44B33",
    "Ridge Logistic Regression" = "#2F7DB2",
    "Elastic Net" = "#E48839",
    "GBM" = "#1C79BA",
    "LightGBM" = "#2AAB5B",
    "XGBoost" = "#479CD4",
    "CatBoost" = "#DBA3DE",
    "KNN" = "#C5A590",
    "SVM" = "#9F7C5E",
    "Random Forest" = "#A05EBA",
    "Decision Tree" = "#43CB7C",
    "Naive Bayes" = "#2CBD9F",
    "AdaBoost" = "#D95F02",
    "MLP" = "#6A5ACD",
    "Extra Trees" = "#008B8B"
  )
  model_linetypes <- setNames(
    rep(c("solid", "22", "42", "13", "longdash"), length.out = length(model_names)),
    model_names
  )

  base_theme <- theme_bw(base_family = "Arial", base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#ECECEC", linewidth = 0.35),
      panel.border = element_blank(),
      axis.line = element_line(color = "#333333", linewidth = 0.5),
      axis.title = element_text(face = "bold", color = "#111111"),
      axis.text = element_text(color = "#111111"),
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.tag = element_text(face = "bold", size = 17),
      plot.margin = margin(8, 10, 8, 8)
    )

  make_roc_plot <- function(dataset_name, include_auc = TRUE, title = NULL, tag = NULL) {
    plot_data <- roc_curves[roc_curves$Dataset == dataset_name, , drop = FALSE]
    auc_data <- performance[performance$Dataset == dataset_name, c("Model", "AUC")]
    auc_lookup <- setNames(auc_data$AUC, auc_data$Model)
    labels <- if (include_auc) {
      setNames(
        sprintf("%s (AUC=%.3f)", model_names, auc_lookup[model_names]),
        model_names
      )
    } else {
      setNames(model_names, model_names)
    }

    ggplot(plot_data, aes(
      x = One_minus_specificity,
      y = Sensitivity,
      color = Model,
      linetype = Model
    )) +
      geom_abline(
        intercept = 0,
        slope = 1,
        color = "#A0A0A0",
        linetype = "dashed",
        linewidth = 0.65
      ) +
      geom_line(linewidth = 0.95) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_color_manual(values = model_colors, breaks = model_names, labels = labels[model_names]) +
      scale_linetype_manual(values = model_linetypes, breaks = model_names, labels = labels[model_names]) +
      labs(
        x = "1-Specificity",
        y = "Sensitivity",
        title = title,
        tag = tag
      ) +
      guides(
        color = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(linewidth = 1.25)),
        linetype = "none"
      ) +
      base_theme +
      theme(legend.position = "bottom")
  }

  make_dca_plot <- function(dataset_name, title = NULL, tag = NULL) {
    plot_data <- dca_curves[dca_curves$Dataset == dataset_name, , drop = FALSE]
    line_order <- c(model_names, "Treat all", "Treat none")
    dca_colors <- c(model_colors, "Treat all" = "#7F7F7F", "Treat none" = "#111111")
    dca_linetypes <- c(model_linetypes, "Treat all" = "dashed", "Treat none" = "dotted")

    ggplot(plot_data, aes(
      x = Threshold,
      y = NetBenefit,
      color = Model,
      linetype = Model
    )) +
      geom_hline(yintercept = 0, color = "#8A8A8A", linewidth = 0.45) +
      geom_line(linewidth = 0.90) +
      coord_cartesian(xlim = c(0, 0.80), ylim = c(-0.10, 0.22), expand = FALSE) +
      scale_color_manual(values = dca_colors, breaks = line_order) +
      scale_linetype_manual(values = dca_linetypes, breaks = line_order) +
      labs(
        x = "Threshold Probability",
        y = "Net Benefit",
        title = title,
        tag = tag
      ) +
      guides(
        color = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(linewidth = 1.25)),
        linetype = "none"
      ) +
      base_theme +
      theme(legend.position = "bottom")
  }

  p_roc_train <- make_roc_plot("Training", TRUE, "Training Set")
  p_roc_test <- make_roc_plot("Testing", TRUE, "Testing Set")
  p_dca_train <- make_dca_plot("Training", "Training Set")
  p_dca_test <- make_dca_plot("Testing", "Testing Set")

  ggsave(
    file.path(output_dir, "Figure_5A_ROC_training.png"),
    p_roc_train,
    width = 13.5,
    height = 9.0,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, "Figure_5B_ROC_testing.png"),
    p_roc_test,
    width = 13.5,
    height = 9.0,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, "Figure_5C_DCA_training.png"),
    p_dca_train,
    width = 13.5,
    height = 9.0,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, "Figure_5D_DCA_testing.png"),
    p_dca_test,
    width = 13.5,
    height = 9.0,
    dpi = 300,
    bg = "white"
  )

  p_roc_train_compact <- make_roc_plot("Training", FALSE, NULL, "A") +
    theme(legend.position = "none")
  p_roc_test_compact <- make_roc_plot("Testing", FALSE, NULL, "B") +
    theme(legend.position = "none")
  p_dca_train_compact <- make_dca_plot("Training", NULL, "C") +
    theme(legend.position = "none")
  p_dca_test_compact <- make_dca_plot("Testing", NULL, "D") +
    theme(legend.position = "none")

  four_panel <- (p_roc_train_compact | p_roc_test_compact) /
    (p_dca_train_compact | p_dca_test_compact)
  ggsave(
    file.path(output_dir, "Figure_5_ROC_DCA.png"),
    four_panel,
    width = 15.5,
    height = 13.5,
    dpi = 300,
    bg = "white"
  )

  make_confusion_plot <- function(dataset_name) {
    plot_data <- confusion_long[confusion_long$Dataset == dataset_name, , drop = FALSE]
    plot_data$Model <- factor(plot_data$Model, levels = model_names)
    plot_data$Actual <- factor(plot_data$Actual, levels = c("Caries", "No caries"))
    plot_data$Predicted <- factor(plot_data$Predicted, levels = c("No caries", "Caries"))
    plot_data$Label <- sprintf("%d\n(%.1f%%)", plot_data$Count, 100 * plot_data$Row_percent)

    ggplot(plot_data, aes(x = Predicted, y = Actual, fill = Row_percent)) +
      geom_tile(color = "white", linewidth = 1.0) +
      geom_text(aes(label = Label), family = "Arial", size = 3.4, color = "#111111") +
      facet_wrap(~Model, ncol = 5) +
      scale_fill_gradient(
        low = "#F7FBFF",
        high = "#2166AC",
        limits = c(0, 1),
        labels = scales::percent
      ) +
      coord_equal() +
      labs(
        x = "Predicted class",
        y = "Actual class",
        fill = "Row %",
        title = paste(dataset_name, "Set")
      ) +
      theme_bw(base_family = "Arial", base_size = 11) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "#F0F0F0", color = NA),
        strip.text = element_text(face = "bold", size = 9.5),
        axis.text = element_text(color = "#111111", size = 9),
        axis.title = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.margin = margin(8, 10, 8, 8)
      )
  }

  confusion_train_figure <- make_confusion_plot("Training")
  confusion_test_figure <- make_confusion_plot("Testing")
  ggsave(
    file.path(output_dir, "Confusion_matrices_training.png"),
    confusion_train_figure,
    width = 16.0,
    height = 10.0,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, "Confusion_matrices_testing.png"),
    confusion_test_figure,
    width = 16.0,
    height = 10.0,
    dpi = 300,
    bg = "white"
  )

  predictor_dictionary <- data.frame(
    Variable = predictor_names,
    Description = c(
      "Supplement vitamin D, ug/day",
      "Age, years",
      "Sex",
      "Race/ethnicity",
      "Poverty income ratio",
      "Smoking status",
      "Alcohol use",
      "Floss/interdental device use, days/week",
      "Last dental visit",
      "Diabetes",
      "Total sugar intake, g/day",
      "Dietary fiber intake, g/day"
    ),
    stringsAsFactors = FALSE
  )

  quality_checks <- data.frame(
    Check = c(
      "Full sample N",
      "Full sample events",
      "Train and test rows cover full sample",
      "No train-test overlap",
      "Exactly 15 models",
      "No missing prediction probabilities",
      "All probabilities within [0,1]",
      "All test thresholds copied from training"
    ),
    Result = c(
      nrow(analysis),
      sum(y == 1L),
      length(unique(c(train_index, test_index))) == nrow(analysis),
      length(intersect(train_index, test_index)) == 0L,
      length(model_names) == 15L,
      !anyNA(predictions[model_names]),
      all(as.matrix(predictions[model_names]) >= 0 & as.matrix(predictions[model_names]) <= 1),
      all(vapply(model_names, function(model) {
        rows <- performance[performance$Model == model, ]
        length(unique(rows$Threshold)) == 1L
      }, logical(1)))
    ),
    stringsAsFactors = FALSE
  )

  write.csv(performance, file.path(output_dir, "model_performance.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(test_rank, file.path(output_dir, "testing_model_ranking.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(confusion_long, file.path(output_dir, "confusion_counts.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(roc_curves, file.path(output_dir, "ROC_coordinates.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(dca_curves, file.path(output_dir, "DCA_coordinates.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(predictions, file.path(output_dir, "individual_predictions.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(split_summary, file.path(output_dir, "split_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(best_parameters, file.path(output_dir, "best_hyperparameters.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(auc_gap, file.path(output_dir, "training_testing_AUC_gap.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(quality_checks, file.path(output_dir, "model_quality_checks.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  workbook_sheets <- list(
    Performance = performance,
    Test_rank = test_rank,
    Confusion_counts = confusion_long,
    Split_summary = split_summary,
    Best_parameters = best_parameters,
    AUC_gap = auc_gap,
    Quality_checks = quality_checks,
    Predictors = predictor_dictionary,
    ROC_coordinates = roc_curves,
    DCA_coordinates = dca_curves,
    Predictions = predictions
  )

  writexl::write_xlsx(
    workbook_sheets,
    path = file.path(output_dir, "model_results.xlsx")
  )

  saveRDS(
    list(
      specification = list(
        data_file = data_file,
        outcome = outcome_name,
        outcome_definition = "At least one untreated coronal-caries tooth surface among 28 teeth",
        predictors = predictor_names,
        model_names = model_names,
        split_ratio = split_ratio,
        split_seed = split_seed,
        cv_seed = cv_seed,
        cv_folds = cv_folds,
        machine_learning_weighted = FALSE,
        threshold_rule = "Training-set Youden index, applied unchanged to the testing set"
      ),
      train_index = train_index,
      test_index = test_index,
      preprocessing = list(
        center = train_center,
        scale = train_scale,
        columns = colnames(design_matrix)
      ),
      performance = performance,
      confusion = confusion_long,
      thresholds = training_thresholds,
      best_parameters = best_parameters,
      predictions = predictions,
      dca = dca_curves,
      quality_checks = quality_checks
    ),
    file.path(output_dir, "model_results.rds")
  )

  capture.output(sessionInfo(), file = file.path(output_dir, "R_session_info.txt"))

  format_metric_line <- function(row) {
    sprintf(
      "%s: AUC %.3f (95%% CI %.3f-%.3f), sensitivity %.3f, specificity %.3f, accuracy %.3f, F1 %.3f.",
      row$Model,
      row$AUC,
      row$AUC_CI_low,
      row$AUC_CI_high,
      row$Sensitivity,
      row$Specificity,
      row$Accuracy,
      row$F1
    )
  }

  top_test_lines <- vapply(seq_len(nrow(test_rank)), function(i) {
    format_metric_line(test_rank[i, ])
  }, character(1))

  readme_lines <- c(
    "Fifteen prediction models with an 80:20 split",
    "",
    sprintf("Data: NHANES 2015-2020 pre-pandemic complete-case sample, n=%d; events=%d.", nrow(analysis), sum(y == 1L)),
    sprintf("Split: stratified 80%% training / 20%% testing; split seed=%d.", split_seed),
    sprintf("Training n=%d, events=%d; testing n=%d, events=%d.", length(y_train), sum(y_train), length(y_test), sum(y_test)),
    "Outcome: latest 28-tooth surface definition; at least one untreated coronal-caries tooth surface.",
    "Predictors: 12 conceptual variables selected at lambda.1se in the latest LASSO analysis.",
    paste("Models:", paste(model_names, collapse = ", "), "."),
    "Tuning: 10-fold cross-validation inside the training set only.",
    "Classification threshold: the training-set Youden threshold for each model, applied unchanged to the testing set.",
    "Machine-learning analyses are unweighted, consistent with the manuscript.",
    "Interpretation: testing-set results are primary; training-set results are shown to evaluate overfitting.",
    "",
    "Testing-set ranking:",
    top_test_lines,
    "",
    sprintf(
      "Largest train-test AUC gap: %s, %.3f.",
      auc_gap$Model[1],
      auc_gap$AUC_gap[1]
    )
  )
  writeLines(readme_lines, file.path(output_dir, "model_analysis_notes.txt"), useBytes = TRUE)

  progress("All 15 models completed and exported")
  cat(paste(readme_lines, collapse = "\n"), "\n")
  cat("\nOutput directory:\n", output_dir, "\n", sep = "")

  lr_context <- list(
    model_name = "Logistic Regression", fit = logistic_fit,
    design_formula = design_formula,
    predictor_names = predictor_names,
    matrix_columns = colnames(design_matrix),
    training_features = analysis[train_index, predictor_names, drop = FALSE],
    testing_features = analysis[test_index, predictor_names, drop = FALSE],
    training_ids = analysis_id[train_index], testing_ids = analysis_id[test_index],
    training_probabilities = prob_train[["Logistic Regression"]],
    testing_probabilities = prob_test[["Logistic Regression"]],
    training_outcome = y_train, testing_outcome = y_test,
    split_seed = split_seed
  )
  saveRDS(lr_context, file.path(output_dir, "logistic_regression_for_SHAP.rds"))
  saveRDS(tuning_details, file.path(output_dir, "tuning_details.rds"))
  saveRDS(list(logistic = logistic_fit, ridge = ridge_fit, elastic_net = elastic_fit,
               gbm = gbm_fit, knn = knn_fit, svm = svm_fit, random_forest = rf_fit,
               decision_tree = dt_fit, naive_bayes = nb_fit, adaboost = adaboost_fit,
               mlp = mlp_fit, extra_trees = extra_fit),
          file.path(output_dir, "fitted_R_models.rds"))
  lightgbm::lgb.save(lgb_fit, file.path(output_dir, "lightgbm_model.txt"))
  xgboost::xgb.save(xgb_fit, file.path(output_dir, "xgboost_model.ubj"))
  catboost::catboost.save_model(cat_fit, file.path(output_dir, "catboost_model.cbm"))
  save_plot_pair(four_panel, file.path(output_dir, "Figure_5_ROC_DCA"), 15.5, 13.5)
  invisible(performance)
}

fit_subgroups <- function(ctx) {
  input <- file.path(ctx$data, "analysis_complete_case.rds")
  prior <- file.path(ctx$tables, "logistic_results_details.csv")
  out <- ctx$subgroups
  require_files(c(input, prior))
  d <- readRDS(input)
  covs <- c('age','sex','race','education','marital','pir','bmi','smoking','alcohol',
            'physical_activity','floss_days','dental_visit','diabetes','hypertension',
            'cvd','ckd','liver_disease','present_teeth_28','energy_kcal','sugar_g',
            'cholesterol_mg','fat_g','fiber_g','protein_g','calcium_mg','phosphorus_mg')
  groups <- c('sex','race','education','marital','smoking','alcohol','physical_activity',
              'dental_visit','diabetes','hypertension','cvd','ckd','liver_disease')
  exposures <- c('dietary_vitd_ug','supplement_vitd_ug')
  y <- 'untreated_surface_caries_28'
  stopifnot(nrow(d)==3089, sum(d[[y]])==592, !anyDuplicated(d$SEQN),
            all(d$age>=60), all(d$present_teeth_28>0),
            setequal(unique(as.character(d$source_period)),c('I','P')),
            !anyNA(d[c(y,covs,exposures,'analysis_weight','survey_psu','survey_strata')]),
            all(d$analysis_weight>0), all(d$present_teeth_28+d$missing_teeth_28==28))
  des <- svydesign(ids=~survey_psu,strata=~survey_strata,weights=~analysis_weight,
                   nest=TRUE,data=d)
  stopifnot(degf(des)==40)
  fit_capture <- function(formula, design) {
    ws <- character()
    fit <- tryCatch(withCallingHandlers(
      svyglm(formula,design=design,family=quasibinomial('logit'),
             control=glm.control(maxit=100,epsilon=1e-9)),
      warning=function(w) {ws <<- c(ws,conditionMessage(w));invokeRestart('muffleWarning')}),
      error=function(e)e)
    list(fit=fit,warnings=paste(unique(ws),collapse='; '))
  }
  fits <- list(); results <- list(); ints <- list(); overall <- list()
  effect <- function(fit, exposure, df) {
    b <- unname(coef(fit)[exposure]); se <- sqrt(unname(vcov(fit)[exposure,exposure]))
    c(beta=b,se=se,OR=exp(b),low=exp(b-qt(.975,df)*se),
      high=exp(b+qt(.975,df)*se),p=2*pt(abs(b/se),df,lower.tail=FALSE))
  }
  for (ex in exposures) {
    f <- fit_capture(reformulate(c(ex,covs),response=y),des)
    overall[[ex]] <- data.frame(exposure=ex,t(effect(f$fit,ex,degf(des))))
    fits[[paste(ex,'overall',sep=':')]]<-f$fit
    for (g in groups) {

      fi <- as.formula(paste(y,'~',paste(c(paste(ex,g,sep='*'),setdiff(covs,g)),collapse='+')))
      f <- fit_capture(fi,des)
      stopifnot(!inherits(f$fit,'error'), isTRUE(f$fit$converged), !anyNA(coef(f$fit)))
      term <- as.formula(paste('~',ex,':',g))
      test <- regTermTest(f$fit,term,method='Wald',df=degf(des))

      jj <- grep(':',names(coef(f$fit)),fixed=TRUE)
      b <- coef(f$fit)[jj]; V<-vcov(f$fit)[jj,jj,drop=FALSE]
      Fcheck<-as.numeric(t(b)%*%solve(V,b))/length(b)
      pcheck<-pf(Fcheck,length(b),degf(des),lower.tail=FALSE)
      stopifnot(abs(pcheck-as.numeric(test$p))<1e-8)
      ints[[length(ints)+1L]] <- data.frame(exposure=ex,group=g,p_interaction=as.numeric(test$p),
                                          numerator_df=length(b),denominator_df=degf(des),
                                          F_statistic=Fcheck,warnings=f$warnings)
      fits[[paste(ex,g,'interaction',sep=':')]]<-f$fit
      for (lev in levels(d[[g]])) {
        keep <- d[[g]]==lev

        ds <- des[keep,]; ds$variables<-droplevels(ds$variables)
        cc <- setdiff(covs,g)
        constant<-cc[vapply(ds$variables[cc],function(z)length(unique(z))<2L,logical(1))]
        cc<-setdiff(cc,constant)
        form<-reformulate(c(ex,cc),response=y)
        f<-fit_capture(form,ds)
        n<-sum(keep); events<-sum(d[[y]][keep]); df<-degf(ds)
        status<-'OK'; vals<-setNames(rep(NA_real_,6),c('beta','se','OR','low','high','p'))
        rank<-NA_integer_; converged<-FALSE; boundary_n<-NA_integer_
        if(inherits(f$fit,'error')) status<-conditionMessage(f$fit) else {
          rank<-f$fit$rank; converged<-isTRUE(f$fit$converged)
          boundary_n<-sum(fitted(f$fit)<1e-8 | fitted(f$fit)>1-1e-8)
          if(!converged) status<-'Nonconvergence'
          else if(anyNA(coef(f$fit))) status<-'Rank deficient'
          else if(grepl('0 or 1|did not converge',f$warnings)) status<-'Separation or nonconvergence'
          else if(df<=0) status<-'No design degrees of freedom'
          else {
            vals<-effect(f$fit,ex,df)
            if(!all(is.finite(vals)) || vals['se']<=0) {status<-'Non-finite estimate';vals[]<-NA_real_}
          }
          fits[[paste(ex,g,lev,sep=':')]]<-f$fit
        }
        results[[length(results)+1L]]<-data.frame(exposure=ex,group=g,level=lev,n=n,
          events=events,nonevents=n-events,design_df=df,rank=rank,
          events_per_coefficient=events/max(1,rank-1),converged=converged,
          fitted_boundary_n=boundary_n,t(vals),p_interaction=as.numeric(test$p),
          status=status,warnings=f$warnings,constant_covariates=paste(constant,collapse=';'),
          formula=paste(deparse(form),collapse=' '))
      }
      cat(ex,g,'completed\n');flush.console()
    }
  }
  res<-do.call(rbind,results); inter<-do.call(rbind,ints); ov<-do.call(rbind,overall)
  stopifnot(nrow(res)==68,nrow(inter)==26, all(aggregate(n~exposure+group,res,sum)$n==3089),
            all(aggregate(events~exposure+group,res,sum)$events==592))

  old<-read.csv(prior,check.names=FALSE)
  old<-old[old$model=='Model 2',]
  for(ex in exposures) {
    a<-ov[ov$exposure==ex,];b<-old[old$variable==ex,]
    stopifnot(abs(a$OR-b$estimate)<1e-7,abs(a$p-b$p_value)<1e-7)
  }
  write.csv(res,file.path(out,'subgroup_results.csv'),row.names=FALSE,fileEncoding='UTF-8')
  write.csv(inter,file.path(out,'interaction_tests.csv'),row.names=FALSE,fileEncoding='UTF-8')
  write.csv(ov,file.path(out,'overall_model2.csv'),row.names=FALSE,fileEncoding='UTF-8')
  saveRDS(fits,file.path(out,'subgroup_models.rds'))
  saveRDS(list(data=d,design=des,covariates=covs,groups=groups,exposures=exposures),file.path(out,'analysis_context.rds'))
  capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
  print(res[,c('exposure','group','level','n','events','OR','low','high','p','status','warnings')],row.names=FALSE)
  print(inter[,c('exposure','group','p_interaction')],row.names=FALSE)
}

make_split_baseline <- function(ctx) {
  work <- ctx$tables
  src <- ctx$data
  ml <- ctx$models
  require_files(c(file.path(src, "analysis_complete_case.rds"),
                  file.path(src, "data_context.rds"),
                  file.path(ml, "individual_predictions.csv")))
  d <- readRDS(file.path(src, 'analysis_complete_case.rds'))
  p <- read.csv(file.path(ml, 'individual_predictions.csv'), check.names=FALSE)
  stopifnot(nrow(d)==3089, nrow(p)==3089, !anyDuplicated(d$SEQN), !anyDuplicated(p$id), setequal(d$SEQN,p$id))
  j <- match(d$SEQN,p$id)
  stopifnot(all(d$untreated_surface_caries_28==p$Outcome[j]))
  d$split <- factor(p$Dataset[j], levels=c('Training','Testing'))
  stopifnot(sum(d$split=='Training')==2470, sum(d$split=='Testing')==619)
  context <- readRDS(file.path(src, "data_context.rds"))
  factor_levels <- context$factor_levels
  continuous_specs <- context$continuous_specs
  categorical_specs <- context$categorical_specs
  needed <- c(continuous_specs$variable,vapply(categorical_specs, `[[`, '', 'variable'))
  stopifnot(all(complete.cases(d[needed])), all(d$analysis_weight>0), all(d$source_period %in% c('I','P')))
  des <- svydesign(ids=~survey_psu,strata=~survey_strata,weights=~analysis_weight,nest=TRUE,data=d)
  sets <- list(Total=des,Training=subset(des,split=='Training'),Testing=subset(des,split=='Testing'))
  rows <- list(); details <- list()
  add <- function(label, vals, p=NA_real_, type, variable, level='') {
   rows[[length(rows)+1]] <<- data.frame(Variables=label,Total=vals[1],Training=vals[2],Testing=vals[3],P=p,row_type=type,variable=variable,level=level)
  }
  for (i in seq_len(nrow(continuous_specs))) {
   s <- continuous_specs[i,]; v <- s$variable
   fit <- svyglm(reformulate('split',response=v),des)
   pv <- unname(coef(summary(fit))[2,4])
   vals <- vapply(names(sets),function(g) {
     z <- svymean(reformulate(v),sets[[g]])
     m <- unname(coef(z)[1]); se <- unname(SE(z)[1])
     details[[length(details)+1]] <<- data.frame(variable=v,level='',group=g,n=nrow(sets[[g]]$variables),estimate=m,se=se)
     sprintf('%.2f +/- %.2f',m,se)
   },'')
   add(s$label, vals, pv, 'continuous', v)
  }
  for (s in categorical_specs) {
   v <- s$variable
   pv <- unname(svychisq(as.formula(paste0('~',v,'+split')),des,statistic='F')$p.value)
   add(s$label,rep('',3),pv,'group',v)
   for (lev in names(s$levels)) {
     vals <- vapply(names(sets),function(g) {
       dd <- sets[[g]]; dd$variables$.ind <- as.numeric(as.character(dd$variables[[v]])==lev)
       z <- svymean(~.ind,dd); pct <- unname(coef(z)[1])*100; n <- sum(dd$variables$.ind)
       details[[length(details)+1]] <<- data.frame(variable=v,level=lev,group=g,n=n,estimate=pct,se=unname(SE(z)[1])*100)
       sprintf('%d (%.2f)',n,pct)
     },'')
     add(unname(s$levels[[lev]]),vals,type='level',variable=v,level=lev)
   }
  }
  out <- do.call(rbind,rows)
  stopifnot(all(is.finite(out$P[out$row_type!='level'])), all(out$P[out$row_type!='level']>=0 & out$P[out$row_type!='level']<=1))
  dt <- do.call(rbind,details)
  for (v in vapply(categorical_specs,`[[`,'','variable')) for (g in names(sets)) {
   zz <- dt[dt$variable==v & dt$group==g,]
   stopifnot(sum(zz$n)==nrow(sets[[g]]$variables), abs(sum(zz$estimate)-100)<1e-8)
  }

  old <- read.csv(file.path(src,'table1_values.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
  stopifnot(nrow(old)==nrow(out), all(gsub('+/-','+/-',old$Total,fixed=TRUE)==out$Total))
  write_json(list(n=list(Total=3089,Training=2470,Testing=619),rows=out),file.path(work,'Table_S1.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
  write.csv(out,file.path(work,'Table_S1.csv'),row.names=FALSE,fileEncoding='UTF-8')
  write.csv(dt,file.path(work,'Table_S1_numeric_details.csv'),row.names=FALSE,fileEncoding='UTF-8')
  write.csv(d[,c('SEQN','split','untreated_surface_caries_28')],file.path(work,'verified_split.csv'),row.names=FALSE)
  audit <- list(source_data=file.path(src,'analysis_complete_case.rds'),source_split=file.path(ml,'individual_predictions.csv'),md5=as.list(tools::md5sum(c(file.path(src,'analysis_complete_case.rds'),file.path(ml,'individual_predictions.csv')))),supplement='DSQTVD: average daily supplement intake in preceding 30 days; same as existing 15-model analysis, not Day1',continuous='Survey-weighted mean and standard error; survey-weighted linear regression Wald t test for split',categorical='Unweighted n and survey-weighted percentage; Rao-Scott F test',weights='I: WTDRD1 * 2/5.2; P: WTDRD1PP * 3.2/5.2; original strata and PSU',split='Saved membership; no resplitting or retraining',disease_definitions='Unchanged from existing model dataset, including self-reported DIQ010 diabetes; no composite diabetes redefinition in this task',counts=as.list(table(d$split)),events=as.list(tapply(d$untreated_surface_caries_28,d$split,sum)),checks='Unique IDs; exact full-ID overlap; identical outcomes; no variable missingness; category n and percentages sum correctly; overall values identical to existing weighted baseline')
  write_json(audit,file.path(work,'Table_S1_audit.json'),auto_unbox=TRUE,pretty=TRUE)
  print(out[,c('Variables','Training','Testing','P')],row.names=FALSE)
  invisible(out)
}

fit_day1_regression <- function(ctx) {
  data_dir <- ctx$rcs
  input<-file.path(data_dir,'full_design_data.rds');before<-tools::md5sum(input)
  all<-readRDS(input)
  all<-all[is.finite(all$analysis_weight)&all$analysis_weight>0&!is.na(all$survey_strata)&!is.na(all$survey_psu),]
  demo<-c('age','sex','race','education','marital','pir')
  covars<-c(demo,'bmi','smoking','alcohol','physical_activity','floss_days','dental_visit','diabetes','hypertension','cvd','ckd','liver_disease','present_teeth_28','energy_kcal','sugar_g','cholesterol_mg','fat_g','fiber_g','protein_g','calcium_mg','phosphorus_mg')
  d<-all[all$included,]
  stopifnot(nrow(d)==2969,sum(d$y)==578,!anyDuplicated(d$SEQN),all(d$age>=60),all(d$present_teeth_28>0),setequal(unique(d$source_period),c('I','P')),all(complete.cases(d[,c(covars,'dietary_vitd_ug','day1_vitd_ug')])) )
  full<-svydesign(ids=~survey_psu,strata=~survey_strata,weights=~analysis_weight,nest=TRUE,data=all)
  design<-subset(full,included)
  df<-degf(design);stopifnot(df==40)

  spec<-list(Crude=character(),'Model 1'=demo,'Model 2'=covars)
  exposures<-c('Dietary vitamin D'='dietary_vitd_ug','Supplement vitamin D'='day1_vitd_ug')
  rows<-list();fits<-list()
  for(label in names(exposures))for(model in names(spec)){
   x<-exposures[[label]];f<-reformulate(c(x,spec[[model]]),'y')
   fit<-svyglm(f,design=design,family=quasibinomial('logit'))
   stopifnot(fit$converged,nobs(fit)==nrow(d),all(is.finite(coef(fit))))
   b<-unname(coef(fit)[x]);se<-sqrt(vcov(fit)[x,x]);t<-b/se
   check<-glm(f,data=d,weights=analysis_weight/mean(analysis_weight),family=quasibinomial('logit'))
   stopifnot(abs(unname(coef(check)[x])-b)<1e-7)

   row<-data.frame(exposure=label,model=model,n=nrow(d),events=sum(d$y),beta=b,se=se,OR=exp(b),CI_low=exp(b-qt(.975,df)*se),CI_high=exp(b+qt(.975,df)*se),P=2*pt(abs(t),df,lower.tail=FALSE),design_df=df,residual_df=fit$df.residual,P_residual_df=2*pt(abs(t),fit$df.residual,lower.tail=FALSE),OR_per10=exp(10*b),formula=paste(deparse(f),collapse=' '))
   rows[[length(rows)+1]]<-row;fits[[paste(label,model)]]<-fit
  }
  res<-do.call(rbind,rows)
  prior<-read_json(file.path(data_dir,'summary.json'),simplifyVector=TRUE)
  z<-res[res$exposure=='Supplement vitamin D'&res$model=='Crude',]
  stopifnot(abs(z$OR-prior$linear_effects$OR[1])<1e-8,abs(z$P-prior$linear_effects$p[1])<1e-8,identical(before,tools::md5sum(input)))
  write.csv(res,file.path(data_dir,'main_logistic_day1_results.csv'),row.names=FALSE,fileEncoding='UTF-8')
  saveRDS(fits,file.path(data_dir,'main_logistic_day1_models.rds'))
  write_json(list(years='2015-March 2020',dietary_exposure='DR1TVD',supplement_exposure='DS1TVD',unit='per 1 microgram on day 1',sample_n=nrow(d),events=sum(d$y),model1=demo,model2=covars,separate_exposures=TRUE,design_df=df,inference='Design-df t test and interval as in prior main regression; residual-df sensitivity also saved',validation='Independent GLM point estimates match; prior day-1 crude OR/P reproduced; sample counts, complete cases, years, ages and input hashes verified',results=res),file.path(data_dir,'main_logistic_day1_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
  print(res[,c('exposure','model','n','events','OR','CI_low','CI_high','P','residual_df','P_residual_df')],row.names=FALSE,digits=9)
}

make_calibration_outputs <- function(ctx) {
  results_file <- file.path(ctx$models, "model_results.rds")
  predictions_file <- file.path(ctx$models, "individual_predictions.csv")
  require_files(c(results_file, predictions_file))
  results <- readRDS(results_file)
  predictions <- read.csv(predictions_file, check.names = FALSE)
  model_names <- results$specification$model_names
  stopifnot(length(model_names) == 15L, all(model_names %in% names(predictions)))
  output_dir <- file.path(ctx$models, "calibration")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  model_labels <- c(
    "Logistic Regression" = "Logistic Regression",
    "Ridge Logistic Regression" = "Ridge Logistic Regression",
    "Elastic Net" = "Elastic Net",
    "GBM" = "GBM",
    "LightGBM" = "LightGBM",
    "XGBoost" = "XGBoost",
    "CatBoost" = "CatBoost",
    "KNN" = "KNN",
    "SVM" = "SVM",
    "Random Forest" = "Random Forest",
    "Decision Tree" = "Decision Tree",
    "Naive Bayes" = "Naive Bayes",
    "AdaBoost" = "AdaBoost",
    "MLP" = "MLP",
    "Extra Trees" = "Extra Trees"
  )

  feature_labels <- c(
    supplement_vitd_ug = "Supplement vitamin D",
    age = "Age",
    sex = "Sex",
    race = "Race/ethnicity",
    pir = "PIR",
    smoking = "Smoking status",
    alcohol = "Alcohol use",
    floss_days = "Floss/interdental use",
    dental_visit = "Last dental visit",
    diabetes = "Diabetes",
    sugar_g = "Total sugar",
    fiber_g = "Dietary fiber"
  )

  theme_publication <- function(base_size = 12) {
    theme_bw(base_size = base_size, base_family = "Arial") +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35),
        strip.background = element_rect(fill = "#EAF2F8", color = "#9AA7B0", linewidth = 0.4),
        strip.text = element_text(face = "bold", color = "#263238", size = base_size - 1),
        plot.title = element_text(face = "bold", size = base_size + 3, color = "#202124"),
        plot.subtitle = element_text(size = base_size, color = "#4D5156"),
        plot.caption = element_text(size = base_size - 2, color = "#5F6368", hjust = 0),
        axis.title = element_text(face = "bold", color = "#303134"),
        axis.text = element_text(color = "#3C4043"),
        legend.position = "top",
        plot.margin = margin(12, 16, 12, 12)
      )
  }

  make_calibration_rows <- function(dataset_name, model_name, probability, outcome, groups = 10L) {
    probability <- pmin(pmax(as.numeric(probability), 1e-6), 1 - 1e-6)
    outcome <- as.integer(outcome)
    bin <- ceiling(rank(probability, ties.method = "first") / length(probability) * groups)
    bin <- pmin(pmax(bin, 1L), groups)

    split_index <- split(seq_along(probability), bin)
    rows <- lapply(names(split_index), function(bin_name) {
      idx <- split_index[[bin_name]]
      data.frame(
        Dataset = dataset_name,
        Model = model_name,
        Bin = as.integer(bin_name),
        N = length(idx),
        Events = sum(outcome[idx]),
        Mean_predicted = mean(probability[idx]),
        Observed_proportion = mean(outcome[idx]),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }

  calculate_calibration_metrics <- function(dataset_name, model_name, probability, outcome) {
    probability <- pmin(pmax(as.numeric(probability), 1e-6), 1 - 1e-6)
    outcome <- as.integer(outcome)
    linear_predictor <- qlogis(probability)
    recalibration <- suppressWarnings(glm(outcome ~ linear_predictor, family = binomial()))
    coefficients <- coef(recalibration)
    data.frame(
      Dataset = dataset_name,
      Model = model_name,
      N = length(outcome),
      Events = sum(outcome),
      Brier_score = mean((outcome - probability)^2),
      Calibration_intercept = unname(coefficients[1]),
      Calibration_slope = unname(coefficients[2]),
      stringsAsFactors = FALSE
    )
  }

  calibration_rows <- list()
  calibration_metrics <- list()
  row_counter <- 0L
  metric_counter <- 0L

  for (dataset_name in c("Training", "Testing")) {
    dataset_rows <- predictions$Dataset == dataset_name
    outcome <- predictions$Outcome[dataset_rows]
    for (model_name in model_names) {
      row_counter <- row_counter + 1L
      metric_counter <- metric_counter + 1L
      probability <- predictions[[model_name]][dataset_rows]
      calibration_rows[[row_counter]] <- make_calibration_rows(
        dataset_name, model_name, probability, outcome, groups = 10L
      )
      calibration_metrics[[metric_counter]] <- calculate_calibration_metrics(
        dataset_name, model_name, probability, outcome
      )
    }
  }

  calibration_data <- do.call(rbind, calibration_rows)
  calibration_metrics <- do.call(rbind, calibration_metrics)

  calibration_metrics$Display_model <- model_labels[calibration_metrics$Model]
  calibration_metrics$Facet <- sprintf(
    "%s\nBrier = %.3f",
    calibration_metrics$Display_model,
    calibration_metrics$Brier_score
  )
  calibration_data <- merge(
    calibration_data,
    calibration_metrics[c("Dataset", "Model", "Facet")],
    by = c("Dataset", "Model"),
    all.x = TRUE,
    sort = FALSE
  )

  make_calibration_plot <- function(dataset_name, display_name) {
    plot_data <- calibration_data[calibration_data$Dataset == dataset_name, ]
    facet_order <- calibration_metrics$Facet[
      calibration_metrics$Dataset == dataset_name &
        match(calibration_metrics$Model, model_names) > 0
    ]
    facet_order <- facet_order[order(match(
      calibration_metrics$Model[calibration_metrics$Dataset == dataset_name],
      model_names
    ))]
    plot_data$Facet <- factor(plot_data$Facet, levels = unique(facet_order))

    ggplot(plot_data, aes(x = Mean_predicted, y = Observed_proportion)) +
      geom_abline(slope = 1, intercept = 0, color = "#4D4D4D", linewidth = 0.65, linetype = "dashed") +
      geom_line(color = "#2F75B5", linewidth = 0.85) +
      geom_point(color = "#C4512A", fill = "white", shape = 21, size = 2.1, stroke = 0.75) +
      facet_wrap(~Facet, ncol = 5, drop = FALSE) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0.01)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0.01)) +
      coord_equal() +
      labs(
        title = paste0("Calibration curves for 15 models: ", display_name),
        subtitle = "Ten equal-frequency groups; dashed line indicates perfect calibration",
        x = "Mean predicted probability",
        y = "Observed event proportion",
        caption = "Outcome: at least one untreated coronal-caries tooth surface. Machine-learning analyses are unweighted."
      ) +
      theme_publication(11) +
      theme(legend.position = "none", panel.spacing = unit(0.9, "lines"))
  }

  p_calibration_training <- make_calibration_plot("Training", "training set")
  p_calibration_testing <- make_calibration_plot("Testing", "testing set")

  ggsave(
    file.path(output_dir, "Figure_Calibration_15_models_training.png"),
    p_calibration_training,
    width = 16, height = 10.2, dpi = 320, bg = "white"
  )
  ggsave(
    file.path(output_dir, "Figure_Calibration_15_models_testing.png"),
    p_calibration_testing,
    width = 16, height = 10.2, dpi = 320, bg = "white"
  )
  ggsave(
    file.path(output_dir, "Figure_Calibration_15_models_training_testing.png"),
    p_calibration_training / p_calibration_testing,
    width = 16, height = 20.5, dpi = 300, bg = "white"
  )

  write.csv(calibration_metrics, file.path(output_dir, "calibration_metrics.csv"), row.names = FALSE)
  write.csv(calibration_data, file.path(output_dir, "calibration_coordinates.csv"), row.names = FALSE)
  save_plot_pair(p_calibration_training,
                 file.path(output_dir, "Calibration_training"), 16, 10.2)
  save_plot_pair(p_calibration_testing,
                 file.path(output_dir, "Calibration_testing"), 16, 10.2)
  invisible(calibration_metrics)
}

predict_logistic_probability <- function(object, newdata) {
  mm <- model.matrix(object$design_formula, data = newdata)
  missing <- setdiff(object$matrix_columns, colnames(mm))
  if (length(missing)) stop("Missing model-matrix columns: ", paste(missing, collapse = ", "))
  mm <- mm[, object$matrix_columns, drop = FALSE]
  as.numeric(predict(object$fit, newdata = data.frame(mm, check.names = TRUE),
                     type = "response"))
}

explain_logistic_regression <- function(ctx) {
  file <- file.path(ctx$models, "logistic_regression_for_SHAP.rds")
  require_files(file)
  object <- readRDS(file)
  stopifnot(identical(object$model_name, "Logistic Regression"),
            inherits(object$fit, "glm"), object$fit$family$family == "binomial",
            object$fit$family$link == "logit")
  background <- object$training_features
  explanation_data <- object$testing_features
  stopifnot(!anyNA(background), !anyNA(explanation_data),
            identical(names(background), names(explanation_data)))
  p_train <- predict_logistic_probability(object, background)
  p_test <- predict_logistic_probability(object, explanation_data)
  stopifnot(all(is.finite(p_train)), all(is.finite(p_test)),
            max(abs(pmin(pmax(p_test, 1e-6), 1 - 1e-6) -
                    object$testing_probabilities)) < 1e-8)
  baseline <- mean(p_train)
  if (ctx$config$shap_nsim < 2L) stop("shap_nsim must be at least 2.")
  set.seed(ctx$config$shap_seed)
  phi <- fastshap::explain(
    object = object, X = background, newdata = explanation_data,
    pred_wrapper = predict_logistic_probability,
    nsim = as.integer(ctx$config$shap_nsim), adjust = TRUE,
    baseline = baseline, shap_only = TRUE, parallel = FALSE
  )
  phi <- as.matrix(phi)
  stopifnot(identical(dim(phi), dim(explanation_data)), all(is.finite(phi)))
  reconstruction_error <- baseline + rowSums(phi) - p_test
  stopifnot(max(abs(reconstruction_error)) < 1e-7)
  importance <- data.frame(feature = colnames(phi),
                           mean_absolute_SHAP = colMeans(abs(phi)),
                           stringsAsFactors = FALSE)
  importance <- importance[order(-importance$mean_absolute_SHAP), ]
  importance$rank <- seq_len(nrow(importance))
  shap_values <- data.frame(SEQN = object$testing_ids,
                           outcome = object$testing_outcome,
                           baseline = baseline, predicted_probability = p_test,
                           reconstruction_error = reconstruction_error,
                           phi, check.names = FALSE)
  write.csv(shap_values, file.path(ctx$shap, "Logistic_Regression_SHAP_values.csv"),
            row.names = FALSE)
  write.csv(importance, file.path(ctx$shap, "Logistic_Regression_SHAP_importance.csv"),
            row.names = FALSE)
  write.csv(data.frame(SEQN = object$testing_ids, explanation_data, check.names = FALSE),
            file.path(ctx$shap, "Logistic_Regression_SHAP_features.csv"), row.names = FALSE)
  display_names <- c(
    supplement_vitd_ug = "Supplement vitamin D", age = "Age", sex = "Gender",
    race = "Race/ethnicity", pir = "PIR", smoking = "Smoking status",
    alcohol = "Alcohol use", floss_days = "Flosser/interdental device use",
    dental_visit = "Last dentist visit", diabetes = "Diabetes",
    sugar_g = "Total sugar", fiber_g = "Dietary fiber"
  )
  plot_phi <- phi
  plot_data <- explanation_data
  colnames(plot_phi) <- unname(display_names[colnames(phi)])
  names(plot_data) <- colnames(plot_phi)
  stopifnot(!anyNA(names(plot_data)))
  sv <- shapviz::shapviz(plot_phi, X = plot_data, baseline = baseline)
  high <- which.max(p_test)
  low <- which.min(p_test)
  p_a <- shapviz::sv_importance(sv, kind = "bar", max_display = ncol(phi)) +
    ggplot2::labs(title = "Logistic Regression", tag = "A",
                  x = "Mean absolute SHAP value")
  p_b <- shapviz::sv_importance(sv, kind = "beeswarm", max_display = ncol(phi)) +
    ggplot2::labs(title = "Logistic Regression", tag = "B",
                  x = "SHAP value (probability scale)")
  p_c <- shapviz::sv_waterfall(sv, row_id = high, max_display = ncol(phi)) +
    ggplot2::labs(title = paste0("Higher predicted probability: SEQN ", object$testing_ids[high]),
                  tag = "C")
  p_d <- shapviz::sv_waterfall(sv, row_id = low, max_display = ncol(phi)) +
    ggplot2::labs(title = paste0("Lower predicted probability: SEQN ", object$testing_ids[low]),
                  tag = "D")
  common <- ggplot2::theme(plot.title = ggplot2::element_text(size = 11),
                            plot.tag = ggplot2::element_text(face = "bold"))
  panels <- patchwork::wrap_plots(p_a + common, p_b + common,
                                  p_c + common, p_d + common, ncol = 2)
  save_plot_pair(panels, file.path(ctx$shap, "Figure_6_Logistic_Regression_SHAP"), 18, 14)
  save_plot_pair(p_a, file.path(ctx$shap, "Figure_6A_SHAP_importance"), 9, 7)
  save_plot_pair(p_b, file.path(ctx$shap, "Figure_6B_SHAP_beeswarm"), 10, 7)
  save_plot_pair(p_c, file.path(ctx$shap, "Figure_6C_SHAP_higher_probability"), 10, 7)
  save_plot_pair(p_d, file.path(ctx$shap, "Figure_6D_SHAP_lower_probability"), 10, 7)
  audit <- list(
    model = "Logistic Regression", family = "binomial", link = "logit",
    explained_output = "Probability of untreated caries",
    background = "Training cohort only", background_n = nrow(background),
    explained_cohort = "Testing", explained_n = nrow(explanation_data),
    method = "Monte Carlo interventional SHAP with local-accuracy adjustment",
    nsim = ctx$config$shap_nsim, seed = ctx$config$shap_seed,
    baseline = baseline, max_additivity_error = max(abs(reconstruction_error)),
    high_probability_id = object$testing_ids[high], low_probability_id = object$testing_ids[low]
  )
  jsonlite::write_json(audit, file.path(ctx$shap, "Logistic_Regression_SHAP_audit.json"),
                       pretty = TRUE, auto_unbox = TRUE, digits = 15)
  saveRDS(list(model = object, SHAP = phi, importance = importance, audit = audit),
          file.path(ctx$shap, "Logistic_Regression_SHAP.rds"))
  invisible(importance)
}

render_dose_response <- function(ctx) {
  require_files(file.path(ctx$rcs, c("summary.json", "rcs_curve.csv",
                                     "segmented_curve_conditional.csv")))
  summary <- jsonlite::read_json(file.path(ctx$rcs, "summary.json"), simplifyVector = TRUE)
  rcs_data <- read.csv(file.path(ctx$rcs, "rcs_curve.csv"))
  segmented_data <- read.csv(file.path(ctx$rcs, "segmented_curve_conditional.csv"))
  fmt_p <- function(p) if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  make_panel <- function(d, title, tag) {
    ggplot2::ggplot(d, ggplot2::aes(x = dose, y = OR)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high),
                          fill = "#64B5F6", alpha = 0.20) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_line(colour = "#1565C0", linewidth = 0.9) +
      ggplot2::coord_cartesian(xlim = c(0, 100)) +
      ggplot2::labs(x = "Day-1 supplement vitamin D (ug/day)",
                    y = "Odds ratio (95% CI)", title = title, tag = tag) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(aspect.ratio = 1, plot.tag = ggplot2::element_text(face = "bold"))
  }
  p_a <- make_panel(rcs_data, "Restricted cubic spline", "A") +
    ggplot2::labs(subtitle = paste0("Overall P ", fmt_p(summary$p_overall),
                                     "; nonlinear P ", fmt_p(summary$p_nonlinear)))
  p_b <- make_panel(segmented_data, "Two-piece logistic regression", "B") +
    ggplot2::geom_vline(xintercept = summary$breakpoint, linetype = "dotted") +
    ggplot2::labs(subtitle = paste0("Estimated breakpoint: ",
                                     signif(summary$breakpoint, 5), " ug/day"))
  save_plot_pair(patchwork::wrap_plots(p_a, p_b, ncol = 1),
                  file.path(ctx$figures, "Figure_2_dose_response"), 8, 15)
}

render_subgroup_forest <- function(ctx) {
  file <- file.path(ctx$subgroups, "subgroup_results.csv")
  require_files(file)
  all <- read.csv(file, stringsAsFactors = FALSE)
  group_names <- c(sex = "Gender", race = "Race/ethnicity", education = "Education",
                   marital = "Marital status", smoking = "Smoking", alcohol = "Alcohol",
                   physical_activity = "Physical activity", dental_visit = "Dentist visit",
                   diabetes = "Diabetes", hypertension = "Hypertension", cvd = "CVD",
                   ckd = "CKD", liver_disease = "Liver disease")
  plots <- list()
  exposures <- c("dietary_vitd_ug", "supplement_vitd_ug")
  labels <- c("Dietary vitamin D", "Supplement vitamin D")
  for (i in seq_along(exposures)) {
    d <- all[all$exposure == exposures[i], ]
    d$row_label <- paste0(unname(group_names[d$group]), ": ", d$level,
                           " (n=", d$n, ")")
    d$row_label <- factor(d$row_label, levels = rev(unique(d$row_label)))
    d$estimable <- d$status == "OK" & is.finite(d$OR) & is.finite(d$low) &
      is.finite(d$high) & d$low > 0 & d$high > 0
    d$display <- ifelse(d$estimable,
                         sprintf("%.2f (%.2f-%.2f)", d$OR, d$low, d$high),
                         "Not estimable")
    p <- ggplot2::ggplot(d, ggplot2::aes(y = row_label)) +
      ggplot2::geom_blank(ggplot2::aes(x = 1)) +
      ggplot2::geom_vline(xintercept = 1, colour = "grey50", linetype = "dashed") +
      ggplot2::geom_segment(data = d[d$estimable, ],
                            ggplot2::aes(x = low, xend = high, yend = row_label),
                            linewidth = 0.5, colour = "#1565C0") +
      ggplot2::geom_point(data = d[d$estimable, ], ggplot2::aes(x = OR),
                          size = 2, colour = "#1565C0") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = labels[i], x = "Adjusted OR per 1 ug/day (95% CI)", y = NULL) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    plots[[i]] <- p
    write.csv(d, file.path(ctx$subgroups, paste0(exposures[i], "_forest_data.csv")),
              row.names = FALSE)
  }
  save_plot_pair(patchwork::wrap_plots(plots, ncol = 2),
                  file.path(ctx$figures, "Figure_S1_subgroups"), 18, 13)
}

draw_flow_diagram <- function(labels, title) {
  n <- length(labels)
  boxes <- data.frame(x = 0, y = rev(seq_len(n)), label = labels)
  arrows <- if (n > 1L) data.frame(x = 0, xend = 0,
                                   y = rev(seq_len(n))[-n] - 0.30,
                                   yend = rev(seq_len(n))[-1] + 0.30) else NULL
  p <- ggplot2::ggplot(boxes) +
    ggplot2::geom_rect(ggplot2::aes(xmin = -1, xmax = 1, ymin = y - 0.29,
                                   ymax = y + 0.29), fill = "#E8F1FA",
                       colour = "#284B63", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = label), size = 3.2) +
    ggplot2::coord_cartesian(xlim = c(-1.04, 1.04), clip = "off") +
    ggplot2::labs(title = title) + ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                    plot.margin = ggplot2::margin(15, 15, 15, 15))
  if (n > 1L) p <- p + ggplot2::geom_segment(
    data = arrows, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.12, "inches")), linewidth = 0.5)
  p
}

render_flows <- function(ctx) {
  file <- file.path(ctx$data, "sample_flow.csv")
  require_files(file)
  flow <- read.csv(file, stringsAsFactors = FALSE)
  labels <- vapply(seq_len(nrow(flow)), function(i) {
    paste0(paste(strwrap(flow$step[i], width = 52), collapse = "\n"),
             "\nn = ", flow$remaining_n[i])
  }, character(1))
  p_a <- draw_flow_diagram(labels, "A. Participant selection")
  n <- tail(flow$remaining_n, 1)
  b_labels <- c(paste0("NHANES 2015-March 2020\nMain analysis: n = ", n),
                "Weighted baseline and logistic regression\nDietary and 30-day supplement intake",
                "Separate day-1 supplement analysis\nRCS and threshold analysis",
                "Spearman correlations and LASSO\nOriginal feature-selection sequence",
                "Stratified training/testing split: 80% / 20%\n15 prediction models",
                "Training-set tuning and thresholds\nROC, DCA and classification metrics",
                "Logistic Regression SHAP analysis\nTraining background; testing explanations")
  p_b <- draw_flow_diagram(b_labels, "B. Analysis workflow")
  save_plot_pair(p_a, file.path(ctx$figures, "Figure_1A_participant_flow"), 8, 12)
  save_plot_pair(p_b, file.path(ctx$figures, "Figure_1B_analysis_workflow"), 8, 11)
}

export_tables <- function(ctx) {
  paths <- c(file.path(ctx$data, "table1_values.csv"),
              file.path(ctx$tables, "logistic_results_wide.csv"),
              file.path(ctx$tables, "Table_3.csv"),
              file.path(ctx$models, "model_performance.csv"))
  require_files(paths)
  table1 <- read.csv(paths[1], check.names = FALSE)
  table2 <- read.csv(paths[2], check.names = FALSE)
  table3 <- read.csv(paths[3], check.names = FALSE)
  performance <- read.csv(paths[4], check.names = FALSE)
  table4 <- performance[performance$Dataset == "Testing", ]
  table_s2 <- performance[performance$Dataset == "Training", ]
  write.csv(table1, file.path(ctx$tables, "Table_1.csv"), row.names = FALSE)
  write.csv(table2, file.path(ctx$tables, "Table_2.csv"), row.names = FALSE)
  write.csv(table4, file.path(ctx$tables, "Table_4.csv"), row.names = FALSE)
  write.csv(table_s2, file.path(ctx$tables, "Table_S2.csv"), row.names = FALSE)
  table_s1 <- make_split_baseline(ctx)
  writexl::write_xlsx(list(Table_1 = table1, Table_2 = table2, Table_3 = table3,
                          Table_4 = table4, Table_S1 = table_s1, Table_S2 = table_s2),
                      file.path(ctx$tables, "All_manuscript_tables.xlsx"))
  invisible(list(Table_1 = table1, Table_2 = table2, Table_3 = table3,
                 Table_4 = table4, Table_S1 = table_s1, Table_S2 = table_s2))
}

run_analysis <- function(config = default_config()) {
  defaults <- default_config()
  unknown <- setdiff(names(config), names(defaults))
  if (length(unknown)) stop("Unknown configuration options: ", paste(unknown, collapse = ", "))
  config <- utils::modifyList(defaults, config)
  stages <- config$stages
  known <- defaults$stages
  if (!length(stages) || any(!stages %in% known)) stop("Invalid analysis stages.")
  if (!is.numeric(config$bootstrap_reps) || config$bootstrap_reps < 0 ||
      config$bootstrap_reps != as.integer(config$bootstrap_reps)) stop("Invalid bootstrap_reps.")
  options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust", scipen = 999,
           timeout = max(config$timeout_seconds, getOption("timeout", 60)))
  if (.Platform$OS.type == "windows") {
    suppressWarnings(try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE))
  }
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  packages <- load_required_packages(stages)
  ctx <- make_context(config)
  config$output_dir <- ctx$root
  jsonlite::write_json(config, file.path(ctx$root, "run_config.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  capture.output(sessionInfo(), file = file.path(ctx$root, "sessionInfo_start.txt"))
  versions <- data.frame(package = packages,
                          version = vapply(packages, function(p) as.character(packageVersion(p)),
                                            character(1)))
  write.csv(versions, file.path(ctx$root, "package_versions.csv"), row.names = FALSE)
  status <- list()
  for (stage in stages) {
    message("Starting stage: ", stage)
    started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    tryCatch({
      switch(stage,
        download = download_inputs(ctx),
        data = prepare_data(ctx),
        regression = fit_main_regression(ctx),
        rcs = { prepare_day1_exposure(ctx); fit_dose_response(ctx) },
        day1_regression = fit_day1_regression(ctx),
        lasso = fit_correlation_lasso(ctx),
        models = fit_prediction_models(ctx),
        calibration = make_calibration_outputs(ctx),
        shap = explain_logistic_regression(ctx),
        subgroups = fit_subgroups(ctx),
        tables = export_tables(ctx),
        figures = { render_flows(ctx); render_dose_response(ctx); render_subgroup_forest(ctx) }
      )
      status[[length(status) + 1L]] <- data.frame(stage = stage, status = "completed",
        started_UTC = started, finished_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
        error = "")
    }, error = function(e) {
      status[[length(status) + 1L]] <<- data.frame(stage = stage, status = "failed",
        started_UTC = started, finished_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
        error = conditionMessage(e))
      write.csv(do.call(rbind, status), file.path(ctx$root, "stage_status.csv"), row.names = FALSE)
      stop("Stage ", stage, " failed: ", conditionMessage(e), call. = FALSE)
    })
    write.csv(do.call(rbind, status), file.path(ctx$root, "stage_status.csv"), row.names = FALSE)
  }
  capture.output(sessionInfo(), file = file.path(ctx$root, "sessionInfo_end.txt"))
  files <- list.files(ctx$root, recursive = TRUE, full.names = TRUE)
  files <- files[!dir.exists(files)]
  root_prefix <- paste0(normalizePath(ctx$root, winslash = "/", mustWork = TRUE), "/")
  full <- normalizePath(files, winslash = "/", mustWork = TRUE)
  manifest <- data.frame(file = substring(full, nchar(root_prefix) + 1L),
                          bytes = file.info(files)$size, md5 = unname(tools::md5sum(files)))
  write.csv(manifest, file.path(ctx$root, "output_manifest.csv"), row.names = FALSE)
  message("Completed requested stages. Output: ", normalizePath(ctx$root, winslash = "/"))
  invisible(ctx)
}

code_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!"--run" %in% args) {
    cat("Code definitions loaded. No analysis was executed.\n",
        "Run deliberately with: Rscript code.R --run\n",
        "Optional arguments: --output=PATH --cache=PATH --stages=download,data,...\n",
        "Or source this file and call run_analysis(default_config()).\n", sep = "")
    return(invisible(NULL))
  }
  config <- default_config()
  recognized <- args == "--run" | grepl("^--(output|cache|stages)=", args)
  if (any(!recognized)) stop("Unknown command-line argument: ", paste(args[!recognized], collapse = ", "))
  for (arg in args) {
    if (startsWith(arg, "--output=")) config$output_dir <- substring(arg, 10)
    if (startsWith(arg, "--cache=")) config$cache_dir <- substring(arg, 9)
    if (startsWith(arg, "--stages=")) config$stages <- strsplit(substring(arg, 10), ",", fixed = TRUE)[[1]]
  }
  run_analysis(config)
}

if (sys.nframe() == 0L) code_main()
