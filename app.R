# Analysis: Kawe: Quantitative and Psychometrics (CFA, Rasch/IRT, & Multilevel Meta-Analysis)
# Principal Investigator: Aditya Kusumawijaya, S.Psi (Universitas Sanata Dharma / Universitas Gadjah Mada)
# Packages: shiny, bslib, lavaan, TAM, WrightMap, metafor, esc, DT, ggplot2, misty, psych, semTools, readxl

library(shiny)
library(bslib)
library(Matrix)
library(MASS)
library(metafor)
library(esc)
library(DT)
library(ggplot2)
library(grid)
library(lavaan)
library(semPlot)
library(psych)
library(semTools)
library(readxl)
library(TAM)
library(WrightMap)
library(misty)
library(parallel)
library(doParallel)
library(foreach)

# Safe CPU core detection for local and WebAssembly (WebR) runtimes
safe_detect_cores <- function() {
  n <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  if (is.null(n) || is.na(n) || length(n) == 0 || n < 1) 1L else as.integer(n)
}
n_sys_cores <- safe_detect_cores()

# Benchmark Dataset Generator (k = 14 RCT Trials in Clinical Psychology)
create_benchmark_data <- function() {
  data.frame(
    Study = c("Beck et al.", "Clark & Wells", "Cuijpers et al.", "Hofmann et al.", 
              "Andersson et al.", "Linehan et al.", "Ehlers et al.", "Barlow et al.", 
              "Butler et al.", "Kessler et al.", "Richards et al.", "Cuijpers & Reynolds",
              "Dimidjian et al.", "Andrews et al."),
    Year = c(2012, 2014, 2015, 2017, 2018, 2019, 2020, 2021, 2022, 2022, 2023, 2023, 2024, 2024),
    Intervention = c("CBT", "CBT", "iCBT", "CBT", "iCBT", "ACT", "CBT", "ACT", "CBT", "iCBT", "BA", "CBT", "BA", "iCBT"),
    Format = c("Individual", "Individual", "Guided Digital", "Group", "Guided Digital", "Individual", "Individual", "Group", "Individual", "Guided Digital", "Individual", "Group", "Individual", "Guided Digital"),
    Control_Type = c("Waitlist", "TAU", "Waitlist", "TAU", "Placebo", "TAU", "Waitlist", "Active Control", "TAU", "Waitlist", "TAU", "Waitlist", "Active Control", "Placebo"),
    Active_Comparator = c("None (Waitlist)", "Treatment as Usual", "None (Waitlist)", "Treatment as Usual", "Attention Placebo", "Treatment as Usual", "None (Waitlist)", "Relaxation Training", "Treatment as Usual", "None (Waitlist)", "Treatment as Usual", "None (Waitlist)", "Supportive Counseling", "Attention Placebo"),
    Population = c("Clinical (MDD)", "Clinical (MDD)", "Subclinical", "Clinical (MDD)", "Subclinical", "Clinical (MDD)", "Clinical (MDD)", "Clinical (Anxiety)", "Clinical (MDD)", "Subclinical", "Clinical (MDD)", "Clinical (MDD)", "Clinical (MDD)", "Subclinical"),
    Age_Group = c("Adults", "Adults", "Adults", "Adults", "Adults", "Adults", "Adults", "Adults", "Adults", "Adolescents", "Adults", "Older Adults", "Adults", "Adults"),
    Therapist_Level = c("Licensed Clinician", "Licensed Clinician", "Automated Digital", "Licensed Clinician", "Automated Digital", "Licensed Clinician", "Licensed Clinician", "Supervised Trainee", "Licensed Clinician", "Automated Digital", "Supervised Trainee", "Licensed Clinician", "Licensed Clinician", "Automated Digital"),
    Risk_of_Bias = c("Low Risk", "Low Risk", "Some Concerns", "Low Risk", "Some Concerns", "Low Risk", "Low Risk", "Some Concerns", "Low Risk", "Some Concerns", "Low Risk", "Some Concerns", "Low Risk", "Low Risk"),
    Duration_Weeks = c(12, 16, 8, 12, 10, 24, 12, 16, 12, 8, 10, 14, 16, 8),
    Total_Sessions = c(12, 16, 8, 12, 8, 24, 12, 14, 12, 8, 10, 14, 16, 8),
    N_Intervention = c(45, 52, 110, 60, 95, 40, 48, 55, 62, 120, 50, 70, 55, 115),
    N_Control = c(42, 50, 105, 58, 90, 38, 46, 52, 60, 115, 48, 68, 52, 110),
    yi = c(0.820, 0.640, 0.280, 0.540, 0.320, 0.880, 0.700, 0.460, 0.650, 0.240, 0.610, 0.480, 0.660, 0.290),
    vi = c(0.0499, 0.0412, 0.0188, 0.0351, 0.0219, 0.0563, 0.0452, 0.0384, 0.0345, 0.0172, 0.0427, 0.0298, 0.0394, 0.0180),
    sei = c(0.2234, 0.2030, 0.1371, 0.1873, 0.1480, 0.2373, 0.2126, 0.1960, 0.1857, 0.1311, 0.2066, 0.1726, 0.1985, 0.1342),
    stringsAsFactors = FALSE
  )
}

# Sample Batch Generator for Calculator (k = 5 — CBT for Depression in Adults with ACEs)
create_calc_sample_batch <- function() {
  data.frame(
    Study = c("Aardal et al", "Ghrun et al", "Krakau et al", "Malkomsen et al", "Wagner et al"),
    Year = c(2026, 2021, 2024, 2025, 2024),
    Intervention = c("CBT", "CBT", "CBT", "CBT", "CBT"),
    Format = c("Interaction", "Main", "Interaction", "Interaction", "Main"),
    Control_Type = c("TAU", "Waitlist", "Active Control", "TAU", "Waitlist"),
    Active_Comparator = c("Treatment as Usual", "None (Waitlist)", "Psychodynamic Therapy", "Treatment as Usual", "None (Waitlist)"),
    Population = c("Non-US", "US", "Non-US", "Non-US", "Non-US"),
    Age_Group = c("Adults", "Adults", "Adults", "Adults", "Adults"),
    Therapist_Level = c("Licensed Clinician", "Licensed Clinician", "Licensed Clinician", "Licensed Clinician", "Licensed Clinician"),
    Risk_of_Bias = c("Low Risk", "Some Concerns", "Low Risk", "Low Risk", "Some Concerns"),
    Duration_Weeks = c(13.5, 12.0, 260.0, 28.0, 6.0),
    Total_Sessions = c(14, 12, 40, 20, 6),
    N_Intervention = c(48, 50, 104, 50, 38),
    N_Control = c(52, 46, 148, 49, 34),
    yi = c(0.5053, -0.0797, -0.3087, 0.2812, 0.8206),
    vi = c(0.2034^2, 0.4099^2, 0.1287^2, 0.2020^2, 0.4807^2),
    sei = c(0.2034, 0.4099, 0.1287, 0.2020, 0.4807),
    stringsAsFactors = FALSE
  )
}

# Dynamic Effect Size Formula Builder for Active Calculator Input
build_active_formula_ui <- function(itype, metric, input_list, res) {
  es_val <- if (!is.null(res) && is.finite(res$es)) res$es else 0.50
  var_val <- if (!is.null(res) && is.finite(res$var)) res$var else 0.05
  se_val <- if (!is.null(res) && is.finite(res$se)) res$se else sqrt(var_val)
  ci_lo <- if (!is.null(res) && is.finite(res$ci.lo)) res$ci.lo else (es_val - 1.96 * se_val)
  ci_hi <- if (!is.null(res) && is.finite(res$ci.hi)) res$ci.hi else (es_val + 1.96 * se_val)
  
  metric_name <- ifelse(metric == "g", "Hedges' g (Unbiased)", "Cohen's d")
  
  if (itype == "mean_sd") {
    n1 <- input_list$grp1n_sd; m1 <- input_list$grp1m_sd; s1 <- input_list$grp1sd_sd
    n2 <- input_list$grp2n_sd; m2 <- input_list$grp2m_sd; s2 <- input_list$grp2sd_sd
    dir_txt <- ifelse(input_list$sd_outcome_dir == "decrease_better", "Skor Lebih Rendah = Efektif (Simtom Menurun)", "Skor Lebih Tinggi = Efektif")
    df_val <- n1 + n2 - 2
    sp_val <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / df_val)
    raw_diff <- if (input_list$sd_outcome_dir == "decrease_better") (m2 - m1) else (m1 - m2)
    d_val <- raw_diff / sp_val
    j_val <- 1 - (3 / (4 * (n1 + n2) - 9))
    g_val <- d_val * j_val
    vd_val <- (n1 + n2) / (n1 * n2) + (d_val^2) / (2 * (n1 + n2))
    vg_val <- if (metric == "g") (j_val^2 * vd_val) else vd_val
    
    tagList(
      div(
        class = "kawe-math-card",
        div(class = "kawe-formula-title", icon("square-root-alt"), " Metode: Post-Test Means & Standard Deviations (Dua Kelompok Independen)"),
        div(class = "kawe-formula-sub", "Rujukan Baku: Cohen (1988); Hedges (1981); Borenstein et al. (2009)."),
        p("Mengestimasi perbedaan rata-rata terstandardisasi antara kelompok intervensi dan kontrol pascates menggunakan deviasi baku gabungan (", tags$em("pooled standard deviation"), ")."),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "1. Landasan Formulasi Matematis:"),
        div(class = "kawe-eq-box",
          HTML("$$\\text{Selisih Rata-rata: } \\Delta\\bar{X} = \\bar{X}_1 - \\bar{X}_2$$"),
          HTML("$$\\text{Deviasi Baku Gabungan: } s_{\\text{pooled}} = \\sqrt{\\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}$$"),
          HTML("$$\\text{Cohen's } d = \\frac{\\Delta\\bar{X}}{s_{\\text{pooled}}}$$"),
          HTML("$$\\text{Faktor Koreksi Hedges: } J = 1 - \\frac{3}{4(n_1 + n_2) - 9}$$"),
          HTML("$$\\text{Hedges' } g = J \\times d$$"),
          HTML("$$\\text{Varians Sampling: } v_d = \\frac{n_1 + n_2}{n_1 n_2} + \\frac{d^2}{2(n_1 + n_2)}, \\quad v_g = J^2 \\times v_d$$"),
          HTML("$$\\text{Standar Error: } SE = \\sqrt{v}, \\quad 95\\% \\text{ CI} = [\\text{ES} \\pm 1.96 \\cdot SE]$$")
        ),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "2. Langkah Perhitungan Terbuka & Substitusi Angka Riil:"),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 1"),
          strong("Ukuran Sampel & Derajat Kebebasan: "),
          sprintf("n1 = %d, n2 = %d, Total N = %d, Derajat Kebebasan (df) = %d + %d - 2 = %d.", n1, n2, n1 + n2, n1, n2, df_val)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 2"),
          strong("Selisih Rata-rata Terarah: "),
          sprintf("Mean1 = %.2f, Mean2 = %.2f | Orientasi: %s -> \\Delta Mean = %.3f.", m1, m2, dir_txt, raw_diff)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 3"),
          strong("Deviasi Baku Gabungan (s_pooled): "),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$s_{\\text{pooled}} = \\sqrt{\\frac{(%d - 1)(%.2f^2) + (%d - 1)(%.2f^2)}{%d}} = \\sqrt{\\frac{%.2f + %.2f}{%d}} = %.4f$$",
                         n1, s1, n2, s2, df_val, (n1-1)*s1^2, (n2-1)*s2^2, df_val, sp_val))
          )
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 4"),
          strong("Standardized Mean Difference (Cohen's d): "),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$d = \\frac{%.3f}{%.4f} = %.4f$$", raw_diff, sp_val, d_val))
          )
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 5"),
          strong("Faktor Pembobot Koreksi Bias Sampel Kecil Hedges (J): "),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$J = 1 - \\frac{3}{4(%d) - 9} = 1 - \\frac{3}{%d} = %.5f$$", n1 + n2, 4 * (n1 + n2) - 9, j_val))
          ),
          sprintf("Estimasi Titik Akhir: %s = %.4f.", metric_name, es_val)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 6"),
          strong("Varians Sampling, Standar Error & 95% Confidence Interval: "),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$v = %.5f, \\quad SE = \\sqrt{%.5f} = %.4f$$", var_val, var_val, se_val)),
            HTML(sprintf("$$95\\%% \\text{ CI} = [%.4f - 1.96(%.4f), \\; %.4f + 1.96(%.4f)] = [%.4f, \\; %.4f]$$",
                         es_val, se_val, es_val, se_val, ci_lo, ci_hi))
          )
        )
      )
    )
  } else if (itype == "mean_se") {
    n1 <- input_list$grp1n_se; m1 <- input_list$grp1m_se; se1 <- input_list$grp1se_se
    n2 <- input_list$grp2n_se; m2 <- input_list$grp2m_se; se2 <- input_list$grp2se_se
    dir_txt <- ifelse(input_list$se_outcome_dir == "decrease_better", "Skor Lebih Rendah = Efektif (Simtom Menurun)", "Skor Lebih Tinggi = Efektif")
    s1_est <- se1 * sqrt(n1)
    s2_est <- se2 * sqrt(n2)
    df_val <- n1 + n2 - 2
    sp_val <- sqrt(((n1 - 1) * s1_est^2 + (n2 - 1) * s2_est^2) / df_val)
    raw_diff <- if (input_list$se_outcome_dir == "decrease_better") (m2 - m1) else (m1 - m2)
    d_val <- raw_diff / sp_val
    j_val <- 1 - (3 / (4 * (n1 + n2) - 9))
    
    tagList(
      div(
        class = "kawe-math-card",
        div(class = "kawe-formula-title", icon("square-root-alt"), " Metode: Means & Standard Errors (Dua Kelompok Independen)"),
        div(class = "kawe-formula-sub", "Rujukan Baku: Higgins et al. (2019); Cochrane Handbook; Borenstein et al. (2009)."),
        p("Ketika deviasi baku (SD) tidak dilaporkan, deviasi baku direkonstruksi matematis dari Galat Baku (SE) dan ukuran sampel n."),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "1. Landasan Formulasi Matematis:"),
        div(class = "kawe-eq-box",
          HTML("$$s_1 = SE_1 \\times \\sqrt{n_1}, \\quad s_2 = SE_2 \\times \\sqrt{n_2}$$"),
          HTML("$$s_{\\text{pooled}} = \\sqrt{\\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}$$"),
          HTML("$$d = \\frac{\\Delta\\bar{X}}{s_{\\text{pooled}}}, \\quad J = 1 - \\frac{3}{4(n_1 + n_2) - 9}, \\quad g = J \\times d$$")
        ),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "2. Langkah Perhitungan Terbuka & Substitusi Angka Riil:"),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 1"),
          strong("Rekonstruksi Deviasi Baku: "),
          sprintf("Kelompok 1: s1 = %.3f * sqrt(%d) = %.4f | Kelompok 2: s2 = %.3f * sqrt(%d) = %.4f.", se1, n1, s1_est, se2, n2, s2_est)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 2"),
          strong("Deviasi Baku Gabungan (s_pooled): "),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$s_{\\text{pooled}} = \\sqrt{\\frac{(%d - 1)(%.4f^2) + (%d - 1)(%.4f^2)}{%d}} = %.4f$$", n1, s1_est, n2, s2_est, df_val, sp_val))
          )
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 3"),
          strong("Estimasi Akhir Besaran Efek: "),
          sprintf("Cohen's d = %.4f | Faktor J = %.5f | %s = %.4f (SE = %.4f, 95%% CI [%.4f, %.4f]).",
                  d_val, j_val, metric_name, es_val, se_val, ci_lo, ci_hi)
        )
      )
    )
  } else if (itype == "mean_gain") {
    gsub <- input_list$gain_sub_type
    r_val <- input_list$gain_r
    dir_txt <- ifelse(input_list$gain_outcome_dir == "decrease_better", "Skor Menurun = Membaik", "Skor Meningkat = Membaik")
    
    tagList(
      div(
        class = "kawe-math-card",
        div(class = "kawe-formula-title", icon("square-root-alt"), " Metode: Mean Gain / Perubahan Pre-Post (Repeated Measures Designs)"),
        div(class = "kawe-formula-sub", "Rujukan Baku: Becker (1988); Morris (2008); Rosenthal (1991); Lipsey & Wilson (2001)."),
        p("Mengestimasi besaran efek dari selisih perubahan skor pre-post antara kelompok intervensi dan kontrol dengan mempertimbangkan korelasi berpasangan (", tags$em("paired correlation r"), ")."),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "1. Landasan Formulasi Matematis:"),
        div(class = "kawe-eq-box",
          HTML("$$\\text{Gain}_1 = \\bar{X}_{post, 1} - \\bar{X}_{pre, 1}, \\quad \\text{Gain}_2 = \\bar{X}_{post, 2} - \\bar{X}_{pre, 2}$$"),
          HTML("$$s_{\\text{change}, i} = \\sqrt{s_{pre, i}^2 + s_{post, i}^2 - 2 r s_{pre, i} s_{post, i}}$$"),
          HTML("$$s_{\\text{pooled, change}} = \\sqrt{\\frac{(n_1 - 1)s_{\\text{change}, 1}^2 + (n_2 - 1)s_{\\text{change}, 2}^2}{n_1 + n_2 - 2}}$$"),
          HTML("$$d_{\\text{gain}} = \\frac{\\text{Gain}_1 - \\text{Gain}_2}{s_{\\text{pooled, change}}}, \\quad g = J \\times d$$"),
          HTML("$$v_d = \\frac{2(1 - r)}{n_1} + \\frac{2(1 - r)}{n_2} + \\frac{d^2}{2(n_1 + n_2)}, \\quad v_g = J^2 \\times v_d$$")
        ),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "2. Langkah Perhitungan Terbuka & Substitusi Angka Riil:"),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Parameter"),
          sprintf("Pola Data: %s | Orientasi Klinis: %s | Asumsi Korelasi Pre-Post (r) = %.2f.",
                  ifelse(gsub == "pre_post_full", "Pre & Post Lengkap", "Perubahan Pre-Post Langsung"), dir_txt, r_val)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Hasil Estimasi"),
          sprintf("Estimasi Titik %s = %.4f | Varians vi = %.5f | SE = %.4f | 95%% CI [%.4f, %.4f].",
                  metric_name, es_val, var_val, se_val, ci_lo, ci_hi)
        )
      )
    )
  } else if (itype == "t_test") {
    t_val <- input_list$t_val
    n1 <- input_list$t_n1; n2 <- input_list$t_n2
    df_val <- n1 + n2 - 2
    d_val <- t_val * sqrt((n1 + n2) / (n1 * n2))
    j_val <- 1 - (3 / (4 * (n1 + n2) - 9))
    g_val <- d_val * j_val
    
    tagList(
      div(
        class = "kawe-math-card",
        div(class = "kawe-formula-title", icon("square-root-alt"), " Metode: Uji-t Dua Kelompok Independen (Independent Student's t-test)"),
        div(class = "kawe-formula-sub", "Rujukan Baku: Dunlap et al. (1996); Rosenthal (1991); Lipsey & Wilson (2001)."),
        p("Mengonversi statistik uji-t independen langsung ke standardized mean difference Cohen's d dan Hedges' g."),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "1. Landasan Formulasi Matematis:"),
        div(class = "kawe-eq-box",
          HTML("$$d = t \\times \\sqrt{\\frac{n_1 + n_2}{n_1 n_2}} = t \\times \\sqrt{\\frac{1}{n_1} + \\frac{1}{n_2}}$$"),
          HTML("$$J = 1 - \\frac{3}{4(n_1 + n_2) - 9}, \\quad g = J \\times d$$"),
          HTML("$$v_d = \\frac{n_1 + n_2}{n_1 n_2} + \\frac{d^2}{2(n_1 + n_2)}, \\quad v_g = J^2 \\times v_d$$")
        ),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "2. Langkah Perhitungan Terbuka & Substitusi Angka Riil:"),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 1"),
          sprintf("Nilai Statistik t = %.3f | Ukuran Kelompok: n1 = %d, n2 = %d (df = %d).", t_val, n1, n2, df_val)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 2"),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$d = %.3f \\times \\sqrt{\\frac{%d + %d}{(%d)(%d)}} = %.3f \\times \\sqrt{%.5f} = %.4f$$",
                         t_val, n1, n2, n1, n2, t_val, (n1+n2)/(n1*n2), d_val))
          )
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 3"),
          sprintf("Faktor Koreksi J = %.5f -> %s = %.4f (SE = %.4f, 95%% CI [%.4f, %.4f]).",
                  j_val, metric_name, es_val, se_val, ci_lo, ci_hi)
        )
      )
    )
  } else if (itype == "f_test") {
    f_val <- input_list$f_val
    n1 <- input_list$f_n1; n2 <- input_list$f_n2
    t_equiv <- sqrt(f_val)
    if (input_list$f_direction == "neg") t_equiv <- -t_equiv
    d_val <- t_equiv * sqrt((n1 + n2) / (n1 * n2))
    j_val <- 1 - (3 / (4 * (n1 + n2) - 9))
    
    tagList(
      div(
        class = "kawe-math-card",
        div(class = "kawe-formula-title", icon("square-root-alt"), " Metode: Uji-F Satu Arah (One-Way ANOVA df1 = 1)"),
        div(class = "kawe-formula-sub", "Rujukan Baku: Lipsey & Wilson (2001); Borenstein et al. (2009)."),
        p("Memanfaatkan kesetaraan statistik antara distribusi F dengan derajat kebebasan numerator df1 = 1 dan distribusi Student's t: t = sqrt(F)."),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "1. Landasan Formulasi Matematis:"),
        div(class = "kawe-eq-box",
          HTML("$$|t| = \\sqrt{F_{(1, df)}}, \\quad d = \\text{sign} \\times \\sqrt{F \\left(\\frac{n_1 + n_2}{n_1 n_2}\\right)}$$"),
          HTML("$$J = 1 - \\frac{3}{4(n_1 + n_2) - 9}, \\quad g = J \\times d$$")
        ),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "2. Langkah Perhitungan Terbuka & Substitusi Angka Riil:"),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 1"),
          sprintf("Nilai F = %.3f -> Ekivalen t = sqrt(%.3f) = %.4f (Arah: %s).",
                  f_val, f_val, abs(t_equiv), ifelse(input_list$f_direction == "pos", "Memihak Intervensi", "Memihak Kontrol"))
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 2"),
          sprintf("Cohen's d = %.4f | Faktor J = %.5f | %s = %.4f (SE = %.4f, 95%% CI [%.4f, %.4f]).",
                  d_val, j_val, metric_name, es_val, se_val, ci_lo, ci_hi)
        )
      )
    )
  } else if (itype == "p_val") {
    p_val <- input_list$p_value
    n1 <- input_list$p_n1; n2 <- input_list$p_n2
    df_val <- n1 + n2 - 2
    t_equiv <- qt(1 - (p_val / 2), df_val)
    if (input_list$p_direction == "neg") t_equiv <- -t_equiv
    d_val <- t_equiv * sqrt((n1 + n2) / (n1 * n2))
    j_val <- 1 - (3 / (4 * (n1 + n2) - 9))
    
    tagList(
      div(
        class = "kawe-math-card",
        div(class = "kawe-formula-title", icon("square-root-alt"), " Metode: Nilai-p Dua Arah & Ukuran Sampel (Inverse Student's t-Distribution)"),
        div(class = "kawe-formula-sub", "Rujukan Baku: Rosenthal (1991); Lipsey & Wilson (2001)."),
        p("Ketika hanya nilai signifikansi statistik (p-value) dan ukuran sampel yang dilaporkan, nilai t direkonstruksi via fungsi inversi kuantil distribusi Student's t."),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "1. Landasan Formulasi Matematis:"),
        div(class = "kawe-eq-box",
          HTML("$$df = n_1 + n_2 - 2, \\quad t = F_{t, df}^{-1}\\left(1 - \\frac{p}{2}\\right)$$"),
          HTML("$$d = \\text{sign} \\times t \\times \\sqrt{\\frac{n_1 + n_2}{n_1 n_2}}, \\quad g = J \\times d$$")
        ),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "2. Langkah Perhitungan Terbuka & Substitusi Angka Riil:"),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 1"),
          sprintf("Nilai p dua arah = %.4f | Sampel: n1 = %d, n2 = %d (df = %d).", p_val, n1, n2, df_val)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 2"),
          sprintf("Rekonstruksi Nilai t: t = qt(1 - %.4f/2, df=%d) = %.4f.", p_val, df_val, abs(t_equiv))
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 3"),
          sprintf("Cohen's d = %.4f | Faktor J = %.5f | %s = %.4f (SE = %.4f, 95%% CI [%.4f, %.4f]).",
                  d_val, j_val, metric_name, es_val, se_val, ci_lo, ci_hi)
        )
      )
    )
  } else if (itype == "table_2x2") {
    a <- input_list$tab_e1; n1 <- input_list$tab_n1; b <- n1 - a
    c <- input_list$tab_e2; n2 <- input_list$tab_n2; d <- n2 - c
    or_val <- (a * d) / (b * c)
    lor_val <- log(or_val)
    v_lor <- 1/a + 1/b + 1/c + 1/d
    d_chinn <- lor_val * (sqrt(3) / pi)
    vd_chinn <- v_lor * (3 / (pi^2))
    rr_val <- (a / n1) / (c / n2)
    rd_val <- (a / n1) - (c / n2)
    
    tagList(
      div(
        class = "kawe-math-card",
        div(class = "kawe-formula-title", icon("square-root-alt"), " Metode: Tabel Kontingensi 2x2 (Data Kategorikal Biner / Proporsi)"),
        div(class = "kawe-formula-sub", "Rujukan Baku: Chinn (2000); Hasselblad & Hedges (1995); Woolf (1955)."),
        p("Mengonversi rasio odds (Odds Ratio) dari luaran biner ke dalam skala standardized mean difference menggunakan metode transformasi logistik kontinu Chinn (2000)."),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "1. Landasan Formulasi Matematis:"),
        div(class = "kawe-eq-box",
          HTML("$$\\text{Odds Ratio: } OR = \\frac{a \\cdot d}{b \\cdot c}, \\quad \\ln(OR) = \\ln\\left(\\frac{a \\cdot d}{b \\cdot c}\\right)$$"),
          HTML("$$\\text{Varians Log OR: } v_{\\ln(OR)} = \\frac{1}{a} + \\frac{1}{b} + \\frac{1}{c} + \\frac{1}{d}$$"),
          HTML("$$\\text{Transformasi Chinn ke d: } d = \\ln(OR) \\times \\frac{\\sqrt{3}}{\\pi} \\approx \\ln(OR) \\times 0.5513$$"),
          HTML("$$\\text{Varians d: } v_d = v_{\\ln(OR)} \\times \\frac{3}{\\pi^2} \\approx v_{\\ln(OR)} \\times 0.3040$$"),
          HTML("$$\\text{Relative Risk: } RR = \\frac{a / n_1}{c / n_2}, \\quad \\text{Risk Difference: } RD = \\frac{a}{n_1} - \\frac{c}{n_2}$$")
        ),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "2. Langkah Perhitungan Terbuka & Substitusi Angka Riil:"),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 1"),
          sprintf("Matriks 2x2: Intervensi [Event a=%d, Non-event b=%d, Total n1=%d] | Kontrol [Event c=%d, Non-event d=%d, Total n2=%d].",
                  a, b, n1, c, d, n2)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 2"),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$OR = \\frac{(%d)(%d)}{(%d)(%d)} = \\frac{%d}{%d} = %.4f, \\quad \\ln(OR) = %.4f$$",
                         a, d, b, c, a*d, b*c, or_val, lor_val)),
            HTML(sprintf("$$v_{\\ln(OR)} = \\frac{1}{%d} + \\frac{1}{%d} + \\frac{1}{%d} + \\frac{1}{%d} = %.5f$$",
                         a, b, c, d, v_lor))
          )
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 3"),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$d = %.4f \\times \\frac{\\sqrt{3}}{\\pi} = %.4f \\times 0.5513 = %.4f$$",
                         lor_val, lor_val, d_chinn)),
            HTML(sprintf("$$v_d = %.5f \\times \\frac{3}{\\pi^2} = %.5f$$", v_lor, vd_chinn))
          ),
          sprintf("Estimasi Titik %s = %.4f (SE = %.4f, 95%% CI [%.4f, %.4f]) | RR = %.3f | RD = %.3f (%.1f%%).",
                  metric_name, es_val, se_val, ci_lo, ci_hi, rr_val, rd_val, rd_val * 100)
        )
      )
    )
  } else if (itype == "correlation") {
    r_val <- input_list$corr_r; n_val <- input_list$corr_n
    d_val <- (2 * r_val) / sqrt(1 - r_val^2)
    j_val <- 1 - (3 / (4 * n_val - 9))
    g_val <- d_val * j_val
    z_val <- 0.5 * log((1 + r_val) / (1 - r_val))
    vz_val <- 1 / (n_val - 3)
    
    tagList(
      div(
        class = "kawe-math-card",
        div(class = "kawe-formula-title", icon("square-root-alt"), " Metode: Koefisien Korelasi (Pearson r / Point-Biserial)"),
        div(class = "kawe-formula-sub", "Rujukan Baku: Borenstein et al. (2009); Jacobs & Viechtbauer (2017); Fisher (1921)."),
        p("Mengonversi koefisien korelasi ke standardized mean difference d dan melaporkan transformasi Fisher's z."),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "1. Landasan Formulasi Matematis:"),
        div(class = "kawe-eq-box",
          HTML("$$d = \\frac{2r}{\\sqrt{1 - r^2}}, \\quad J = 1 - \\frac{3}{4N - 9}, \\quad g = J \\times d$$"),
          HTML("$$\\text{Fisher's } z = \\frac{1}{2}\\ln\\left(\\frac{1 + r}{1 - r}\\right) = \\text{arctanh}(r), \\quad v_z = \\frac{1}{N - 3}, \\quad SE_z = \\frac{1}{\\sqrt{N - 3}}$$")
        ),
        
        div(class = "fw-bold text-primary mt-3 mb-1", "2. Langkah Perhitungan Terbuka & Substitusi Angka Riil:"),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 1"),
          sprintf("Koefisien Korelasi r = %.3f | Ukuran Sampel Total N = %d.", r_val, n_val)
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 2"),
          div(class = "kawe-eq-box",
            HTML(sprintf("$$d = \\frac{2(%.3f)}{\\sqrt{1 - %.3f^2}} = \\frac{%.3f}{\\sqrt{%.4f}} = %.4f$$",
                         r_val, r_val, 2*r_val, 1 - r_val^2, d_val))
          )
        ),
        div(
          class = "kawe-formula-step",
          tags$span(class = "kawe-step-num", "Langkah 3"),
          sprintf("Faktor J = %.5f -> %s = %.4f (SE = %.4f, 95%% CI [%.4f, %.4f]) | Fisher's z = %.4f (SE_z = %.4f).",
                  j_val, metric_name, es_val, se_val, ci_lo, ci_hi, z_val, sqrt(vz_val))
        )
      )
    )
  }
}

# ==============================================================================
# CONFIRMATORY FACTOR ANALYSIS (CFA) & RASCH (TAM) ENGINES & HELPERS
# ==============================================================================

# Workaround for S7 / ggplot2 4.0 binary operator dispatch if needed
if (!exists("+")) {
  `+` <- function(e1, e2) {
    if (missing(e2)) return(.Primitive("+")(e1))
    if (inherits(e1, "ggplot") || inherits(e1, "gg")) {
      return(ggplot2:::add_ggplot(e1, e2, deparse(substitute(e2))))
    }
    .Primitive("+")(e1, e2)
  }
}

# --- 1. CFA Simulation Function ---
simulate_cfa_dataset <- function(n_obs = 300, n_items = 12, n_factors = 3, factor_cor = 0.35, lambda_min = 0.55, lambda_max = 0.85, n_likert = 5) {
  set.seed(as.integer(Sys.time()) %% 100000 + sample(1:1000, 1))
  Sigma <- matrix(factor_cor, nrow = n_factors, ncol = n_factors)
  diag(Sigma) <- 1.0
  eigen_min <- min(eigen(Sigma)$values)
  if (eigen_min <= 0) {
    diag(Sigma) <- diag(Sigma) - eigen_min + 0.05
  }
  L <- chol(Sigma)
  F_scores <- matrix(rnorm(n_obs * n_factors), nrow = n_obs) %*% L
  
  items_per_fac <- ceiling(n_items / n_factors)
  fac_assignment <- rep(1:n_factors, each = items_per_fac)[1:n_items]
  
  lambdas <- runif(n_items, lambda_min, lambda_max)
  y_cont <- matrix(0, nrow = n_obs, ncol = n_items)
  for (i in 1:n_items) {
    f_idx <- fac_assignment[i]
    eps <- rnorm(n_obs, mean = 0, sd = sqrt(max(0.15, 1 - lambdas[i]^2)))
    y_cont[, i] <- lambdas[i] * F_scores[, f_idx] + eps
  }
  
  cuts <- qnorm(seq(1 / (n_likert + 1), n_likert / (n_likert + 1), length.out = n_likert - 1))
  y_likert <- apply(y_cont, 2, function(col) {
    findInterval(col, cuts) + 1
  })
  colnames(y_likert) <- paste0("Item.", 1:n_items)
  as.data.frame(y_likert)
}

# --- 2. CFA Descriptive Statistics Calculator ---
cfa_descriptive_stats <- function(df) {
  df_num <- df[, sapply(df, is.numeric), drop = FALSE]
  if (ncol(df_num) == 0) return(data.frame())
  desc <- as.data.frame(psych::describe(df_num))
  n_vec <- pmax(4, as.numeric(desc$n))
  se_skew <- sqrt((6 * n_vec * (n_vec - 1)) / ((n_vec - 2) * (n_vec + 1) * (n_vec + 3)))
  se_kurt <- 2 * se_skew * sqrt((n_vec^2 - 1) / ((n_vec - 3) * (n_vec + 5)))
  
  out <- data.frame(
    Butir = rownames(desc),
    N = as.integer(desc$n),
    Mean = round(desc$mean, 3),
    SD = round(desc$sd, 3),
    Median = round(desc$median, 3),
    Min = desc$min,
    Max = desc$max,
    Skewness = round(desc$skew, 3),
    Kurtosis = round(desc$kurtosis, 3),
    SE_Skew = round(se_skew, 3),
    SE_Kurt = round(se_kurt, 3),
    Status_Normalitas = ifelse(abs(desc$skew) <= 2.0 & abs(desc$kurtosis) <= 7.0, 
                               "Normal Univariat (Kline, 2016)", 
                               "Deviasi Non-Normal"),
    stringsAsFactors = FALSE
  )
  out
}

# --- 3. Reliability Calculators ---
calc_cronbach_alpha <- function(df) {
  df <- na.omit(df)
  k <- ncol(df)
  if (k < 2) return(NA_real_)
  item_vars <- apply(df, 2, var)
  total_var <- var(rowSums(df))
  if (is.na(total_var) || total_var <= 0) return(NA_real_)
  alpha <- (k / (k - 1)) * (1 - sum(item_vars) / total_var)
  return(round(alpha, 3))
}

calc_composite_reliability <- function(loadings) {
  loadings <- loadings[!is.na(loadings)]
  if (length(loadings) < 2) return(list(cr = NA_real_, ave = NA_real_))
  sum_l <- sum(loadings)
  sum_theta <- sum(pmax(0.01, 1 - loadings^2))
  cr <- (sum_l^2) / (sum_l^2 + sum_theta)
  ave <- mean(loadings^2)
  list(cr = round(cr, 3), ave = round(ave, 3))
}

# --- 4. Rasch Simulation Function (TAM) ---
simulate_rasch_dataset <- function(n_persons = 300, n_items = 12, n_dims = 1, item_type = "dichotomous", dim_cor = 0.40, n_categories = 4) {
  set.seed(as.integer(Sys.time()) %% 100000 + sample(1:1000, 1))
  
  if (n_dims == 1) {
    thetas <- matrix(rnorm(n_persons, mean = 0, sd = 1), ncol = 1)
  } else {
    Sigma <- matrix(dim_cor, nrow = n_dims, ncol = n_dims)
    diag(Sigma) <- 1.0
    L <- chol(Sigma)
    thetas <- matrix(rnorm(n_persons * n_dims), nrow = n_persons) %*% L
  }
  
  items_per_dim <- ceiling(n_items / n_dims)
  dim_map <- rep(1:n_dims, each = items_per_dim)[1:n_items]
  betas <- seq(-1.8, 1.8, length.out = n_items)
  
  resp_mat <- matrix(0, nrow = n_persons, ncol = n_items)
  
  if (item_type == "dichotomous") {
    for (i in 1:n_items) {
      d <- dim_map[i]
      p_prob <- 1 / (1 + exp(-(thetas[, d] - betas[i])))
      resp_mat[, i] <- rbinom(n_persons, 1, p_prob)
    }
  } else {
    K <- n_categories - 1
    taus <- seq(-1.2, 1.2, length.out = K)
    for (i in 1:n_items) {
      d <- dim_map[i]
      th <- thetas[, d]
      b <- betas[i]
      numerators <- matrix(0, nrow = n_persons, ncol = n_categories)
      for (k in 1:K) {
        step_thresholds <- b + taus[1:k]
        numerators[, k + 1] <- k * th - sum(step_thresholds)
      }
      max_num <- apply(numerators, 1, max)
      exp_num <- exp(numerators - max_num)
      probs <- exp_num / rowSums(exp_num)
      for (p in 1:n_persons) {
        resp_mat[p, i] <- sample(0:K, size = 1, prob = probs[p, ])
      }
    }
  }
  
  colnames(resp_mat) <- paste0("Item.", 1:n_items)
  as.data.frame(resp_mat)
}


# UI Definition
ui <- page_navbar(
  id = "main_navbar",
  window_title = "Kawe: Quantitative and Psychometrics",
  fillable = FALSE,
  title = tags$div(
    class = "kawe-brand-horizontal"
  ),
  theme = bs_theme(
    version = 5,
    bg = "#FFFFFF",
    fg = "#231618",
    primary = "#650915",
    secondary = "#C5A059",
    success = "#650915",
    info = "#C5A059",
    warning = "#800020",
    danger = "#4A050D"
  ),
  header = tagList(
    tags$head(
      tags$link(rel = "shortcut icon", href = "favicon.svg", type = "image/svg+xml"),
      tags$style(HTML("
        /* Google Fonts */
        @import url('https://fonts.googleapis.com/css2?family=Cinzel:wght@500;600;700&family=Cormorant+Garamond:ital,wght@0,400;0,500;0,600;0,700;1,400;1,600&family=EB+Garamond:ital,wght@0,400;0,500;0,600;1,400&family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

        :root {
          /* Palet Warna KAWE Academic Press */
          --kawe-maroon: #650915;
          --kawe-maroon-deep: #4A050D;
          --kawe-maroon-light: #7E0D1B;
          --kawe-gold: #C5A059;
          --kawe-gold-light: #E5D3B3;
          --kawe-gold-trans: rgba(197, 160, 89, 0.35);
          --kawe-bg-white: #FFFFFF;
          --kawe-bg-paper: #FAF8F5;
          --kawe-text-main: #231618;
          --kawe-text-muted: #5A4B4C;

          /* Tipografi */
          --font-logo: 'Cinzel', 'Cormorant Garamond', Georgia, serif;
          --font-serif: 'Cormorant Garamond', 'EB Garamond', Georgia, serif;
          --font-sans: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
          --font-mono: 'JetBrains Mono', 'Fira Code', monospace;
        }

        body {
          background-color: var(--kawe-bg-white);
          color: var(--kawe-text-main);
          font-family: var(--font-sans);
          line-height: 1.65;
          scrollbar-width: thin;
          scrollbar-color: var(--kawe-gold) var(--kawe-bg-paper);
        }

        /* Headings */
        h1, h2, h3, .brand-title {
          font-family: var(--font-logo) !important;
          color: var(--kawe-maroon) !important;
          text-transform: uppercase;
          letter-spacing: 0.18em;
          font-weight: 600;
        }

        h4, h5, h6 {
          font-family: var(--font-serif) !important;
          color: var(--kawe-maroon) !important;
          letter-spacing: 0.06em;
          font-weight: 600;
        }

        .brand-subtitle {
          font-family: var(--font-serif) !important;
          color: var(--kawe-gold) !important;
          text-transform: uppercase;
          letter-spacing: 0.15em;
          font-size: 0.88rem;
          font-weight: 600;
        }

        /* Dividers */
        hr, .gold-divider {
          border: 0 !important;
          border-top: 1px solid var(--kawe-gold) !important;
          opacity: 0.45;
          margin: 1.25rem 0;
        }

        /* Navbar Layout & Styling (Full Flat Maroon Header) */
        .navbar {
          position: sticky !important;
          top: 0 !important;
          z-index: 1050 !important;
          background-color: var(--kawe-maroon) !important;
          border-bottom: 2px solid var(--kawe-gold) !important;
          box-shadow: 0 4px 14px rgba(43, 24, 26, 0.25);
          padding: 0.35rem 1.15rem !important;
        }

        /* Navbar Brand Horizontal (Icon Left, Divider, Text Centered Right) */
        .navbar-brand {
          display: flex !important;
          align-items: center !important;
          padding: 0.15rem 0 !important;
          margin-right: 1.25rem !important;
        }

        .kawe-brand-horizontal {
          display: inline-flex;
          align-items: center;
          justify-content: center;
          gap: 0.85rem; /* Jarak proporsional antara ikon, garis emas, dan grup teks */
          padding: 0.10rem 0.25rem;
          background-color: transparent;
        }

        .kawe-icon {
          height: 38px; /* Ukuran proporsional di tengah vertikal header navbar */
          width: auto;
          object-fit: contain;
          display: block;
          margin: 0 auto;
          filter: drop-shadow(0 1px 3px rgba(0, 0, 0, 0.40));
        }

        .kawe-divider {
          height: 30px;
          width: 1px;
          background-color: #C5A059; /* Antique Gold */
          opacity: 0.70;
          align-self: center;
        }

        .kawe-text-group {
          display: flex;
          flex-direction: column;
          justify-content: center;
          align-items: center;
          text-align: center;
        }

        .kawe-title {
          font-family: 'Cinzel', 'Cormorant Garamond', serif;
          font-size: 1.28rem;
          font-weight: 600;
          color: #FAF8F5 !important; /* Crisp Ivory White pada Maroon */
          letter-spacing: 0.36em;
          padding-left: 0.36em; /* Kompensasi trailing letter-spacing agar KAWE tepat di tengah */
          line-height: 1.0;
          margin-bottom: 0.12rem;
          text-align: center;
          display: block;
          width: 100%;
        }

        .kawe-subtitle {
          font-family: 'Cormorant Garamond', serif;
          font-size: 0.58rem;
          font-weight: 600;
          color: #E5D3B3 !important; /* Soft Antique Gold */
          letter-spacing: 0.12em;
          padding-left: 0.12em;
          line-height: 1.1;
          text-transform: uppercase;
          white-space: nowrap;
          text-align: center;
          display: block;
          width: 100%;
        }

        /* Hero Card Brand Override (Light Background) */
        .kawe-hero-card .kawe-brand-horizontal {
          display: flex;
          justify-content: center;
          align-items: center;
          width: fit-content;
          margin: 0 auto;
        }

        .kawe-hero-card .kawe-title {
          color: var(--kawe-maroon) !important;
        }

        .kawe-hero-card .kawe-subtitle {
          color: var(--kawe-text-main) !important;
        }

        .navbar-nav {
          flex-wrap: wrap !important;
          overflow: visible !important;
        }

        .navbar .nav-link, .navbar-nav > li > a {
          font-family: var(--font-sans);
          font-size: 0.82rem !important;
          font-weight: 500;
          color: #EAE0D8 !important; /* Soft ivory text pada latar marun */
          padding: 0.50rem 0.75rem !important;
          border-bottom: 2px solid transparent;
          white-space: nowrap !important;
          transition: all 0.2s ease;
        }

        .navbar .nav-link:hover, .navbar-nav > li > a:hover {
          color: #FFFFFF !important;
          border-bottom: 2px solid var(--kawe-gold) !important;
          background-color: rgba(255, 255, 255, 0.08) !important;
        }

        .navbar .nav-link.active, .navbar-nav > li.active > a {
          color: #FFFFFF !important;
          font-weight: 600;
          border-bottom: 2px solid var(--kawe-gold) !important;
          background-color: rgba(0, 0, 0, 0.20) !important;
          border-radius: 3px 3px 0 0;
        }

        /* Navbar Dropdown Menu (Tab Meta-Analisis ke bawah) */
        .navbar .dropdown-menu {
          background-color: var(--kawe-maroon) !important;
          border: 1px solid var(--kawe-gold) !important;
          border-radius: 4px !important;
          box-shadow: 0 8px 24px rgba(35, 22, 24, 0.35) !important;
          padding: 0.35rem 0 !important;
          margin-top: 0.25rem !important;
        }

        .navbar .dropdown-menu .dropdown-item {
          font-family: var(--font-sans) !important;
          font-size: 0.82rem !important;
          font-weight: 500 !important;
          color: #EAE0D8 !important;
          padding: 0.55rem 1.15rem !important;
          display: flex !important;
          align-items: center !important;
          gap: 0.55rem !important;
          transition: all 0.18s ease;
        }

        .navbar .dropdown-menu .dropdown-item:hover,
        .navbar .dropdown-menu .dropdown-item:focus {
          background-color: var(--kawe-maroon-deep) !important;
          color: #FFFFFF !important;
          border-left: 3px solid var(--kawe-gold) !important;
        }

        .navbar .dropdown-menu .dropdown-item.active {
          background-color: rgba(197, 160, 89, 0.22) !important;
          color: #FFFFFF !important;
          font-weight: 600 !important;
          border-left: 3px solid var(--kawe-gold) !important;
        }

        /* Sub-tabs inside Cards (navset_card_tab) */
        .card-header .nav-link, .nav-tabs .nav-link, .nav-pills .nav-link {
          font-family: var(--font-sans) !important;
          font-size: 0.82rem !important;
          font-weight: 500;
          color: var(--kawe-text-muted) !important;
          padding: 0.45rem 0.85rem !important;
          border-bottom: 2px solid transparent !important;
          transition: all 0.2s ease;
        }

        .card-header .nav-link:hover, .nav-tabs .nav-link:hover {
          color: var(--kawe-maroon) !important;
          background-color: rgba(101, 9, 21, 0.04) !important;
        }

        .card-header .nav-link.active, .nav-tabs .nav-link.active, .nav-pills .nav-link.active {
          color: var(--kawe-maroon) !important;
          font-weight: 600 !important;
          border-bottom: 2px solid var(--kawe-maroon) !important;
          background-color: rgba(101, 9, 21, 0.07) !important;
          border-radius: 4px 4px 0 0;
        }

        /* Cards & Containers */
        .card {
          background-color: var(--kawe-bg-paper) !important;
          border: 1px solid var(--kawe-gold-trans) !important;
          border-radius: 4px !important;
          box-shadow: none !important;
        }

        .card-header {
          background-color: #FAF8F5 !important;
          border-bottom: 1px solid var(--kawe-gold-trans) !important;
          color: var(--kawe-maroon) !important;
          font-family: var(--font-serif);
          font-weight: 600;
          letter-spacing: 0.04em;
          padding: 0.65rem 1rem !important;
        }

        /* Value Boxes Academic Press Theming */
        .bslib-value-box, .value-box {
          border: 1px solid var(--kawe-gold-trans) !important;
          border-radius: 4px !important;
          box-shadow: none !important;
        }

        .bslib-value-box.bg-primary {
          background-color: var(--kawe-maroon) !important;
          color: #FAF8F5 !important;
        }

        .bslib-value-box.bg-primary .value-box-title {
          color: #E5D3B3 !important;
          font-family: var(--font-serif) !important;
          letter-spacing: 0.08em;
          text-transform: uppercase;
          font-size: 0.72rem;
        }

        .bslib-value-box.bg-primary .value-box-value {
          color: #FAF8F5 !important;
          font-family: var(--font-logo) !important;
          font-weight: 600;
        }

        .bslib-value-box.bg-secondary, .bslib-value-box.bg-info, .bslib-value-box.bg-warning, .bslib-value-box.bg-success {
          background-color: #FFFFFF !important;
          border: 1px solid var(--kawe-gold-trans) !important;
          color: var(--kawe-maroon) !important;
        }

        .bslib-value-box.bg-secondary .value-box-title,
        .bslib-value-box.bg-info .value-box-title,
        .bslib-value-box.bg-warning .value-box-title,
        .bslib-value-box.bg-success .value-box-title {
          color: var(--kawe-gold) !important;
          font-family: var(--font-serif) !important;
          letter-spacing: 0.08em;
          text-transform: uppercase;
          font-size: 0.72rem;
          font-weight: 600;
        }

        .bslib-value-box.bg-secondary .value-box-value,
        .bslib-value-box.bg-info .value-box-value,
        .bslib-value-box.bg-warning .value-box-value,
        .bslib-value-box.bg-success .value-box-value {
          color: var(--kawe-maroon) !important;
          font-family: var(--font-logo) !important;
          font-weight: 600;
        }

        .value-box-showcase i {
          color: rgba(197, 160, 89, 0.65) !important;
        }

        /* Form Controls */
        .form-control, .form-select {
          border: 1px solid #D6C2A0 !important;
          border-radius: 4px;
          font-size: 0.88rem;
          color: var(--kawe-text-main);
        }

        .form-control:focus, .form-select:focus {
          border-color: var(--kawe-maroon) !important;
          box-shadow: 0 0 0 0.2rem rgba(101, 9, 21, 0.15) !important;
        }

        /* Buttons */
        .btn-primary {
          background-color: var(--kawe-maroon) !important;
          border: 1px solid var(--kawe-gold) !important;
          color: #FFFFFF !important;
          font-weight: 500;
          letter-spacing: 0.03em;
          border-radius: 4px;
          transition: all 0.2s ease;
        }

        .btn-primary:hover, .btn-primary:focus {
          background-color: var(--kawe-maroon-deep) !important;
          border-color: #FFFFFF !important;
          box-shadow: 0 0 0 0.15rem rgba(197, 160, 89, 0.35);
        }

        .btn-outline-primary {
          color: var(--kawe-maroon) !important;
          border-color: var(--kawe-maroon) !important;
          background-color: transparent !important;
          border-radius: 4px;
        }

        .btn-outline-primary:hover {
          background-color: var(--kawe-maroon) !important;
          color: #FFFFFF !important;
        }

        .btn-outline-secondary {
          color: var(--kawe-text-muted) !important;
          border-color: var(--kawe-gold) !important;
          background-color: transparent !important;
        }

        .btn-outline-secondary:hover {
          background-color: var(--kawe-gold) !important;
          color: #FFFFFF !important;
        }

        /* Tables & DataTables */
        table.dataTable thead th {
          background-color: var(--kawe-maroon) !important;
          color: #FFFFFF !important;
          font-family: var(--font-serif) !important;
          letter-spacing: 0.04em;
          border-bottom: 2px solid var(--kawe-gold) !important;
          font-size: 0.86rem;
        }

        table.dataTable tbody tr:hover {
          background-color: #F8F4EE !important;
        }

        /* Hero & Pillar Cards */
        .kawe-hero-card {
          background: linear-gradient(180deg, #FFFFFF 0%, #FAF8F5 100%);
          border: 1px solid var(--kawe-gold-trans);
          border-radius: 4px;
          padding: 2.2rem 1.8rem;
          text-align: center;
        }

        .kawe-pillar-card {
          background-color: #FFFFFF;
          border: 1px solid var(--kawe-gold-trans);
          border-top: 3px solid var(--kawe-maroon);
          border-radius: 4px;
          padding: 1.35rem;
          height: 100%;
          transition: all 0.25s ease;
        }

        .kawe-pillar-card:hover {
          border-top-color: var(--kawe-gold);
          box-shadow: 0 3px 12px rgba(101, 9, 21, 0.08);
        }

        .active-data-pill {
          background-color: #4A050D;
          border: 1px solid #C5A059;
          color: #FAF8F5;
          padding: 0.25rem 0.70rem;
          border-radius: 3px;
          font-size: 0.72rem;
          letter-spacing: 0.04em;
          display: inline-flex;
          align-items: center;
          gap: 0.4rem;
        }

        /* Global Active Dataset & Model Status Bar */
        .kawe-global-bar {
          background-color: #FAF8F5;
          border-bottom: 1px solid var(--kawe-gold-trans);
          padding: 0.40rem 1.25rem;
          font-size: 0.78rem;
          display: flex;
          align-items: center;
          justify-content: space-between;
          flex-wrap: wrap;
          gap: 0.5rem;
        }

        .kawe-clinical-pill {
          background-color: #FFFFFF;
          border: 1px solid var(--kawe-gold-trans);
          border-left: 3px solid var(--kawe-maroon);
          padding: 0.45rem 0.75rem;
          border-radius: 4px;
          font-size: 0.80rem;
        }

        .kawe-step-footer {
          background-color: #FAF8F5;
          border: 1px solid var(--kawe-gold-trans);
          border-radius: 4px;
          padding: 0.65rem 1rem;
          margin-top: 0.85rem;
          display: flex;
          justify-content: space-between;
          align-items: center;
          flex-wrap: wrap;
          gap: 0.5rem;
        }

        /* MathJax & Formula Cards */
        .kawe-math-card {
          background-color: #FFFFFF;
          border: 1px solid var(--kawe-gold-trans);
          border-left: 4px solid var(--kawe-maroon);
          border-radius: 4px;
          padding: 1.15rem 1.35rem;
          margin-bottom: 1.15rem;
          box-shadow: 0 1px 4px rgba(0,0,0,0.03);
        }

        .kawe-formula-step {
          background-color: #FAF8F5;
          border: 1px solid rgba(197, 160, 89, 0.25);
          border-radius: 4px;
          padding: 0.95rem 1.15rem;
          margin-bottom: 0.85rem;
        }

        .kawe-formula-title {
          font-family: var(--font-serif);
          color: var(--kawe-maroon);
          font-weight: 700;
          font-size: 1.05rem;
          margin-bottom: 0.35rem;
          letter-spacing: 0.04em;
        }

        .kawe-formula-sub {
          color: var(--kawe-text-muted);
          font-size: 0.84rem;
          margin-bottom: 0.65rem;
        }

        .kawe-eq-box {
          background-color: #FFFFFF;
          border-radius: 4px;
          padding: 0.75rem 1rem;
          margin: 0.5rem 0;
          border: 1px dashed var(--kawe-gold-trans);
          font-size: 0.96rem;
          overflow-x: auto;
        }

        .kawe-step-num {
          display: inline-block;
          background-color: var(--kawe-maroon);
          color: #FAF8F5;
          font-size: 0.72rem;
          font-weight: 600;
          padding: 0.15rem 0.55rem;
          border-radius: 3px;
          margin-right: 0.45rem;
        }
      ")),
      withMathJax(),
      tags$script(type = "text/x-mathjax-config", HTML("
        MathJax.Hub.Config({
          tex2jax: {
            inlineMath: [['$','$'], ['\\\\(','\\\\)']],
            displayMath: [['$$','$$'], ['\\\\[','\\\\]']],
            processEscapes: true
          }
        });
      ")),
      tags$script(HTML("
        Shiny.addCustomMessageHandler('copy_to_clipboard', function(message) {
          if (navigator.clipboard && message.text) {
            navigator.clipboard.writeText(message.text);
          }
        });
        if ('scrollRestoration' in history) {
          history.scrollRestoration = 'manual';
        }
        function keepScrollAtTop() {
          window.scrollTo(0, 0);
          document.documentElement.scrollTop = 0;
          document.body.scrollTop = 0;
        }
        window.addEventListener('beforeunload', keepScrollAtTop);
        document.addEventListener('DOMContentLoaded', keepScrollAtTop);
        $(document).ready(keepScrollAtTop);
        $(document).on('shiny:connected', keepScrollAtTop);
        $(document).on('shown.bs.tab', function(e) {
          keepScrollAtTop();
        });
        $(document).on('shiny:value', function(event) {
          if (window.MathJax) {
            MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
          }
        });
      "))
    ),
    uiOutput("global_status_banner")
  ),
  
  # Tab 1: Research Hub & Academic Portal
  nav_panel(
    value = "tab_hub",
    icon = icon("university"),
    title = "Home",
    withMathJax(),
    
    # Profil
    div(
      class = "kawe-hero-card w-100 my-3",
      style = "width: 100%; max-width: 100%;",
      div(
        class = "row align-items-center w-100 m-0", # Sebaris, sejajar tengah, & lebar penuh
        
        # Bagian Kiri (Lebar 3 dari 12 kolom = 25% agar ruang kanan maksimal)
        div(
          class = "col-md-3 kawe-brand-horizontal mb-3 mb-md-0",
          style = "gap: 1.5rem; justify-content: center;",
          img(
            src = "KW.png",
            alt = "KAWE Pic",
            style = "max-height: 420px; width: 100%; object-fit: contain; filter: drop-shadow(0 2px 8px rgba(101, 9, 21, 0.12));"
          )
        ),
        
        # Bagian Kanan (Lebar 9 dari 12 kolom = 75%)
        div(
          class = "col-md-9",
          
          # Judul & Sub-judul di sebelah kanan (lebar penuh rata kiri)
          div(
            class = "kawe-text-group w-100 mb-3",
            style = "display: flex; flex-direction: column; align-items: flex-start;",
            span(class = "kawe-title", style = "font-size: 1.8rem;", "KAWE : QUANTITATIVE AND PSYCHOMETRICS"),
            span(class = "kawe-subtitle mt-1", style = "color: #C5A059; font-weight: 600; letter-spacing: 0.18em;", "INTEGRATED PLATFORM FOR LATENT MODELING, PSYCHOMETRICS, & EVIDENCE SYNTHESIS")
          ),
          
          # Garis emas merentang 100% penuh
          div(class = "gold-divider w-100 my-3"),
          
          # Paragraf rata kanan-kiri (justify) & tanpa batas max-width
          p(
            class = "lead text-dark w-100 mb-0",
            style = "width: 100%; max-width: none; font-size: 1.02rem; line-height: 1.75; text-align: justify;",
            "The KAWE emblem unifies three foundational pillars of our discipline: the classical Greek letter Psi representing psychology, the Gaussian normal distribution curve symbolizing statistical inference, and the logistic S-curve of Item Response Theory (IRT) reflecting modern psychometric precision. Rendered in Deep Maroon and Antique Gold, the dual-toned geometry embodies our core mission—bridging observable behavioral data (manifest variables) with underlying human attributes (latent traits) through rigorous, transparent, and gold-standard quantitative methodology.",
          )
        )
      )
    ),
    
    # Modular Academic Research Workflow Portal
    div(
      class = "kawe-hero-card my-4 p-4",
      style = "background: #FFFFFF; border: 1px solid #E5D3B3; border-radius: 14px; box-shadow: 0 4px 20px rgba(101, 9, 21, 0.05);",
      div(
        class = "row g-4 mt-1",
        
        # 1. Analisis Faktor Konfirmatori (CFA)
        div(
          class = "col-lg-4 col-md-12",
          div(
            class = "card h-100 border p-3 d-flex flex-column",
            style = "border-radius: 12px; border-top: 5px solid #650915 !important; background-color: #FAF8F5;",
            div(
              class = "d-flex justify-content-between align-items-start mb-3 w-100",
              div(
                class = "d-flex align-items-center gap-3",
                div(
                  style = "width: 46px; height: 46px; border-radius: 10px; background: #650915; color: #FFFFFF; display: flex; align-items: center; justify-content: center; font-size: 1.25rem; flex-shrink: 0;",
                  icon("project-diagram")
                ),
                div(class = "fw-bold text-dark", style = "font-size: 1.05rem;", "Analisis Faktor (CFA)")
              ),
              span(class = "badge text-white shadow-sm", style = "background-color: #650915; font-size: 0.75rem; padding: 6px 10px; border-radius: 6px;", "Lavaan Engine")
            ),
            p(
              class = "small text-muted mb-3",
              style = "line-height: 1.6; min-height: 58px;",
              "Evaluasi dimensionalitas instrumen, muatan faktor terstandardisasi (λ), eliminasi butir redundan, dan komparasi kelayakan model (Uni, Multi, Second-Order, Bi-Factor)."
            ),
            div(
              class = "d-flex flex-wrap justify-content-center gap-1 mb-3 text-center",
              span(class = "badge bg-white text-secondary border", "Simulasi Data"),
              span(class = "badge bg-white text-secondary border", "Eliminasi Butir"),
              span(class = "badge bg-white text-secondary border", "Path Diagram"),
              span(class = "badge bg-white text-secondary border", "Fit Indices")
            ),
            div(
              class = "d-flex flex-column gap-2 mt-auto",
              actionButton("hero_btn_cfa_model", "Mulai Analisis CFA", class = "btn text-white fw-bold py-2 shadow-sm", style = "background-color: #650915; border-color: #650915;", icon = icon("arrow-circle-right")),
              actionButton("hero_btn_cfa_comp", "Perbandingan 4 Model Fit", class = "btn btn-outline-dark btn-sm py-2", icon = icon("balance-scale"))
            )
          )
        ),
        
        # 2. Pemodelan Rasch & IRT (TAM)
        div(
          class = "col-lg-4 col-md-12",
          div(
            class = "card h-100 border p-3 d-flex flex-column",
            style = "border-radius: 12px; border-top: 5px solid #C5A059 !important; background-color: #FAF8F5;",
            div(
              class = "d-flex justify-content-between align-items-start mb-3 w-100",
              div(
                class = "d-flex align-items-center gap-3",
                div(
                  style = "width: 46px; height: 46px; border-radius: 10px; background: #C5A059; color: #FFFFFF; display: flex; align-items: center; justify-content: center; font-size: 1.25rem; flex-shrink: 0;",
                  icon("ruler-combined")
                ),
                div(class = "fw-bold text-dark", style = "font-size: 1.05rem;", "Pemodelan Rasch / IRT")
              ),
              span(class = "badge text-dark shadow-sm border", style = "background-color: #E5D3B3; font-size: 0.75rem; padding: 6px 10px; border-radius: 6px;", "TAM Engine")
            ),
            p(
              class = "small text-muted mb-3",
              style = "line-height: 1.6; min-height: 58px;",
              "Kalibrasi objektif tingkat kesulitan butir (\\delta) & kemampuan person (\\theta), evaluasi Infit/Outfit MNSQ, Thurstonian Thresholds, dan Wright Map multi-panel resolusi tinggi."
            ),
            div(
              class = "d-flex flex-wrap justify-content-center gap-1 mb-3 text-center",
              span(class = "badge bg-white text-secondary border", "Infit / Outfit"),
              span(class = "badge bg-white text-secondary border", "Thurstonian"),
              span(class = "badge bg-white text-secondary border", "Wright Map 4K"),
              span(class = "badge bg-white text-secondary border", "RSM vs PCM")
            ),
            div(
              class = "d-flex flex-column gap-2 mt-auto",
              actionButton("hero_btn_rasch_model", "Mulai Pemodelan Rasch", class = "btn text-white fw-bold py-2 shadow-sm", style = "background-color: #8B6B23; border-color: #8B6B23;", icon = icon("arrow-circle-right")),
              actionButton("hero_btn_rasch_comp", "Komparasi Fit RSM, PCM, MRSM, MPCM", class = "btn btn-outline-dark btn-sm py-2", icon = icon("balance-scale"))
            )
          )
        ),
        
        # 3. Meta-Analisis Klinis
        div(
          class = "col-lg-4 col-md-12",
          div(
            class = "card h-100 border p-3 d-flex flex-column",
            style = "border-radius: 12px; border-top: 5px solid #231618 !important; background-color: #FAF8F5;",
            div(
              class = "d-flex justify-content-between align-items-start mb-3 w-100",
              div(
                class = "d-flex align-items-center gap-3",
                div(
                  style = "width: 46px; height: 46px; border-radius: 10px; background: #231618; color: #FFFFFF; display: flex; align-items: center; justify-content: center; font-size: 1.25rem; flex-shrink: 0;",
                  icon("chart-pie")
                ),
                div(class = "fw-bold text-dark", style = "font-size: 1.05rem;", "Meta-Analisis Klinis")
              ),
              span(class = "badge text-white shadow-sm", style = "background-color: #231618; font-size: 0.75rem; padding: 6px 10px; border-radius: 6px;", "Metafor & ESC")
            ),
            p(
              class = "small text-muted mb-3",
              style = "line-height: 1.6; min-height: 58px;",
              "Kalkulasi besaran efek terstandardisasi (g/d/r/OR), sintesis bukti random-effects, analisis heterogenitas (τ², I²), forest plot interaktif, meta-regresi, dan diagnostik bias publikasi."
            ),
            div(
              class = "d-flex flex-wrap justify-content-center gap-1 mb-3 text-center",
              span(class = "badge bg-white text-secondary border", "Effect Size"),
              span(class = "badge bg-white text-secondary border", "Forest Plot"),
              span(class = "badge bg-white text-secondary border", "Meta-Regresi"),
              span(class = "badge bg-white text-secondary border", "Funnel & Fail-Safe")
            ),
            div(
              class = "d-flex flex-column gap-2 mt-auto",
              actionButton("hero_btn_calc", "Kalkulator Effect Size", class = "btn text-white fw-bold py-2 shadow-sm", style = "background-color: #231618; border-color: #231618;", icon = icon("calculator")),
              actionButton("hero_btn_data", "Model Meta-Analisis & Forest Plot", class = "btn btn-outline-dark btn-sm py-2", icon = icon("chart-pie"))
            )
          )
        )
      )
    ),
    # Academic Press Footer
    div(
      class = "text-center py-4 my-2 border-top",
      style = "border-color: #E5D3B3 !important;",
      div(
        class = "brand-title text-muted mb-1",
        style = "font-size: 0.90rem; letter-spacing: 0.18em;",
        "KAWE | QUANTITATIVE & PSYCHOMETRIC RESEARCH"
      ),
      div(
        class = "small text-muted",
        "Aditya Kusumawijaya · Yogyakarta, Indonesia · Versi 0.11"
      )
    )
  ),
  
  # ============================================================================
  # DROPDOWN MENU: CONFIRMATORY FACTOR ANALYSIS (CFA)
  # ============================================================================
  nav_menu(
    title = "CFA",
    icon = icon("project-diagram"),
    
    # 2.1 Model & Simulasi CFA
    nav_panel(
      title = "Model & Simulasi",
      value = "tab_cfa_model",
      icon = icon("sliders-h"),
      layout_sidebar(
        fillable = FALSE,
        sidebar = sidebar(
          width = 360,
          title = "Konfigurasi CFA & Simulasi",
          h6(class = "text-primary fw-bold mb-2", icon("database"), " Sumber Data Kuesioner:"),
          radioButtons(
            "cfa_source_type",
            NULL,
            choices = c(
              "Simulasi Data CFA" = "sim",
              "Dataset Bawaan Politomi Likert (4 Dimensi, k=20)" = "builtin",
              "Unggah Berkas (CSV / Excel)" = "upload"
            ),
            selected = "sim"
          ),
          conditionalPanel(
            condition = "input.cfa_source_type == 'sim'",
            div(
              class = "p-2 bg-white rounded border mb-2",
              sliderInput("cfa_sim_n", "Jumlah Responden (N):", min = 100, max = 1500, value = 300, step = 50),
              sliderInput("cfa_sim_k", "Total Butir Kuesioner (k):", min = 6, max = 30, value = 12, step = 1),
              sliderInput("cfa_sim_factors", "Jumlah Dimensi Laten:", min = 1, max = 5, value = 3, step = 1),
              sliderInput("cfa_sim_cor", "Korelasi Antardimensi:", min = 0.0, max = 0.70, value = 0.35, step = 0.05),
              sliderInput("cfa_sim_likert", "Kategori Skala Likert (1 - m):", min = 4, max = 7, value = 5, step = 1),
              actionButton("btn_cfa_gen_sim", "Bangkitkan Data Simulasi", class = "btn-primary w-100", icon = icon("dice"))
            )
          ),
          conditionalPanel(
            condition = "input.cfa_source_type == 'upload'",
            div(
              class = "p-2 bg-white rounded border mb-2",
              fileInput("cfa_upload_file", "Pilih Berkas (.csv / .xlsx):", accept = c(".csv", ".xlsx", ".xls"))
            )
          ),
          hr(class = "my-2"),
          h6(class = "text-primary fw-bold mb-2", icon("project-diagram"), " Spesifikasi Model CFA:"),
          selectInput(
            "cfa_model_structure",
            "Pilih Struktur Model:",
            choices = c(
              "Unidimensionalitas (1 Faktor Umum)" = "unidim",
              "Multidimensionalitas (Kustom Dimensi)" = "multidim",
              "Second-Order (Faktor Orde Kedua)" = "second_order",
              "Bi-Factor (General + Specific Factors)" = "bifactor"
            ),
            selected = "multidim"
          ),
          conditionalPanel(
            condition = "input.cfa_model_structure == 'multidim' || input.cfa_model_structure == 'second_order' || input.cfa_model_structure == 'bifactor'",
            radioButtons(
              "cfa_dim_assign_mode",
              "Metode Pengisian Dimensi:",
              choices = c(
                "Alokasi Butir Otomatis Rata" = "auto",
                "Pilih Butir per Dimensi Sendiri" = "custom_assign",
                "Ketik Sintaks Lavaan Bebas" = "syntax"
              ),
              selected = "auto"
            ),
            conditionalPanel(
              condition = "input.cfa_dim_assign_mode == 'custom_assign'",
              sliderInput("cfa_n_dims_custom", "Tentukan Jumlah Dimensi:", min = 2, max = 5, value = 3, step = 1),
              uiOutput("cfa_custom_dim_inputs_ui")
            ),
            conditionalPanel(
              condition = "input.cfa_dim_assign_mode == 'syntax'",
              textAreaInput("cfa_custom_syntax_text", "Sintaks Model Lavaan:", value = "F1 =~ Item.1 + Item.2 + Item.3 + Item.4\nF2 =~ Item.5 + Item.6 + Item.7 + Item.8\nF3 =~ Item.9 + Item.10 + Item.11 + Item.12", rows = 4)
            )
          ),
          selectInput(
            "cfa_estimator",
            "Estimator Parameter (Lavaan):",
            choices = list(
              "Family: Maximum Likelihood (ML)" = c(
                "ML - Standard Maximum Likelihood" = "ML",
                "MLR - Robust ML (Huber-White / Yuan-Bentler)" = "MLR",
                "MLM - Satorra-Bentler Scaled (Robust SE)" = "MLM",
                "MLMV - Mean & Variance-Adjusted ML" = "MLMV",
                "MLMVS - Scaled & Shifted ML" = "MLMVS",
                "MLF - First-Order Derivatives Standard Errors" = "MLF"
              ),
              "Family: Weighted Least Squares (WLS / DWLS)" = c(
                "WLSMV - Robust DWLS Mean- & Variance-Adjusted (Gold Standard Ordinal/Likert)" = "WLSMV",
                "WLSM - Robust DWLS Mean-Adjusted" = "WLSM",
                "WLSMVS - Scaled & Shifted DWLS" = "WLSMVS",
                "DWLS - Diagonally Weighted Least Squares" = "DWLS",
                "WLS - Weighted Least Squares / ADF (Browne, 1984)" = "WLS"
              ),
              "Family: Unweighted Least Squares (ULS)" = c(
                "ULSMV - Robust ULS Mean- & Variance-Adjusted" = "ULSMV",
                "ULSM - Robust ULS Mean-Adjusted" = "ULSM",
                "ULSMVS - Scaled & Shifted ULS" = "ULSMVS",
                "ULS - Standard Unweighted Least Squares" = "ULS"
              ),
              "Family: Generalized Least Squares & Other" = c(
                "GLS - Generalized Least Squares" = "GLS",
                "PML - Pairwise Maximum Likelihood (Categorical)" = "PML"
              )
            ),
            selected = "MLR"
          ),
          actionButton("btn_run_cfa", "Jalankan Estimasi CFA", class = "btn-success w-100 mt-2", icon = icon("play-circle"))
        ),
        div(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Status Model",
              value = textOutput("vb_cfa_status"),
              showcase = icon("check-circle"),
              theme = "primary"
            ),
            value_box(
              title = "Total Butir (k)",
              value = textOutput("vb_cfa_k"),
              showcase = icon("list-ol"),
              theme = "secondary"
            ),
            value_box(
              title = "Ukuran Sampel (N)",
              value = textOutput("vb_cfa_n"),
              showcase = icon("users"),
              theme = "info"
            ),
            value_box(
              title = "Estimator Terpilih",
              value = textOutput("vb_cfa_estimator"),
              showcase = icon("cogs"),
              theme = "warning"
            )
          ),
          hr(class = "my-3"),
          card(
            class = "mb-3",
            card_header(class = "bg-light d-flex justify-content-between align-items-center",
              span(class = "fw-bold text-primary", icon("code"), " Formulasi Sintaks Model Lavaan yang Digunakan"),
              span(class = "badge bg-primary", "Sintaks Otomatis")
            ),
            card_body(
              verbatimTextOutput("cfa_syntax_preview")
            )
          ),
          card(
            full_screen = TRUE,
            card_header(class = "bg-light d-flex justify-content-between align-items-center",
              span(class = "fw-bold text-primary", icon("table"), " Matriks Respons Kuesioner"),
              uiOutput("cfa_data_badge_ui")
            ),
            card_body(
              style = "overflow: visible; padding-bottom: 25px;",
              div(
                style = "clear: both; overflow-x: auto; width: 100%;",
                DTOutput("table_cfa_data_preview")
              )
            )
          )
        )
      )
    ),
    
    # 2.2 Statistik Deskriptif Butir
    nav_panel(
      title = "Statistik Deskriptif",
      value = "tab_cfa_desc",
      icon = icon("calculator"),
      card(
        full_screen = TRUE,
        card_header(
          class = "bg-light d-flex justify-content-between align-items-center",
          span(class = "fw-bold text-primary", icon("chart-bar"), " Statistik Deskriptif & Pengujian Normalitas Univariat"),
          downloadButton("btn_dl_cfa_desc_csv", "Unduh CSV Deskriptif", class = "btn-sm btn-outline-success")
        ),
        card_body(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Rata-rata Skor",
              value = textOutput("vb_cfa_desc_mean"),
              showcase = icon("calculator"),
              theme = "primary"
            ),
            value_box(
              title = "Rata-rata Skewness",
              value = textOutput("vb_cfa_desc_skew"),
              showcase = icon("balance-scale-left"),
              theme = "info"
            ),
            value_box(
              title = "Rata-rata Kurtosis",
              value = textOutput("vb_cfa_desc_kurt"),
              showcase = icon("chart-area"),
              theme = "warning"
            ),
            value_box(
              title = "Normal Univariat (%)",
              value = textOutput("vb_cfa_desc_norm_pct"),
              showcase = icon("shield-alt"),
              theme = "success"
            )
          ),
          hr(class = "my-3"),
          div(
            class = "alert alert-light border my-2 py-2 small",
            icon("info-circle", class = "text-primary me-2"),
            span("Kriteria Standar Emas Normalitas Univariat (Kline, 2016; West et al., 1995): Nilai absolut |Skewness| < 2.0 dan |Kurtosis| < 7.0 mengindikasikan distribusi butir masih memenuhi asumsi normalitas univariat untuk estimator berbasis Maximum Likelihood.")
          ),
          div(
            style = "clear: both; overflow-x: auto; width: 100%;",
            DTOutput("table_cfa_descriptive")
          )
        )
      ),
      card(
        class = "mt-3",
        card_header(class = "bg-light fw-bold text-primary", icon("eye"), " Eksplorasi Sebaran Empiris Respons per Butir"),
        card_body(
          layout_columns(
            col_widths = c(4, 8),
            div(
              h6(class = "fw-bold text-primary mb-2", "Pilih Butir Kuesioner:"),
              selectInput("cfa_desc_item_picker", label = NULL, choices = NULL, selectize = FALSE, width = "100%"),
              p(class = "small text-muted mt-2", "Grafik di samping menampilkan sebaran frekuensi empiris pilihan respons responden pada skala Likert.")
            ),
            plotOutput("cfa_desc_item_plot", height = "320px")
          )
        )
      )
    ),
    
    # 2.3 Matriks Korelasi (Pearson & Polikorik)
    nav_panel(
      title = "Matriks Korelasi",
      value = "tab_cfa_cor",
      icon = icon("th"),
      navset_card_tab(
        id = "cfa_cor_subtabs",
        nav_panel(
          title = "Korelasi Pearson (Linear Kontinu)",
          value = "subtab_cor_pearson",
          icon = icon("chart-line"),
          div(
            class = "d-flex justify-content-between align-items-center mb-3",
            div(
              h5(class = "fw-bold text-primary mb-1", icon("th"), " Matriks Korelasi Pearson Bivariat (r)"),
              p(class = "small text-muted mb-0", "Mengukur hubungan linear bivariat antar-butir kontinu. Digunakan sebagai matriks input estimasi Maximum Likelihood (ML) konvensional.")
            ),
            downloadButton("btn_dl_cfa_pearson_csv", "Unduh CSV Matriks Pearson", class = "btn-sm btn-outline-success")
          ),
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Determinan Matriks (|R|)",
              value = textOutput("vb_cfa_pearson_det"),
              showcase = icon("cube"),
              theme = "primary",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_pearson_det_sub"))
            ),
            value_box(
              title = "KMO Sampling Adequacy",
              value = textOutput("vb_cfa_pearson_kmo"),
              showcase = icon("check-circle"),
              theme = "info",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_pearson_kmo_sub"))
            ),
            value_box(
              title = "Bartlett's Sphericity (p)",
              value = textOutput("vb_cfa_pearson_bartlett"),
              showcase = icon("project-diagram"),
              theme = "success",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_pearson_bartlett_sub"))
            ),
            value_box(
              title = "Rata-rata Korelasi (r)",
              value = textOutput("vb_cfa_pearson_mean_r"),
              showcase = icon("calculator"),
              theme = "warning",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_pearson_range_r"))
            )
          ),
          hr(class = "my-3"),
          div(
            class = "alert alert-light border my-2 py-2 small",
            icon("info-circle", class = "text-primary me-2"),
            span("Kriteria Diagnostik Matriks R: Nilai Determinan |R| > 0.00001 mengindikasikan matriks definit positif dan bebas dari singularitas/multikolinearitas ekstrem. KMO MSA >= 0.70 (Kaiser, 1974) dan Bartlett's Test of Sphericity (p < 0.001) mengonfirmasi kecukupan interkorelasi butir untuk diestimasi dalam model faktor.")
          ),
          layout_columns(
            col_widths = c(7, 5),
            card(
              card_header(class = "bg-light fw-bold text-primary", icon("table"), " Tabel Matriks Korelasi Pearson (r)"),
              card_body(
                div(
                  style = "clear: both; overflow-x: auto; width: 100%;",
                  DTOutput("table_cfa_pearson_cor")
                )
              )
            ),
            card(
              card_header(
                class = "bg-light d-flex justify-content-between align-items-center",
                span(class = "fw-bold text-primary", icon("image"), " Visualisasi Heatmap Pearson"),
                div(
                  class = "d-flex gap-2 align-items-center",
                  checkboxInput("cfa_pearson_show_vals", "Tampilkan Angka", value = TRUE, width = "auto"),
                  selectInput("cfa_pearson_plot_type", NULL, choices = c("Segitiga Atas" = "upper", "Segitiga Bawah" = "lower", "Matriks Penuh" = "full"), selected = "upper", width = "140px")
                )
              ),
              card_body(
                plotOutput("cfa_pearson_heatmap", height = "520px")
              )
            )
          )
        ),
        nav_panel(
          title = "Korelasi Polikorik (IFA / DWLS / WLSMV)",
          value = "subtab_cor_poly",
          icon = icon("layer-group"),
          div(
            class = "d-flex justify-content-between align-items-center mb-3",
            div(
              h5(class = "fw-bold text-primary mb-1", icon("layer-group"), " Matriks Korelasi Polikorik (\u03c1)"),
              p(class = "small text-muted mb-0", "Mengestimasi korelasi laten kontinum yang mendasari respons skala ordinal/kategorikal terurut (Likert). Input wajib untuk Item Factor Analysis (IFA) dan estimator Robust DWLS / WLSMV.")
            ),
            downloadButton("btn_dl_cfa_poly_csv", "Unduh CSV Matriks Polikorik", class = "btn-sm btn-outline-success")
          ),
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Determinan Matriks (|\u03c1|)",
              value = textOutput("vb_cfa_poly_det"),
              showcase = icon("cube"),
              theme = "primary",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_poly_det_sub"))
            ),
            value_box(
              title = "KMO Polikorik MSA",
              value = textOutput("vb_cfa_poly_kmo"),
              showcase = icon("check-circle"),
              theme = "info",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_poly_kmo_sub"))
            ),
            value_box(
              title = "Bartlett Polikorik (p)",
              value = textOutput("vb_cfa_poly_bartlett"),
              showcase = icon("project-diagram"),
              theme = "success",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_poly_bartlett_sub"))
            ),
            value_box(
              title = "Rata-rata Korelasi (\u03c1)",
              value = textOutput("vb_cfa_poly_mean_rho"),
              showcase = icon("calculator"),
              theme = "warning",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_poly_range_rho"))
            )
          ),
          hr(class = "my-3"),
          div(
            class = "alert alert-light border my-2 py-2 small",
            icon("info-circle", class = "text-primary me-2"),
            span("Landasan Metodologis IFA / DWLS: Penggunaan korelasi Pearson pada data ordinal/Likert mengatenuasi besaran korelasi (Holgado-Tello et al., 2010; Muth\u00e9n, 1984). Korelasi polikorik memulihkan parameter asosiasi laten tanpa bias restriksi kategorikal, menjadikannya input fundamental bagi estimator DWLS / WLSMV di lavaan.")
          ),
          layout_columns(
            col_widths = c(7, 5),
            card(
              card_header(class = "bg-light fw-bold text-primary", icon("table"), " Tabel Matriks Korelasi Polikorik (\u03c1)"),
              card_body(
                div(
                  style = "clear: both; overflow-x: auto; width: 100%;",
                  DTOutput("table_cfa_poly_cor")
                )
              )
            ),
            card(
              card_header(
                class = "bg-light d-flex justify-content-between align-items-center",
                span(class = "fw-bold text-primary", icon("image"), " Visualisasi Heatmap Polikorik"),
                div(
                  class = "d-flex gap-2 align-items-center",
                  checkboxInput("cfa_poly_show_vals", "Tampilkan Angka", value = TRUE, width = "auto"),
                  selectInput("cfa_poly_plot_type", NULL, choices = c("Segitiga Atas" = "upper", "Segitiga Bawah" = "lower", "Matriks Penuh" = "full"), selected = "upper", width = "140px")
                )
              ),
              card_body(
                plotOutput("cfa_poly_heatmap", height = "520px")
              )
            )
          )
        ),
        
        # Sub-tab 3: Delta Atenuasi (Di sebelah Tab Polikorik)
        nav_panel(
          title = "Delta Atenuasi (\\Delta = \\rho - r)",
          value = "subtab_cor_delta",
          icon = icon("balance-scale"),
          div(
            class = "d-flex justify-content-between align-items-center mb-3",
            div(
              h5(class = "fw-bold text-danger mb-1", icon("balance-scale"), " Matriks Delta Atenuasi (\\Delta = \\rho - r)"),
              p(class = "small text-muted mb-0", "Mengkuantifikasi derajat pemulihan parameter asosiasi laten dan atenuasi koefisien Pearson akibat restriksi respons kategorikal (Likert).")
            ),
            downloadButton("btn_dl_cfa_delta_matrix_csv", "Unduh CSV Matriks Delta (\\Delta)", class = "btn-sm btn-outline-danger")
          ),
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Rata-rata Atenuasi (\\Delta)",
              value = textOutput("vb_cfa_cor_diff_mean"),
              showcase = icon("arrow-up"),
              theme = "danger",
              p(class = "small text-muted mb-0", "Rata-rata selisih \u03c1 - r")
            ),
            value_box(
              title = "Atenuasi Terbesar (\\Delta_max)",
              value = textOutput("vb_cfa_delta_max"),
              showcase = icon("exclamation-triangle"),
              theme = "warning",
              p(class = "small text-muted mb-0", textOutput("vb_cfa_delta_max_pair"))
            ),
            value_box(
              title = "Rata-rata Atenuasi (%)",
              value = textOutput("vb_cfa_delta_pct_mean"),
              showcase = icon("percentage"),
              theme = "info",
              p(class = "small text-muted mb-0", "Relatif terhadap \u03c1")
            ),
            value_box(
              title = "Pasangan Teratenuasi",
              value = textOutput("vb_cfa_delta_pair_count"),
              showcase = icon("check-double"),
              theme = "primary",
              p(class = "small text-muted mb-0", "Pasangan dengan \u03c1 > r")
            )
          ),
          hr(class = "my-3"),
          div(
            class = "alert alert-light border my-2 py-2 small",
            icon("info-circle", class = "text-danger me-2"),
            span("Kriteria Evaluasi Delta Atenuasi: Nilai \\Delta > 0 membuktikan bahwa korelasi Pearson mengalami underestimasi terhadap asosiasi laten yang sebenarnya. Semakin tinggi \\Delta, semakin rentan model CFA mengalami atenuasi muatan faktor (\\lambda) serta inflasi indeks kecocokan model semu jika memaksakan estimator Maximum Likelihood konvensional.")
          ),
          layout_columns(
            col_widths = c(7, 5),
            card(
              card_header(class = "bg-light fw-bold text-danger", icon("table"), " Tabel Matriks Delta Atenuasi (\\Delta = \\rho - r)"),
              card_body(
                div(
                  style = "clear: both; overflow-x: auto; width: 100%;",
                  DTOutput("table_cfa_delta_matrix")
                )
              )
            ),
            card(
              card_header(
                class = "bg-light d-flex justify-content-between align-items-center",
                span(class = "fw-bold text-danger", icon("image"), " Visualisasi Heatmap Delta (\\Delta)"),
                div(
                  class = "d-flex gap-2 align-items-center",
                  checkboxInput("cfa_delta_show_vals", "Tampilkan Angka", value = TRUE, width = "auto"),
                  selectInput("cfa_delta_plot_type", NULL, choices = c("Segitiga Atas" = "upper", "Segitiga Bawah" = "lower", "Matriks Penuh" = "full"), selected = "upper", width = "140px")
                )
              ),
              card_body(
                plotOutput("cfa_delta_heatmap", height = "520px")
              )
            )
          ),
          card(
            class = "mt-3",
            card_header(
              class = "bg-light d-flex justify-content-between align-items-center",
              span(class = "fw-bold text-primary", icon("list-ol"), " Pemeringkatan Pasangan Butir Berdasarkan Derajat Atenuasi Terbesar"),
              downloadButton("btn_dl_cfa_cor_diff_csv", "Unduh CSV Ranking Pasangan", class = "btn-sm btn-outline-secondary")
            ),
            card_body(
              p(class = "small text-muted mb-2", "Daftar pasangan butir diurutkan berdasarkan nilai absolut selisih (|\\Delta|) secara menurun untuk mengidentifikasi butir-butir yang mengalami distorsi kategorikal paling ekstrem."),
              div(
                style = "clear: both; overflow-x: auto; width: 100%;",
                DTOutput("table_cfa_cor_diff")
              )
            )
          )
        )
      )
    ),
    
    # 2.4 Eliminasi Butir Cepat (Fast Item Elimination)
    nav_panel(
      title = "Eliminasi Butir Cepat",
      value = "tab_cfa_elim",
      icon = icon("broom"),
      layout_sidebar(
        sidebar = sidebar(
          width = 360,
          title = "Parameter Eliminasi Cepat",
          sliderInput("cfa_elim_max_cor", "Batas Korelasi Item Maks (Redundansi):", min = 0.60, max = 0.95, value = 0.80, step = 0.05),
          p(class = "small text-muted mb-2", "Butir yang saling berkorelasi r > ambang batas (misal r > 0.80) ditandai sebagai redundan/kolinear."),
          sliderInput("cfa_elim_min_load", "Batas Muatan Faktor Minimum (λ):", min = 0.30, max = 0.65, value = 0.40, step = 0.05),
          p(class = "small text-muted mb-2", "Butir dengan muatan faktor terstandardisasi λ < batas akan dieliminasi."),
          hr(class = "my-2"),
          checkboxInput("cfa_elim_apply_cap", "Batasi Jumlah Butir per Dimensi", value = FALSE),
          conditionalPanel(
            condition = "input.cfa_elim_apply_cap == true",
            sliderInput("cfa_elim_cap_val", "Maksimal Butir per Dimensi:", min = 2, max = 8, value = 4, step = 1),
            p(class = "small text-muted", "Hanya mempertahankan top butir terbaik dengan loading tertinggi pada setiap dimensi.")
          ),
          actionButton("btn_run_elim", "Analisis & Seleksi Butir", class = "btn-primary w-100 mb-2", icon = icon("filter")),
          actionButton("btn_apply_elim", "Terapkan Butir Terpilih ke CFA", class = "btn-success w-100", icon = icon("check-double"))
        ),
        div(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Butir Awal (k_awal)",
              value = textOutput("vb_elim_k_orig"),
              showcase = icon("list"),
              theme = "secondary"
            ),
            value_box(
              title = "Dieliminasi",
              value = textOutput("vb_elim_k_dropped"),
              showcase = icon("trash-alt"),
              theme = "warning"
            ),
            value_box(
              title = "Dipertahankan (k_akhir)",
              value = textOutput("vb_elim_k_kept"),
              showcase = icon("check-circle"),
              theme = "primary"
            ),
            value_box(
              title = "Pasangan r > Batas",
              value = textOutput("vb_elim_high_cor_count"),
              showcase = icon("exclamation-triangle"),
              theme = "info"
            )
          ),
          hr(class = "my-3"),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("shield-alt"), " Perbandingan Reliabilitas Skala: Sebelum vs Setelah Seleksi Butir (Paket misty)"),
            card_body(
              tableOutput("table_cfa_elim_rel_comparison"),
              div(
                class = "alert alert-light border my-2 py-2 small",
                icon("info-circle", class = "text-primary me-2"),
                span("Reliabilitas diestimasi menggunakan paket ", tags$strong("misty"), " (Yanagida, 2026) dengan koefisien ", tags$strong("Cronbach's Alpha (α)"), " dan ", tags$strong("McDonald's Omega (ω)"), " beserta 95% Confidence Interval. McDonald's Omega mengevaluasi konsistensi internal bebas dari asumsi tau-ekuivalen esensial yang kaku.")
              )
            )
          ),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("chart-pie"), " Perbandingan Model Fit CFA: Sebelum vs Setelah Eliminasi Butir"),
            card_body(
              layout_columns(
                col_widths = c(3, 3, 3, 3),
                value_box(
                  title = "Perubahan ΔCFI",
                  value = textOutput("vb_elim_delta_cfi"),
                  showcase = icon("chart-line"),
                  theme = "success"
                ),
                value_box(
                  title = "Perubahan ΔTLI",
                  value = textOutput("vb_elim_delta_tli"),
                  showcase = icon("chart-bar"),
                  theme = "primary"
                ),
                value_box(
                  title = "Perubahan ΔRMSEA",
                  value = textOutput("vb_elim_delta_rmsea"),
                  showcase = icon("compress-arrows-alt"),
                  theme = "info"
                ),
                value_box(
                  title = "Perubahan ΔAIC",
                  value = textOutput("vb_elim_delta_aic"),
                  showcase = icon("balance-scale"),
                  theme = "secondary"
                )
              ),
              hr(class = "my-2"),
              tableOutput("table_cfa_elim_fit_comparison"),
              div(
                class = "alert alert-light border my-2 py-2 small",
                icon("info-circle", class = "text-primary me-2"),
                span("Evaluasi Model Fit Butir Tereliminasi: Pangkasan butir yang redundan (r > 0.80) atau memiliki muatan faktor rendah (λ < 0.40) mengoptimalkan parsimoni skala dan mengurangi deviasi residual matriks kovarians, yang secara langsung meningkatkan indeks kesesuaian komparatif (CFI/TLI) dan mereduksi kesalahan aproksimasi (RMSEA).")
              )
            )
          ),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("list-check"), " Status Seleksi Setiap Butir Kuesioner"),
            card_body(
              DTOutput("table_cfa_elim_item_details")
            )
          ),
          card(
            card_header(class = "bg-light fw-bold text-primary", icon("th"), " Matriks Korelasi Antar-Butir (Heatmap Redundansi)"),
            card_body(
              plotOutput("cfa_elim_heatmap_plot", height = "520px")
            )
          )
        )
      )
    ),
    
    # 2.4 Model Fit & Muatan Faktor
    nav_panel(
      title = "Model Fit & Loadings",
      value = "tab_cfa_fit",
      icon = icon("check-double"),
      card(
        card_header(class = "bg-light fw-bold text-primary", icon("balance-scale"), " Evaluasi Kelayakan Model (CFA Global Fit Indices)"),
        card_body(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "CFI (Comparative Fit)",
              value = textOutput("vb_cfa_cfi"),
              showcase = icon("chart-line"),
              theme = "primary"
            ),
            value_box(
              title = "TLI (Tucker-Lewis)",
              value = textOutput("vb_cfa_tli"),
              showcase = icon("chart-bar"),
              theme = "info"
            ),
            value_box(
              title = "RMSEA [90% CI]",
              value = textOutput("vb_cfa_rmsea"),
              showcase = icon("compress-arrows-alt"),
              theme = "warning"
            ),
            value_box(
              title = "SRMR (Residuals)",
              value = textOutput("vb_cfa_srmr"),
              showcase = icon("wave-square"),
              theme = "secondary"
            )
          ),
          hr(class = "my-3"),
          tableOutput("table_cfa_fit_indices"),
          div(
            class = "alert alert-light border my-2 py-2 small",
            icon("info-circle", class = "text-primary me-2"),
            span("Pedoman Hu & Bentler (1999): Model dinyatakan memiliki kesesuaian yang baik (Good Fit) apabila CFI >= 0.95, TLI >= 0.95, RMSEA <= 0.06, dan SRMR <= 0.08.")
          ),
          hr(class = "my-3"),
          h6(class = "fw-bold text-primary", icon("list-ol"), " Muatan Faktor Terstandardisasi (Standardized Factor Loadings λ):"),
          DTOutput("table_cfa_loadings"),
          hr(class = "my-3"),
          h6(class = "fw-bold text-primary", icon("project-diagram"), " Matriks Kovarians & Korelasi Antar-Faktor Laten:"),
          tableOutput("table_cfa_factor_cors"),
          hr(class = "my-3"),
          h6(class = "fw-bold text-primary", icon("shield-alt"), " Estimasi Reliabilitas Konstruk Skala (Paket misty: Cronbach's Alpha & McDonald's Omega):"),
          tableOutput("table_cfa_fit_reliability"),
          div(
            class = "alert alert-light border my-2 py-2 small",
            icon("info-circle", class = "text-primary me-2"),
            span("Estimasi reliabilitas skala dan masing-masing dimensi laten menggunakan fungsi ", tags$code("misty::item.alpha()"), " dan ", tags$code("misty::item.omega()"), " (Yanagida, 2026). Nilai koefisien >= 0.70 merefleksikan konsistensi internal yang memadai untuk instrumen psikologis.")
          )
        )
      )
    ),
    
    # 2.5 Diagram Jalur (Path Diagram CFA)
    nav_panel(
      title = "Diagram Jalur (Path)",
      value = "tab_cfa_path",
      icon = icon("bezier-curve"),
      layout_sidebar(
        sidebar = sidebar(
          width = 340,
          title = "Konfigurasi Diagram Jalur",
          selectInput(
            "cfa_path_model_struct",
            "Struktur Model Diagram:",
            choices = c(
              "Multidimensional (Faktor Berkorelasi)" = "multi",
              "Unidimensional (1 Faktor)" = "unidim",
              "Second-Order (Hierarkis Orde Kedua)" = "second_order",
              "Bi-Factor (Umum & Spesifik Ortogonal)" = "bifactor"
            ),
            selected = "multi"
          ),
          selectInput(
            "cfa_path_style",
            "Gaya Visualisasi Diagram:",
            choices = c(
              "Academic Kawe (Eksklusif Maroon & Gold)" = "kawe",
              "Path Diagram Biasa (Standar Klasik)" = "biasa"
            ),
            selected = "kawe"
          ),
          selectInput(
            "cfa_path_layout",
            "Pola Tata Letak (Layout):",
            choices = c(
              "Tree Horizontal (Kiri ke Kanan - Anti Tabrakan)" = "tree_h",
              "Tree Vertikal (Atas ke Bawah)" = "tree_v",
              "Circle (Melingkar)" = "circle",
              "Spring (Organik / Gaya Bebas)" = "spring"
            ),
            selected = "tree_h"
          ),
          sliderInput("cfa_path_cex", "Ukuran Simbol & Label:", min = 0.5, max = 1.3, value = 0.85, step = 0.05),
          sliderInput("cfa_path_edge_width", "Ketebalan Garis Panah (APA 7th Standar):", min = 0.8, max = 2.0, value = 1.2, step = 0.1),
          sliderInput("cfa_path_height", "Tinggi Kanvas Diagram (px):", min = 600, max = 1300, value = 820, step = 40),
          checkboxInput("cfa_path_show_std", "Tampilkan Nilai Loading Terstandardisasi", value = TRUE),
          checkboxInput("cfa_path_show_resid", "Tampilkan Varians Residual", value = FALSE),
          hr(class = "my-2"),
          downloadButton("btn_dl_cfa_path_png", "Unduh Gambar (PNG 300 DPI)", class = "btn-primary w-100 mb-2"),
          downloadButton("btn_dl_cfa_path_pdf", "Unduh Format Vektor (PDF)", class = "btn-secondary w-100")
        ),
        card(
          card_header(
            class = "bg-light d-flex justify-content-between align-items-center flex-wrap gap-2",
            div(class = "fw-bold text-primary", icon("bezier-curve"), " Visualisasi Diagram Jalur Analisis Faktor Konfirmatori"),
            uiOutput("cfa_path_badges")
          ),
          card_body(
            class = "p-2 d-flex justify-content-center align-items-center",
            style = "background-color: #FFFFFF; min-height: 760px; overflow-x: auto;",
            uiOutput("cfa_path_plot_ui")
          )
        )
      )
    ),
    
    # 2.6 Perbandingan Model Fit (CFA Model Comparison)
    nav_panel(
      title = "Perbandingan Model Fit",
      value = "tab_cfa_comp",
      icon = icon("balance-scale"),
      layout_sidebar(
        sidebar = sidebar(
          width = 340,
          title = "Komparasi Struktur Model CFA",
          p(class = "small text-muted mb-3", "Mengevaluasi kesesuaian empiris antara model bersaing: Unidimensi, Multidimensi Orde-1, Orde-2 (Second-Order), dan Bi-Factor."),
          actionButton("btn_run_cfa_comp", "Bandingkan 4 Model CFA", class = "btn-primary w-100 mb-3", icon = icon("calculator")),
          div(
            class = "p-2 bg-light border rounded small",
            strong("Struktur Model yang Dibandingkan:"),
            tags$ul(
              class = "ps-3 mb-0 mt-1",
              tags$li(strong("1. Unidimensi:"), " 1 faktor umum menjelaskan seluruh butir."),
              tags$li(strong("2. Multidimensi:"), " Faktor-faktor laten orde-1 saling berkorelasi."),
              tags$li(strong("3. Second-Order:"), " Faktor umum hierarkis menjelaskan faktor-faktor orde-1."),
              tags$li(strong("4. Bi-Factor:"), " Faktor umum ortogonal + faktor spesifik ortogonal.")
            )
          ),
          hr(class = "my-2"),
          downloadButton("btn_dl_cfa_comp_csv", "Unduh CSV Tabel Komparasi", class = "btn-outline-success w-100")
        ),
        div(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Model Terbaik (AIC)",
              value = textOutput("vb_cfa_comp_best_aic"),
              showcase = icon("trophy"),
              theme = "primary"
            ),
            value_box(
              title = "Model Terbaik (BIC)",
              value = textOutput("vb_cfa_comp_best_bic"),
              showcase = icon("award"),
              theme = "info"
            ),
            value_box(
              title = "Model Terbaik (CFI/RMSEA)",
              value = textOutput("vb_cfa_comp_best_cfi"),
              showcase = icon("check-circle"),
              theme = "success"
            ),
            value_box(
              title = "Total Sampel & Butir",
              value = textOutput("vb_cfa_comp_n_k"),
              showcase = icon("users"),
              theme = "secondary"
            )
          ),
          hr(class = "my-3"),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("table"), " Tabel Evaluasi & Komparasi Kesesuaian Model CFA (Global Fit Indices)"),
            card_body(
              DTOutput("table_cfa_comp_summary")
            )
          ),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("vial"), " Uji Rasio Kemiripan Model Bersarang (Likelihood Ratio Test / Δχ²)"),
            card_body(
              tableOutput("table_cfa_comp_lrt")
            )
          ),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("chart-bar"), " Visualisasi Diagnostik Perbandingan Indeks Fit (CFI, TLI, RMSEA, SRMR)"),
            card_body(
              plotOutput("cfa_comp_bar_plot", height = "420px")
            )
          ),
          card(
            card_header(class = "bg-light fw-bold text-primary", icon("graduation-cap"), " Catatan Metodologis & Evaluasi Parsimoni Model"),
            card_body(
              div(
                class = "alert alert-light border mb-0 small",
                style = "line-height: 1.7;",
                p(strong("Pedoman Penilaian Model CFA (Hu & Bentler, 1999; Kline, 2016):"), 
                  " Pemilihan model pengukuran terbaik tidak semata-mata bertumpu pada minimisasi nilai Chi-Square (χ²), melainkan menyeimbangkan kecukupan aproksimasi kovarians empiris (CFI/TLI >= 0.95, RMSEA <= 0.06, SRMR <= 0.08) dengan prinsip parsimoni model yang direfleksikan oleh indeks informasi AIC dan BIC."),
                p(strong("Dinamika Model Bi-Factor vs Multidimensi:"),
                  " Pada konstruksi skala psikologi, model bi-factor sering kali menghasilkan peningkatan nilai fit mekanis akibat relaksasi ortogonalitas. Namun demikian, peneliti perlu memeriksa apakah muatan faktor spesifik tetap substansial (λ >= 0.30) atau mengalami atenuasi tajam setelah varians umum diekstraksi oleh faktor general.")
              )
            )
          )
        )
      )
    )
  ),

  # ============================================================================
  # DROPDOWN MENU: RASCH ANALYSIS (TAM ENGINE)
  # ============================================================================
  nav_menu(
    title = "Rasch",
    icon = icon("ruler-combined"),
    
    # 3.1 Model & Simulasi Rasch TAM
    nav_panel(
      title = "Model & Simulasi TAM",
      value = "tab_rasch_model",
      icon = icon("sliders-h"),
      layout_sidebar(
        fillable = FALSE,
        sidebar = sidebar(
          width = 360,
          title = "Konfigurasi Model Rasch (TAM)",
          h6(class = "text-primary fw-bold mb-2", icon("database"), " Sumber Data Respons Butir:"),
          radioButtons(
            "rasch_source_type",
            NULL,
            choices = c(
              "Simulasi Rasch (TAM)" = "sim",
              "Dataset Bawaan Politomi (4 Dimensi, k=20)" = "builtin_qq",
              "Dataset Bawaan Dikotomi (NeoUni k=8)" = "builtin",
              "Unggah Berkas (CSV / Excel)" = "upload"
            ),
            selected = "sim"
          ),
          conditionalPanel(
            condition = "input.rasch_source_type == 'sim'",
            div(
              class = "p-2 bg-white rounded border mb-2",
              sliderInput("rasch_sim_n", "Jumlah Responden (Person N):", min = 100, max = 1500, value = 300, step = 50),
              sliderInput("rasch_sim_k", "Total Butir (k):", min = 5, max = 30, value = 10, step = 1),
              radioButtons(
                "rasch_sim_format",
                "Format Respons:",
                choices = c("Dikotomi (0/1)" = "dichotomous", "Politomi Likert (Rating Scale/PCM)" = "polytomous"),
                selected = "dichotomous"
              ),
              radioButtons(
                "rasch_sim_dim_count",
                "Dimensionalitas Laten:",
                choices = c("Unidimensi (1 Dimensi)" = "uni", "Multidimensi (2 Dimensi)" = "multi_2", "Multidimensi (3 Dimensi)" = "multi_3"),
                selected = "uni"
              ),
              conditionalPanel(
                condition = "input.rasch_sim_dim_count != 'uni'",
                sliderInput("rasch_sim_cor", "Korelasi Antardimensi:", min = 0.1, max = 0.8, value = 0.40, step = 0.05)
              ),
              conditionalPanel(
                condition = "input.rasch_sim_format == 'polytomous'",
                sliderInput("rasch_sim_n_cats", "Kategori Respons (0 s/d m):", min = 3, max = 5, value = 4, step = 1)
              ),
              actionButton("btn_rasch_gen_sim", "Bangkitkan Data Simulasi Rasch", class = "btn-primary w-100", icon = icon("dice"))
            )
          ),
          conditionalPanel(
            condition = "input.rasch_source_type == 'upload'",
            div(
              class = "p-2 bg-white rounded border mb-2",
              fileInput("rasch_upload_file", "Pilih Berkas Respons (.csv / .xlsx):", accept = c(".csv", ".xlsx", ".xls"))
            )
          ),
          hr(class = "my-2"),
          h6(class = "text-primary fw-bold mb-2", icon("cogs"), " Spesifikasi Pemodelan TAM:"),
          radioButtons(
            "rasch_model_dim_spec",
            "Spesifikasi Dimensi TAM:",
            choices = c("Unidimensi" = "uni", "Multidimensi (Q-Matrix)" = "multi"),
            selected = "uni"
          ),
          conditionalPanel(
            condition = "input.rasch_model_dim_spec == 'multi'",
            div(
              class = "p-2 bg-light rounded border mb-2",
              sliderInput("rasch_multi_q_dims", "Jumlah Dimensi (D):", min = 2, max = 5, value = 4, step = 1),
              radioButtons(
                "rasch_qmatrix_mode",
                "Metode Desain Q-Matrix:",
                choices = c(
                  "Pemetaan Dimensi Interaktif" = "interactive",
                  "Sintaks R Studio / Biner" = "syntax",
                  "Otomatis (Seimbang)" = "auto"
                ),
                selected = "interactive"
              ),
              conditionalPanel(
                condition = "input.rasch_qmatrix_mode == 'interactive'",
                actionButton("btn_load_std_qmatrix", "Muat Desain Q-Matrix Standar (4 Dimensi)", class = "btn-outline-primary btn-sm w-100 mb-2", icon = icon("layer-group")),
                uiOutput("rasch_qmatrix_interactive_ui")
              ),
              conditionalPanel(
                condition = "input.rasch_qmatrix_mode == 'syntax'",
                textAreaInput(
                  "rasch_qmatrix_syntax_text",
                  "Sintaks Matriks Q (R Studio):",
                  value = "Q <- matrix(c(\n  1, 0, 0,\n  0, 1, 0,\n  0, 0, 1\n), nrow = 3, ncol = 3, byrow = TRUE)",
                  rows = 6
                ),
                p(class = "small text-muted mb-0", "Definisikan matriks biner (0/1) k baris x D kolom atau sintaksis R matrix(...).")
              )
            )
          ),
          selectInput(
            "rasch_irt_model_type",
            "Tipe Model & Estimator IRT (TAM):",
            choices = list(
              "Rasch Family (Fixed Slopes / MML)" = c(
                "1PL - One-Parameter Logistic / Rasch Model (Dikotomi)" = "1PL",
                "PCM - Partial Credit Model (Masters, 1982)" = "PCM",
                "PCM2 - ConQuest Parameterization (item + item:step)" = "PCM2",
                "RSM - Rating Scale Model (Andrich, 1978)" = "RSM"
              ),
              "Slope Estimation / 2PL Family (tam.mml.2pl)" = c(
                "2PL - Two-Parameter Logistic Model (Birnbaum, 1968)" = "2PL",
                "GPCM - Generalized Partial Credit Model (Muraki, 1992)" = "GPCM",
                "2PL.groups - Constrained 2PL Slopes" = "2PL.groups",
                "GPCM.groups - Constrained GPCM Slopes" = "GPCM.groups"
              ),
              "Guessing / 3PL Family (tam.mml.3pl)" = c(
                "3PL - Three-Parameter Logistic Model (Guessing)" = "3PL"
              ),
              "Joint Maximum Likelihood (JML)" = c(
                "JML - Joint Maximum Likelihood (tam.jml)" = "JML"
              )
            ),
            selected = "1PL"
          ),
          hr(class = "my-2"),
          h6(class = "text-primary fw-bold mb-2", icon("microchip"), " Akselerasi Komputasi Multi-Core:"),
          checkboxInput(
            "rasch_use_parallel",
            tags$strong("Aktifkan Paralelisasi (Multi-Core CPU)"),
            value = TRUE
          ),
          conditionalPanel(
            condition = "input.rasch_use_parallel == true",
            div(
              class = "p-2 bg-light rounded border mb-2",
              sliderInput(
                "rasch_n_cores",
                "Alokasi Core CPU (doParallel):",
                min = 1,
                max = max(1, n_sys_cores),
                value = max(1, min(4, n_sys_cores)),
                step = 1
              ),
              div(
                class = "small text-muted d-flex align-items-center justify-content-between",
                span(icon("server"), sprintf(" Sistem: %d Cores Terdeteksi", n_sys_cores)),
                span(class = "badge bg-success", "doParallel Siap")
              )
            )
          ),
          actionButton("btn_run_rasch", "Estimasi Model Rasch (TAM)", class = "btn-success w-100 mt-2", icon = icon("play-circle"))
        ),
        div(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Deviance (-2LL)",
              value = textOutput("vb_rasch_dev"),
              showcase = icon("chart-pie"),
              theme = "primary"
            ),
            value_box(
              title = "AIC / BIC",
              value = textOutput("vb_rasch_aic_bic"),
              showcase = icon("balance-scale"),
              theme = "secondary"
            ),
            value_box(
              title = "EAP / WLE Reliability",
              value = textOutput("vb_rasch_rel"),
              showcase = icon("shield-alt"),
              theme = "info"
            ),
            value_box(
              title = "Dimensi Laten",
              value = textOutput("vb_rasch_dims"),
              showcase = icon("project-diagram"),
              theme = "warning"
            )
          ),
          hr(class = "my-3"),
          conditionalPanel(
            condition = "input.rasch_model_dim_spec == 'multi'",
            card(
              class = "mb-3",
              card_header(
                class = "bg-light d-flex justify-content-between align-items-center",
                span(class = "fw-bold text-primary", icon("th"), " Pratinjau Matriks Desain Q (Q-Matrix Inspection)"),
                span(class = "badge bg-secondary", textOutput("rasch_qmatrix_badge_dims", inline = TRUE))
              ),
              card_body(
                p(class = "small text-muted mb-2", "Matriks Q memetakan butir manifest (baris) ke dimensi laten (kolom). Nilai 1 menandakan butir mengukur dimensi tersebut."),
                div(style = "overflow-x: auto; max-height: 380px;", tableOutput("rasch_qmatrix_preview_table"))
              )
            )
          ),
          card(
            fill = FALSE,
            class = "mb-3 shadow-sm border",
            card_header(class = "bg-light fw-bold text-primary", icon("project-diagram"), " Matriks Varians-Kovarians & Korelasi Dimensi Laten (TAM)"),
            card_body(
              h6(class = "fw-bold text-primary small", icon("th-large"), " Matriks Varians-Kovarians Dimensi Laten:"),
              tableOutput("table_rasch_latent_cov"),
              hr(class = "my-2"),
              h6(class = "fw-bold text-primary small", icon("link"), " Matriks Korelasi Antar-Dimensi Laten:"),
              tableOutput("table_rasch_latent_cor"),
              hr(class = "my-2"),
              h6(class = "fw-bold text-primary small", icon("shield-alt"), " Reliabilitas Pengukuran per Dimensi (WLE & EAP Reliability):"),
              tableOutput("table_rasch_reliability_dims")
            )
          ),
          card(
            fill = FALSE,
            class = "mb-3 shadow-sm border",
            card_header(class = "bg-light fw-bold text-primary", icon("terminal"), " Output Konsol Estimasi Paket TAM (Robitzsch, 2025)"),
            card_body(
              verbatimTextOutput("rasch_tam_console_summary")
            )
          ),
          card(
            fill = FALSE,
            class = "mb-3 shadow-sm border",
            card_header(class = "bg-light d-flex justify-content-between align-items-center",
              span(class = "fw-bold text-primary", icon("table"), " Matriks Data Respons Responden (TAM Input)"),
              uiOutput("rasch_data_badge_ui")
            ),
            card_body(
              style = "overflow: visible; clear: both;",
              div(
                style = "overflow-x: auto; width: 100%; clear: both;",
                DTOutput("table_rasch_resp_preview")
              )
            )
          )
        )
      )
    ),
    
    # 3.2 Item Fit (Infit & Outfit Butir)
    nav_panel(
      title = "Item Fit (Infit/Outfit)",
      value = "tab_rasch_itemfit",
      icon = icon("check-double"),
      div(
        class = "d-flex flex-column gap-3",
        # Kartu 1: Value Boxes & Tabel Item Fit
        card(
          class = "mb-2 shadow-sm border",
          card_header(
            class = "bg-light d-flex justify-content-between align-items-center flex-wrap gap-2",
            span(class = "fw-bold text-primary", icon("check-double"), " Estimasi Parameter Kesulitan & Kesesuaian Butir (Item Infit & Outfit)"),
            div(
              class = "d-flex gap-2 align-items-center",
              actionButton("btn_goto_rasch_reanalysis_from_item", "Analisis Ulang Butir Misfit", class = "btn-sm btn-primary", icon = icon("filter")),
              downloadButton("btn_dl_rasch_itemfit_csv", "Unduh CSV Item Fit", class = "btn-sm btn-outline-success")
            )
          ),
          card_body(
            layout_columns(
              col_widths = c(3, 3, 3, 3),
              value_box(
                title = "Butir Produktif / Fit",
                value = textOutput("vb_rasch_item_fit_count"),
                showcase = icon("check-circle"),
                theme = "success"
              ),
              value_box(
                title = "Underfit / Noise (MNSQ > 1.5)",
                value = textOutput("vb_rasch_item_underfit_count"),
                showcase = icon("exclamation-circle"),
                theme = "warning"
              ),
              value_box(
                title = "Overfit / Redundan (MNSQ < 0.5)",
                value = textOutput("vb_rasch_item_overfit_count"),
                showcase = icon("compress"),
                theme = "info"
              ),
              value_box(
                title = "Rata-rata Item Difficulty (\\delta)",
                value = textOutput("vb_rasch_item_mean_beta"),
                showcase = icon("ruler-horizontal"),
                theme = "secondary"
              )
            ),
            hr(class = "my-3"),
            div(
              class = "alert alert-light border my-2 py-2 small",
              icon("info-circle", class = "text-primary me-2"),
              span("Kriteria Linacre (2002): Rentang nilai Mean Square (MNSQ) antara 0.50 hingga 1.50 dikategorikan produktif untuk pengukuran. Nilai MNSQ > 1.50 mengindikasikan adanya noise / variasi acak (underfit), sedangkan MNSQ < 0.50 mengindikasikan pola respons yang terlalu terduga / redundan (overfit).")
            ),
            div(
              style = "width: 100%; overflow-x: auto; clear: both; margin-top: 15px; margin-bottom: 15px;",
              DTOutput("table_rasch_itemfit_full")
            )
          )
        ),
        
        # Kartu 2: Scatterplot Diagnostik
        card(
          class = "mb-2 shadow-sm border",
          card_header(class = "bg-light fw-bold text-primary", icon("chart-line"), " Scatterplot Diagnostik Item Fit (Infit MNSQ vs Outfit MNSQ)"),
          card_body(
            plotOutput("rasch_itemfit_scatterplot", height = "460px")
          )
        ),
        
        # Kartu 3: Matriks Thurstonian Thresholds
        card(
          class = "mb-2 shadow-sm border",
          card_header(class = "bg-light fw-bold text-primary", icon("layer-group"), " Matriks Nilai Ambang Batas Kategori (Thurstonian Thresholds Matrix)"),
          card_body(
            div(style = "width: 100%; overflow-x: auto;", tableOutput("table_rasch_thurstonian_thresholds")),
            div(
              class = "alert alert-light border my-2 py-2 small",
              icon("info-circle", class = "text-primary me-2"),
              span("Thurstonian Thresholds mengindikasikan titik kontinum laten (θ) di mana probabilitas kumulatif untuk merespons kategori k atau lebih tinggi mencapai 50%. Nilai ambang batas yang monoton meningkat mencerminkan keterurutan kategori respons (ordered response categories).")
            )
          )
        )
      )
    ),
    
    # 3.3 Person Fit (Infit & Outfit Responden)
    nav_panel(
      title = "Person Fit (Infit/Outfit)",
      value = "tab_rasch_personfit",
      icon = icon("user-check"),
      div(
        class = "d-flex flex-column gap-3",
        # Kartu 1: Value Boxes & Dimension-wise Person Fit
        card(
          class = "mb-2 shadow-sm border",
          card_header(
            class = "bg-light d-flex justify-content-between align-items-center flex-wrap gap-2",
            span(class = "fw-bold text-primary", icon("users"), " Parameter Kemampuan & Kesesuaian Responden (Person Measures & Fit)"),
            div(
              class = "d-flex gap-2 align-items-center",
              actionButton("btn_goto_rasch_reanalysis_from_person", "Analisis Ulang Responden Aberan", class = "btn-sm btn-primary", icon = icon("user-minus")),
              downloadButton("btn_dl_rasch_personfit_csv", "Unduh CSV Person Fit", class = "btn-sm btn-outline-success")
            )
          ),
          card_body(
            layout_columns(
              col_widths = c(3, 3, 3, 3),
              value_box(
                title = "Rata-rata Kemampuan (θ)",
                value = textOutput("vb_rasch_person_mean_theta"),
                showcase = icon("user"),
                theme = "primary"
              ),
              value_box(
                title = "Standar Deviasi θ",
                value = textOutput("vb_rasch_person_sd_theta"),
                showcase = icon("arrows-alt-h"),
                theme = "secondary"
              ),
              value_box(
                title = "Responden Fit Normal (%)",
                value = textOutput("vb_rasch_person_fit_pct"),
                showcase = icon("user-shield"),
                theme = "success"
              ),
              value_box(
                title = "Pola Aberan / Tebakan (%)",
                value = textOutput("vb_rasch_person_misfit_pct"),
                showcase = icon("user-times"),
                theme = "warning"
              )
            ),
            hr(class = "my-3"),
            div(
              class = "alert alert-light border my-2 py-2 small",
              icon("info-circle", class = "text-primary me-2"),
              span("Kesesuaian Responden (Person-Fit): Responden dengan Person Infit/Outfit MNSQ > 1.50 menunjukkan kemungkinan perilaku menjawab secara aberan, tebakan acak, atau kelelahan kognitif saat mengerjakan instrumen.")
            ),
            h6(class = "fw-bold text-primary small mt-3 mb-2", icon("list-check"), " Evaluasi Kesesuaian Responden per Dimensi (Dimension-wise Person Fit Summary):"),
            div(style = "width: 100%; overflow-x: auto;", tableOutput("table_rasch_person_fit_dim_summary"))
          )
        ),
        
        # Kartu 2: Tabel Lengkap Responden
        card(
          class = "mb-2 shadow-sm border",
          card_header(class = "bg-light fw-bold text-primary", icon("table"), " Tabel Lengkap Kemampuan & Kesesuaian Responden (Person Fit Full Data)"),
          card_body(
            div(
              style = "width: 100%; overflow-x: auto; clear: both; margin-top: 15px; margin-bottom: 15px;",
              DTOutput("table_rasch_personfit_full")
            )
          )
        ),
        
        # Kartu 3: Distribusi Theta
        card(
          class = "mb-2 shadow-sm border",
          card_header(class = "bg-light fw-bold text-primary", icon("chart-area"), " Distribusi Kemampuan Responden (Latent Trait θ Distribution)"),
          card_body(
            plotOutput("rasch_person_theta_dist_plot", height = "420px")
          )
        )
      )
    ),
    
    # 3.4 Analisis Ulang & Pemurnian (Re-Analysis)
    nav_panel(
      title = "Analisis Ulang (Re-Analysis)",
      value = "tab_rasch_reanalysis",
      icon = icon("recycle"),
      layout_sidebar(
        sidebar = sidebar(
          width = 380,
          title = "Konfigurasi Re-Analisis & Pemurnian",
          div(
            class = "alert alert-light border py-2 px-3 small mb-3",
            style = "line-height: 1.5;",
            icon("info-circle", class = "text-primary me-1"),
            span("Lakukan kalibrasi ulang (*re-estimation*) setelah memurnikan instrumen dari butir misfit (underfit/overfit) dan/atau responden aberan (*satisficing* / tebakan acak).")
          ),
          
          # Card A: Filter Butir Misfit
          div(
            class = "p-2 bg-white rounded border mb-2",
            checkboxInput("rasch_re_filter_items", strong("1. Eliminasi Butir Misfit"), value = TRUE),
            conditionalPanel(
              condition = "input.rasch_re_filter_items == true",
              radioButtons(
                "rasch_re_item_criteria",
                "Kriteria Deteksi Misfit Butir:",
                choices = c(
                  "Infit MNSQ di luar batas" = "infit",
                  "Outfit MNSQ di luar batas" = "outfit",
                  "Infit ATAU Outfit di luar batas" = "either"
                ),
                selected = "either"
              ),
              sliderInput("rasch_re_item_upper_mnsq", "Batas Atas MNSQ (Underfit / Noise):", min = 1.20, max = 2.00, value = 1.50, step = 0.05),
              sliderInput("rasch_re_item_lower_mnsq", "Batas Bawah MNSQ (Overfit / Redundan):", min = 0.30, max = 0.80, value = 0.50, step = 0.05),
              selectizeInput(
                "rasch_re_manual_drop_items",
                "Pilihan Manual Butir Tambahan:",
                choices = NULL,
                multiple = TRUE,
                options = list(plugins = list('remove_button'), placeholder = "Pilih butir untuk dikeluarkan...")
              )
            )
          ),
          
          # Card B: Filter Responden Aberan
          div(
            class = "p-2 bg-white rounded border mb-3",
            checkboxInput("rasch_re_filter_persons", strong("2. Eliminasi Responden Aberan"), value = TRUE),
            conditionalPanel(
              condition = "input.rasch_re_filter_persons == true",
              radioButtons(
                "rasch_re_person_criteria",
                "Kriteria Person Misfit:",
                choices = c(
                  "Person Infit MNSQ > Batas" = "infit",
                  "Person Outfit MNSQ > Batas" = "outfit",
                  "Person Infit ATAU Outfit > Batas" = "either"
                ),
                selected = "either"
              ),
              sliderInput("rasch_re_person_mnsq_thresh", "Batas Maksimal Person MNSQ:", min = 1.20, max = 2.50, value = 1.50, step = 0.05),
              p(class = "small text-muted mb-0", "Responden dengan Person MNSQ > ambang batas mengindikasikan deviasi pola respons (tebakan acak, kepatuhan parsial, atau careless responding).")
            )
          ),
          
          # Action Buttons
          actionButton("btn_run_rasch_reanalysis", "Jalankan Analisis Ulang Rasch", class = "btn-primary w-100 mb-2 py-2 fw-bold shadow-sm", icon = icon("calculator")),
          actionButton("btn_apply_rasch_reanalysis", "Terapkan Model Pemurnian ke Seluruh Tab", class = "btn-success w-100 mb-2 py-2 fw-bold shadow-sm", icon = icon("check-double")),
          actionButton("btn_reset_rasch_reanalysis", "Reset ke Data Lengkap Awal", class = "btn-outline-danger w-100 mb-2 py-1", icon = icon("undo")),
          hr(class = "my-2"),
          downloadButton("btn_dl_rasch_reanalysis_comp_csv", "Unduh CSV Tabel Komparasi", class = "btn-outline-secondary w-100 btn-sm")
        ),
        
        # Main Display
        div(
          # Row 1: Value Boxes
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Butir Terpilih / Aktif",
              value = textOutput("vb_rasch_re_items"),
              showcase = icon("list-check"),
              theme = "primary"
            ),
            value_box(
              title = "Responden Lolos Filter",
              value = textOutput("vb_rasch_re_persons"),
              showcase = icon("users"),
              theme = "info"
            ),
            value_box(
              title = "Reliabilitas EAP Pasca-Pemurnian",
              value = textOutput("vb_rasch_re_eap"),
              showcase = icon("shield-alt"),
              theme = "success"
            ),
            value_box(
              title = "Proporsi Responden Fit",
              value = textOutput("vb_rasch_re_fit_pct"),
              showcase = icon("user-shield"),
              theme = "secondary"
            )
          ),
          hr(class = "my-3"),
          
          # Row 2: Tabel Komparasi Model Awal vs Setelah Pemurnian
          card(
            class = "mb-3 shadow-sm border",
            card_header(
              class = "bg-light d-flex justify-content-between align-items-center",
              span(class = "fw-bold text-primary", icon("balance-scale"), " Komparasi Metrik Psikometris: Model Awal vs Model Setelah Pemurnian"),
              span(class = "badge bg-primary text-white", "Evaluasi Kalibrasi Ulang TAM")
            ),
            card_body(
              div(style = "width: 100%; overflow-x: auto;", tableOutput("table_rasch_reanalysis_comp")),
              div(
                class = "alert alert-light border my-2 py-2 small",
                icon("info-circle", class = "text-primary me-2"),
                span("Sesuai pedoman psikometrika Linacre (2002) dan Bond & Fox (2015), pemurnian instrumen melalui pengeluaran butir underfit/overfit dan responden aberan memitigasi distorsi estimasi kesulitan butir (item difficulty) serta mengeliminasi inflasi varians error.")
              )
            )
          ),
          
          # Row 3: Rincian Butir & Responden yang Dieliminasi & Pergeseran Butir
          card(
            class = "mb-3 shadow-sm border",
            card_header(class = "bg-light fw-bold text-primary", icon("exchange-alt"), " Rincian Elemen yang Dieliminasi & Pergeseran Parameter Butir (Item Difficulty & Fit Shifts)"),
            card_body(
              navset_card_tab(
                nav_panel(
                  title = "Pergeseran Parameter Butir (Δ Kesulitan δ & Fit)",
                  icon = icon("table"),
                  div(
                    class = "alert alert-light border my-2 py-2 small",
                    icon("info-circle", class = "text-primary me-2"),
                    span("Tabel ini menyajikan perubahan tingkat kesulitan butir (Item Difficulty δ) serta pergeseran nilai Infit dan Outfit MNSQ sebelum vs sesudah pemurnian. Penghapusan responden aberan memitigasi distorsi pada estimasi tingkat kesulitan dan mengembalikan nilai MNSQ menuju nilai ideal 1.00.")
                  ),
                  div(style = "width: 100%; overflow-x: auto; margin-top: 10px;", DTOutput("table_rasch_re_items_comp"))
                ),
                nav_panel(
                  title = "Butir Misfit Tereliminasi",
                  icon = icon("list"),
                  div(style = "width: 100%; overflow-x: auto; margin-top: 10px;", DTOutput("table_rasch_re_dropped_items"))
                ),
                nav_panel(
                  title = "Responden Aberan Tereliminasi",
                  icon = icon("user-times"),
                  div(style = "width: 100%; overflow-x: auto; margin-top: 10px;", DTOutput("table_rasch_re_dropped_persons"))
                )
              )
            )
          ),
          
          # Row 4: Wright Map Model Pemurnian
          card(
            class = "mb-3 shadow-sm border",
            card_header(
              class = "bg-light d-flex justify-content-between align-items-center",
              span(class = "fw-bold text-primary", icon("map"), " Wright Map Model Pasca-Pemurnian (Peta Responden-Butir Bersih)"),
              downloadButton("btn_dl_rasch_re_wright_png", "Unduh Wright Map Bersih (PNG)", class = "btn-sm btn-outline-primary")
            ),
            card_body(
              plotOutput("plot_rasch_reanalysis_wright_map", height = "560px")
            )
          )
        )
      )
    ),
    
    # 3.5 Wright Map (Uni & Multi)
    nav_panel(
      title = "Wright Map (Uni & Multi)",
      value = "tab_rasch_wright",
      icon = icon("map"),
      card(
        fill = FALSE,
        card_header(
          class = "bg-light d-flex justify-content-between align-items-center flex-wrap gap-2",
          span(class = "fw-bold text-primary", icon("map"), " Peta Responden-Butir (Person-Item Wright Map)"),
          div(
            class = "d-flex gap-2 align-items-center",
            downloadButton("btn_dl_wright_png", "Unduh Gambar (PNG 300 DPI)", class = "btn-sm btn-primary"),
            downloadButton("btn_dl_wright_pdf", "Unduh Format Vektor (PDF)", class = "btn-sm btn-secondary")
          )
        ),
        card_body(
          layout_columns(
            col_widths = c(4, 4, 4),
            radioButtons(
              "rasch_wm_view_mode",
              "Pilihan Skala Wright Map:",
              choices = c("Unidimensi" = "uni", "Multidimensi" = "multi"),
              selected = "uni"
            ),
            radioButtons(
              "rasch_wm_style",
              "Gaya Visualisasi Wright Map:",
              choices = c("Gaya Akademis KAWE (High Contrast)" = "kawe", "Format Standar WrightMap" = "standard"),
              selected = "kawe"
            ),
            sliderInput("rasch_wm_cex_scale", "Skala Ukuran Teks:", min = 0.6, max = 1.4, value = 0.95, step = 0.05)
          ),
          hr(class = "my-2"),
          plotOutput("rasch_wright_map_display", height = "760px", width = "100%"),
          hr(class = "my-3"),
          div(
            class = "alert alert-light border my-2 py-2 small",
            tags$strong("Prinsip Wright Map (Wilson, 2005; Linacre, 2002): "),
            "Wright Map menjajarkan sebaran kemampuan responden (θ) di panel kiri dengan tingkat kesulitan butir / batas kategori (β) di panel kanan pada kontinum skala logit yang sama. Pada pemodelan multidimensi, distribusi kemampuan untuk setiap dimensi ditampilkan berdampingan guna mengevaluasi targeting alat ukur terhadap atribut laten partisipan."
          )
        )
      )
    ),
    
    # 3.5 Perbandingan Model Fit (RSM, PCM, MRSM, MPCM)
    nav_panel(
      title = "Perbandingan Fit (RSM/PCM/MRSM/MPCM)",
      value = "tab_rasch_comp",
      icon = icon("balance-scale"),
      layout_sidebar(
        sidebar = sidebar(
          width = 340,
          title = "Komparasi Model Rasch (TAM)",
          p(class = "small text-muted mb-3", "Membandingkan kesesuaian empiris antara 4 model politomi IRT: RSM vs PCM (Unidimensi) dan MRSM vs MPCM (Multidimensi)."),
          actionButton("btn_run_rasch_comp", "Bandingkan 4 Model Rasch", class = "btn-primary w-100 mb-3", icon = icon("calculator")),
          div(
            class = "p-2 bg-light border rounded small",
            strong("Klasifikasi Model IRT yang Diuji:"),
            tags$ul(
              class = "ps-3 mb-0 mt-1",
              tags$li(strong("1. RSM (Rating Scale):"), " Asumsi jarak ambang batas kategori identik antar-butir (unidimensi)."),
              tags$li(strong("2. PCM (Partial Credit):"), " Jarak ambang batas diestimasi secara bebas per-butir (unidimensi)."),
              tags$li(strong("3. MRSM (Multidimensional RSM):"), " Asumsi ambang batas invarian pada struktur dimensi jamak."),
              tags$li(strong("4. MPCM (Multidimensional PCM):"), " Parameter ambang batas bebas per-butir pada struktur dimensi jamak.")
            )
          ),
          hr(class = "my-2"),
          downloadButton("btn_dl_rasch_comp_csv", "Unduh CSV Tabel Komparasi", class = "btn-outline-success w-100")
        ),
        div(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Model Terbaik (AIC)",
              value = textOutput("vb_rasch_comp_best_aic"),
              showcase = icon("trophy"),
              theme = "primary"
            ),
            value_box(
              title = "Model Terbaik (BIC)",
              value = textOutput("vb_rasch_comp_best_bic"),
              showcase = icon("award"),
              theme = "info"
            ),
            value_box(
              title = "Deviance Terendah (-2LL)",
              value = textOutput("vb_rasch_comp_best_dev"),
              showcase = icon("chart-line"),
              theme = "success"
            ),
            value_box(
              title = "Evaluasi Parsimoni",
              value = textOutput("vb_rasch_comp_parsimony"),
              showcase = icon("balance-scale"),
              theme = "secondary"
            )
          ),
          hr(class = "my-3"),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("table"), " Tabel Komparasi Kesesuaian Model IRT / Rasch (Global Fit & Information Criteria)"),
            card_body(
              DTOutput("table_rasch_comp_summary")
            )
          ),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("vial"), " Uji Rasio Kemiripan Model (Likelihood Ratio Test / ANOVA: RSM vs PCM & MRSM vs MPCM)"),
            card_body(
              tableOutput("table_rasch_comp_anova")
            )
          ),
          card(
            class = "mb-3",
            card_header(class = "bg-light fw-bold text-primary", icon("chart-bar"), " Visualisasi Komparasi Kriteria Informasi (Kriteria AIC & BIC)"),
            card_body(
              plotOutput("rasch_comp_bar_plot", height = "420px")
            )
          ),
          card(
            card_header(class = "bg-light fw-bold text-primary", icon("graduation-cap"), " Catatan Metodologis & Evaluasi Psikometrika Akademis"),
            card_body(
              div(
                class = "alert alert-light border mb-0 small",
                style = "line-height: 1.7;",
                p(strong("Dilema RSM vs PCM (Invariansi Ambang Batas):"),
                  " Rating Scale Model (RSM) merupakan model yang paling hemat parameter (parsimoni tertinggi) karena mengasumsikan parameter transisi antarkategori respons (step parameters) konstan untuk seluruh butir instrumen. Sebaliknya, Partial Credit Model (PCM) merelaksasi asumsi tersebut dengan mengestimasi ambang batas secara spesifik pada masing-masing butir. Penurunan nilai Deviance yang signifikan pada uji LRT mengindikasikan bahwa asumsi jarak kategori yang seragam tidak didukung secara empiris oleh data respons."),
                p(strong("Implikasi Multidimensionalitas (MRSM & MPCM):"),
                  " Pemodelan multidimensi memperhitungkan heterogenitas atribut laten partisipan secara simultan melalui matriks varians-kovarians. Pengabaian dimensi laten yang terbukti terpisah secara empiris dapat memicu atenuasi reliabilitas dan distorsi estimasi kemampuan individu (θ).")
              )
            )
          )
        )
      )
    )
  ),

  # Dropdown Menu: Meta-Analisis
  nav_menu(
    title = "Meta-Analisis",
    icon = icon("chart-pie"),
    
    # Tab 2: Kalkulator Effect Size
    nav_panel(
      value = "tab_calc",
      icon = icon("calculator"),
      title = "Effect Size",
      layout_columns(
        col_widths = c(6, 6),
        
        # Left Column: Input Form (Structured 2-card layout)
        div(
          card(
            class = "mb-3",
            card_header(
              span(icon("calculator"), " HASIL EKSTRAKSI DATA STUDI PRIMER")
            ),
            card_body(
              selectInput(
                "input_type",
                "Format Pelaporan Data Hasil Studi Primer:",
                choices = c(
                  "Perubahan Skor Prates-Pascates (Mean Difference / Gain Score)" = "mean_gain",
                  "Skor Pascates Dua Kelompok (Between-Group Means & SDs)" = "mean_sd",
                  "Rata-Rata & Standar Error (Means & SE)" = "mean_se",
                  "Statistik Uji-t Dua Kelompok Independen (Independent t-test)" = "t_test",
                  "Signifikansi Nilai p Dua Kelompok (Two-tailed p-value & n)" = "p_val",
                  "Statistik ANOVA Satu Arah Dua Kelompok (One-Way F-test)" = "f_test",
                  "Tabel Kontingensi 2×2 (Proporsi Biner / Kejadian Dikotomi)" = "table_2x2",
                  "Koefisien Korelasi Bivariat (Pearson's r)" = "correlation"
                ),
                selected = "mean_gain"
              ),
              
              # Sub-panel for esc_mean_gain
              conditionalPanel(
                condition = "input.input_type == 'mean_gain'",
                layout_columns(
                  col_widths = c(6, 6),
                  radioButtons(
                    "gain_sub_type",
                    "Format Pelaporan Data Pengukuran:",
                    choices = c(
                      "Prates & Pascates Lengkap (Means, SDs, n)" = "pre_post_full",
                      "Selisih Perubahan Langsung (Mean Gain / Change Score)" = "gain_direct"
                    ),
                    selected = "pre_post_full"
                  ),
                  selectInput(
                    "gain_outcome_dir",
                    "Arah Efektivitas Skala Pengukuran:",
                    choices = c(
                      "Skor Pascates Menurun = Pemulihan (misal: Depresi, Kecemasan)" = "decrease_better",
                      "Skor Pascates Meningkat = Pemulihan (misal: Resiliensi, Efikasi)" = "increase_better"
                    )
                  )
                ),
                
                conditionalPanel(
                  condition = "input.gain_sub_type == 'pre_post_full'",
                  h6("Kelompok Eksperimen / Intervensi (Treatment Group):", class = "text-primary fw-bold mt-2"),
                  layout_columns(
                    col_widths = c(4),
                    numericInput("pre1_n", "Ukuran Sampel Eksperimen (n1):", value = 40, min = 2, step = 1)
                  ),
                  layout_columns(
                    col_widths = c(6, 6),
                    numericInput("pre1_m", "Rata-Rata Prates (Pre Mean):", value = 28.4, step = 0.1),
                    numericInput("pre1_sd", "Standar Deviasi Prates (Pre SD):", value = 5.2, min = 0.01, step = 0.1)
                  ),
                  layout_columns(
                    col_widths = c(6, 6),
                    numericInput("post1_m", "Rata-Rata Pascates (Post Mean):", value = 16.2, step = 0.1),
                    numericInput("post1_sd", "Standar Deviasi Pascates (Post SD):", value = 4.8, min = 0.01, step = 0.1)
                  ),
                  h6("Kelompok Pembanding / Kontrol (Control Group):", class = "text-primary fw-bold mt-2"),
                  layout_columns(
                    col_widths = c(6),
                    numericInput("pre2_n", "Ukuran Sampel Kontrol (n2):", value = 40, min = 2, step = 1)
                  ),
                  layout_columns(
                    col_widths = c(6, 6),
                    numericInput("pre2_m", "Rata-Rata Prates (Pre Mean):", value = 27.9, step = 0.1),
                    numericInput("pre2_sd", "Standar Deviasi Prates (Pre SD):", value = 5.0, min = 0.01, step = 0.1)
                  ),
                  layout_columns(
                    col_widths = c(6, 6),
                    numericInput("post2_m", "Rata-Rata Pascates (Post Mean):", value = 24.1, step = 0.1),
                    numericInput("post2_sd", "Standar Deviasi Pascates (Post SD):", value = 5.1, min = 0.01, step = 0.1)
                  )
                ),
                
                conditionalPanel(
                  condition = "input.gain_sub_type == 'gain_direct'",
                  h6("Kelompok Eksperimen / Intervensi (Treatment Group):", class = "text-primary fw-bold mt-2"),
                  layout_columns(
                    col_widths = c(4, 4, 4),
                    numericInput("gain1_m", "Rata-Rata Selisih (Gain Mean):", value = 12.2, step = 0.1),
                    numericInput("gain1_sd", "Standar Deviasi Selisih (Gain SD):", value = 4.5, min = 0.01, step = 0.1),
                    numericInput("gain1_n", "Ukuran Sampel (n1):", value = 40, min = 2, step = 1)
                  ),
                  h6("Kelompok Pembanding / Kontrol (Control Group):", class = "text-secondary fw-bold mt-2"),
                  layout_columns(
                    col_widths = c(4, 4, 4),
                    numericInput("gain2_m", "Rata-Rata Selisih (Gain Mean):", value = 3.8, step = 0.1),
                    numericInput("gain2_sd", "Standar Deviasi Selisih (Gain SD):", value = 4.2, min = 0.01, step = 0.1),
                    numericInput("gain2_n", "Ukuran Sampel (n2):", value = 40, min = 2, step = 1)
                  )
                ),
                
                layout_columns(
                  col_widths = c(6),
                  numericInput("gain_r", "Korelasi Prates-Pascates (r):", value = 0.50, min = -0.99, max = 0.99, step = 0.05)
                ),
                
                layout_columns(
                  col_widths = c(6),
                  div(class = "small text-muted pt-9", "Asumsi Rosenthal (1991) r = 0.50 bila korelasi prates-pascates tidak dilaporkan dalam artikel primer.")
                )
              ),
              
              # Means & SDs
              conditionalPanel(
                condition = "input.input_type == 'mean_sd'",
                selectInput(
                  "sd_outcome_dir",
                  "Arah Efektivitas Skala Pengukuran:",
                  choices = c(
                    "Skor Pascates Menurun = Pemulihan (misal: Depresi, Kecemasan)" = "decrease_better",
                    "Skor Pascates Meningkat = Pemulihan (misal: Resiliensi, Kesejahteraan)" = "increase_better"
                  ),
                  selected = "decrease_better"
                ),
                h6("Kelompok Eksperimen / Intervensi (Treatment Group):", class = "text-primary fw-bold"),
                layout_columns(
                  col_widths = c(4, 4, 4),
                  numericInput("grp1m_sd", "Rata-Rata (Mean):", value = 14.5, step = 0.1),
                  numericInput("grp1sd_sd", "Standar Deviasi (SD):", value = 3.8, min = 0.01, step = 0.1),
                  numericInput("grp1n_sd", "Ukuran Sampel (n1):", value = 50, min = 2, step = 1)
                ),
                h6("Kelompok Pembanding / Kontrol (Control Group):", class = "text-secondary fw-bold mt-2"),
                layout_columns(
                  col_widths = c(4, 4, 4),
                  numericInput("grp2m_sd", "Rata-Rata (Mean):", value = 18.2, step = 0.1),
                  numericInput("grp2sd_sd", "Standar Deviasi (SD):", value = 4.1, min = 0.01, step = 0.1),
                  numericInput("grp2n_sd", "Ukuran Sampel (n2):", value = 48, min = 2, step = 1)
                )
              ),
              
              # Means & SE
              conditionalPanel(
                condition = "input.input_type == 'mean_se'",
                selectInput(
                  "se_outcome_dir",
                  "Arah Efektivitas Skala Pengukuran:",
                  choices = c(
                    "Skor Pascates Menurun = Pemulihan (misal: Depresi, Kecemasan)" = "decrease_better",
                    "Skor Pascates Meningkat = Pemulihan (misal: Resiliensi, Kesejahteraan)" = "increase_better"
                  ),
                  selected = "decrease_better"
                ),
                h6("Kelompok Eksperimen / Intervensi (Treatment Group):", class = "text-primary fw-bold"),
                layout_columns(
                  col_widths = c(4, 4, 4),
                  numericInput("grp1m_se", "Rata-Rata (Mean):", value = 14.5, step = 0.1),
                  numericInput("grp1se_se", "Standar Error (SE):", value = 0.54, min = 0.001, step = 0.01),
                  numericInput("grp1n_se", "Ukuran Sampel (n1):", value = 50, min = 2, step = 1)
                ),
                h6("Kelompok Pembanding / Kontrol (Control Group):", class = "text-secondary fw-bold mt-2"),
                layout_columns(
                  col_widths = c(4, 4, 4),
                  numericInput("grp2m_se", "Rata-Rata (Mean):", value = 18.2, step = 0.1),
                  numericInput("grp2se_se", "Standar Error (SE):", value = 0.59, min = 0.001, step = 0.01),
                  numericInput("grp2n_se", "Ukuran Sampel (n2):", value = 48, min = 2, step = 1)
                )
              ),
              
              # t-test
              conditionalPanel(
                condition = "input.input_type == 't_test'",
                numericInput("t_val", "Nilai Statistik Uji-t (Independent t):", value = 2.85, step = 0.01),
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput("t_n1", "Ukuran Sampel Intervensi (n1):", value = 45, min = 2, step = 1),
                  numericInput("t_n2", "Ukuran Sampel Kontrol (n2):", value = 45, min = 2, step = 1)
                )
              ),
              
              # p-value
              conditionalPanel(
                condition = "input.input_type == 'p_val'",
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput("p_value", "Nilai Signifikansi p (Two-tailed p):", value = 0.015, min = 0.0001, max = 0.9999, step = 0.005),
                  selectInput("p_direction", "Arah Efek Intervensi:", choices = c("Memihak Kelompok Intervensi" = "pos", "Memihak Kelompok Kontrol" = "neg"))
                ),
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput("p_n1", "Ukuran Sampel Intervensi (n1):", value = 50, min = 2, step = 1),
                  numericInput("p_n2", "Ukuran Sampel Kontrol (n2):", value = 50, min = 2, step = 1)
                )
              ),
              
              # F-test
              conditionalPanel(
                condition = "input.input_type == 'f_test'",
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput("f_val", "Nilai Statistik Uji-F (df1 = 1):", value = 7.42, min = 0.001, step = 0.1),
                  selectInput("f_direction", "Arah Efek Intervensi:", choices = c("Memihak Kelompok Intervensi" = "pos", "Memihak Kelompok Kontrol" = "neg"))
                ),
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput("f_n1", "Ukuran Sampel Intervensi (n1):", value = 50, min = 2, step = 1),
                  numericInput("f_n2", "Ukuran Sampel Kontrol (n2):", value = 50, min = 2, step = 1)
                )
              ),
              
              # 2x2 table
              conditionalPanel(
                condition = "input.input_type == 'table_2x2'",
                h6("Kelompok Eksperimen / Intervensi (Treatment Group):", class = "text-primary fw-bold"),
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput("tab_e1", "Jumlah Kejadian / Respons (Event):", value = 35, min = 0, step = 1),
                  numericInput("tab_n1", "Total Ukuran Sampel (n1):", value = 50, min = 1, step = 1)
                ),
                h6("Kelompok Pembanding / Kontrol (Control Group):", class = "text-secondary fw-bold mt-2"),
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput("tab_e2", "Jumlah Kejadian / Respons (Event):", value = 18, min = 0, step = 1),
                  numericInput("tab_n2", "Total Ukuran Sampel (n2):", value = 50, min = 1, step = 1)
                )
              ),
              
              # Correlation
              conditionalPanel(
                condition = "input.input_type == 'correlation'",
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput("corr_r", "Koefisien Korelasi Pearson (r):", value = 0.35, min = -0.999, max = 0.999, step = 0.01),
                  numericInput("corr_n", "Total Ukuran Sampel (N):", value = 90, min = 4, step = 1)
                )
              )
            )
          ),
          
          card(
            card_header(
              span(icon("tag"), " IDENTIFIKASI STUDI PRIMER & KODE MODERATOR")
            ),
            card_body(
              layout_columns(
                col_widths = c(6, 6),
                textInput("study_name", "Label Studi (Peneliti & Tahun / Author et al.):", value = "Studi Baru et al."),
                numericInput("study_year", "Tahun Publikasi (Publication Year):", value = 2024, min = 1950, max = 2035, step = 1)
              ),
              layout_columns(
                col_widths = c(6, 6),
                selectInput(
                  "study_interv",
                  "Modalitas Intervensi (Intervention):",
                  choices = c(
                    "CBT (Cognitive Behavioral Therapy)" = "CBT",
                    "ACT (Acceptance and Commitment Therapy)" = "ACT",
                    "BA (Behavioral Activation)" = "BA",
                    "iCBT (Internet-based CBT)" = "iCBT",
                    "MBI (Mindfulness-Based Intervention)" = "Mindfulness",
                    "IPT (Interpersonal Psychotherapy)" = "IPT",
                    "Psychodynamic Therapy" = "Psychodynamic",
                    "Lainnya (Other Interventions)" = "Lainnya"
                  )
                ),
                selectInput(
                  "study_format",
                  "Format Pengantaran (Delivery Format):",
                  choices = c(
                    "Individual (Face-to-Face)" = "Individual",
                    "Kelompok (Group Face-to-Face)" = "Group",
                    "Digital Terpandu (Guided Digital / iCBT)" = "Guided Digital",
                    "Mandiri Tanpa Panduan (Unguided Self-Help)" = "Unguided Self-Help"
                  )
                )
              ),
              layout_columns(
                col_widths = c(6, 6),
                selectInput(
                  "study_ctrl",
                  "Kondisi Pembanding / Kontrol (Control Condition):",
                  choices = c(
                    "Daftar Tunggu (Waitlist / WL)" = "Waitlist",
                    "Perawatan Rutin (Treatment as Usual / TAU)" = "TAU",
                    "Plasebo / Kontrol Atensi (Placebo / Attention)" = "Placebo",
                    "Komparator Terapi Aktif (Active Comparator / Control)" = "Active Control"
                  )
                ),
                conditionalPanel(
                  condition = "input.study_ctrl == 'Active Control'",
                  textInput("study_active_comparator", "Spesifikasi Komparator Aktif:", value = "Supportive Counseling", placeholder = "Misal: Supportive Counseling, Relaxation, Psikoedukasi")
                ),
                conditionalPanel(
                  condition = "input.study_ctrl != 'Active Control'",
                  p(class = "small text-muted pt-4 mb-0", "Komparator standar tanpa komponen perlakuan psikologis aktif.")
                )
              ),
              layout_columns(
                col_widths = c(6, 6),
                selectInput(
                  "study_pop",
                  "Karakteristik Populasi Klinis (Population):",
                  choices = c(
                    "Klinis: Gangguan Depresi Mayor (Clinical MDD)" = "Clinical (MDD)",
                    "Klinis: Gangguan Kecemasan (Clinical Anxiety)" = "Clinical (Anxiety)",
                    "Klinis: Gangguan Stres Pascatrauma (Clinical PTSD)" = "Clinical (PTSD / Trauma)",
                    "Subklinis (Subclinical Elevated Symptoms)" = "Subclinical",
                    "Populasi Umum / Komunitas (Community / Non-clinical)" = "Community"
                  )
                ),
                selectInput(
                  "study_age",
                  "Kohort Usia Partisipan (Age Cohort):",
                  choices = c(
                    "Dewasa (Adults, 18–64 tahun)" = "Adults",
                    "Remaja (Adolescents, 12–17 tahun)" = "Adolescents",
                    "Lansia (Older Adults, ≥ 65 tahun)" = "Older Adults",
                    "Anak-anak (Children, < 12 tahun)" = "Children"
                  ),
                  selected = "Adults"
                )
              ),
              layout_columns(
                col_widths = c(4, 4, 4),
                selectInput(
                  "study_therapist",
                  "Kualifikasi Terapis / Fasilitator:",
                  choices = c(
                    "Klinisi Berlisensi (Licensed Clinician)" = "Licensed Clinician",
                    "Terapis Tersupervisi (Supervised Trainee)" = "Supervised Trainee",
                    "Sistem Digital Otomatis (Automated Digital)" = "Automated Digital"
                  )
                ),
                selectInput(
                  "study_rob",
                  "Risiko Bias Metodologis (RoB 2):",
                  choices = c(
                    "Risiko Rendah (Low Risk of Bias)" = "Low Risk",
                    "Kekhawatiran Terbatas (Some Concerns)" = "Some Concerns",
                    "Risiko Tinggi (High Risk of Bias)" = "High Risk"
                  )
                ),
                numericInput("study_dur", "Durasi Intervensi (Minggu):", value = 12, min = 1, step = 1)
              ),
              layout_columns(
                col_widths = c(6, 6),
                numericInput("study_sessions", "Jumlah Sesi Terapi (Total Sessions):", value = 12, min = 1, step = 1),
                radioButtons(
                  "es_metric",
                  "Metrik Standardized Mean Difference (SMD):",
                  choices = c(
                    "Hedges' g (Terkoreksi Bias / Unbiased)" = "g",
                    "Cohen's d (Estimasi Bias Sampel Kecil)" = "d"
                  ),
                  selected = "g",
                  inline = TRUE
                )
              ),
              div(
                class = "d-flex gap-2 mt-2",
                actionButton("btn_calculate", "Hitung Parameter Efek", class = "btn-outline-primary flex-grow-1", icon = icon("calculator")),
                actionButton("btn_append_study", "Simpan ke Matriks Studi", class = "btn-primary flex-grow-1", icon = icon("plus-circle")),
                actionButton("btn_reset_calc_form", "Atur Ulang Formulir", class = "btn-outline-secondary", icon = icon("redo"))
              )
            )
          )
        ),
        
        # Right Column: Output & Realtime History
        div(
          card(
            card_header(
              span(icon("chart-line"), " BESARAN EFEK/EFFECT SIZE")
            ),
            card_body(
              navset_tab(
                id = "calc_results_navset",
                
                # Panel 1: Hasil Estimasi & Ringkasan
                nav_panel(
                  title = span(icon("poll"), " Ringkasan Estimasi"),
                  value = "panel_calc_results",
                  layout_columns(
                    col_widths = c(4, 4, 4),
                    div(
                      class = "kawe-kpi-card text-center p-3 rounded border",
                      style = "background: #fdfcfb; border: 1px solid #d9d2c7 !important; border-radius: 8px;",
                      div(
                        class = "d-flex justify-content-center align-items-center mb-1 gap-2",
                        span(icon("balance-scale"), style = "font-size: 1.15rem; color: #6b1d2f;"),
                        span(class = "fw-bold text-uppercase small", style = "letter-spacing: 0.5px; color: #495057;", "Besaran Efek / Effect Size")
                      ),
                      div(
                        class = "fw-bold",
                        style = "font-size: 1.9rem; line-height: 1.2; color: #6b1d2f; font-family: 'Cinzel', serif, Georgia;",
                        textOutput("vb_calc_es", inline = TRUE)
                      )
                    ),
                    div(
                      class = "kawe-kpi-card text-center p-3 rounded border",
                      style = "background: #fdfcfb; border: 1px solid #d9d2c7 !important; border-radius: 8px;",
                      div(
                        class = "d-flex justify-content-center align-items-center mb-1 gap-2",
                        span(icon("wave-square"), style = "font-size: 1.15rem; color: #8c6d3b;"),
                        span(class = "fw-bold text-uppercase small", style = "letter-spacing: 0.5px; color: #495057;", "Varians Sampling (vi)")
                      ),
                      div(
                        class = "fw-bold",
                        style = "font-size: 1.9rem; line-height: 1.2; color: #6b1d2f; font-family: 'Cinzel', serif, Georgia;",
                        textOutput("vb_calc_var", inline = TRUE)
                      )
                    ),
                    div(
                      class = "kawe-kpi-card text-center p-3 rounded border",
                      style = "background: #fdfcfb; border: 1px solid #d9d2c7 !important; border-radius: 8px;",
                      div(
                        class = "d-flex justify-content-center align-items-center mb-1 gap-2",
                        span(icon("arrows-alt-h"), style = "font-size: 1.15rem; color: #31708f;"),
                        span(class = "fw-bold text-uppercase small", style = "letter-spacing: 0.5px; color: #495057;", "Standar Error (SE)")
                      ),
                      div(
                        class = "fw-bold",
                        style = "font-size: 1.9rem; line-height: 1.2; color: #6b1d2f; font-family: 'Cinzel', serif, Georgia;",
                        textOutput("vb_calc_se", inline = TRUE)
                      )
                    )
                  ),
                  
                  hr(class = "my-2"),
                  
                  h6("Visualisasi Effect Size dengan 95% CI:", class = "fw-bold text-dark mb-1"),
                  plotOutput("calc_visual_preview_plot", height = "135px"),
                  
                  h6("Hasil Effect Size:", class = "fw-bold text-dark mt-2"),
                  tableOutput("calc_summary_table"),
                  
                  # Quick callout for active formula
                  uiOutput("calc_quick_formula_callout"),
                  
                  div(
                    class = "mt-3",
                    actionButton("btn_view_calc_compendium", "Buka Kumpulan Rumus dan Panduan Perhitungan Lengkap →", class = "btn-sm btn-outline-secondary w-100", icon = icon("book-open"))
                  )
                ),
                
                # Panel 2: Kompendium Seluruh Rumus Effect Size
                nav_panel(
                  title = span(icon("book-open"), " Kumpulan Rumus dan Panduan Perhitungan"),
                  value = "panel_formula_catalog",
                  div(
                    class = "pt-2",
                    div(
                      class = "mb-2",
                      actionButton("btn_back_to_calc_results2", "← Kembali ke Ringkasan", class = "btn-sm btn-outline-secondary", icon = icon("arrow-left"))
                    ),
                    accordion(
                      id = "acc_all_formulas",
                      open = FALSE,
                      
                      # 1. Post-Test Means & SDs
                      accordion_panel(
                        title = span(icon("chart-bar"), " 1. Standardized Mean Difference: Rata-rata Pascates & Deviasi Baku (Post-Test Means & SDs)"),
                        value = "acc_post_means",
                        div(
                          p("Desain eksperimen dua kelompok independen (treatment vs control) pascates mengevaluasi perbedaan rata-rata terstandardisasi. Formulasi baku didasarkan pada Cohen (1988) dan Hedges (1981):"),
                          div(class = "kawe-eq-box",
                            HTML("$$\\text{Selisih Rata-rata: } \\Delta\\bar{X} = \\bar{X}_1 - \\bar{X}_2$$"),
                            HTML("$$\\text{Deviasi Baku Gabungan: } s_{\\text{pooled}} = \\sqrt{\\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}$$"),
                            HTML("$$\\text{Cohen's } d = \\frac{\\bar{X}_1 - \\bar{X}_2}{s_{\\text{pooled}}}$$"),
                            HTML("$$\\text{Faktor Koreksi Bias Hedges: } J = 1 - \\frac{3}{4(n_1 + n_2) - 9} = 1 - \\frac{3}{4df - 1}$$"),
                            HTML("$$\\text{Hedges' } g = J \\times d$$"),
                            HTML("$$\\text{Varians Sampling: } v_d = \\frac{n_1 + n_2}{n_1 n_2} + \\frac{d^2}{2(n_1 + n_2)}, \\quad v_g = J^2 \\times v_d$$"),
                            HTML("$$\\text{Standar Error: } SE_g = \\sqrt{v_g}, \\quad 95\\% \\text{ CI} = [g - 1.96 \\cdot SE_g, \\; g + 1.96 \\cdot SE_g]$$")
                          ),
                          p(class = "small text-muted mt-2",
                            strong("Asumsi & Konvensi: "),
                            "Diasumsikan varians homogen antar-kelompok (homoscedasticity). Pada ukuran sampel kecil (n < 50), penggunaan Hedges' g wajib menggantikan Cohen's d guna mengoreksi upward positive bias (Hedges, 1981; Cochrane Handbook 2023)."
                          )
                        )
                      ),
                      
                      # 2. Mean Gain Scores
                      accordion_panel(
                        title = span(icon("exchange-alt"), " 2. Mean Gain / Perubahan Pre-Post (Repeated Measures Designs)"),
                        value = "acc_gain",
                        div(
                          p("Pada studi klinis dengan desain pretest-posttest-control group (PPCR), besaran efek dihitung dari selisih perubahan skor (gain) antara kelompok intervensi dan kontrol (Becker, 1988; Morris, 2008):"),
                          div(class = "kawe-eq-box",
                            HTML("$$\\text{Gain Intervensi: } \\text{Gain}_1 = \\bar{X}_{post,1} - \\bar{X}_{pre,1}, \\quad \\text{Gain Kontrol: } \\text{Gain}_2 = \\bar{X}_{post,2} - \\bar{X}_{pre,2}$$"),
                            HTML("$$\\text{Deviasi Baku Perubahan: } s_{\\text{change}, i} = \\sqrt{s_{pre, i}^2 + s_{post, i}^2 - 2 r s_{pre, i} s_{post, i}}$$"),
                            HTML("$$\\text{Deviasi Baku Gabungan: } s_{\\text{pooled, change}} = \\sqrt{\\frac{(n_1 - 1)s_{\\text{change},1}^2 + (n_2 - 1)s_{\\text{change},2}^2}{n_1 + n_2 - 2}}$$"),
                            HTML("$$d_{\\text{gain}} = \\frac{\\text{Gain}_1 - \\text{Gain}_2}{s_{\\text{pooled, change}}}, \\quad g_{\\text{gain}} = J \\times d_{\\text{gain}}$$"),
                            HTML("$$v_d = \\frac{2(1 - r_1)}{n_1} + \\frac{2(1 - r_2)}{n_2} + \\frac{d^2}{2(n_1 + n_2)}, \\quad v_g = J^2 \\times v_d$$")
                          ),
                          p(class = "small text-muted mt-2",
                            strong("Asumsi Korelasi Pre-Post (r): "),
                            "Bila artikel primer tidak melaporkan r pre-post secara eksplisit, konvensi Rosenthal (1991) menetapkan nilai default konservatif r = 0.50 (Morris, 2008)."
                          )
                        )
                      ),
                      
                      # 3. Means & SE
                      accordion_panel(
                        title = span(icon("arrows-alt-h"), " 3. Rekonstruksi dari Mean & Standar Error (SE)"),
                        value = "acc_mean_se",
                        div(
                          p("Jika studi primer hanya melaporkan Standar Error (SE) dan ukuran sampel (n), deviasi baku kelompok (s) direkonstruksi matematis:"),
                          div(class = "kawe-eq-box",
                            HTML("$$s_1 = SE_1 \\times \\sqrt{n_1}, \\quad s_2 = SE_2 \\times \\sqrt{n_2}$$"),
                            HTML("$$s_{\\text{pooled}} = \\sqrt{\\frac{(n_1 - 1)(SE_1^2 \\cdot n_1) + (n_2 - 1)(SE_2^2 \\cdot n_2)}{n_1 + n_2 - 2}}$$"),
                            HTML("$$d = \\frac{\\bar{X}_1 - \\bar{X}_2}{s_{\\text{pooled}}}, \\quad g = J \\times d$$")
                          )
                        )
                      ),
                      
                      # 4. Inferential t, F, p
                      accordion_panel(
                        title = span(icon("balance-scale"), " 4. Konversi dari Statistik Inferensial Parametrik (Uji-t, Uji-F, Nilai-p)"),
                        value = "acc_inferential",
                        div(
                          p("Ketika mean dan deviasi baku tidak dilaporkan, besaran efek direkonstruksi dari statistik uji parametrik (Dunlap et al., 1996; Lipsey & Wilson, 2001):"),
                          div(class = "kawe-eq-box",
                            HTML("$$\\text{Uji-t Independen: } d = t \\times \\sqrt{\\frac{n_1 + n_2}{n_1 n_2}}, \\quad v_d = \\frac{n_1 + n_2}{n_1 n_2} + \\frac{d^2}{2(n_1 + n_2)}$$"),
                            HTML("$$\\text{Uji-F Satu Arah (} df_1=1 \\text{): } |t| = \\sqrt{F}, \\quad d = \\pm \\sqrt{F \\left(\\frac{n_1 + n_2}{n_1 n_2}\\right)}$$"),
                            HTML("$$\\text{Nilai-p Dua Arah: } t = F_{t, df}^{-1}(1 - p/2), \\quad df = n_1 + n_2 - 2$$"),
                            HTML("$$g = J \\times d, \\quad v_g = J^2 \\times v_d$$")
                          ),
                          p(class = "small text-muted mt-2",
                            strong("Arah Efek: "),
                            "Statistik F dan nilai-p tidak memuat informasi tanda matematis (+/-); penentuan arah keunggulan intervensi wajib diverifikasi dari laporan deskriptif naratif studi primer."
                          )
                        )
                      ),
                      
                      # 5. Binary & 2x2 Table
                      accordion_panel(
                        title = span(icon("th-large"), " 5. Luaran Dikotomis / Tabel Kontingensi 2x2 (Odds Ratio ke SMD)"),
                        value = "acc_binary",
                        div(
                          p("Konversi rasio odds (OR) dari luaran biner (misalnya respons klinis, remisi) ke skala standardized mean difference d (Chinn, 2000; Hasselblad & Hedges, 1995):"),
                          div(class = "kawe-eq-box",
                            HTML("$$\\text{Odds Ratio: } OR = \\frac{a \\cdot d}{b \\cdot c}, \\quad \\ln(OR) = \\ln(a) - \\ln(b) - \\ln(c) + \\ln(d)$$"),
                            HTML("$$\\text{Varians Woolf: } v_{\\ln(OR)} = \\frac{1}{a} + \\frac{1}{b} + \\frac{1}{c} + \\frac{1}{d}$$"),
                            HTML("$$\\text{Konversi ke } d: d = \\ln(OR) \\times \\frac{\\sqrt{3}}{\\pi} \\approx 0.5513 \\times \\ln(OR)$$"),
                            HTML("$$\\text{Varians } d: v_d = v_{\\ln(OR)} \\times \\frac{3}{\\pi^2} \\approx 0.3040 \\times v_{\\ln(OR)}, \\quad g = J \\times d$$")
                          ),
                          p(class = "small text-muted mt-2",
                            strong("Landasan Teoretis: "),
                            "Transformasi ini didasarkan pada asumsi variabel laten logistik kontinu yang mendasari dikotomisasi ambang batas luaran klinis (Chinn, 2000)."
                          )
                        )
                      ),
                      
                      # 6. Correlation r
                      accordion_panel(
                        title = span(icon("project-diagram"), " 6. Koefisien Asosiasi & Korelasi (r & Fisher z)"),
                        value = "acc_corr",
                        div(
                          p("Transformasi korelasi Pearson atau point-biserial r ke skala standardized mean difference d (Borenstein et al., 2009; Jacobs & Viechtbauer, 2017):"),
                          div(class = "kawe-eq-box",
                            HTML("$$d = \\frac{2r}{\\sqrt{1 - r^2}}, \\quad v_d = \\frac{4 v_r}{(1 - r^2)^3}, \\quad v_r = \\frac{(1 - r^2)^2}{N - 1}$$"),
                            HTML("$$\\text{Transformasi Fisher's } z: z = \\frac{1}{2}\\ln\\left(\\frac{1 + r}{1 - r}\\right) = \\text{arctanh}(r), \\quad v_z = \\frac{1}{N - 3}$$"),
                            HTML("$$g = J \\times d, \\quad J = 1 - \\frac{3}{4N - 9}$$")
                          )
                        )
                      ),
                      
                      # 7. J-Factor & Varians
                      accordion_panel(
                        title = span(icon("sliders-h"), " 7. Faktor Koreksi Bias Sampel Kecil Hedges (J) & Formulasi Varians Asimtotik"),
                        value = "acc_j_var",
                        div(
                          p("Formulasi eksak dan aproksimasi koreksi bias sampel kecil Hedges (1981):"),
                          div(class = "kawe-eq-box",
                            HTML("$$\\text{Faktor Koreksi Eksak (Gamma): } J_{\\text{exact}} = \\frac{\\Gamma(df/2)}{\\sqrt{df/2} \\cdot \\Gamma((df - 1)/2)}, \\quad df = n_1 + n_2 - 2$$"),
                            HTML("$$\\text{Aproksimasi Polinomial Hedges (1981): } J = 1 - \\frac{3}{4(n_1 + n_2) - 9} = 1 - \\frac{3}{4df - 1}$$"),
                            HTML("$$\\text{Varians Sampling Asimtotik: } v_d = \\frac{n_1 + n_2}{n_1 n_2} + \\frac{d^2}{2(n_1 + n_2)}$$"),
                            HTML("$$\\text{Varians Hedges' g: } v_g = J^2 \\times v_d, \\quad \\text{Bobot Studi (Weight): } w_i = \\frac{1}{v_g}$$")
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          ),
          
          card(
            class = "mt-3",
            card_header(
              class = "d-flex justify-content-between align-items-center",
              span(icon("history"), " STUDI TERSIMPAN"),
              uiOutput("calc_saved_count_badge")
            ),
            card_body(
              p(class = "small text-muted mb-2", "Berikut adalah daftar studi yang telah terkumpul dari kalkulator effect size"),
              tableOutput("table_calc_recent_preview"),
              div(
                class = "d-flex gap-2 mt-2",
                actionButton("btn_goto_compiled", "Buka Matriks Gabungan →", class = "btn-sm btn-outline-primary flex-grow-1", icon = icon("layer-group")),
                actionButton("btn_calc_direct_model", "Langsung ke Pemodelan Sintesis →", class = "btn-sm btn-primary flex-grow-1", icon = icon("chart-pie"))
              )
            )
          )
        )
      )
    ),
    # Tab 3: Kumpulan Effect Size
    nav_panel(
      value = "tab_compiled",
      icon = icon("layer-group"),
      title = "Kumpulan Effect Size",
      card(
        card_header(
          class = "bg-light d-flex justify-content-between align-items-center flex-wrap gap-2",
          span(class = "fw-bold text-primary", icon("layer-group"), " Matriks Studi Hasil Kompilasi Kalkulator"),
          div(
            class = "d-flex gap-1 flex-wrap",
            actionButton("btn_use_calc_as_active", "Jadikan Dataset Analisis Utama", class = "btn-sm btn-primary", icon = icon("check-circle")),
            actionButton("btn_calc_to_model", "Analisis di Pemodelan Meta-Analisis →", class = "btn-sm btn-primary", icon = icon("chart-pie")),
            actionButton("btn_append_calc_to_active", "Gabungkan ke Data Aktif", class = "btn-sm btn-outline-primary", icon = icon("plus")),
            actionButton("btn_load_calc_sample", "Muat Contoh Batch (k=5)", class = "btn-sm btn-outline-secondary", icon = icon("magic")),
            actionButton("btn_delete_calc_row", "Hapus Baris Terpilih", class = "btn-sm btn-outline-danger", icon = icon("minus-circle")),
            actionButton("btn_clear_calc_compiled", "Kosongkan Semua", class = "btn-sm btn-outline-danger", icon = icon("trash")),
            downloadButton("btn_download_calc_csv", "Unduh CSV", class = "btn-sm btn-outline-success")
          )
        ),
        card_body(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title = "Studi Terkompilasi",
              value = textOutput("vb_calc_comp_k"),
              showcase = icon("database"),
              theme = "primary"
            ),
            value_box(
              title = "Rata-rata Effect Size",
              value = textOutput("vb_calc_comp_mean_g"),
              showcase = icon("chart-bar"),
              theme = "info"
            ),
            value_box(
              title = "Rentang Besaran Efek",
              value = textOutput("vb_calc_comp_range_g"),
              showcase = icon("arrows-alt-h"),
              theme = "secondary"
            ),
            value_box(
              title = "Total Partisipan",
              value = textOutput("vb_calc_comp_n"),
              showcase = icon("users"),
              theme = "success"
            )
          ),
          
          uiOutput("calc_compiled_status_alert"),
          
          DTOutput("table_calc_compiled")
        )
      )
    ),
    
    # Tab 4: Data Sintesis
    nav_panel(
      value = "tab_data",
      icon = icon("database"),
      title = "Data Sintesis",
      card(
        card_header(
          class = "bg-light d-flex justify-content-between align-items-center flex-wrap gap-2",
          span(class = "fw-bold text-primary", icon("database"), " DATASET ANALISIS"),
          div(
            class = "d-flex gap-1 flex-wrap",
            actionButton("btn_load_benchmark", "Muat Data Benchmark (k=14)", class = "btn-sm btn-outline-primary", icon = icon("database")),
            actionButton("btn_delete_active_row", "Hapus Studi Terpilih", class = "btn-sm btn-outline-danger", icon = icon("minus-circle")),
            actionButton("btn_clear_data", "Reset Dataset", class = "btn-sm btn-outline-danger", icon = icon("trash")),
            downloadButton("btn_download_csv", "Unduh Data Aktif (.csv)", class = "btn-sm btn-outline-success"),
            downloadButton("btn_download_template_csv", "Unduh Template CSV", class = "btn-sm btn-outline-secondary"),
            actionButton("btn_data_to_model", "Lanjut ke Pemodelan Sintesis →", class = "btn-sm btn-primary", icon = icon("chart-pie"))
          )
        ),
        card_body(
          layout_columns(
            col_widths = c(4, 4, 4),
            div(
              h6("1. Pilih Sumber Dataset Aktif:", class = "fw-bold text-primary mb-2"),
              radioButtons(
                "dataset_source_choice",
                NULL,
                choices = c(
                  "Dataset Aplikasi (14 Uji Klinis RCT Psikologi)" = "benchmark",
                  "Kompilasi Hasil Kalkulator Effect Size" = "calculator",
                  "Unggah Berkas Mandiri (.csv / .xlsx)" = "upload"
                ),
                selected = "benchmark"
              ),
              conditionalPanel(
                condition = "input.dataset_source_choice == 'upload'",
                fileInput("upload_file", "Pilih Berkas (.csv / .xlsx):", accept = c(".csv", ".xlsx", ".xls"))
              )
            ),
            div(
              h6("2. Penyaring Cepat / Sub-Analisis Sensitivitas:", class = "fw-bold text-primary mb-2"),
              selectInput(
                "quick_sensitivity_filter",
                "Filter Inklusi Studi Efektivitas:",
                choices = c(
                  "Seluruh Studi dalam Dataset (Tanpa Filter)" = "all",
                  "Hanya Studi Risiko Bias Rendah (Low Risk of Bias)" = "low_rob",
                  "Kecualikan Kontrol Daftar Tunggu (Tanpa Waitlist)" = "no_waitlist",
                  "Hanya Populasi Klinis Terdiagnosis (Clinical MDD/Anxiety)" = "clinical_only"
                ),
                selected = "all"
              ),
              p(class = "small text-muted mb-0", "Klik dua kali pada sel tabel di bawah untuk menyunting nilai secara langsung, atau pilih baris dan klik 'Hapus Studi Terpilih'.")
            ),
            div(
              class = "p-3 bg-light rounded border",
              strong("Ketentuan Kolom Data Sintesis Bukti:"),
              tags$ul(
                class = "mb-0 ps-3 small text-muted",
                tags$li("Kolom besaran efek wajib: ", tags$code("yi"), " (atau ", tags$code("es"), ", ", tags$code("hedges_g"), ")"),
                tags$li("Kolom presisi wajib: ", tags$code("vi"), " (varians) atau ", tags$code("sei"), " (standar error)"),
                tags$li("Kolom moderator psikologi: ", tags$code("Control_Type"), ", ", tags$code("Format"), ", ", tags$code("Intervention"), ", ", tags$code("Population"), ", ", tags$code("Duration_Weeks"))
              )
            )
          ),
          
          hr(class = "my-3"),
          
          DTOutput("table_active_dataset")
        )
      )
    ),
    
    # Tab 5: Model Meta Analisis
    nav_panel(
      value = "tab_model",
      icon = icon("chart-pie"),
      title = "Model Meta Analisis",
      layout_sidebar(
        sidebar = sidebar(
          width = 330,
          title = "Spesifikasi Model Sintesis",
          
          selectInput(
            "meta_estimator",
            "Estimator Varians (tau^2):",
            choices = c(
              "Restricted Maximum Likelihood (REML)" = "REML",
              "DerSimonian-Laird (DL)" = "DL",
              "Maximum Likelihood (ML)" = "ML",
              "Paule-Mandel (PM)" = "PM",
              "Empirical Bayes (EB)" = "EB",
              "Sidik-Jonkman (SJ)" = "SJ",
              "Fixed-Effect (FE / Equal-Effects)" = "FE"
            ),
            selected = "REML"
          ),
          
          selectInput(
            "meta_test_type",
            "Metode Inferensi & Interval Kepercayaan:",
            choices = c(
              "Distribusi Normal Standar (z-test)" = "z",
              "Koreksi Hartung-Knapp / HKSJ (t-test)" = "knha"
            ),
            selected = "z"
          ),
          
          selectInput(
            "conf_level",
            "Tingkat Kepercayaan (CI):",
            choices = c("95%" = 0.95, "99%" = 0.99),
            selected = 0.95
          ),
          
          hr(),
          
          h6("Meta-Regresi & Klaster Multilevel:", class = "fw-bold text-dark"),
          selectInput(
            "meta_continuous_covar",
            "Kovariat Kontinu (Meta-Regresi):",
            choices = c("Tanpa Kovariat Kontinu" = "none")
          ),
          selectInput(
            "ml_cluster_var",
            "Variabel Klaster Level-3 (rma.mv):",
            choices = c(
              "Pendekatan Intervensi (Intervention)" = "Intervention",
              "Populasi Klinis (Population)" = "Population",
              "Tipe Kelompok Kontrol (Control_Type)" = "Control_Type",
              "Format Penghantaran (Format)" = "Format",
              "Identitas Studi (Study)" = "Study"
            ),
            selected = "Intervention"
          )
        ),
        
        card(
          card_header(
            class = "bg-light d-flex justify-content-between align-items-center flex-wrap gap-2",
            span(class = "fw-bold text-primary", icon("chart-pie"), " Rangkuman Metrik Sintesis & Heterogenitas"),
            div(
              class = "d-flex gap-1",
              actionButton("btn_model_to_subgroup", "Analisis Subkelompok →", class = "btn-sm btn-outline-primary", icon = icon("object-group")),
              actionButton("btn_model_to_forest", "Forest Plot →", class = "btn-sm btn-outline-primary", icon = icon("tree")),
              actionButton("btn_model_to_bias", "Diagnostik Bias →", class = "btn-sm btn-outline-secondary", icon = icon("shield-alt"))
            )
          ),
          card_body(
            layout_columns(
              col_widths = c(3, 3, 3, 3),
              value_box(
                title = "Pooled Effect Size (Hedges' g)",
                value = textOutput("vb_pooled_es"),
                showcase = icon("bullseye"),
                theme = "primary"
              ),
              value_box(
                title = "Heterogenitas (I^2)",
                value = textOutput("vb_hetero_i2"),
                showcase = icon("percent"),
                theme = "warning"
              ),
              value_box(
                title = "Varians Antarstudi (tau^2)",
                value = textOutput("vb_hetero_tau2"),
                showcase = icon("project-diagram"),
                theme = "info"
              ),
              value_box(
                title = "95% Prediction Interval",
                value = textOutput("vb_pred_interval"),
                showcase = icon("arrows-alt-h"),
                theme = "secondary"
              )
            ),
            
            hr(class = "my-3"),
            
            navset_card_tab(
              nav_panel(
                title = "Estimasi Model Formal (metafor)",
                div(
                  h6("Luaran Estimasi Formal Model Sintesis (metafor::rma):", class = "fw-bold text-primary mb-2"),
                  verbatimTextOutput("meta_model_summary")
                )
              ),
              nav_panel(
                title = "Meta-Regresi Kovariat Kontinu",
                uiOutput("meta_continuous_covar_ui")
              ),
              nav_panel(
                title = "Sensitivitas (Leave-One-Out) & Influensial",
                layout_columns(
                  col_widths = c(6, 6),
                  div(
                    h6("Diagnostik Studi Influensial & Outlier (Baujat Plot):", class = "fw-bold text-primary"),
                    plotOutput("meta_baujat_plot", height = "360px"),
                    p(class = "small text-muted mt-1", "Studi di kuadran kanan atas menyumbang heterogenitas Cochran Q terbesar sekaligus memberikan pengaruh paling kuat terhadap pergeseran efek gabungan.")
                  ),
                  div(
                    h6("Analisis Sensitivitas Leave-One-Out (Eksklusi Sekuensial):", class = "fw-bold text-primary"),
                    div(style = "max-height: 390px; overflow-y: auto;", tableOutput("table_leave1out"))
                  )
                )
              ),
              nav_panel(
                title = "Pemodelan Multilevel (3-Level rma.mv)",
                layout_columns(
                  col_widths = c(6, 6),
                  div(
                    h6("Dekomposisi Varians Tiga Aras (Cheung, 2014):", class = "fw-bold text-primary"),
                    tableOutput("table_multilevel_decomp"),
                    hr(class = "my-2"),
                    h6("Ringkasan Model Multilevel (rma.mv):", class = "fw-bold text-primary"),
                    verbatimTextOutput("multilevel_model_summary")
                  ),
                  div(
                    h6("Proporsi Varians Multilevel (Level 1 vs Level 2 vs Level 3):", class = "fw-bold text-primary"),
                    plotOutput("multilevel_variance_plot", height = "360px")
                  )
                )
              )
            )
          )
        )
      )
    ),
    
    # Tab 6: Forest Plot
    nav_panel(
      value = "tab_forest",
      icon = icon("tree"),
      title = "Forest Plot",
      layout_sidebar(
        sidebar = sidebar(
          width = 340,
          title = "Konfigurasi Visual Forest Plot",
          
          radioButtons(
            "forest_style",
            "Format Layout Publikasi:",
            choices = c(
              "Gaya Standar (meta::forest)" = "meta_style",
              "Gaya Klasik (metafor::forest)" = "metafor_style"
            ),
            selected = "meta_style"
          ),
          
          selectInput(
            "forest_sort",
            "Pengurutan Studi:",
            choices = c(
              "Sesuai Urutan Berkas Data" = "none",
              "Besaran Efek (Menaik / Ascending)" = "yi_asc",
              "Besaran Efek (Menurun / Descending)" = "yi_desc",
              "Tahun Publikasi" = "year",
              "Tipe Kelompok Kontrol (Control_Type)" = "control",
              "Format Penghantaran (Format)" = "format",
              "Pendekatan Terapi (Intervention)" = "interv",
              "Risiko Bias (Risk_of_Bias)" = "rob"
            ),
            selected = "none"
          ),
          
          checkboxInput("forest_show_pi", "Tampilkan Prediction Interval (Bar Merah)", value = TRUE),
          checkboxInput("forest_show_weights", "Tampilkan Kolom Bobot Studi (Weight %)", value = TRUE),
          
          conditionalPanel(
            condition = "input.forest_style == 'metafor_style'",
            textInput("forest_label_left", "Label Sumbu Kiri (Kontrol Unggul):", value = "Control Favorable"),
            textInput("forest_label_right", "Label Sumbu Kanan (Intervensi Unggul):", value = "Intervention Favorable")
          ),
          
          sliderInput("forest_cex", "Skala Ukuran Teks / Titik:", min = 0.75, max = 1.30, value = 0.95, step = 0.05),
          
          hr(),
          
          downloadButton("btn_dl_forest_png", "Unduh PNG", class = "btn-primary w-100 mb-2"),
          downloadButton("btn_dl_forest_pdf", "Unduh PDF", class = "btn-secondary w-100")
        ),
        
        card(
          card_header(
            class = "bg-light fw-bold text-primary",
            icon("tree"), " FOREST PLOT "
          ),
          card_body(
            plotOutput("forest_plot_display", height = "680px")
          )
        )
      )
    ),
    
    # Tab 7: Bias Diagnostic
    nav_panel(
      value = "tab_bias",
      icon = icon("shield-alt"),
      title = "Bias Diagnostic",
      layout_sidebar(
        sidebar = sidebar(
          width = 330,
          title = "Pengaturan Plot Corong & Bias",
          
          radioButtons(
            "funnel_type",
            "Jenis Funnel Plot:",
            choices = c(
              "Contour-Enhanced (Zona Signifikansi)" = "contour",
              "Standard Funnel Plot (Pseudo 95% CI)" = "standard",
              "Trim and Fill Funnel Plot (Studi Imputasi)" = "trimfill"
            ),
            selected = "contour"
          ),
          
          numericInput(
            "bias_k_min",
            "Paksa Minimal Studi (k.min Egger & Fail-Safe N):",
            value = 5,
            min = 2,
            max = 10,
            step = 1
          ),
          
          numericInput(
            "fsn_orwin_target",
            "Target Besaran Efek Orwin Fail-Safe N (g):",
            value = 0.20,
            min = 0.05,
            max = 1.00,
            step = 0.05
          ),
          
          div(
            class = "small text-muted mb-3",
            tags$em("Contour-enhanced funnel plot"), " menampilkan zona signifikansi statistik (p < 0.10, p < 0.05, p < 0.01). Parameter ", tags$code("k.min = 5"), " memaksa komputasi uji asimetri Egger dan indeks ", tags$em("Fail-Safe N"), " (Rosenthal, Orwin, Rosenberg) tetap menghasilkan nilai meskipun jumlah studi terbatas (k >= 5)."
          ),
          
          hr(),
          
          downloadButton("btn_dl_funnel_png", "Unduh Funnel Plot PNG", class = "btn-primary w-100 mb-2"),
          downloadButton("btn_dl_funnel_pdf", "Unduh Funnel Plot PDF", class = "btn-secondary w-100")
        ),
        
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header(class = "bg-light fw-bold text-primary", icon("filter"), " FUNNEL PLOT"),
            card_body(
              plotOutput("funnel_plot_display", height = "560px")
            )
          ),
          card(
            card_header(class = "bg-light fw-bold text-primary", icon("shield-alt"), " UJI ASIMETRI & KETAHANAN PUBLIKASI"),
            card_body(
              h6("1. Uji Regresi Linear Egger (Egger's Test, k.min = 5):", class = "fw-bold text-dark"),
              tableOutput("table_egger_test"),
              
              hr(class = "my-2"),
              
              h6("2. Uji Korelasi Peringkat Begg & Mazumdar:", class = "fw-bold text-dark"),
              tableOutput("table_begg_test"),
              
              hr(class = "my-2"),
              
              h6("3. Analisis Imputasi Trim and Fill (Duval & Tweedie):", class = "fw-bold text-dark"),
              tableOutput("table_trimfill_summary"),
              
              hr(class = "my-2"),
              
              h6("4. Indeks Ketahanan Bukti Fail-Safe N (Rosenthal, Orwin & Rosenberg — Paksa k.min = 5):", class = "fw-bold text-dark"),
              tableOutput("table_failsafe_n"),
              
              div(
                class = "alert alert-light border mt-2 small",
                "Metode ", tags$em("Trim and Fill"), " mengestimasi jumlah studi yang hilang akibat bias publikasi, sedangkan ", tags$em("Fail-Safe N"), " (Rosenthal, Orwin, Rosenberg) mengevaluasi ketahanan sintesis bukti terhadap efek laci arsip (", tags$em("file-drawer problem"), ")."
              )
            )
          )
        )
      )
    ),
    
    # Tab 8: Sub Group & Multi Level Analysis
    nav_panel(
      value = "tab_subgroup",
      icon = icon("object-group"),
      title = "Sub Group & Multi Level Analysis",
      layout_sidebar(
        sidebar = sidebar(
          width = 330,
          title = "Konfigurasi Subkelompok",
          
          radioButtons(
            "subgroup_forest_style",
            "Format Layout Subgroup Forest:",
            choices = c(
              "Gaya Standar meta::forest (Author, Year, g, SE + Red PI)" = "meta_style",
              "Gaya Klasik metafor::forest (KAWE Academic Press)" = "metafor_style"
            ),
            selected = "meta_style"
          ),
          
          selectInput(
            "subgroup_var_select",
            "Variabel Subkelompok Klinis:",
            choices = c(
              "Tipe Kelompok Kontrol (Control_Type)" = "Control_Type",
              "Format Penghantaran (Format)" = "Format",
              "Pendekatan Terapi (Intervention)" = "Intervention",
              "Risiko Bias (Risk_of_Bias)" = "Risk_of_Bias",
              "Populasi Sasaran (Population)" = "Population",
              "Kualifikasi Terapis (Therapist_Level)" = "Therapist_Level",
              "Kelompok Usia (Age_Group)" = "Age_Group"
            ),
            selected = "Control_Type"
          ),
          
          sliderInput("subgroup_forest_cex", "Skala Ukuran Teks / Titik:", min = 0.7, max = 1.3, value = 0.90, step = 0.05),
          
          hr(),
          
          downloadButton("btn_dl_subgroup_png", "Unduh PNG Subgroup", class = "btn-primary w-100 mb-2"),
          downloadButton("btn_dl_subgroup_pdf", "Unduh PDF Subgroup", class = "btn-secondary w-100 mb-2"),
          downloadButton("btn_dl_subgroup_csv", "Unduh Tabel Komparasi (.csv)", class = "btn-outline-primary w-100")
        ),
        
        card(
          card_header(
            class = "bg-light fw-bold text-primary",
            icon("object-group"), " Sintesis Bukti Komparatif Antar-Subkelompok"
          ),
          card_body(
            layout_columns(
              col_widths = c(3, 3, 3, 3),
              value_box(
                title = "Uji Beda Subkelompok (QB)",
                value = textOutput("vb_sub_qb"),
                showcase = icon("balance-scale"),
                theme = "primary"
              ),
              value_box(
                title = "Heterogenitas Residual (QE)",
                value = textOutput("vb_sub_qe"),
                showcase = icon("wave-square"),
                theme = "info"
              ),
              value_box(
                title = "Varians Dijelaskan (R^2)",
                value = textOutput("vb_sub_r2"),
                showcase = icon("chart-pie"),
                theme = "warning"
              ),
              value_box(
                title = "Jumlah Subkelompok",
                value = textOutput("vb_sub_k_cats"),
                showcase = icon("layer-group"),
                theme = "secondary"
              )
            ),
            
            hr(class = "my-3"),
            
            navset_card_tab(
              nav_panel(
                title = "Subgroup Forest Plot",
                plotOutput("subgroup_forest_plot_display", height = "720px")
              ),
              nav_panel(
                title = "Tabel Komparasi Subkelompok",
                tableOutput("table_subgroup_comprehensive")
              )
            )
          )
        )
      )
    )
  ),
  
  # Tab 9: Pedoman Metodologis
  nav_panel(
    value = "tab_guide",
    icon = icon("book"),
    title = "Pedoman Metodologis",
    div(
      class = "container-fluid py-3",
      style = "max-width: 1400px; margin: 0 auto;",
      
      # Header Card
      card(
        class = "mb-4 border-0 shadow-sm",
        style = "background: linear-gradient(135deg, #FAF8F5 0%, #FFFFFF 100%); border-left: 5px solid #650915 !important;",
        card_body(
          class = "py-3 px-4",
          div(
            class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
            div(
              h4(class = "fw-bold mb-1", style = "color: #650915; font-family: 'Cinzel', serif; letter-spacing: 0.05em;",
                 icon("book-open"), " KOMPENDIUM METODOLOGIS: PSIKOMETRIKA TERAPAN & SINTESIS BUKTI"),
              p(class = "text-muted small mb-0", "Pedoman konseptual, landasan matematis, dan standar baku pelaporan untuk Meta-Analisis Klinis, Confirmatory Factor Analysis (CFA), dan Teori Respons Butir Modern (Rasch/TAM).")
            ),
            div(
              span(class = "badge", style = "background-color: #650915; color: #FAF8F5; font-size: 0.8rem; padding: 0.45rem 0.75rem;", "Gold Standard Psychometrics & EBP")
            )
          )
        )
      ),
      
      # Sub-Tabs Pedoman Metodologis
      navset_tab(
        id = "guide_subtabs",
        
        # ======================================================================
        # SUB-TAB 1: SINTESIS BUKTI KLINIS & META-ANALISIS
        # ======================================================================
        nav_panel(
          title = span(icon("layer-group"), " 1. Sintesis Bukti & Meta-Analisis"),
          value = "subtab_guide_meta",
          div(
            class = "pt-3",
            
            # Bagian 1: Estimasi Besaran Efek (d vs g)
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("balance-scale"), " 1. ESTIMASI BESARAN EFEK TERSTANDARDISASI (COHEN'S d vs HEDGES' g)")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Standardized Mean Difference: Koreksi Bias Sampel Kecil", class = "text-primary fw-bold mb-2"),
              p(
                "Pada uji klinis psikologi dengan ukuran sampel per kelompok kecil (", tags$code("n < 50"), "), ",
                "estimator standar ", tags$em("Cohen's d"), " secara matematis mengalami bias positif (", tags$em("upward positive bias"), "). ",
                "Kondisi ini menyebabkan estimasi efektivitas intervensi tampak lebih tinggi daripada parameter populasi sebenarnya."
              ),
              p(
                "Hedges (1981) memperkenalkan faktor pembobot koreksi ", tags$em("J"), " berbasis fungsi Gamma dan derajat kebebasan (", tags$code("df = n1 + n2 - 2"), "). ",
                "Faktor ini menghasilkan ", tags$em("unbiased standardized mean difference"), " (", tags$em("Hedges' g"), "). ",
                "Berdasarkan pedoman Cochrane Handbook dan PRISMA 2020, pelaporan ", tags$em("Hedges' g"), " merupakan standar emas yang wajib diutamakan pada seluruh sintesis bukti kuantitatif."
              )
            ),
            div(
              h6("Formulasi Matematis & Notasi:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$d = \\frac{\\bar{X}_1 - \\bar{X}_2}{s_{\\text{pooled}}}, \\quad s_{\\text{pooled}} = \\sqrt{\\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}$$"),
                HTML("$$J = 1 - \\frac{3}{4(n_1 + n_2) - 9}, \\quad g = J \\times d, \\quad v_g = J^2 \\left[ \\frac{n_1 + n_2}{n_1 n_2} + \\frac{d^2}{2(n_1 + n_2)} \\right]$$")
              ),
              h6("Pedoman Magnitudo Efek (Cohen, 1988):", class = "fw-bold text-dark mt-2 mb-1"),
              tags$ul(
                class = "small text-muted mb-0",
                tags$li(strong("|g| < 0.20: "), "Efek dapat diabaikan (negligible effect)"),
                tags$li(strong("0.20 <= |g| < 0.50: "), "Besaran efek kecil (small effect)"),
                tags$li(strong("0.50 <= |g| < 0.80: "), "Besaran efek sedang (medium effect)"),
                tags$li(strong("|g| >= 0.80: "), "Besaran efek besar (large clinical effect)")
              )
            )
          )
        )
      ),
      
      # Bagian 2: Spesifikasi Model Sintesis (FE vs RE)
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("cubes"), " 2. SPESIFIKASI MODEL SINTESIS (FIXED-EFFECT vs RANDOM-EFFECTS)")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Model Efek Tetap (Fixed-Effect / FE) vs Efek Acak (Random-Effects / RE)", class = "text-primary fw-bold mb-2"),
              p(
                strong("Model Efek Tetap (Fixed-Effect / FE): "),
                "Mengasumsikan bahwa seluruh studi membagikan satu besaran efek sejati tunggal (", tags$code("theta"), "). ",
                "Seluruh variabilitas yang teramati antarstudi diasumsikan murni akibat kesalahan pengambilan sampel (", tags$em("sampling error"), "). ",
                "Asumsi ini hampir selalu tidak realistis pada uji intervensi klinis psikologis di mana karakteristik subjek, kualifikasi klinisi, dan format terapi bervariasi antarstudi."
              ),
              p(
                strong("Model Efek Acak (Random-Effects / RE): "),
                "Mengasumsikan bahwa besaran efek sejati bervariasi antarpopulasi dan berdistribusi normal di sekitar rata-rata populasi (", tags$code("mu"), ") dengan varians antarstudi sejati ",
                tags$code("tau^2"), ". Inferensi RE memungkinkan generalisasi eksternal melampaui studi-studi primer yang teramati."
              )
            ),
            div(
              h6("Formulasi & Struktur Pembobotan:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$\\text{Model FE: } y_i = \\theta + e_i, \\quad w_i = \\frac{1}{v_i}, \\quad SE(\\hat{\\theta}) = \\sqrt{\\frac{1}{\\sum w_i}}$$"),
                HTML("$$\\text{Model RE: } y_i = \\mu + u_i + e_i, \\quad w_i^* = \\frac{1}{v_i + \\hat{\\tau}^2}, \\quad SE(\\hat{\\mu}) = \\sqrt{\\frac{1}{\\sum w_i^*}}$$")
              ),
              div(
                class = "alert alert-light border small mb-0 py-2",
                strong("Rekomendasi Metodologis: "),
                "Gunakan Model Efek Acak (RE) dengan estimator REML (Restricted Maximum Likelihood) atau Paule-Mandel, disertai koreksi Knapp-Hartung (Hartung-Knapp/Sidik-Jonkman) untuk mengontrol tingkat kesalahan Tipe I bila jumlah studi primer k terbatas."
              )
            )
          )
        )
      ),
      
      # Bagian 3: Heterogenitas & Prediction Interval
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("chart-area"), " 3. EVALUASI HETEROGENITAS & PREDICTION INTERVAL")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Demistifikasi Heterogenitas: tau^2, I^2 & Interval Prediksi (PI)", class = "text-primary fw-bold mb-2"),
              p(
                "Statistik ", tags$em("I^2"), " bukan merupakan ukuran absolut dari besaran varians antarstudi, melainkan rasio variabilitas sejati terhadap variabilitas total. ",
                "Nilai ", tags$em("I^2"), " sangat dipengaruhi oleh presisi sampling studi primer; peningkatan ukuran sampel N secara otomatis menaikkan ", tags$em("I^2"), " kendati dispersi sejati tetap konstan (Borenstein et al., 2017)."
              ),
              p(
                "Oleh karena itu, varians antarstudi sejati ", tags$code("tau^2"), " (atau deviasi baku antarstudi ", tags$code("tau"), ") dan ",
                tags$strong("Interval Prediksi 95% (Prediction Interval / PI)"), " wajib dilaporkan. ",
                "Interval Prediksi menggambarkan rentang di mana besaran efek intervensi yang sesungguhnya diharapkan berada pada populasi klinis baru di masa depan."
              )
            ),
            div(
              h6("Formulasi Statistik Heterogenitas:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$Q = \\sum_{i=1}^k w_i(y_i - \\hat{\\theta}_{\\text{FE}})^2 \\sim \\chi^2_{k-1}, \\quad I^2 = \\max\\left(0, \\frac{Q - (k-1)}{Q}\\right) \\times 100\\%$$"),
                HTML("$$\\text{Prediction Interval 95\\%: } \\hat{\\mu} \\pm t_{k-2, 0.975} \\sqrt{\\hat{\\tau}^2 + SE(\\hat{\\mu})^2}$$")
              ),
              h6("Klasifikasi Heterogenitas Higgins et al. (2003):", class = "fw-bold text-dark mt-2 mb-1"),
              tags$ul(
                class = "small text-muted mb-0",
                tags$li(strong("I^2 < 25%: "), "Heterogenitas rendah (low)"),
                tags$li(strong("25% <= I^2 < 50%: "), "Heterogenitas moderat (moderate)"),
                tags$li(strong("50% <= I^2 < 75%: "), "Heterogenitas substansial (substantial)"),
                tags$li(strong("I^2 >= 75%: "), "Heterogenitas tinggi (considerable heterogeneity)")
              )
            )
          )
        )
      ),
      
      # Bagian 4: Bias Publikasi & Egger
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("shield-alt"), " 4. EVALUASI BIAS PUBLIKASI & ANALISIS SENSITIVITAS")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Evaluasi Bias Publikasi & Efek Studi Kecil (Small-Study Effects)", class = "text-primary fw-bold mb-2"),
              p(
                "Asimetri pada ", tags$em("funnel plot"), " tidak selalu identik dengan bias publikasi; asimetri dapat mencerminkan heterogenitas metodologis atau efek studi kecil (", tags$em("small-study effects"), "). ",
                tags$em("Contour-enhanced funnel plot"), " memfasilitasi pembuktian apakah studi-studi kecil yang hilang terkonsentrasi secara selektif pada area ketidaksignifikanan statistik (", tags$code("p > 0.05"), ") versus area signifikan (Peters et al., 2008)."
              ),
              p(
                "Uji formal Egger mengevaluasi asimetri corong secara matematis melalui model regresi terbobot. ",
                "Metode ", tags$em("Trim and Fill"), " mengestimasi jumlah studi yang hilang serta memodelkan efek gabungan terkoreksi, sementara ",
                tags$em("Fail-Safe N"), " (Rosenthal, Orwin, Rosenberg) mengevaluasi ketahanan terhadap ", tags$em("file-drawer problem"), "."
              )
            ),
            div(
              h6("Formulasi Regresi Uji Asimetri Egger:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$\\frac{y_i}{SE_i} = \\beta_0 + \\beta_1 \\left(\\frac{1}{SE_i}\\right) + \\epsilon_i, \\quad H_0: \\beta_0 = 0$$")
              ),
              h6("Ambang Batas Ketahanan Rosenthal Fail-Safe N:", class = "fw-bold text-dark mt-2 mb-1"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$N_{\\text{fs}} = \\frac{\\left(\\sum z_i\\right)^2}{2.706} - k, \\quad \\text{Kriteria Kuat: } N_{\\text{fs}} > 5k + 10$$")
              ),
              p(class = "small text-muted mb-0", "Bila N_fs melampaui ambang 5k + 10, kesimpulan sintesis bukti terbukti tahan terhadap potensi studi yang tidak terpublikasi.")
            )
          )
        )
      ),
      
      # Bagian 5: Taksonomi Rumus Effect Size (Ditulis Mengalir Kebawah Penuh)
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("square-root-alt"), " 5. TAKSONOMI FORMULASI BESARAN EFEK (EFFECT SIZES)")
        ),
        card_body(
          p(class = "text-muted mb-3",
            "Berikut merupakan kompendium formulasi matematis baku untuk konversi berbagai format data hasil studi primer ke dalam metrik Standardized Mean Difference (SMD: Cohen's d & Hedges' g) dan varians samplingnya (vi), sesuai rujukan baku Doing Meta-Analysis in R (Harrer et al., 2021; Bab 3)."),
          
          # 5.1 Post-Test Means & SDs
          div(
            class = "p-3 mb-3 rounded border",
            style = "background-color: #FAF8F5; border-color: rgba(197, 160, 89, 0.3) !important;",
            layout_columns(
              col_widths = c(6, 6),
              div(
                h6(class = "fw-bold mb-1", style = "color: #650915;", "5.1 Skor Pascates Dua Kelompok Independen (Between-Group Means & SDs)"),
                p(class = "small mb-1", "Diterapkan pada uji acak terkontrol (RCT) standar yang melaporkan rerata (M), deviasi baku (SD), dan ukuran sampel (n) per kelompok pascates (Cohen, 1988; Hedges, 1981):"),
                tags$ul(
                  class = "small text-muted mb-0",
                  tags$li("Membagi selisih rerata terarah dengan deviasi baku gabungan (s_pooled)."),
                  tags$li("Koreksi Hedges J wajib diterapkan untuk menghilangkan bias positif sampel kecil."),
                  tags$li("Varians sampling vg dihitung untuk pembobotan presisi sintesis.")
                )
              ),
              div(
                div(class = "kawe-eq-box my-0",
                  HTML("$$s_{\\text{pooled}} = \\sqrt{\\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}, \\quad d = \\frac{\\bar{X}_1 - \\bar{X}_2}{s_{\\text{pooled}}}$$"),
                  HTML("$$J = 1 - \\frac{3}{4(n_1 + n_2) - 9}, \\quad g = J \\cdot d, \\quad v_g = J^2 \\left[ \\frac{n_1 + n_2}{n_1 n_2} + \\frac{d^2}{2(n_1 + n_2)} \\right]$$")
                )
              )
            )
          ),
          
          # 5.2 Pre-Post Mean Gain
          div(
            class = "p-3 mb-3 rounded border",
            style = "background-color: #FAF8F5; border-color: rgba(197, 160, 89, 0.3) !important;",
            layout_columns(
              col_widths = c(6, 6),
              div(
                h6(class = "fw-bold mb-1", style = "color: #650915;", "5.2 Skor Perubahan Prates-Pascates (Repeated Measures PPCR / Gain Score)"),
                p(class = "small mb-1", "Diterapkan ketika studi melaporkan perubahan skor pre-post (Gain) antara kelompok intervensi dan kontrol (Becker, 1988; Morris, 2008):"),
                tags$ul(
                  class = "small text-muted mb-0",
                  tags$li("Deviasi baku perubahan mengintegrasikan korelasi berpasangan (r)."),
                  tags$li("Asumsi Rosenthal r = 0.50 digunakan bila r tidak dilaporkan dalam artikel primer."),
                  tags$li("Varians memperhitungkan pengurangan sampling error akibat desain pengukuran berulang.")
                )
              ),
              div(
                div(class = "kawe-eq-box my-0",
                  HTML("$$s_{\\text{change}, i} = \\sqrt{s_{pre, i}^2 + s_{post, i}^2 - 2 r s_{pre, i} s_{post, i}}$$"),
                  HTML("$$d_{\\text{gain}} = \\frac{\\Delta\\text{Gain}}{s_{\\text{pooled, change}}}, \\quad v_d = \\frac{2(1 - r_1)}{n_1} + \\frac{2(1 - r_2)}{n_2} + \\frac{d^2}{2(n_1 + n_2)}$$")
                )
              )
            )
          ),
          
          # 5.3 Parametrik t, F, p
          div(
            class = "p-3 mb-3 rounded border",
            style = "background-color: #FAF8F5; border-color: rgba(197, 160, 89, 0.3) !important;",
            layout_columns(
              col_widths = c(6, 6),
              div(
                h6(class = "fw-bold mb-1", style = "color: #650915;", "5.3 Rekonstruksi dari Nilai Statistik Inferensial (Parametrik t, F, p)"),
                p(class = "small mb-1", "Digunakan ketika artikel primer tidak melaporkan mean dan SD secara lengkap, melainkan hanya nilai uji statistik (Dunlap et al., 1996; Lipsey & Wilson, 2001):"),
                tags$ul(
                  class = "small text-muted mb-0",
                  tags$li("Uji-t independen: Konversi langsung berbasis nilai t dan ukuran sampel."),
                  tags$li("Uji-F satu arah (df1 = 1): Memanfaatkan relasi matematis t = sqrt(F)."),
                  tags$li("Nilai-p dua arah: Rekonstruksi nilai t via inversi kuantil Student's t distribution.")
                )
              ),
              div(
                div(class = "kawe-eq-box my-0",
                  HTML("$$\\text{Uji-t: } d = t \\times \\sqrt{\\frac{n_1 + n_2}{n_1 n_2}}, \\quad \\text{Uji-F: } d = \\pm \\sqrt{F \\left(\\frac{n_1 + n_2}{n_1 n_2}\\right)}$$"),
                  HTML("$$\\text{Nilai-p: } t = F_{t, df}^{-1}(1 - p/2), \\quad g = J \\times d$$")
                )
              )
            )
          ),
          
          # 5.4 Biner & Tabel 2x2
          div(
            class = "p-3 mb-3 rounded border",
            style = "background-color: #FAF8F5; border-color: rgba(197, 160, 89, 0.3) !important;",
            layout_columns(
              col_widths = c(6, 6),
              div(
                h6(class = "fw-bold mb-1", style = "color: #650915;", "5.4 Transformasi Luaran Kategorikal Biner ke Skala d (Tabel Kontingensi 2×2 & Odds Ratio)"),
                p(class = "small mb-1", "Mengonversi rasio odds (Odds Ratio) dari data dikotomi kejadian klinis (event / non-event) ke skala continuous standardized mean difference (Chinn, 2000; Hasselblad & Hedges, 1995):"),
                tags$ul(
                  class = "small text-muted mb-0",
                  tags$li("Berdasarkan asumsi variabel laten berdistribusi logistik kontinu."),
                  tags$li("Faktor pengali sqrt(3)/pi ~ 0.5513 menyelaraskan varians logistik ke normal baku."),
                  tags$li("Varians sampling d ditransformasikan secara proporsional dari varians Woolf.")
                )
              ),
              div(
                div(class = "kawe-eq-box my-0",
                  HTML("$$OR = \\frac{a \\cdot d}{b \\cdot c}, \\quad v_{\\ln(OR)} = \\frac{1}{a} + \\frac{1}{b} + \\frac{1}{c} + \\frac{1}{d}$$"),
                  HTML("$$d = \\ln(OR) \\times \\frac{\\sqrt{3}}{\\pi} \\approx 0.5513 \\ln(OR), \\quad v_d = v_{\\ln(OR)} \\times \\frac{3}{\\pi^2}$$")
                )
              )
            )
          ),
          
          # 5.5 Korelasi r & Fisher z
          div(
            class = "p-3 mb-0 rounded border",
            style = "background-color: #FAF8F5; border-color: rgba(197, 160, 89, 0.3) !important;",
            layout_columns(
              col_widths = c(6, 6),
              div(
                h6(class = "fw-bold mb-1", style = "color: #650915;", "5.5 Koefisien Asosiasi & Korelasi (Pearson r & Fisher's z)"),
                p(class = "small mb-1", "Mengintegrasikan korelasi point-biserial atau Pearson r ke dalam skala standardized mean difference (Borenstein et al., 2009; Jacobs & Viechtbauer, 2017):"),
                tags$ul(
                  class = "small text-muted mb-0",
                  tags$li("Korelasi point-biserial r dikonversi ke d dengan formula Jacobs."),
                  tags$li("Transformasi Fisher's z menstabilkan varians sampling independen dari nilai r populasi."),
                  tags$li("Varians sampling Fisher z: vz = 1 / (N - 3).")
                )
              ),
              div(
                div(class = "kawe-eq-box my-0",
                  HTML("$$d = \\frac{2r}{\\sqrt{1 - r^2}}, \\quad J = 1 - \\frac{3}{4N - 9}, \\quad g = J \\times d$$"),
                  HTML("$$z = \\frac{1}{2}\\ln\\left(\\frac{1 + r}{1 - r}\\right) = \\text{arctanh}(r), \\quad v_z = \\frac{1}{N - 3}$$")
                )
              )
            )
          )
        )
      )
    ) # End of pt-3 div in Sub-Tab 1
  ), # End of Sub-Tab 1 (Meta-Analisis)

  # ============================================================================
  # SUB-TAB 2: CONFIRMATORY FACTOR ANALYSIS (CFA & SEM LATEN)
  # ============================================================================
  nav_panel(
    title = span(icon("project-diagram"), " 2. Confirmatory Factor Analysis (CFA)"),
    value = "subtab_guide_cfa",
    div(
      class = "pt-3",
      
      # 2.1 Paradigma Pengukuran Laten & Asumsi
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("sitemap"), " 1. PARADIGMA PENGUKURAN LATEN & ASUMSI DASAR CFA")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Model Reflektif vs Formatif & Independensi Lokal", class = "text-primary fw-bold mb-2"),
              p(
                "Dalam pemodelan faktor konfirmatori (CFA), indikator empiris diasumsikan beroperasi di bawah ", tags$strong("paradigma reflektif"), ". ",
                "Artinya, atribut laten (faktor) dianggap sebagai entitas psikologis yang mendasari munculnya variasi skor pada butir-butir instrumen. ",
                "Arah kausalitas teoretis bergerak dari konstruk laten menuju butir manifestasi (Bollen, 1989)."
              ),
              p(
                tags$strong("Asumsi Independensi Lokal (Local Independence): "),
                "Mensyaratkan bahwa kovarians antar-butir instrumen harus dijelaskan secara tuntas oleh variabel laten yang dimodelkan. ",
                "Setelah mengontrol skor faktor laten, sisa korelasi antara galat pengukuran (", tags$em("residual errors"), ") dari dua butir mana pun harus mendekati nol. ",
                "Pelanggaran asumsi ini (misal akibat kemiripan redaksional) menghasilkan kovarians galat yang memerlukan justifikasi konseptual khusus."
              ),
              h6("Taksonomi Struktur Model Faktor (Brown, 2015):", class = "fw-bold text-dark mt-2 mb-1"),
              tags$ul(
                class = "small text-muted mb-0",
                tags$li(strong("Model Unidimensi: "), "Seluruh butir memuat pada satu faktor tunggal umum."),
                tags$li(strong("Model Multidimensi Berkorelasi: "), "Setiap kelompok butir memuat pada dimensi unik dengan estimasi kovarians bebas antar-faktor."),
                tags$li(strong("Model Hierarkis Orde-2: "), "Faktor umum orde-tinggi menjelaskan matriks korelasi antar-faktor orde-pertama."),
                tags$li(strong("Model Bifaktor: "), "Setiap butir secara simultan memuat pada faktor umum (General) dan satu faktor spesifik (Group/Specific) yang saling ortogonal.")
              )
            ),
            div(
              h6("Formulasi Matematis Model Pengukuran CFA:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$\\mathbf{y} = \\boldsymbol{\\nu} + \\boldsymbol{\\Lambda} \\boldsymbol{\\eta} + \\boldsymbol{\\epsilon}$$"),
                HTML("$$\\boldsymbol{\\Sigma}(\\boldsymbol{\\theta}) = \\boldsymbol{\\Lambda} \\boldsymbol{\\Phi} \\boldsymbol{\\Lambda}^T + \\boldsymbol{\\Theta}_\\epsilon$$")
              ),
              p(class = "small text-muted mb-2",
                "Di mana y adalah vektor indikator teramati, v adalah intersep, Lambda adalah matriks muatan faktor (factor loadings), eta adalah vektor variabel laten, dan epsilon adalah galat pengukuran. Sigma(theta) adalah matriks kovarians model yang diestimasikan."
              ),
              div(
                class = "alert alert-light border small mb-0 py-2",
                strong("Pemberian Skala Laten (Identification Constraint): "),
                "Variabel laten tidak memiliki metrik bawaan. Identifikasi model wajib dilakukan melalui penetapan butir penanda (Marker Variable, lambda_1 = 1) atau pembakuan varians laten (Fixed Factor Variance, phi_11 = 1)."
              )
            )
          )
        )
      ),
      
      # 2.2 Pemilihan Estimator Berdasarkan Karakteristik Skala Data
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("calculator"), " 2. PEMILIHAN ESTIMATOR BERDASARKAN SKALA DATA & DISTRIBUSI")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Estimasi Kontinu vs Kategorikal Ordinal (Finney & DiStefano, 2013)", class = "text-primary fw-bold mb-2"),
              p(
                "Ketepatan estimasi parameter muatan faktor dan galat baku sangat ditentukan oleh kesesuaian antara tipe skala data empiris dan fungsi estimasi (", tags$em("discrepancy function"), ") yang digunakan:"
              ),
              tags$ul(
                class = "small mb-2",
                tags$li(strong("Maximum Likelihood (ML): "), "Mengasumsikan normalitas multivariat pada data berskala kontinu interval/rasio. Rentan menggelembungkan statistik uji chi-kuadrat dan meremehkan galat baku bila asumsi normalitas dilanggar."),
                tags$li(strong("Robust ML (MLM / Satorra-Bentler): "), "Menerapkan koreksi Satorra-Bentler scaled chi-square dan galat baku robust tipe sandwich untuk menangani ketidaknormalan distribusi kontinu."),
                tags$li(strong("Huber-White Sandwich (MLR): "), "Menggunakan estimator kuasi-ML dengan galat baku robust dan penskalaan Yuan-Bentler. Sangat tangguh pada sampel sedang dan mendukung data hilang (FIML)."),
                tags$li(strong("Diagonally Weighted Least Squares (WLSMV): "), "Standar emas untuk data ordinal/kategorikal (skala Likert 3–7 pilihan) yang dianalisis melalui matriks korelasi polikorik (polychoric correlation). Menghindari distorsi atenuasi muatan faktor akibat pelanggaran skala kontinu.")
              )
            ),
            div(
              h6("Matriks Keputusan Operasional Pemilihan Estimator:", class = "fw-bold text-dark mb-2"),
              div(
                class = "table-responsive",
                tags$table(
                  class = "table table-sm table-bordered small mb-2",
                  tags$thead(
                    class = "table-light text-center",
                    tags$tr(tags$th("Estimator"), tags$th("Skala Indikator"), tags$th("Normalitas"), tags$th("Koreksi SE & Chi-Square"))
                  ),
                  tags$tbody(
                    tags$tr(tags$td(strong("ML")), tags$td("Kontinu (>= 7 opsi)"), tags$td("Multivariate Normal"), tags$td("Standar Asimtotik")),
                    tags$tr(tags$td(strong("MLM")), tags$td("Kontinu"), tags$td("Non-Normal (Moderate/Severe)"), tags$td("Satorra-Bentler Scaled")),
                    tags$tr(tags$td(strong("MLR")), tags$td("Kontinu / Kompleks"), tags$td("Non-Normal / Missing Data"), tags$td("Huber-White Sandwich")),
                    tags$tr(tags$td(strong("WLSMV")), tags$td("Ordinal / Likert (2-6 opsi)"), tags$td("Asumsi Laten Polikorik"), tags$td("Mean & Variance Adjusted"))
                  )
                )
              ),
              p(class = "small text-muted mb-0",
                "Rekomendasi Praktis: Pada skala psikologi Likert 4-5 opsi responss, gunakan WLSMV atau MLR guna menjamin kestabilan estimasi parameter."
              )
            )
          )
        )
      ),
      
      # 2.3 Kriteria Kesesuaian Model (Goodness-of-Fit Criteria)
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("check-circle"), " 3. EVALUASI KESESUAIAN MODEL (GOODNESS-OF-FIT CRITERIA)")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Demistifikasi Goodness-of-Fit: Absolut vs Komparatif", class = "text-primary fw-bold mb-2"),
              p(
                "Evaluasi kesesuaian model CFA tidak boleh bertumpu pada satu indeks tunggal. ",
                "Pendekatan multi-indeks standar Hu & Bentler (1999) dan Kline (2016) mengombinasikan uji fit absolut, komparatif, dan berbasis residu:"
              ),
              tags$ul(
                class = "small mb-2",
                tags$li(strong("Uji Chi-Kuadrat (chi^2, df, p-value): "), "Menguji hipotesis nol kesesuaian sempurna antara matriks kovarians sampel dan matriks kovarians model. Catatan Kritis: Chi-kuadrat sangat sensitif terhadap ukuran sampel N; pada sampel N > 200, penolakan H0 (p < 0.05) merupakan keniscayaan matematis akibat tingginya kekuatan uji (statistical power), bukan bukti penolakan substantif model."),
                tags$li(strong("Comparative Fit Index (CFI) & Tucker-Lewis Index (TLI): "), "Mengevaluasi peningkatan kesesuaian model dibandingkan model dasar (baseline/null model). Ambang batas: >= 0.90 (dapat diterima) dan >= 0.95 (kesesuaian sangat baik)."),
                tags$li(strong("Root Mean Square Error of Approximation (RMSEA): "), "Mengukur kesalahan aproksimasi kovarians per derajat kebebasan pada tingkat populasi (dengan penalti kompleksitas parameter). Ambang batas: <= 0.05 (close fit), <= 0.08 (reasonable fit), disertai batas atas 90% CI < 0.10."),
                tags$li(strong("Standardized Root Mean Square Residual (SRMR): "), "Rata-rata selisih terstandardisasi antara korelasi observasi dan korelasi model. Ambang batas: <= 0.08.")
              )
            ),
            div(
              h6("Formulasi Matematis Indeks Fit Utama:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$\\text{CFI} = 1 - \\frac{\\max(0, \\chi_T^2 - df_T)}{\\max(0, \\chi_B^2 - df_B, \\chi_T^2 - df_T)}$$"),
                HTML("$$\\text{RMSEA} = \\sqrt{\\max\\left(0, \\frac{\\chi_T^2 - df_T}{df_T(N - 1)}\\right)}$$"),
                HTML("$$\\text{SRMR} = \\sqrt{\\frac{2 \\sum_{i \\le j} \\left(\\frac{s_{ij} - \\hat{\\sigma}_{ij}}{s_{ii} s_{jj}}\\right)^2}{p(p + 1)}}$$")
              ),
              div(
                class = "alert alert-light border small mb-0 py-2",
                strong("Peringatan Metodologis (Marsh et al., 2004): "),
                "Ambang 'golden rules' Hu & Bentler (1999) tidak boleh dijadikan kriteria mutlak mekanistik tanpa mempertimbangkan kompleksitas model, jumlah butir, dan ukuran sampel empiris."
              )
            )
          )
        )
      ),
      
      # 2.4 Parameter Muatan Faktor, Validitas Konvergen & Reliabilitas
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("shield-alt"), " 4. MUATAN FAKTOR, VALIDITAS KONVERGEN & RELIABILITAS KOMPOSIT")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Muatan Faktor & Demistifikasi Cronbach's Alpha vs McDonald's Omega", class = "text-primary fw-bold mb-2"),
              p(
                tags$strong("Standar Retensi Muatan Faktor (Standardized Factor Loading / lambda): "),
                "Muatan faktor terstandardisasi merefleksikan korelasi antara indikator dan variabel laten. ",
                "Nilai ", tags$code("lambda >= 0.40"), " merupakan batas retensi minimal (Hair et al., 2019), sedangkan ",
                tags$code("lambda >= 0.50"), " (menjelaskan >= 25% varians butir) atau ",
                tags$code("lambda >= 0.70"), " (menjelaskan >= 50% varians butir) mencerminkan kualitas indikator yang sangat baik."
              ),
              p(
                tags$strong("Kelemahan Mendasar Cronbach's Alpha (alpha): "),
                "Cronbach's alpha bertumpu pada asumsi ", tags$em("essential tau-equivalence"), " yang mengasumsikan seluruh butir memiliki muatan faktor yang persis sama (lambda_1 = lambda_2 = ... = lambda_k). ",
                "Pada instrumen psikologi empiris, asumsi ini hampir selalu dilanggar, menyebabkan koefisien alpha meremehkan reliabilitas sejati instrumen (atenuasi reliabilitas; Sijtsma, 2009; McNeish, 2018)."
              ),
              p(
                tags$strong("Keunggulan McDonald's Omega (omega) via Paket misty (Yanagida, 2026): "),
                "McDonald's omega memodelkan reliabilitas berbasis estimasi parameter muatan faktor dan varians galat unik setiap butir (congeneric model), disertai estimasi 95% Confidence Interval untuk pelaporan ilmiah yang presisi."
              )
            ),
            div(
              h6("Formulasi Validitas Konvergen & Reliabilitas:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$\\text{Average Variance Extracted (AVE)} = \\frac{\\sum_{i=1}^k \\lambda_i^2}{\\sum_{i=1}^k \\lambda_i^2 + \\sum_{i=1}^k \\theta_{ii}} \\ge 0.50$$"),
                HTML("$$\\text{McDonald's Omega (\\omega)} = \\frac{\\left(\\sum_{i=1}^k \\lambda_i\\right)^2}{\\left(\\sum_{i=1}^k \\lambda_i\\right)^2 + \\sum_{i=1}^k \\theta_{ii}}$$"),
                HTML("$$\\text{Cronbach's Alpha (\\alpha)} = \\frac{k}{k-1} \\left(1 - \\frac{\\sum s_i^2}{s_{\\text{total}}^2}\\right)$$")
              ),
              p(class = "small text-muted mb-0",
                "Kriteria Fornell-Larcker: Skala menunjukkan validitas konvergen memadai bila nilai AVE >= 0.50, mengindikasikan bahwa varians yang dijelaskan oleh konstruk laten melampaui varians galat pengukuran."
              )
            )
          )
        )
      ),
      
      # 2.5 Protokol Seleksi Butir & Modifikasi Model
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("broom"), " 5. PROTOKOL SELEKSI BUTIR CEPAT & MODIFIKASI MODEL")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Deteksi Redundansi & Protokol Pemangkasan Skala", class = "text-primary fw-bold mb-2"),
              p(
                "Dalam penyusunan instrumen psikometrik ringkas (", tags$em("short-form inventories"), "), seleksi butir harus dilakukan secara terukur dengan memadukan kriteria statistik dan integritas konseptual:"
              ),
              tags$ul(
                class = "small mb-0",
                tags$li(strong("Deteksi Kolinearitas & Redundansi: "), "Pasangan butir dengan korelasi bivariat r > 0.80 mengindikasikan duplikasi konten (content overlap). Salah satu butir perlu dieliminasi untuk mencegah multikolinearitas dan inflasi korelasi galat."),
                tags$li(strong("Eliminasi Muatan Rendah: "), "Butir dengan muatan faktor lambda < 0.40 dievaluasi untuk dipangkas karena lebih banyak memuat varians galat (uniqueness) daripada varians konstruk laten."),
                tags$li(strong("Pembatasan Kapasitas per Dimensi: "), "Mempertahankan 3-5 butir terbaik dengan muatan faktor tertinggi pada setiap dimensi memastikan parsimoni model tanpa mengorbankan stabilitas identifikasi.")
              )
            ),
            div(
              h6("Etika & Kaidah Modification Indices (MI):", class = "fw-bold text-dark mb-2"),
              p(class = "small text-muted mb-2",
                "Modification Indices (MI) memperkirakan penurunan statistik chi-kuadrat bila suatu parameter bebas ditambahkan ke dalam model. Namun demikian, modifikasi model tunduk pada rambu-rambu ketat:"
              ),
              tags$ul(
                class = "small text-muted mb-2",
                tags$li(strong("Pantangan Data Dredging: "), "Dilarang membebaskan kovarians antar-galat (correlated errors) semata-mata karena nilai MI tinggi tanpa justifikasi teoretis yang jelas."),
                tags$li(strong("Justifikasi Kovarians Galat: "), "Pembebasan kovarians theta_ij hanya dapat dibenarkan apabila terdapat kesamaan format butir yang mencolok, kata kunci identik, atau efek pembalikan skor (reverse-coded items)."),
                tags$li(strong("Kebutuhan Replikasi: "), "Model hasil modifikasi wajib diuji ulang pada sampel validasi independen (cross-validation) guna membuktikan bahwa perbaikan fit bukan murni artefak idiosyncrasy sampel.")
              )
            )
          )
        )
      )
    ) # End of pt-3 div in Sub-Tab 2
  ), # End of Sub-Tab 2 (CFA)

  # ============================================================================
  # SUB-TAB 3: TEORI RESPONS BUTIR MODERN & MODEL RASCH (TAM)
  # ============================================================================
  nav_panel(
    title = span(icon("chart-bar"), " 3. Pemodelan Rasch & IRT Modern"),
    value = "subtab_guide_rasch",
    div(
      class = "pt-3",
      
      # 3.1 Paradigma Pengukuran Rasch vs CTT
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("compass"), " 1. PARADIGMA PENGUKURAN RASCH vs TEORI TES KLASIK (CTT)")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Transformasi Ordinal ke Skala Interval Invarian (Rasch, 1960)", class = "text-primary fw-bold mb-2"),
              p(
                "Teori Tes Klasik (CTT) memiliki kelemahan mendasar: skor mentah total (", tags$em("raw sum score"), ") berskala ordinal dan memiliki spasi antarkategori yang tidak setara. ",
                "Selain itu, indeks kesulitan butir CTT bergantung pada kemampuan sampel (", tags$em("sample-dependent"), "), sedangkan estimasi kemampuan responden bergantung pada paket tes tertentu (", tags$em("test-dependent"), ")."
              ),
              p(
                tags$strong("Prinsip Objektivitas Khusus (Specific Objectivity): "),
                "Model Rasch mengatasi keterbatasan tersebut dengan memisahkan parameter kemampuan subjek (", tags$code("theta"), ") dan tingkat kesulitan butir (", tags$code("beta"), ") secara independen. ",
                "Melalui transformasi log-odds (logit), model Rasch mengonversi data respons ordinal diskrit ke dalam kontinum linear berskala interval yang invarian lintas kelompok partisipan (Bond & Fox, 2015)."
              )
            ),
            div(
              h6("Formulasi Model Rasch Dikotomi:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$P(X_{ni} = 1 | \\theta_n, \\beta_i) = \\frac{\\exp(\\theta_n - \\beta_i)}{1 + \\exp(\\theta_n - \\beta_i)}$$"),
                HTML("$$\\ln\\left(\\frac{P_{ni1}}{1 - P_{ni1}}\\right) = \\theta_n - \\beta_i$$")
              ),
              h6("Karakteristik Pengukuran Fundamental Rasch:", class = "fw-bold text-dark mt-2 mb-1"),
              tags$ul(
                class = "small text-muted mb-0",
                tags$li(strong("Sifat Konjungtif: "), "Bila theta_n = beta_i, probabilitas responden menjawab benar atau menyetujui butir tepat 0.50."),
                tags$li(strong("Bebas Skala Mentah: "), "Skor mentah merupakan statistik cukup (sufficient statistic) untuk estimasi parameter logit."),
                tags$li(strong("Invariansi Skala: "), "Hierarki kesulitan butir tetap stabil terlepas dari tingkat kemampuan kelompok responden yang diuji.")
              )
            )
          )
        )
      ),
      
      # 3.2 Taksonomi Model IRT Politomi
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("layer-group"), " 2. TAKSONOMI MODEL IRT POLITOMI (RSM, PCM, MRSM, MPCM)")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Perbandingan RSM vs PCM: Seragam vs Unik (Masters, 1982)", class = "text-primary fw-bold mb-2"),
              p(
                "Pada instrumen skala Likert politomi multi-poin, pemodelan Rasch diekstensikan menjadi dua keluarga model utama berbasis struktur parameter ambang kategori (", tags$em("step thresholds"), "):"
              ),
              tags$ul(
                class = "small mb-2",
                tags$li(strong("Rating Scale Model (RSM - Andrich, 1978): "), "Mengasumsikan bahwa jarak relatif antar-ambang kategori (threshold tau_k) berjarak seragam di seluruh butir instrumen. Model ini memiliki parsimoni tinggi (hemat parameter) dan cocok untuk kuesioner dengan opsi responss identik dan format terstandardisasi ketat."),
                tags$li(strong("Partial Credit Model (PCM - Masters, 1982): "), "Membebaskan asumsi jarak seragam, memungkinkan setiap butir memiliki parameter ambang batas langkah (step thresholds beta_ik) yang unik dan independen. Model ini lebih fleksibel dalam menangkap variasi psikometrik antar-butir."),
                tags$li(strong("Model Multidimensi (MRSM & MPCM via Engine TAM): "), "Memperluas RSM dan PCM ke dalam ruang vektor laten multidimensi (theta_1, theta_2, ..., theta_D). Kovarians antar-dimensi diestimasi secara simultan dengan integrasi kuadratur Gauss-Hermite (Adams et al., 1997; Wu et al., 2012).")
              )
            ),
            div(
              h6("Formulasi Model Politomi (PCM & RSM):", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$\\text{PCM: } P(X_{ni} = x) = \\frac{\\exp\\left(\\sum_{j=0}^x (\\theta_n - \\beta_{ij})\\right)}{\\sum_{k=0}^{m_i} \\exp\\left(\\sum_{j=0}^k (\\theta_n - \\beta_{ij})\\right)}$$"),
                HTML("$$\\text{RSM: } P(X_{ni} = x) = \\frac{\\exp\\left(x(\\theta_n - \\beta_i) - \\sum_{j=0}^x \\tau_j\\right)}{\\sum_{k=0}^m \\exp\\left(k(\\theta_n - \\beta_i) - \\sum_{j=0}^k \\tau_j\\right)}$$")
              ),
              div(
                class = "alert alert-light border small mb-0 py-2",
                strong("Kriteria Seleksi Model (AIC / BIC): "),
                "Pilihlah model dengan nilai Akaike Information Criterion (AIC) dan Bayesian Information Criterion (BIC) terendah. Likelihood Ratio Test (LRT) digunakan untuk menguji signifikansi perbaikan fit antara model tersarang (nested models)."
              )
            )
          )
        )
      ),
      
      # 3.3 Evaluasi Kesesuaian Butir (Item Fit) & Responden (Person Fit)
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("check-double"), " 3. EVALUASI KESESUAIAN BUTIR (ITEM FIT) & RESPONDEN (PERSON FIT)")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Statistik Infit & Outfit Mean Square (MNSQ)", class = "text-primary fw-bold mb-2"),
              p(
                "Kesesuaian data dengan model dievaluasi melalui statistik residu kuadrat rata-rata (", tags$em("Mean Square / MNSQ"), ") terstandardisasi. Nilai ekspektasi model adalah persis ", tags$strong("1.00"), ":"
              ),
              tags$ul(
                class = "small mb-2",
                tags$li(strong("Infit MNSQ (Information-Weighted): "), "Residu dibobot oleh varians informasi butir, menjadikannya sensitif terhadap pola responss yang tidak terduga pada butir yang sesuai dengan tingkat kemampuan responden (respons terarah)."),
                tags$li(strong("Outfit MNSQ (Unweighted): "), "Residu unweighted yang sangat sensitif terhadap respons ekstrem atau anomali acak (outliers) pada butir yang jauh lebih sulit atau jauh lebih mudah daripada kemampuan responden."),
                tags$li(strong("Rentang Kesesuaian Ideal (Wright & Linacre, 1994): "), "Untuk survei psikologi, kuesioner klinis, dan skala kepribadian: rentang kesesuaian adalah ", tags$strong("0.50 hingga 1.50"), ". Untuk tes prestasi presisi tinggi: rentang ideal 0.70 hingga 1.30.")
              ),
              p(class = "small text-muted mb-0",
                "Z-Standardized (ZSTD): Digunakan bersama MNSQ; nilai |ZSTD| > 2.0 mengindikasikan signifikansi statistik ketidaksesuaian pada taraf alpha = 0.05 bila ukuran sampel N memadai."
              )
            ),
            div(
              h6("Diagnostik Distorsi & Pola Respons Responden:", class = "fw-bold text-dark mb-2"),
              div(
                class = "p-3 rounded border mb-2",
                style = "background-color: #FAF8F5;",
                h6(class = "fw-bold small text-danger mb-1", "Underfitting (MNSQ > 1.50):"),
                p(class = "small mb-2 text-muted", "Menunjukkan adanya variasi acak (noise) berlebih, kebingungan responden terhadap stimulus butir, atau adanya multidimensionalitas yang tidak tertangkap oleh model tunggal."),
                h6(class = "fw-bold small text-warning mb-1", "Overfitting (MNSQ < 0.50):"),
                p(class = "small mb-0 text-muted", "Pola respons terlalu dapat diprediksi (Guttman-like pattern), mengindikasikan redundansi butir atau ketergantungan lokal (local item dependency).")
              ),
              div(
                class = "alert alert-light border small mb-0 py-2",
                strong("Analisis Person Fit: "),
                "Responden dengan Person Fit MNSQ tinggi mengindikasikan perilaku respons atipikal (aberrant responses), seperti satisficing, kecerobohan (careless responding), tebakan acak, atau pemalsuan respons (faking)."
              )
            )
          )
        )
      ),
      
      # 3.4 Peta Responden-Butir (Wright Map) & Penargetan
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("map"), " 4. PETA RESPONDEN-BUTIR (WRIGHT MAP) & PENARGETAN INSTRUMEN")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Prinsip Pemetaan Wilson (2005) & Targeting Skala", class = "text-primary fw-bold mb-2"),
              p(
                tags$strong("Peta Responden-Butir (Wright Map): "),
                "Menjajarkan distribusi parameter kemampuan responden (", tags$code("theta"), ") di panel kiri dan tingkat kesulitan butir / batas kategori (", tags$code("beta / Thurstonian thresholds"), ") di panel kanan pada kontinum logit yang sama. ",
                "Visualisasi ini memungkinkan evaluasi langsung terhadap kesesuaian jangkauan pengukuran (", tags$em("targeting"), ") instrumen terhadap populasi yang diukur."
              ),
              tags$ul(
                class = "small text-muted mb-0",
                tags$li(strong("Targeting Sempurna: "), "Rerata sebaran kemampuan responden (theta) berimpit dekat dengan rerata tingkat kesulitan butir (beta = 0 logit)."),
                tags$li(strong("Floor Effect: "), "Instrumen terlalu sulit; sebaran theta responden berada jauh di bawah butir termudah."),
                tags$li(strong("Ceiling Effect: "), "Instrumen terlalu mudah; sebaran theta responden berada jauh di atas butir tersulit."),
                tags$li(strong("Measurement Gaps: "), "Kekosongan butir pada rentang logit tertentu mengindikasikan rendahnya presisi pengukuran pada stratum kemampuan tersebut.")
              )
            ),
            div(
              h6("Fitur Wright Map Gaya Akademis Berkontras Tinggi:", class = "fw-bold text-dark mb-2"),
              tags$ul(
                class = "small text-muted mb-3",
                tags$li(strong("Titik Lokasi Butir (\\delta): "), "Titik oranye tegas berbingkai gelap merefleksikan tingkat kesulitan rata-rata butir."),
                tags$li(strong("Garis Ambang Kategori (T1, T2, T3, ...): "), "Kode warna yang berbeda menunjukkan batas transisi logit di mana responden memiliki peluang 50:50 untuk melangkah ke kategori respons yang lebih tinggi."),
                tags$li(strong("Dukungan Multidimensi: "), "Pada pemodelan multidimensi (MRSM/MPCM), histogram distribusi kemampuan responden ditampilkan per subdimensi secara berdampingan tanpa mengorbankan keterbacaan sumbu butir.")
              ),
              div(
                class = "alert alert-light border small mb-0 py-2",
                strong("Rekomendasi Revisi Butir: "),
                "Butir yang saling bertumpuk pada posisi logit yang persis sama menunjukkan redundansi pengukuran; salah satu butir dapat dieliminasi guna mempersingkat durasi tes."
              )
            )
          )
        )
      ),
      
      # 3.5 Indeks Separasi & Reliabilitas Rasch
      card(
        class = "mb-4 shadow-sm",
        card_header(
          class = "bg-light d-flex align-items-center",
          span(class = "fw-bold text-primary", icon("chart-line"), " 5. INDEKS SEPARASI & RELIABILITAS RASCH")
        ),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            div(
              h6("Person Reliability vs Item Reliability (Linacre, 2012)", class = "text-primary fw-bold mb-2"),
              p(
                "Berbeda dengan CTT yang hanya melaporkan satu koefisien reliabilitas tunggal, pemodelan Rasch membedakan secara tegas antara reliabilitas subjek dan reliabilitas butir instrumen:"
              ),
              tags$ul(
                class = "small mb-0",
                tags$li(strong("Person Separation Reliability (Rp): "), "Mengukur proporsi varians sejati kemampuan responden terhadap varians teramati. Berkorespondensi erat dengan reliabilitas WLE (Warm's Weighted Likelihood Estimates) pada TAM. Ambang batas: >= 0.70 (dapat diterima), >= 0.80 (baik), >= 0.90 (sangat baik untuk keputusan klinis individual)."),
                tags$li(strong("Item Separation Reliability (Ri): "), "Mengukur replikabilitas hierarki kesulitan butir bila instrumen diujikan pada sampel responden lain dengan karakteristik setara. Nilai Ri umumnya tinggi (>= 0.90) pada ukuran sampel yang memadai.")
              )
            ),
            div(
              h6("Formulasi Indeks Separasi & Strata:", class = "fw-bold text-dark mb-2"),
              div(class = "kawe-eq-box mb-2",
                HTML("$$R_p = \\frac{\\text{Var}(\\hat{\\theta}) - \\text{MSE}(\\hat{\\theta})}{\\text{Var}(\\hat{\\theta})}, \\quad G_p = \\sqrt{\\frac{R_p}{1 - R_p}}$$"),
                HTML("$$\\text{Person Strata } (H) = \\frac{4G_p + 1}{3}$$")
              ),
              div(
                class = "alert alert-light border small mb-0 py-2",
                strong("Interpretasi Strata Kemampuan: "),
                "Nilai strata H mengindikasikan berapa banyak tingkatan kemampuan laten (misal: rendah, sedang, tinggi) yang dapat dibedakan secara andal dan bebas galat oleh skala pengukuran."
              )
            )
          )
        )
      )
    ) # End of pt-3 div in Sub-Tab 3
  ) # End of Sub-Tab 3 (Rasch)
) # End of navset_tab (guide_subtabs)
    ) # End of container-fluid div
  ) # End of nav_panel tab_guide
) # End of page_navbar

# Server Logic
server <- function(input, output, session) {

  # ============================================================================
  # 1. CFA SERVER LOGIC (SIMULASI, DESKRIPTIF, ELIMINASI, FIT, PATH)
  # ============================================================================
  
  # Reactive Data Storage for CFA
  cfa_data_store <- reactiveVal(simulate_cfa_dataset(n_obs = 300, n_items = 12, n_factors = 3, factor_cor = 0.35, lambda_min = 0.55, lambda_max = 0.85, n_likert = 5))
  cfa_data_label <- reactiveVal("Simulasi Default CFA (N = 300, k = 12, 3 Dimensi)")
  cfa_active_items <- reactiveVal(paste0("Item.", 1:12))
  
  # Handle CFA Simulation Button
  observeEvent(input$btn_cfa_gen_sim, {
    req(input$cfa_sim_n, input$cfa_sim_k, input$cfa_sim_factors)
    sim_df <- simulate_cfa_dataset(
      n_obs = input$cfa_sim_n,
      n_items = input$cfa_sim_k,
      n_factors = input$cfa_sim_factors,
      factor_cor = input$cfa_sim_cor,
      n_likert = input$cfa_sim_likert
    )
    cfa_data_store(sim_df)
    cfa_active_items(colnames(sim_df))
    cfa_data_label(sprintf("Simulasi CFA Baru (N = %d, k = %d, %d Dimensi, Likert 1-%d)", 
                           input$cfa_sim_n, input$cfa_sim_k, input$cfa_sim_factors, input$cfa_sim_likert))
    showNotification("Data simulasi CFA berhasil dibangkitkan!", type = "message", duration = 4)
  })
  
  # Handle CFA Built-in Benchmark Data (N=999, k=20, 4 Dimensi, Likert 1-7)
  observeEvent(input$cfa_source_type, {
    if (input$cfa_source_type == "builtin") {
      poly_path <- "02_Working_Data/Polytomous_Benchmark_Data.csv"
      if (file.exists(poly_path)) {
        df_poly <- read.csv(poly_path, check.names = FALSE)
      } else {
        set.seed(2026)
        df_poly <- as.data.frame(matrix(sample(1:7, 999 * 20, replace = TRUE), nrow = 999, ncol = 20))
        colnames(df_poly) <- paste0("Item.", 1:20)
      }
      cfa_data_store(df_poly)
      cfa_active_items(colnames(df_poly))
      cfa_data_label(sprintf("Dataset Bawaan Politomi Likert (N = %d, k = %d, 4 Dimensi, Skala 1-7)", nrow(df_poly), ncol(df_poly)))
      
      updateSliderInput(session, "cfa_n_dims_custom", value = 4)
      cfa_syntax_default <- paste(
        "# Model CFA: Model Pengukuran Politomi (4 Dimensi, k = 20)",
        "Faktor_1 =~ Item.1 + Item.2 + Item.3 + Item.4 + Item.5",
        "Faktor_2 =~ Item.6 + Item.7 + Item.8 + Item.9 + Item.10",
        "Faktor_3 =~ Item.11 + Item.12 + Item.13 + Item.14 + Item.15",
        "Faktor_4 =~ Item.16 + Item.17 + Item.18 + Item.19 + Item.20",
        sep = "\n"
      )
      updateTextAreaInput(session, "cfa_custom_syntax_text", value = cfa_syntax_default)
      showNotification("Dataset bawaan Politomi (4 Dimensi Laten, k = 20, N = 999) berhasil dimuat.", type = "message")
    }
  })
  
  # Handle CFA File Upload
  observeEvent(input$cfa_upload_file, {
    req(input$cfa_upload_file)
    fpath <- input$cfa_upload_file$datapath
    fname <- input$cfa_upload_file$name
    ext <- tools::file_ext(fname)
    tryCatch({
      if (ext == "csv") {
        df_up <- read.csv(fpath, check.names = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        df_up <- as.data.frame(readxl::read_excel(fpath))
      } else {
        stop("Format berkas harus .csv atau .xlsx")
      }
      num_cols <- names(df_up)[sapply(df_up, is.numeric)]
      if (length(num_cols) < 3) stop("Berkas harus memuat minimal 3 kolom butir numerik.")
      df_up <- df_up[, num_cols, drop = FALSE]
      cfa_data_store(df_up)
      cfa_active_items(colnames(df_up))
      cfa_data_label(sprintf("Unggahan Mandiri: %s (N = %d, k = %d)", fname, nrow(df_up), ncol(df_up)))
      showNotification(sprintf("Berhasil memuat %d butir dari berkas %s", ncol(df_up), fname), type = "message")
    }, error = function(e) {
      showNotification(paste("Gagal memuat berkas:", e$message), type = "error")
    })
  })
  
  # Dynamic UI for Custom Dimension Inputs
  output$cfa_custom_dim_inputs_ui <- renderUI({
    df <- cfa_data_store()
    items <- colnames(df)
    n_dims <- input$cfa_n_dims_custom
    req(n_dims, length(items) >= n_dims)
    is_builtin_poly <- (input$cfa_source_type == "builtin") || all(paste0("Item.", 1:20) %in% items)
    
    if (is_builtin_poly && n_dims == 4) {
      poly_dims <- list(
        list(name = "Faktor_1", desc = "Dimensi 1 (5 Butir: 1 s/d 5)", items = paste0("Item.", 1:5)),
        list(name = "Faktor_2", desc = "Dimensi 2 (5 Butir: 6 s/d 10)", items = paste0("Item.", 6:10)),
        list(name = "Faktor_3", desc = "Dimensi 3 (5 Butir: 11 s/d 15)", items = paste0("Item.", 11:15)),
        list(name = "Faktor_4", desc = "Dimensi 4 (5 Butir: 16 s/d 20)", items = paste0("Item.", 16:20))
      )
      dim_cards <- lapply(1:4, function(d) {
        div(
          class = "p-2 mb-2 bg-light rounded border",
          textInput(paste0("cfa_dim_name_", d), sprintf("Nama Dimensi %d:", d), value = poly_dims[[d]]$name),
          span(class = "badge bg-secondary mb-1", poly_dims[[d]]$desc),
          selectizeInput(
            paste0("cfa_dim_items_", d),
            sprintf("Pilih Butir Dimensi %d:", d),
            choices = items,
            selected = intersect(poly_dims[[d]]$items, items),
            multiple = TRUE,
            options = list(plugins = list('remove_button'))
          )
        )
      })
      tagList(dim_cards)
    } else {
      items_per_dim <- ceiling(length(items) / n_dims)
      dim_cards <- lapply(1:n_dims, function(d) {
        start_i <- (d - 1) * items_per_dim + 1
        end_i <- min(d * items_per_dim, length(items))
        def_selected <- if (start_i <= length(items)) items[start_i:end_i] else items[1]
        
        div(
          class = "p-2 mb-2 bg-light rounded border",
          textInput(paste0("cfa_dim_name_", d), sprintf("Nama Dimensi %d:", d), value = paste0("Faktor_", d)),
          selectizeInput(
            paste0("cfa_dim_items_", d),
            sprintf("Pilih Butir Dimensi %d:", d),
            choices = items,
            selected = def_selected,
            multiple = TRUE,
            options = list(plugins = list('remove_button'))
          )
        )
      })
      tagList(dim_cards)
    }
  })
  
  # Reactive CFA Syntax Generator
  cfa_syntax_reactive <- reactive({
    df <- cfa_data_store()
    items <- cfa_active_items()
    req(length(items) >= 3)
    
    st_type <- input$cfa_model_structure
    is_builtin_poly <- (input$cfa_source_type == "builtin") || all(paste0("Item.", 1:20) %in% items)
    
    if (st_type == "unidim") {
      model_name <- if (is_builtin_poly) "General" else "F1"
      paste0(model_name, " =~ ", paste(items, collapse = " + "))
    } else if (input$cfa_dim_assign_mode == "syntax") {
      req(input$cfa_custom_syntax_text)
      input$cfa_custom_syntax_text
    } else if (input$cfa_dim_assign_mode == "custom_assign") {
      n_dims <- input$cfa_n_dims_custom
      req(n_dims)
      lines <- c()
      fac_names <- c()
      for (d in 1:n_dims) {
        fn <- input[[paste0("cfa_dim_name_", d)]]
        if (is.null(fn) || fn == "") fn <- paste0("Faktor_", d)
        itms <- input[[paste0("cfa_dim_items_", d)]]
        itms <- intersect(itms, items)
        if (length(itms) >= 2) {
          lines <- c(lines, paste0(fn, " =~ ", paste(itms, collapse = " + ")))
          fac_names <- c(fac_names, fn)
        }
      }
      if (length(lines) == 0) {
        if (is_builtin_poly) {
          d1 <- intersect(paste0("Item.", 1:5), items)
          d2 <- intersect(paste0("Item.", 6:10), items)
          d3 <- intersect(paste0("Item.", 11:15), items)
          d4 <- intersect(paste0("Item.", 16:20), items)
          lines <- c(
            paste0("Faktor_1 =~ ", paste(d1, collapse = " + ")),
            paste0("Faktor_2 =~ ", paste(d2, collapse = " + ")),
            paste0("Faktor_3 =~ ", paste(d3, collapse = " + ")),
            paste0("Faktor_4 =~ ", paste(d4, collapse = " + "))
          )
          fac_names <- c("Faktor_1", "Faktor_2", "Faktor_3", "Faktor_4")
        } else {
          k <- length(items)
          n_d <- min(3, max(2, ceiling(k / 4)))
          spl <- split(items, cut(seq_along(items), n_d, labels = FALSE))
          lines <- sapply(seq_along(spl), function(i) paste0("F", i, " =~ ", paste(spl[[i]], collapse = " + ")))
          fac_names <- paste0("F", seq_along(spl))
        }
      }
      if (st_type == "second_order" && length(fac_names) >= 2) {
        sec_name <- if (is_builtin_poly) "General" else "General"
        lines <- c(lines, paste0(sec_name, " =~ ", paste(fac_names, collapse = " + ")))
      } else if (st_type == "bifactor") {
        gen_name <- if (is_builtin_poly) "General" else "General"
        lines <- c(paste0(gen_name, " =~ ", paste(items, collapse = " + ")), lines)
      }
      paste(lines, collapse = "\n")
    } else {
      # Auto split
      if (is_builtin_poly) {
        d1 <- intersect(paste0("Item.", 1:5), items)
        d2 <- intersect(paste0("Item.", 6:10), items)
        d3 <- intersect(paste0("Item.", 11:15), items)
        d4 <- intersect(paste0("Item.", 16:20), items)
        fac_lines <- c(
          paste0("Faktor_1 =~ ", paste(d1, collapse = " + ")),
          paste0("Faktor_2 =~ ", paste(d2, collapse = " + ")),
          paste0("Faktor_3 =~ ", paste(d3, collapse = " + ")),
          paste0("Faktor_4 =~ ", paste(d4, collapse = " + "))
        )
        if (st_type == "second_order") {
          fac_lines <- c(fac_lines, "General =~ Faktor_1 + Faktor_2 + Faktor_3 + Faktor_4")
        } else if (st_type == "bifactor") {
          fac_lines <- c(paste0("General =~ ", paste(items, collapse = " + ")), fac_lines)
        }
        paste(fac_lines, collapse = "\n")
      } else {
        k <- length(items)
        n_d <- min(4, max(2, ceiling(k / 4)))
        spl <- split(items, cut(seq_along(items), n_d, labels = FALSE))
        fac_lines <- sapply(seq_along(spl), function(i) paste0("F", i, " =~ ", paste(spl[[i]], collapse = " + ")))
        if (st_type == "second_order" && length(spl) >= 3) {
          fac_lines <- c(fac_lines, paste0("General =~ ", paste(paste0("F", seq_along(spl)), collapse = " + ")))
        } else if (st_type == "bifactor") {
          fac_lines <- c(paste0("General =~ ", paste(items, collapse = " + ")), fac_lines)
        }
        paste(fac_lines, collapse = "\n")
      }
    }
  })
  
  output$cfa_syntax_preview <- renderPrint({
    cat(cfa_syntax_reactive())
  })
  
  # Notification when CFA estimation is executed (tanpa berpindah tab otomatis)
  observeEvent(input$btn_run_cfa, {
    showNotification(
      sprintf("Estimasi CFA berhasil dijalankan dengan estimator %s.", input$cfa_estimator),
      type = "message",
      duration = 4
    )
  })
  
  cfa_fit_trigger <- reactiveVal(1)
  
  # Reactive Fitted CFA Model
  cfa_model_fit <- reactive({
    cfa_fit_trigger()
    input$btn_run_cfa
    df <- cfa_data_store()
    syntax <- cfa_syntax_reactive()
    est <- input$cfa_estimator
    st_type <- input$cfa_model_structure
    
    req(nrow(df) >= 30, nchar(syntax) > 5)
    
    is_ortho <- (st_type == "bifactor")
    use_std_lv <- (st_type != "second_order")
    is_cat_est <- est %in% c("WLSMV", "WLSM", "WLSMVS", "DWLS", "ULSMV", "ULSM", "ULSMVS", "PML")
    
    fit <- tryCatch({
      if (is_cat_est) {
        lavaan::cfa(
          model = syntax,
          data = df,
          estimator = est,
          ordered = TRUE,
          orthogonal = is_ortho,
          std.lv = use_std_lv
        )
      } else {
        lavaan::cfa(
          model = syntax,
          data = df,
          estimator = est,
          orthogonal = is_ortho,
          std.lv = use_std_lv
        )
      }
    }, error = function(e) {
      # Fallback ke skala kontinu jika data bukan kategori diskrit
      tryCatch({
        lavaan::cfa(
          model = syntax,
          data = df,
          estimator = est,
          ordered = FALSE,
          orthogonal = is_ortho,
          std.lv = use_std_lv
        )
      }, error = function(e2) {
        NULL
      })
    })
    fit
  })
  
  # CFA Model Value Boxes
  output$vb_cfa_status <- renderText({
    fit <- cfa_model_fit()
    if (is.null(fit)) return("Belum Diestimasi / Gagal Konvergensi")
    if (lavInspect(fit, "converged")) "Konvergen (Tepat)" else "Tidak Konvergen"
  })
  output$vb_cfa_k <- renderText({
    paste(length(cfa_active_items()), "Butir")
  })
  output$vb_cfa_n <- renderText({
    df <- cfa_data_store()
    paste(nrow(df), "Responden")
  })
  output$vb_cfa_estimator <- renderText({
    input$cfa_estimator
  })
  
  output$cfa_data_badge_ui <- renderUI({
    span(class = "badge bg-primary", cfa_data_label())
  })
  
  output$table_cfa_data_preview <- renderDT({
    df <- cfa_data_store()
    datatable(
      df,
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 20, 25, 50, 100),
        scrollX = TRUE,
        autoWidth = FALSE,
        dom = "lfrtip",
        columnDefs = list(
          list(className = "dt-center", targets = "_all")
        )
      ),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })
  
  # Update CFA Descriptive Item Picker
  observe({
    items <- cfa_active_items()
    updateSelectInput(session, "cfa_desc_item_picker", choices = items, selected = items[1])
  })
  
  # CFA Descriptive Table
  cfa_desc_df_reactive <- reactive({
    df <- cfa_data_store()
    items <- cfa_active_items()
    df_sub <- df[, intersect(items, colnames(df)), drop = FALSE]
    cfa_descriptive_stats(df_sub)
  })
  
  output$table_cfa_descriptive <- renderDT({
    d <- cfa_desc_df_reactive()
    if (nrow(d) == 0) return(NULL)
    datatable(
      d,
      options = list(
        paging = FALSE,
        searching = TRUE,
        info = FALSE,
        scrollX = TRUE
      ),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })
  
  output$vb_cfa_desc_mean <- renderText({
    d <- cfa_desc_df_reactive()
    if (nrow(d) == 0) return("-")
    sprintf("%.3f", mean(d$Mean, na.rm = TRUE))
  })
  output$vb_cfa_desc_skew <- renderText({
    d <- cfa_desc_df_reactive()
    if (nrow(d) == 0) return("-")
    sprintf("%.3f", mean(d$Skewness, na.rm = TRUE))
  })
  output$vb_cfa_desc_kurt <- renderText({
    d <- cfa_desc_df_reactive()
    if (nrow(d) == 0) return("-")
    sprintf("%.3f", mean(d$Kurtosis, na.rm = TRUE))
  })
  output$vb_cfa_desc_norm_pct <- renderText({
    d <- cfa_desc_df_reactive()
    if (nrow(d) == 0) return("-")
    n_norm <- sum(d$Status_Normalitas == "Normal Univariat (Kline, 2016)")
    sprintf("%.1f%% (%d/%d)", (n_norm / nrow(d)) * 100, n_norm, nrow(d))
  })
  
  output$btn_dl_cfa_desc_csv <- downloadHandler(
    filename = function() paste0("CFA_Item_Descriptive_", Sys.Date(), ".csv"),
    content = function(file) write.csv(cfa_desc_df_reactive(), file, row.names = FALSE)
  )
  
  # CFA Item Distribution Plot
  output$cfa_desc_item_plot <- renderPlot({
    req(input$cfa_desc_item_picker)
    df <- cfa_data_store()
    itm <- input$cfa_desc_item_picker
    req(itm %in% colnames(df))
    
    vals <- na.omit(df[[itm]])
    par(mar = c(4.5, 4.5, 2.5, 1.5), bg = "#FFFFFF")
    
    h <- hist(vals, breaks = seq(min(vals) - 0.5, max(vals) + 0.5, by = 1), plot = FALSE)
    max_cnt <- max(h$counts, na.rm = TRUE)
    bp <- barplot(
      h$counts,
      names.arg = h$mids,
      col = "#650915",
      border = "#C5A059",
      xlab = paste("Pilihan Respons Butir:", itm),
      ylab = "Frekuensi Responden",
      main = paste("Sebaran Empiris Respons Butir:", itm),
      las = 1,
      col.main = "#650915",
      ylim = c(0, max_cnt * 1.15)
    )
    grid(nx = NA, ny = NULL, col = "#E5D3B3")
    text(
      x = bp,
      y = h$counts,
      labels = h$counts,
      pos = 3,
      cex = 0.95,
      font = 2,
      col = "#650915"
    )
  })
  
  # --- CFA Correlation Matrix Engine (Pearson & Polychoric) ---
  cfa_calc_matrix_diagnostics <- function(R, n_obs) {
    if (is.null(R) || nrow(R) < 2) return(list(det = NA, kmo = NA, bart_chi = NA, bart_p = NA))
    det_val <- tryCatch(det(R), error = function(e) NA)
    
    kmo_val <- tryCatch({
      res <- psych::KMO(R)
      res$MSA
    }, error = function(e) {
      tryCatch({
        R_sm <- psych::cor.smooth(R)
        psych::KMO(R_sm)$MSA
      }, error = function(e2) NA)
    })
    
    bart_res <- tryCatch({
      res <- psych::cortest.bartlett(R, n = n_obs)
      list(chisq = res$chisq, p.value = res$p.value)
    }, error = function(e) {
      list(chisq = NA, p.value = NA)
    })
    
    list(det = det_val, kmo = kmo_val, bart_chi = bart_res$chisq, bart_p = bart_res$p.value)
  }

  cfa_pearson_cor_reactive <- reactive({
    df <- cfa_data_store()
    items <- cfa_active_items()
    req(df, items)
    df_sub <- df[, intersect(items, colnames(df)), drop = FALSE]
    req(ncol(df_sub) >= 2)
    df_num <- as.data.frame(lapply(df_sub, function(col) as.numeric(as.character(col))))
    cor(df_num, use = "pairwise.complete.obs")
  })

  cfa_poly_cor_reactive <- reactive({
    df <- cfa_data_store()
    items <- cfa_active_items()
    req(df, items)
    df_sub <- df[, intersect(items, colnames(df)), drop = FALSE]
    req(ncol(df_sub) >= 2)
    df_num <- as.data.frame(lapply(df_sub, function(col) as.numeric(as.character(col))))
    tryCatch({
      poly_res <- psych::polychoric(df_num, max.cat = 15, progress = FALSE)
      poly_res$rho
    }, error = function(e) {
      cor(df_num, use = "pairwise.complete.obs")
    })
  })

  cfa_pearson_diag_reactive <- reactive({
    R <- cfa_pearson_cor_reactive()
    df <- cfa_data_store()
    req(R, df)
    cfa_calc_matrix_diagnostics(R, nrow(df))
  })

  cfa_poly_diag_reactive <- reactive({
    rho <- cfa_poly_cor_reactive()
    df <- cfa_data_store()
    req(rho, df)
    cfa_calc_matrix_diagnostics(rho, nrow(df))
  })

  # Pearson Value Boxes
  output$vb_cfa_pearson_det <- renderText({
    d <- cfa_pearson_diag_reactive()
    if (is.na(d$det)) return("-")
    if (abs(d$det) < 0.0001) sprintf("%.3e", d$det) else sprintf("%.5f", d$det)
  })
  output$vb_cfa_pearson_det_sub <- renderText({
    d <- cfa_pearson_diag_reactive()
    if (is.na(d$det)) return("-")
    if (d$det > 0.00001) "Definit Positif (|R| > 0.00001)" else "Peringatan: Potensi Singularitas"
  })

  output$vb_cfa_pearson_kmo <- renderText({
    d <- cfa_pearson_diag_reactive()
    if (is.na(d$kmo)) return("-")
    sprintf("%.3f", d$kmo)
  })
  output$vb_cfa_pearson_kmo_sub <- renderText({
    d <- cfa_pearson_diag_reactive()
    if (is.na(d$kmo)) return("-")
    if (d$kmo >= 0.80) "Sangat Baik (Meritorious)"
    else if (d$kmo >= 0.70) "Layak Difaktorkan (Middling)"
    else if (d$kmo >= 0.60) "Cukup Lemah (Mediocre)"
    else "Tidak Memadai (Unacceptable)"
  })

  output$vb_cfa_pearson_bartlett <- renderText({
    d <- cfa_pearson_diag_reactive()
    if (is.na(d$bart_p)) return("-")
    if (d$bart_p < 0.001) "p < 0.001" else sprintf("p = %.3f", d$bart_p)
  })
  output$vb_cfa_pearson_bartlett_sub <- renderText({
    d <- cfa_pearson_diag_reactive()
    if (is.na(d$bart_chi)) return("-")
    sprintf("Chi-sq = %.2f (Signifikan)", d$bart_chi)
  })

  output$vb_cfa_pearson_mean_r <- renderText({
    R <- cfa_pearson_cor_reactive()
    req(R)
    r_off <- R[lower.tri(R)]
    sprintf("%.3f", mean(r_off, na.rm = TRUE))
  })
  output$vb_cfa_pearson_range_r <- renderText({
    R <- cfa_pearson_cor_reactive()
    req(R)
    r_off <- R[lower.tri(R)]
    sprintf("Rentang: [%.3f, %.3f]", min(r_off, na.rm = TRUE), max(r_off, na.rm = TRUE))
  })

  # Polychoric Value Boxes
  output$vb_cfa_poly_det <- renderText({
    d <- cfa_poly_diag_reactive()
    if (is.na(d$det)) return("-")
    if (abs(d$det) < 0.0001) sprintf("%.3e", d$det) else sprintf("%.5f", d$det)
  })
  output$vb_cfa_poly_det_sub <- renderText({
    d <- cfa_poly_diag_reactive()
    if (is.na(d$det)) return("-")
    if (d$det > 0.00001) "Definit Positif (|\u03c1| > 0.00001)" else "Peringatan: Potensi Singularitas"
  })

  output$vb_cfa_poly_kmo <- renderText({
    d <- cfa_poly_diag_reactive()
    if (is.na(d$kmo)) return("-")
    sprintf("%.3f", d$kmo)
  })
  output$vb_cfa_poly_kmo_sub <- renderText({
    d <- cfa_poly_diag_reactive()
    if (is.na(d$kmo)) return("-")
    if (d$kmo >= 0.80) "Sangat Baik (Meritorious)"
    else if (d$kmo >= 0.70) "Layak Difaktorkan (Middling)"
    else if (d$kmo >= 0.60) "Cukup Lemah (Mediocre)"
    else "Tidak Memadai (Unacceptable)"
  })

  output$vb_cfa_poly_bartlett <- renderText({
    d <- cfa_poly_diag_reactive()
    if (is.na(d$bart_p)) return("-")
    if (d$bart_p < 0.001) "p < 0.001" else sprintf("p = %.3f", d$bart_p)
  })
  output$vb_cfa_poly_bartlett_sub <- renderText({
    d <- cfa_poly_diag_reactive()
    if (is.na(d$bart_chi)) return("-")
    sprintf("Chi-sq = %.2f (Signifikan)", d$bart_chi)
  })

  output$vb_cfa_poly_mean_rho <- renderText({
    rho <- cfa_poly_cor_reactive()
    req(rho)
    rho_off <- rho[lower.tri(rho)]
    sprintf("%.3f", mean(rho_off, na.rm = TRUE))
  })
  output$vb_cfa_poly_range_rho <- renderText({
    rho <- cfa_poly_cor_reactive()
    req(rho)
    rho_off <- rho[lower.tri(rho)]
    sprintf("Rentang: [%.3f, %.3f]", min(rho_off, na.rm = TRUE), max(rho_off, na.rm = TRUE))
  })

  # Pearson Table & Plot
  output$table_cfa_pearson_cor <- renderDT({
    R <- cfa_pearson_cor_reactive()
    req(R)
    d_mat <- cbind(data.frame(Butir = rownames(R)), as.data.frame(round(R, 3)))
    datatable(
      d_mat,
      rownames = FALSE,
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = "_all"),
          list(className = "dt-center font-monospace fw-bold text-primary", targets = 0)
        )
      ),
      class = "display compact cell-border stripe hover"
    )
  })

  output$cfa_pearson_heatmap <- renderPlot({
    R <- cfa_pearson_cor_reactive()
    req(R)
    type_sel <- input$cfa_pearson_plot_type
    if (is.null(type_sel)) type_sel <- "upper"
    show_num <- isTRUE(input$cfa_pearson_show_vals)
    
    col_pal <- colorRampPalette(c("#1D4ED8", "#F8FAFC", "#650915"))(200)
    
    par(mar = c(1, 1, 3, 1), bg = "#FFFFFF")
    corrplot::corrplot(
      R,
      method = "color",
      type = type_sel,
      col = col_pal,
      addCoef.col = if (show_num && ncol(R) <= 30) "#1E293B" else NULL,
      number.cex = if (ncol(R) <= 10) 0.85 else if (ncol(R) <= 18) 0.70 else 0.55,
      tl.col = "#650915",
      tl.srt = 45,
      tl.cex = if (ncol(R) <= 12) 0.90 else 0.75,
      diag = TRUE,
      mar = c(0, 0, 2, 0)
    )
    title(main = "Heatmap Matriks Korelasi Pearson (r)", col.main = "#650915", font.main = 2, cex.main = 1.1)
  })

  # Polychoric Table & Plot
  output$table_cfa_poly_cor <- renderDT({
    rho <- cfa_poly_cor_reactive()
    req(rho)
    d_mat <- cbind(data.frame(Butir = rownames(rho)), as.data.frame(round(rho, 3)))
    datatable(
      d_mat,
      rownames = FALSE,
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = "_all"),
          list(className = "dt-center font-monospace fw-bold text-primary", targets = 0)
        )
      ),
      class = "display compact cell-border stripe hover"
    )
  })

  output$cfa_poly_heatmap <- renderPlot({
    rho <- cfa_poly_cor_reactive()
    req(rho)
    type_sel <- input$cfa_poly_plot_type
    if (is.null(type_sel)) type_sel <- "upper"
    show_num <- isTRUE(input$cfa_poly_show_vals)
    
    col_pal <- colorRampPalette(c("#1D4ED8", "#F8FAFC", "#650915"))(200)
    
    par(mar = c(1, 1, 3, 1), bg = "#FFFFFF")
    corrplot::corrplot(
      rho,
      method = "color",
      type = type_sel,
      col = col_pal,
      addCoef.col = if (show_num && ncol(rho) <= 30) "#1E293B" else NULL,
      number.cex = if (ncol(rho) <= 10) 0.85 else if (ncol(rho) <= 18) 0.70 else 0.55,
      tl.col = "#650915",
      tl.srt = 45,
      tl.cex = if (ncol(rho) <= 12) 0.90 else 0.75,
      diag = TRUE,
      mar = c(0, 0, 2, 0)
    )
    title(main = "Heatmap Matriks Korelasi Polikorik (\u03c1)", col.main = "#650915", font.main = 2, cex.main = 1.1)
  })

  # Delta Attenuation Matrix Reactive, DT, Heatmap & Download
  cfa_delta_matrix_reactive <- reactive({
    R <- cfa_pearson_cor_reactive()
    rho <- cfa_poly_cor_reactive()
    req(R, rho)
    common_items <- intersect(rownames(R), rownames(rho))
    req(length(common_items) >= 2)
    R_sub <- R[common_items, common_items]
    rho_sub <- rho[common_items, common_items]
    
    delta_mat <- rho_sub - R_sub
    diag(delta_mat) <- 0
    delta_mat
  })

  output$table_cfa_delta_matrix <- renderDT({
    delta_mat <- cfa_delta_matrix_reactive()
    req(delta_mat)
    d_fmt <- as.data.frame(apply(delta_mat, 2, function(col) sprintf("%+.3f", col)))
    diag(d_fmt) <- "-"
    d_out <- cbind(data.frame(Butir = rownames(delta_mat)), d_fmt)
    datatable(
      d_out,
      rownames = FALSE,
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = "_all"),
          list(className = "dt-center font-monospace fw-bold text-danger", targets = 0)
        )
      ),
      class = "display compact cell-border stripe hover"
    )
  })

  output$cfa_delta_heatmap <- renderPlot({
    delta_mat <- cfa_delta_matrix_reactive()
    req(delta_mat)
    type_sel <- input$cfa_delta_plot_type
    if (is.null(type_sel)) type_sel <- "upper"
    show_num <- isTRUE(input$cfa_delta_show_vals)
    
    col_pal <- colorRampPalette(c("#1D4ED8", "#FFFFFF", "#650915"))(200)
    
    par(mar = c(1, 1, 3, 1), bg = "#FFFFFF")
    corrplot::corrplot(
      delta_mat,
      is.corr = FALSE,
      method = "color",
      type = type_sel,
      col = col_pal,
      addCoef.col = if (show_num && ncol(delta_mat) <= 30) "#1E293B" else NULL,
      number.cex = if (ncol(delta_mat) <= 10) 0.85 else if (ncol(delta_mat) <= 18) 0.70 else 0.55,
      tl.col = "#650915",
      tl.srt = 45,
      tl.cex = if (ncol(delta_mat) <= 12) 0.90 else 0.75,
      diag = FALSE,
      mar = c(0, 0, 2, 0)
    )
    title(main = expression(paste("Heatmap Matriks Delta Atenuasi (", Delta, " = ", rho, " - r)")), col.main = "#650915", font.main = 2, cex.main = 1.05)
  })

  output$btn_dl_cfa_delta_matrix_csv <- downloadHandler(
    filename = function() paste0("CFA_Delta_Attenuation_Matrix_", Sys.Date(), ".csv"),
    content = function(file) {
      delta_mat <- cfa_delta_matrix_reactive()
      req(delta_mat)
      d_out <- cbind(data.frame(Butir = rownames(delta_mat)), as.data.frame(round(delta_mat, 4)))
      write.csv(d_out, file, row.names = FALSE)
    }
  )

  # Discrepancy (Polychoric vs Pearson) Reactive, DT, Value Box
  cfa_cor_diff_df_reactive <- reactive({
    R <- cfa_pearson_cor_reactive()
    rho <- cfa_poly_cor_reactive()
    req(R, rho)
    common_items <- intersect(rownames(R), rownames(rho))
    req(length(common_items) >= 2)
    R_sub <- R[common_items, common_items]
    rho_sub <- rho[common_items, common_items]
    
    pairs <- which(lower.tri(R_sub), arr.ind = TRUE)
    diff_df <- data.frame(
      Butir_A = colnames(R_sub)[pairs[, 2]],
      Butir_B = rownames(R_sub)[pairs[, 1]],
      r_Pearson = round(R_sub[lower.tri(R_sub)], 3),
      rho_Polikorik = round(rho_sub[lower.tri(rho_sub)], 3),
      Selisih_Atenuasi = round(rho_sub[lower.tri(rho_sub)] - R_sub[lower.tri(R_sub)], 3),
      stringsAsFactors = FALSE
    )
    diff_df <- diff_df[order(abs(diff_df$Selisih_Atenuasi), decreasing = TRUE), ]
    diff_df
  })

  output$vb_cfa_cor_diff_mean <- renderText({
    d <- cfa_cor_diff_df_reactive()
    if (is.null(d) || nrow(d) == 0) return("-")
    sprintf("%+.3f", mean(d$Selisih_Atenuasi, na.rm = TRUE))
  })

  output$vb_cfa_delta_max <- renderText({
    d <- cfa_cor_diff_df_reactive()
    if (is.null(d) || nrow(d) == 0) return("-")
    sprintf("%+.3f", max(d$Selisih_Atenuasi, na.rm = TRUE))
  })

  output$vb_cfa_delta_max_pair <- renderText({
    d <- cfa_cor_diff_df_reactive()
    if (is.null(d) || nrow(d) == 0) return("-")
    max_row <- d[which.max(d$Selisih_Atenuasi), ]
    paste0("Pasangan: ", max_row$Butir_A, " - ", max_row$Butir_B)
  })

  output$vb_cfa_delta_pct_mean <- renderText({
    d <- cfa_cor_diff_df_reactive()
    if (is.null(d) || nrow(d) == 0) return("-")
    pcts <- (d$Selisih_Atenuasi / abs(d$rho_Polikorik)) * 100
    pcts <- pcts[is.finite(pcts)]
    sprintf("%.1f%%", mean(pcts, na.rm = TRUE))
  })

  output$vb_cfa_delta_pair_count <- renderText({
    d <- cfa_cor_diff_df_reactive()
    if (is.null(d) || nrow(d) == 0) return("-")
    n_pos <- sum(d$Selisih_Atenuasi > 0, na.rm = TRUE)
    sprintf("%d/%d Butir (%.1f%% Positif)", n_pos, nrow(d), (n_pos / nrow(d)) * 100)
  })

  output$table_cfa_cor_diff <- renderDT({
    d <- cfa_cor_diff_df_reactive()
    req(d)
    datatable(
      d,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        lengthMenu = list(c(10, 15, 25, 50, -1), c("10", "15", "25", "50", "Semua Pasangan")),
        scrollX = TRUE,
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      class = "display compact cell-border stripe hover"
    )
  })

  # CSV Download Handlers
  output$btn_dl_cfa_pearson_csv <- downloadHandler(
    filename = function() paste0("CFA_Pearson_Correlation_Matrix_", Sys.Date(), ".csv"),
    content = function(file) {
      R <- cfa_pearson_cor_reactive()
      req(R)
      d_out <- cbind(data.frame(Butir = rownames(R)), as.data.frame(round(R, 4)))
      write.csv(d_out, file, row.names = FALSE)
    }
  )

  output$btn_dl_cfa_poly_csv <- downloadHandler(
    filename = function() paste0("CFA_Polychoric_Correlation_Matrix_", Sys.Date(), ".csv"),
    content = function(file) {
      rho <- cfa_poly_cor_reactive()
      req(rho)
      d_out <- cbind(data.frame(Butir = rownames(rho)), as.data.frame(round(rho, 4)))
      write.csv(d_out, file, row.names = FALSE)
    }
  )

  output$btn_dl_cfa_cor_diff_csv <- downloadHandler(
    filename = function() paste0("CFA_Correlation_Attenuation_Discrepancy_", Sys.Date(), ".csv"),
    content = function(file) {
      d_diff <- cfa_cor_diff_df_reactive()
      req(d_diff)
      write.csv(d_diff, file, row.names = FALSE)
    }
  )

  # --- CFA Fast Item Elimination Engine ---
  cfa_elim_results_store <- reactiveVal(NULL)
  
  observeEvent(input$btn_run_elim, {
    df <- cfa_data_store()
    items <- cfa_active_items()
    df_sub <- df[, intersect(items, colnames(df)), drop = FALSE]
    req(ncol(df_sub) >= 4)
    
    cor_thresh <- input$cfa_elim_max_cor
    min_load <- input$cfa_elim_min_load
    apply_cap <- input$cfa_elim_apply_cap
    cap_val <- input$cfa_elim_cap_val
    
    # 1. Correlation Matrix
    cor_mat <- cor(df_sub, use = "pairwise.complete.obs")
    k <- ncol(df_sub)
    is_builtin_poly <- (input$cfa_source_type == "builtin") || all(paste0("Item.", 1:20) %in% items)
    
    # Helper ekstrak struktur faktor dari sintaks CFA aktif agar mendukung model Uni maupun Multidimensi
    syn_active <- tryCatch(isolate(cfa_syntax_reactive()), error = function(e) NULL)
    build_subset_fac_map <- function(syn_str, allowed_items, full_cols) {
      f_map <- list()
      if (is.null(syn_str) || !nzchar(syn_str)) return(f_map)
      lines <- unlist(strsplit(syn_str, "\n"))
      for (ln in lines[grepl("=~", lines)]) {
        parts <- strsplit(ln, "=~")[[1]]
        if (length(parts) != 2) next
        lhs <- trimws(parts[1])
        if (grepl("^#", lhs)) next
        rhs_tokens <- trimws(unlist(strsplit(parts[2], "\\+")))
        rhs_clean <- sub(".*\\*", "", rhs_tokens)
        if (all(rhs_clean %in% full_cols)) {
          kept_rhs <- intersect(rhs_clean, allowed_items)
          if (length(kept_rhs) >= 2) f_map[[lhs]] <- kept_rhs
        }
      }
      if (length(f_map) >= 3) {
        lens <- sapply(f_map, length)
        gen_cands <- names(f_map)[lens == max(lens) & lens > min(lens) * 1.5]
        if (length(gen_cands) == 1 && grepl("general|umum", gen_cands[1], ignore.case = TRUE)) {
          f_map[[gen_cands[1]]] <- NULL
        }
      }
      f_map
    }
    
    pre_fac_map <- build_subset_fac_map(syn_active, colnames(df_sub), colnames(df_sub))
    if (length(pre_fac_map) == 0 && is_builtin_poly) {
      pre_fac_map <- list(
        Faktor_1 = intersect(paste0("Item.", 1:5), colnames(df_sub)),
        Faktor_2 = intersect(paste0("Item.", 6:10), colnames(df_sub)),
        Faktor_3 = intersect(paste0("Item.", 11:15), colnames(df_sub)),
        Faktor_4 = intersect(paste0("Item.", 16:20), colnames(df_sub))
      )
      pre_fac_map <- pre_fac_map[sapply(pre_fac_map, length) >= 2]
    }
    syn_pre <- if (length(pre_fac_map) >= 1) {
      paste(sapply(names(pre_fac_map), function(fn) paste0(fn, " =~ ", paste(pre_fac_map[[fn]], collapse = " + "))), collapse = "\n")
    } else {
      paste0("F1 =~ ", paste(colnames(df_sub), collapse = " + "))
    }
    
    # 2. Estimate factor loadings sesuai spesifikasi dimensi aktif
    loadings_vec <- tryCatch({
      fit_pre <- lavaan::cfa(syn_pre, data = df_sub, std.lv = TRUE)
      std_sol <- lavaan::standardizedSolution(fit_pre)
      loads <- std_sol[std_sol$op == "=~" & std_sol$rhs %in% colnames(df_sub), ]
      setNames(loads$est.std, loads$rhs)
    }, error = function(e) {
      apply(df_sub, 2, function(x) abs(cor(x, rowSums(df_sub) - x)))
    })
    
    # High correlation pairs detection (|r| > cor_thresh)
    high_cor_pairs <- which(abs(cor_mat) > cor_thresh & lower.tri(cor_mat), arr.ind = TRUE)
    dropped_due_cor <- c()
    cor_partners <- list()
    
    if (nrow(high_cor_pairs) > 0) {
      for (idx in seq_len(nrow(high_cor_pairs))) {
        r_item <- rownames(cor_mat)[high_cor_pairs[idx, 1]]
        c_item <- colnames(cor_mat)[high_cor_pairs[idx, 2]]
        
        # Keep track of partner
        cor_partners[[r_item]] <- c(cor_partners[[r_item]], paste0(c_item, " (r=", round(cor_mat[r_item, c_item], 2), ")"))
        cor_partners[[c_item]] <- c(cor_partners[[c_item]], paste0(r_item, " (r=", round(cor_mat[r_item, c_item], 2), ")"))
        
        # Drop the one with lower loading
        l1 <- if (!is.na(loadings_vec[r_item])) loadings_vec[r_item] else 0.5
        l2 <- if (!is.na(loadings_vec[c_item])) loadings_vec[c_item] else 0.5
        drop_target <- if (l1 <= l2) r_item else c_item
        dropped_due_cor <- unique(c(dropped_due_cor, drop_target))
      }
    }
    
    # Evaluate each item
    item_status_list <- lapply(colnames(df_sub), function(it) {
      ld <- if (!is.na(loadings_vec[it])) round(loadings_vec[it], 3) else 0.50
      max_r <- max(abs(cor_mat[it, setdiff(colnames(df_sub), it)]))
      partner_txt <- if (!is.null(cor_partners[[it]])) paste(unique(cor_partners[[it]]), collapse = ", ") else "-"
      
      status <- "Dipertahankan"
      reason <- "Memenuhi kriteria fit & daya diskriminasi"
      
      if (it %in% dropped_due_cor) {
        status <- "Dieliminasi"
        reason <- paste("Redundansi tinggi (r >", cor_thresh, ")")
      } else if (ld < min_load) {
        status <- "Dieliminasi"
        reason <- paste0("Muatan faktor rendah (λ = ", ld, " < ", min_load, ")")
      }
      
      data.frame(
        Butir = it,
        Muatan_Faktor_λ = ld,
        Korelasi_Maks_r = round(max_r, 3),
        Pasangan_Redundan = partner_txt,
        Status = status,
        Alasan_Seleksi = reason,
        stringsAsFactors = FALSE
      )
    })
    
    details_df <- do.call(rbind, item_status_list)
    
    # Apply dimension cap if active (per-dimensi jika model multidimensi)
    if (isTRUE(apply_cap)) {
      if (length(pre_fac_map) >= 2) {
        for (d_name in names(pre_fac_map)) {
          d_items <- intersect(pre_fac_map[[d_name]], details_df$Butir)
          d_kept <- d_items[details_df$Status[details_df$Butir %in% d_items] == "Dipertahankan"]
          if (length(d_kept) > cap_val) {
            sub_ld <- details_df$Muatan_Faktor_λ[details_df$Butir %in% d_kept]
            order_ld <- order(sub_ld, decreasing = TRUE)
            excess <- d_kept[order_ld][(cap_val + 1):length(d_kept)]
            details_df$Status[details_df$Butir %in% excess] <- "Dieliminasi"
            details_df$Alasan_Seleksi[details_df$Butir %in% excess] <- paste0("Melebihi kuota dimensi ", d_name, " (", cap_val, " butir terbaik)")
          }
        }
      } else {
        retained_candidates <- details_df$Butir[details_df$Status == "Dipertahankan"]
        if (length(retained_candidates) > cap_val) {
          order_ld <- order(details_df$Muatan_Faktor_λ[details_df$Butir %in% retained_candidates], decreasing = TRUE)
          excess_items <- retained_candidates[order_ld][(cap_val + 1):length(retained_candidates)]
          details_df$Status[details_df$Butir %in% excess_items] <- "Dieliminasi"
          details_df$Alasan_Seleksi[details_df$Butir %in% excess_items] <- paste("Melebihi kuota maksimal (", cap_val, "butir terbaik)")
        }
      }
    }
    
    retained_final <- details_df$Butir[details_df$Status == "Dipertahankan"]
    if (length(retained_final) < 3) {
      # Safeguard
      retained_final <- details_df$Butir[order(details_df$Muatan_Faktor_λ, decreasing = TRUE)[1:3]]
      details_df$Status[details_df$Butir %in% retained_final] <- "Dipertahankan (Safeguard min 3)"
    }
    
    # Reliability comparison strictly via package misty (Yanagida, 2026): Alpha and Omega only
    a_orig_obj <- tryCatch(misty::item.alpha(df_sub, output = FALSE)$result$alpha, error = function(e) NULL)
    a_sel_obj  <- tryCatch(misty::item.alpha(df_sub[, retained_final, drop = FALSE], output = FALSE)$result$alpha, error = function(e) NULL)
    
    o_orig_obj <- tryCatch(misty::item.omega(df_sub, output = FALSE)$result$omega, error = function(e) NULL)
    o_sel_obj  <- tryCatch(misty::item.omega(df_sub[, retained_final, drop = FALSE], output = FALSE)$result$omega, error = function(e) NULL)
    
    a_orig_val <- if (!is.null(a_orig_obj)) sprintf("%.3f [95%% CI: %.3f, %.3f]", a_orig_obj$alpha, a_orig_obj$low, a_orig_obj$upp) else "-"
    a_sel_val  <- if (!is.null(a_sel_obj)) sprintf("%.3f [95%% CI: %.3f, %.3f]", a_sel_obj$alpha, a_sel_obj$low, a_sel_obj$upp) else "-"
    
    o_orig_val <- if (!is.null(o_orig_obj)) sprintf("%.3f [95%% CI: %.3f, %.3f]", o_orig_obj$omega, o_orig_obj$low, o_orig_obj$upp) else "-"
    o_sel_val  <- if (!is.null(o_sel_obj)) sprintf("%.3f [95%% CI: %.3f, %.3f]", o_sel_obj$omega, o_sel_obj$low, o_sel_obj$upp) else "-"
    
    d_alpha_str <- if (!is.null(a_orig_obj) && !is.null(a_sel_obj)) sprintf("%+.3f", a_sel_obj$alpha - a_orig_obj$alpha) else "-"
    d_omega_str <- if (!is.null(o_orig_obj) && !is.null(o_sel_obj)) sprintf("%+.3f", o_sel_obj$omega - o_orig_obj$omega) else "-"
    
    comp_df <- data.frame(
      Metrik_Reliabilitas = c(
        "Cronbach's Alpha (α) [misty::item.alpha]",
        "McDonald's Omega (ω) [misty::item.omega]",
        "Jumlah Butir Instrumen (k)"
      ),
      Sebelum_Seleksi = c(a_orig_val, o_orig_val, paste(k, "Butir")),
      Setelah_Seleksi = c(a_sel_val, o_sel_val, paste(length(retained_final), "Butir")),
      `Perubahan (\\Delta)` = c(d_alpha_str, d_omega_str, sprintf("%+d Butir", length(retained_final) - k)),
      Evaluasi_Psikometrika = c(
        if (!is.null(a_sel_obj) && a_sel_obj$alpha >= 0.70) "Konsistensi Internal Sangat Baik (>= 0.70)" else "Dapat Diterima / Butuh Penguatan",
        if (!is.null(o_sel_obj) && o_sel_obj$omega >= 0.70) "Reliabilitas Konstruk Presisi (Bebas Asumsi Tau-Ekuivalen)" else "Dapat Diterima",
        "Optimal & Hemat Parameter"
      ),
      stringsAsFactors = FALSE
    )
    
    # Estimasi Perbandingan Model Fit CFA (Sebelum vs Setelah Eliminasi Butir)
    fit_full_elim <- tryCatch({
      lavaan::cfa(syn_pre, data = df_sub, estimator = "MLR", std.lv = TRUE)
    }, error = function(e) NULL)
    
    sel_fac_map <- build_subset_fac_map(syn_pre, retained_final, colnames(df_sub))
    syn_sel <- if (length(sel_fac_map) >= 2) {
      paste(sapply(names(sel_fac_map), function(fn) paste0(fn, " =~ ", paste(sel_fac_map[[fn]], collapse = " + "))), collapse = "\n")
    } else {
      paste0("F1 =~ ", paste(retained_final, collapse = " + "))
    }
    
    fit_sel_elim <- tryCatch({
      lavaan::cfa(syn_sel, data = df_sub[, retained_final, drop = FALSE], estimator = "MLR", std.lv = TRUE)
    }, error = function(e) NULL)
    
    fit_comp_table <- NULL
    delta_fit_metrics <- list(d_cfi = 0, d_tli = 0, d_rmsea = 0, d_aic = 0)
    
    if (!is.null(fit_full_elim) && !is.null(fit_sel_elim)) {
      fm_f <- tryCatch(lavaan::fitMeasures(fit_full_elim), error = function(e) NULL)
      fm_s <- tryCatch(lavaan::fitMeasures(fit_sel_elim), error = function(e) NULL)
      if (!is.null(fm_f) && !is.null(fm_s)) {
        get_val <- function(fm, scaled_key, unscaled_key) {
          if (!is.null(fm[scaled_key]) && !is.na(fm[scaled_key])) fm[scaled_key] else fm[unscaled_key]
        }
        chisq_f <- get_val(fm_f, "chisq.scaled", "chisq")
        df_f <- get_val(fm_f, "df.scaled", "df")
        p_f <- get_val(fm_f, "pvalue.scaled", "pvalue")
        cfi_f <- get_val(fm_f, "cfi.scaled", "cfi")
        tli_f <- get_val(fm_f, "tli.scaled", "tli")
        rmsea_f <- get_val(fm_f, "rmsea.scaled", "rmsea")
        srmr_f <- fm_f["srmr"]
        aic_f <- fm_f["aic"]
        bic_f <- fm_f["bic"]
        
        chisq_s <- get_val(fm_s, "chisq.scaled", "chisq")
        df_s <- get_val(fm_s, "df.scaled", "df")
        p_s <- get_val(fm_s, "pvalue.scaled", "pvalue")
        cfi_s <- get_val(fm_s, "cfi.scaled", "cfi")
        tli_s <- get_val(fm_s, "tli.scaled", "tli")
        rmsea_s <- get_val(fm_s, "rmsea.scaled", "rmsea")
        srmr_s <- fm_s["srmr"]
        aic_s <- fm_s["aic"]
        bic_s <- fm_s["bic"]
        
        delta_fit_metrics <- list(
          d_cfi = cfi_s - cfi_f,
          d_tli = tli_s - tli_f,
          d_rmsea = rmsea_s - rmsea_f,
          d_aic = aic_s - aic_f
        )
        
        fit_comp_table <- data.frame(
          Indeks_Fit = c("Chi-Square (χ²)", "Derajat Bebas (df)", "Nilai p (p-value)", "CFI (Comparative Fit)", "TLI (Tucker-Lewis)", "RMSEA [Error]", "SRMR [Residual]", "AIC (Akaike IC)", "BIC (Bayesian IC)"),
          Model_Awal = c(sprintf("%.3f", chisq_f), sprintf("%d", as.integer(df_f)), sprintf("%.3f", p_f), sprintf("%.3f", cfi_f), sprintf("%.3f", tli_f), sprintf("%.3f", rmsea_f), sprintf("%.3f", srmr_f), sprintf("%.3f", aic_f), sprintf("%.3f", bic_f)),
          Model_Setelah_Eliminasi = c(sprintf("%.3f", chisq_s), sprintf("%d", as.integer(df_s)), sprintf("%.3f", p_s), sprintf("%.3f", cfi_s), sprintf("%.3f", tli_s), sprintf("%.3f", rmsea_s), sprintf("%.3f", srmr_s), sprintf("%.3f", aic_s), sprintf("%.3f", bic_s)),
          `Perubahan (\\Delta)` = c(sprintf("%+.3f", chisq_s - chisq_f), sprintf("%+d", as.integer(df_s - df_f)), sprintf("%+.3f", p_s - p_f), sprintf("%+.3f", cfi_s - cfi_f), sprintf("%+.3f", tli_s - tli_f), sprintf("%+.3f", rmsea_s - rmsea_f), sprintf("%+.3f", srmr_s - srmr_f), sprintf("%+.3f", aic_s - aic_f), sprintf("%+.3f", bic_s - bic_f)),
          Evaluasi_Fit = c(
            "Penurunan deviasi kovarians empiris",
            "Peningkatan derajat kebebasan model",
            "Signifikansi statistik",
            ifelse(cfi_s >= cfi_f, "Peningkatan Kesesuaian Komparatif (Fit Naik)", "Atenuasi minor"),
            ifelse(tli_s >= tli_f, "Peningkatan Kesesuaian Relatif (Fit Naik)", "Atenuasi minor"),
            ifelse(rmsea_s <= rmsea_f, "Penurunan Error Aproksimasi (Optimal)", "Perubahan toleransi error"),
            ifelse(srmr_s <= srmr_f, "Penurunan Residual Standar (Optimal)", "Perubahan residual"),
            "Penurunan AIC (Model Lebih Hemat/Efisien)",
            "Penurunan BIC (Penalti Parameter Efektif)"
          ),
          stringsAsFactors = FALSE
        )
      }
    }
    
    res <- list(
      details = details_df,
      comparison = comp_df,
      fit_comparison = fit_comp_table,
      delta_fit = delta_fit_metrics,
      retained = retained_final,
      dropped = setdiff(details_df$Butir, retained_final),
      cor_mat = cor_mat,
      high_cor_n = nrow(high_cor_pairs)
    )
    cfa_elim_results_store(res)
    # Otomatis mutakhirkan butir aktif dan picu re-estimasi CFA menyeluruh tanpa perlu kembali ke Tab 1
    cfa_active_items(retained_final)
    cfa_fit_trigger(cfa_fit_trigger() + 1)
    showNotification(
      sprintf("Seleksi butir selesai: %d butir dipertahankan (%d dieliminasi). Seluruh tabel CFA, reliabilitas, dan diagram jalur otomatis diperbarui!", 
              length(retained_final), length(details_df$Butir) - length(retained_final)), 
      type = "message", 
      duration = 6
    )
  })
  
  # Value Boxes in Elimination Tab
  output$vb_elim_k_orig <- renderText({
    res <- cfa_elim_results_store()
    if (is.null(res)) paste(length(cfa_active_items()), "Butir") else paste(nrow(res$details), "Butir")
  })
  output$vb_elim_k_dropped <- renderText({
    res <- cfa_elim_results_store()
    if (is.null(res)) "-" else paste(length(res$dropped), "Butir")
  })
  output$vb_elim_k_kept <- renderText({
    res <- cfa_elim_results_store()
    if (is.null(res)) "-" else paste(length(res$retained), "Butir")
  })
  output$vb_elim_high_cor_count <- renderText({
    res <- cfa_elim_results_store()
    if (is.null(res)) "-" else paste(res$high_cor_n, "Pasangan")
  })
  
  output$table_cfa_elim_rel_comparison <- renderTable({
    res <- cfa_elim_results_store()
    if (is.null(res)) return(data.frame(Status = "Silakan klik 'Analisis & Seleksi Butir' untuk melihat komparasi reliabilitas."))
    res$comparison
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$table_cfa_elim_fit_comparison <- renderTable({
    res <- cfa_elim_results_store()
    if (is.null(res) || is.null(res$fit_comparison)) {
      return(data.frame(Status = "Silakan klik 'Analisis & Seleksi Butir' untuk melihat perbandingan indeks fit model sebelum vs sesudah eliminasi."))
    }
    res$fit_comparison
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$vb_elim_delta_cfi <- renderText({
    res <- cfa_elim_results_store()
    if (is.null(res) || is.null(res$delta_fit)) return("-")
    sprintf("%+.3f", res$delta_fit$d_cfi)
  })
  output$vb_elim_delta_tli <- renderText({
    res <- cfa_elim_results_store()
    if (is.null(res) || is.null(res$delta_fit)) return("-")
    sprintf("%+.3f", res$delta_fit$d_tli)
  })
  output$vb_elim_delta_rmsea <- renderText({
    res <- cfa_elim_results_store()
    if (is.null(res) || is.null(res$delta_fit)) return("-")
    sprintf("%+.3f", res$delta_fit$d_rmsea)
  })
  output$vb_elim_delta_aic <- renderText({
    res <- cfa_elim_results_store()
    if (is.null(res) || is.null(res$delta_fit)) return("-")
    sprintf("%+.3f", res$delta_fit$d_aic)
  })
  
  output$table_cfa_elim_item_details <- renderDT({
    res <- cfa_elim_results_store()
    if (is.null(res)) return(NULL)
    datatable(
      res$details,
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE
      ),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })
  
  output$cfa_elim_heatmap_plot <- renderPlot({
    res <- cfa_elim_results_store()
    if (is.null(res)) {
      df <- cfa_data_store()
      items <- cfa_active_items()
      cor_mat <- cor(df[, intersect(items, colnames(df))], use = "pairwise.complete.obs")
    } else {
      cor_mat <- res$cor_mat
    }
    
    par(mar = c(5, 5, 3, 2), bg = "#FFFFFF")
    k <- ncol(cor_mat)
    image(1:k, 1:k, t(cor_mat[k:1, ]), col = colorRampPalette(c("#FFFFFF", "#FAF8F5", "#C5A059", "#650915"))(50),
          axes = FALSE, xlab = "", ylab = "", main = "Inter-Item Correlation Heatmap (Korelasi > 0.80 Ditandai Gelap)")
    axis(1, at = 1:k, labels = colnames(cor_mat), las = 2, cex.axis = 0.8)
    axis(2, at = 1:k, labels = rev(colnames(cor_mat)), las = 2, cex.axis = 0.8)
    # Add values
    for (i in 1:k) {
      for (j in 1:k) {
        val <- cor_mat[i, j]
        is_high <- (abs(val) > input$cfa_elim_max_cor && i != j)
        text(i, k - j + 1, sprintf("%.2f", val), cex = 0.75, font = ifelse(is_high, 2, 1),
             col = ifelse(is_high, "#FFFFFF", ifelse(val > 0.5, "#231618", "#5A4B4C")))
      }
    }
  })
  
  # Apply Eliminated Items to CFA Model
  observeEvent(input$btn_apply_elim, {
    res <- cfa_elim_results_store()
    req(res, length(res$retained) >= 3)
    cfa_active_items(res$retained)
    cfa_fit_trigger(cfa_fit_trigger() + 1)
    showNotification(sprintf("Berhasil! Model CFA kini dikonfigurasi dengan %d butir terpilih tanpa perlu kembali ke Tab 1.", length(res$retained)), type = "message", duration = 5)
  })
  
  # --- CFA Model Fit & Parameter Tables ---
  cfa_fit_indices_table <- reactive({
    fit <- cfa_model_fit()
    req(fit)
    
    fm <- fitMeasures(fit)
    get_fm <- function(sc, un) {
      if (!is.null(fm[sc]) && !is.na(fm[sc])) unname(fm[sc]) else unname(fm[un])
    }
    
    chisq_val <- round(get_fm("chisq.scaled", "chisq"), 3)
    df_val <- as.integer(get_fm("df.scaled", "df"))
    pval_val <- round(get_fm("pvalue.scaled", "pvalue"), 3)
    cfi_val <- round(get_fm("cfi.scaled", "cfi"), 3)
    tli_val <- round(get_fm("tli.scaled", "tli"), 3)
    rmsea_val <- round(get_fm("rmsea.scaled", "rmsea"), 3)
    rmsea_ci <- round(c(get_fm("rmsea.ci.lower.scaled", "rmsea.ci.lower"), get_fm("rmsea.ci.upper.scaled", "rmsea.ci.upper")), 3)
    srmr_val <- round(unname(fm["srmr"]), 3)
    aic_val <- if (!is.na(fm["aic"])) round(unname(fm["aic"]), 3) else NA_real_
    bic_val <- if (!is.na(fm["bic"])) round(unname(fm["bic"]), 3) else NA_real_
    
    data.frame(
      Indeks_Fit = c("Chi-Square (χ²)", "CFI (Comparative Fit Index)", "TLI (Tucker-Lewis Index)", 
                     "RMSEA (Root Mean Square Error)", "SRMR (Standardized RMR)", "AIC", "BIC"),
      Estimasi_Model = c(
        sprintf("%.3f (df = %d, p = %.3f)", chisq_val, df_val, pval_val),
        sprintf("%.3f", cfi_val),
        sprintf("%.3f", tli_val),
        if (!is.na(rmsea_ci[1]) && !is.na(rmsea_ci[2])) sprintf("%.3f [90%% CI: %.3f, %.3f]", rmsea_val, rmsea_ci[1], rmsea_ci[2]) else sprintf("%.3f", rmsea_val),
        sprintf("%.3f", srmr_val),
        if (!is.na(aic_val)) sprintf("%.3f", aic_val) else "N/A (Non-ML)",
        if (!is.na(bic_val)) sprintf("%.3f", bic_val) else "N/A (Non-ML)"
      ),
      Kriteria_Cutoff_Baku = c(
        "p > 0.05 (Sensitif Sampel Besar)",
        ">= 0.95 (Good Fit), >= 0.90 (Acceptable)",
        ">= 0.95 (Good Fit), >= 0.90 (Acceptable)",
        "<= 0.06 (Good Fit), <= 0.08 (Acceptable)",
        "<= 0.08 (Good Fit)",
        "Nilai Lebih Rendah Lebih Baik",
        "Nilai Lebih Rendah Lebih Baik"
      ),
      Evaluasi_Fit = c(
        ifelse(!is.na(pval_val) && pval_val > 0.05, "Exact Fit", "Sensitif terhadap N besar"),
        ifelse(!is.na(cfi_val) && cfi_val >= 0.95, "Good Fit", ifelse(!is.na(cfi_val) && cfi_val >= 0.90, "Acceptable Fit", "Poor Fit")),
        ifelse(!is.na(tli_val) && tli_val >= 0.95, "Good Fit", ifelse(!is.na(tli_val) && tli_val >= 0.90, "Acceptable Fit", "Poor Fit")),
        ifelse(!is.na(rmsea_val) && rmsea_val <= 0.06, "Good Fit", ifelse(!is.na(rmsea_val) && rmsea_val <= 0.08, "Acceptable Fit", "Poor Fit")),
        ifelse(!is.na(srmr_val) && srmr_val <= 0.08, "Good Fit", "Poor Fit"),
        if (!is.na(aic_val)) "Kriteria Informasi (Komparasi Model)" else "Tidak Tersedia (Estimator WLS/DWLS/ULS)",
        if (!is.na(bic_val)) "Kriteria Informasi (Komparasi Model)" else "Tidak Tersedia (Estimator WLS/DWLS/ULS)"
      ),
      stringsAsFactors = FALSE
    )
  })
  
  output$table_cfa_fit_indices <- renderTable({
    t <- tryCatch(cfa_fit_indices_table(), error = function(e) NULL)
    if (is.null(t)) return(data.frame(Status = "Model CFA belum diestimasi atau mengalami singularitas data."))
    t
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$vb_cfa_cfi <- renderText({
    fit <- cfa_model_fit()
    if (is.null(fit)) return("-")
    fm <- fitMeasures(fit)
    val <- if (!is.null(fm["cfi.scaled"]) && !is.na(fm["cfi.scaled"])) fm["cfi.scaled"] else fm["cfi"]
    sprintf("%.3f", unname(val))
  })
  output$vb_cfa_tli <- renderText({
    fit <- cfa_model_fit()
    if (is.null(fit)) return("-")
    fm <- fitMeasures(fit)
    val <- if (!is.null(fm["tli.scaled"]) && !is.na(fm["tli.scaled"])) fm["tli.scaled"] else fm["tli"]
    sprintf("%.3f", unname(val))
  })
  output$vb_cfa_rmsea <- renderText({
    fit <- cfa_model_fit()
    if (is.null(fit)) return("-")
    fm <- fitMeasures(fit)
    r_val <- if (!is.null(fm["rmsea.scaled"]) && !is.na(fm["rmsea.scaled"])) fm["rmsea.scaled"] else fm["rmsea"]
    r_lo <- if (!is.null(fm["rmsea.ci.lower.scaled"]) && !is.na(fm["rmsea.ci.lower.scaled"])) fm["rmsea.ci.lower.scaled"] else fm["rmsea.ci.lower"]
    r_hi <- if (!is.null(fm["rmsea.ci.upper.scaled"]) && !is.na(fm["rmsea.ci.upper.scaled"])) fm["rmsea.ci.upper.scaled"] else fm["rmsea.ci.upper"]
    sprintf("%.3f [%.3f, %.3f]", unname(r_val), unname(r_lo), unname(r_hi))
  })
  output$vb_cfa_srmr <- renderText({
    fit <- cfa_model_fit()
    if (is.null(fit)) return("-")
    sprintf("%.3f", unname(fitMeasures(fit, "srmr")))
  })
  
  output$table_cfa_loadings <- renderDT({
    fit <- cfa_model_fit()
    req(fit)
    std_sol <- lavaan::standardizedSolution(fit)
    loads <- std_sol[std_sol$op == "=~", ]
    
    out <- data.frame(
      Faktor_Laten = loads$lhs,
      Butir_Indikator = loads$rhs,
      Loading_Baku_λ = round(loads$est.std, 3),
      SE = round(loads$se, 3),
      z_value = round(loads$z, 3),
      p_value = round(loads$pvalue, 3),
      Komunalitas_R2 = round(loads$est.std^2, 3),
      Evaluasi = ifelse(loads$est.std >= 0.70, "Ideal (>= 0.70)",
                        ifelse(loads$est.std >= 0.50, "Cukup (>= 0.50)", "Rendah (< 0.50)")),
      stringsAsFactors = FALSE
    )
    datatable(
      out,
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE
      ),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })
  
  output$table_cfa_factor_cors <- renderTable({
    fit <- cfa_model_fit()
    req(fit)
    std_sol <- lavaan::standardizedSolution(fit)
    lv_names <- lavaan::lavNames(fit, "lv")
    cors <- std_sol[std_sol$op == "~~" & std_sol$lhs != std_sol$rhs & std_sol$lhs %in% lv_names & std_sol$rhs %in% lv_names, ]
    if (nrow(cors) == 0) return(data.frame(Keterangan = "Tidak ada korelasi antar-faktor (ortogonal atau model unidimensional)"))
    data.frame(
      Faktor_1 = cors$lhs,
      Faktor_2 = cors$rhs,
      Korelasi_r = round(cors$est.std, 3),
      SE = round(cors$se, 3),
      z_value = round(cors$z, 3),
      p_value = round(cors$pvalue, 3),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Estimasi Reliabilitas Konstruk CFA (misty::item.alpha & misty::item.omega)
  output$table_cfa_fit_reliability <- renderTable({
    df <- cfa_data_store()
    items <- cfa_active_items()
    df_sub <- df[, intersect(items, colnames(df)), drop = FALSE]
    req(ncol(df_sub) >= 2)
    
    calc_misty_row <- function(sub_df, name) {
      if (ncol(sub_df) < 2) {
        return(data.frame(Dimensi_Skala = name, Jumlah_Butir = ncol(sub_df), Cronbach_Alpha = "-", McDonald_Omega = "-", Evaluasi = "Butir < 2", stringsAsFactors = FALSE))
      }
      a <- tryCatch(misty::item.alpha(sub_df, output = FALSE)$result$alpha, error = function(e) NULL)
      o <- tryCatch(misty::item.omega(sub_df, output = FALSE)$result$omega, error = function(e) NULL)
      eval_txt <- if (!is.null(o) && !is.na(o$omega)) {
        if (o$omega >= 0.70) "Reliabilitas Konstruk Sangat Baik (>= 0.70)" else if (o$omega >= 0.60) "Cukup (>= 0.60)" else "Kurang Memadai (< 0.60)"
      } else if (!is.null(a) && !is.na(a$alpha)) {
        if (a$alpha >= 0.70) "Alpha Baik (>= 0.70; Omega butuh >= 3 butir)" else if (a$alpha >= 0.60) "Cukup (>= 0.60)" else "Kurang Memadai (< 0.60)"
      } else {
        "Tidak Dapat Diestimasi"
      }
      data.frame(
        Dimensi_Skala = name,
        Jumlah_Butir = ncol(sub_df),
        Cronbach_Alpha = if (!is.null(a)) sprintf("%.3f [95%% CI: %.3f, %.3f]", a$alpha, a$low, a$upp) else "-",
        McDonald_Omega = if (!is.null(o)) sprintf("%.3f [95%% CI: %.3f, %.3f]", o$omega, o$low, o$upp) else "-",
        Evaluasi = eval_txt,
        stringsAsFactors = FALSE
      )
    }
    
    fac_map <- list()
    fit <- tryCatch(cfa_model_fit(), error = function(e) NULL)
    if (!is.null(fit)) {
      pt <- lavaan::parTable(fit)
      pt_items <- pt[pt$op == "=~" & pt$rhs %in% colnames(df_sub), , drop = FALSE]
      if (nrow(pt_items) > 0) {
        fac_names <- unique(pt_items$lhs)
        for (fn in fac_names) {
          f_its <- unique(pt_items$rhs[pt_items$lhs == fn])
          if (length(f_its) >= 1) fac_map[[fn]] <- f_its
        }
        # Pada Bi-Factor, keluarkan faktor General yang memuat seluruh butir jika ada faktor spesifik
        if (length(fac_map) >= 3) {
          lens <- sapply(fac_map, length)
          gen_cands <- names(fac_map)[lens == max(lens) & lens > min(lens) * 1.5]
          if (length(gen_cands) == 1 && grepl("general|umum", gen_cands[1], ignore.case = TRUE)) {
            fac_map[[gen_cands[1]]] <- NULL
          }
        }
      }
    }
    
    if (length(fac_map) < 2 && ((input$cfa_source_type == "builtin") || all(paste0("Item.", 1:20) %in% colnames(df_sub)))) {
      d1 <- intersect(c("Item.1", "Item.4", "Item.8", "Item.11", "Item.13", "Item.16"), colnames(df_sub))
      d2 <- intersect(c("Item.2", "Item.6", "Item.7", "Item.9", "Item.12", "Item.15", "Item.17"), colnames(df_sub))
      d3 <- intersect(c("Item.3", "Item.5", "Item.10", "Item.14"), colnames(df_sub))
      fac_map <- list()
      if (length(d1) >= 2) fac_map[["D1 (Pembatasan Peran)"]] <- d1
      if (length(d2) >= 2) fac_map[["D2 (Pelepasan Keterikatan)"]] <- d2
      if (length(d3) >= 2) fac_map[["D3 (Ketidaksepadanan Usaha)"]] <- d3
    }
    
    if (length(fac_map) >= 2) {
      rows_list <- lapply(names(fac_map), function(fn) {
        calc_misty_row(df_sub[, intersect(fac_map[[fn]], colnames(df_sub)), drop = FALSE], fn)
      })
      rows_list[[length(rows_list) + 1]] <- calc_misty_row(df_sub, "Skala Keseluruhan (Total Scale)")
      do.call(rbind, rows_list)
    } else {
      calc_misty_row(df_sub, "Skala Instrumen Lengkap")
    }
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # --- CFA Path Diagram Model Engine ---
  cfa_path_model_fit <- reactive({
    cfa_fit_trigger()
    struct_type <- input$cfa_path_model_struct
    df <- cfa_data_store()
    items <- cfa_active_items()
    df_sub <- df[, intersect(items, colnames(df)), drop = FALSE]
    req(ncol(df_sub) >= 3, nrow(df_sub) >= 20)
    
    is_builtin_poly <- (input$cfa_source_type == "builtin") || all(paste0("Item.", 1:20) %in% items)
    est_meth <- isolate(input$cfa_estimator)
    if (is.null(est_meth)) est_meth <- "MLR"
    
    d1_it <- intersect(paste0("Item.", 1:5), colnames(df_sub))
    d2_it <- intersect(paste0("Item.", 6:10), colnames(df_sub))
    d3_it <- intersect(paste0("Item.", 11:15), colnames(df_sub))
    d4_it <- intersect(paste0("Item.", 16:20), colnames(df_sub))
    
    if (!is_builtin_poly || length(d1_it) < 2 || length(d2_it) < 2 || length(d3_it) < 2 || length(d4_it) < 2) {
      k_tot <- ncol(df_sub)
      n_d <- min(4, max(2, ceiling(k_tot / 4)))
      spl <- split(colnames(df_sub), cut(seq_len(k_tot), n_d, labels = FALSE))
      d1_it <- spl[[1]]
      d2_it <- if (n_d >= 2) spl[[2]] else character(0)
      d3_it <- if (n_d >= 3) spl[[3]] else character(0)
      d4_it <- if (n_d >= 4) spl[[4]] else character(0)
    }
    
    syn <- switch(struct_type,
      "unidim" = {
        model_name <- if (is_builtin_poly) "General" else "F1"
        paste0(model_name, " =~ ", paste(colnames(df_sub), collapse = " + "))
      },
      "multi" = {
        lines <- c()
        if (length(d1_it) >= 2) lines <- c(lines, paste0("Faktor_1 =~ ", paste(d1_it, collapse = " + ")))
        if (length(d2_it) >= 2) lines <- c(lines, paste0("Faktor_2 =~ ", paste(d2_it, collapse = " + ")))
        if (length(d3_it) >= 2) lines <- c(lines, paste0("Faktor_3 =~ ", paste(d3_it, collapse = " + ")))
        if (length(d4_it) >= 2) lines <- c(lines, paste0("Faktor_4 =~ ", paste(d4_it, collapse = " + ")))
        if (length(lines) < 2) paste0("F1 =~ ", paste(colnames(df_sub), collapse = " + ")) else paste(lines, collapse = "\n")
      },
      "second_order" = {
        lines <- c()
        fac_names <- c()
        if (length(d1_it) >= 2) { lines <- c(lines, paste0("Faktor_1 =~ ", paste(d1_it, collapse = " + "))); fac_names <- c(fac_names, "Faktor_1") }
        if (length(d2_it) >= 2) { lines <- c(lines, paste0("Faktor_2 =~ ", paste(d2_it, collapse = " + "))); fac_names <- c(fac_names, "Faktor_2") }
        if (length(d3_it) >= 2) { lines <- c(lines, paste0("Faktor_3 =~ ", paste(d3_it, collapse = " + "))); fac_names <- c(fac_names, "Faktor_3") }
        if (length(d4_it) >= 2) { lines <- c(lines, paste0("Faktor_4 =~ ", paste(d4_it, collapse = " + "))); fac_names <- c(fac_names, "Faktor_4") }
        if (length(fac_names) >= 2) {
          sec_name <- if (is_builtin_poly) "General" else "General"
          lines <- c(lines, paste0(sec_name, " =~ ", paste(fac_names, collapse = " + ")))
          paste(lines, collapse = "\n")
        } else {
          paste0("F1 =~ ", paste(colnames(df_sub), collapse = " + "))
        }
      },
      "bifactor" = {
        gen_name <- if (is_builtin_poly) "General" else "General"
        lines <- c(paste0(gen_name, " =~ ", paste(colnames(df_sub), collapse = " + ")))
        fac_names <- c()
        if (length(d1_it) >= 2) { lines <- c(lines, paste0("Faktor_1 =~ ", paste(d1_it, collapse = " + "))); fac_names <- c(fac_names, "Faktor_1") }
        if (length(d2_it) >= 2) { lines <- c(lines, paste0("Faktor_2 =~ ", paste(d2_it, collapse = " + "))); fac_names <- c(fac_names, "Faktor_2") }
        if (length(d3_it) >= 2) { lines <- c(lines, paste0("Faktor_3 =~ ", paste(d3_it, collapse = " + "))); fac_names <- c(fac_names, "Faktor_3") }
        if (length(d4_it) >= 2) { lines <- c(lines, paste0("Faktor_4 =~ ", paste(d4_it, collapse = " + "))); fac_names <- c(fac_names, "Faktor_4") }
        for (fn in fac_names) {
          lines <- c(lines, paste0(gen_name, " ~~ 0*", fn))
        }
        if (length(fac_names) >= 2) {
          for (i in 1:(length(fac_names)-1)) {
            for (j in (i+1):length(fac_names)) {
              lines <- c(lines, paste0(fac_names[i], " ~~ 0*", fac_names[j]))
            }
          }
        }
        paste(lines, collapse = "\n")
      }
    )
    
    is_bi <- (struct_type == "bifactor")
    use_std_lv <- (struct_type != "second_order")
    is_cat_est <- est_meth %in% c("WLSMV", "WLSM", "WLSMVS", "DWLS", "ULSMV", "ULSM", "ULSMVS", "PML")
    fit <- tryCatch({
      if (is_cat_est) {
        lavaan::cfa(syn, data = df_sub, estimator = est_meth, ordered = TRUE, orthogonal = is_bi, std.lv = use_std_lv)
      } else {
        lavaan::cfa(syn, data = df_sub, estimator = est_meth, orthogonal = is_bi, std.lv = use_std_lv)
      }
    }, error = function(e) {
      tryCatch({
        lavaan::cfa(syn, data = df_sub, estimator = est_meth, ordered = FALSE, orthogonal = is_bi, std.lv = use_std_lv)
      }, error = function(e2) {
        tryCatch({
          lavaan::cfa(syn, data = df_sub, estimator = "ML", orthogonal = is_bi, std.lv = use_std_lv)
        }, error = function(e3) NULL)
      })
    })
    
    list(fit = fit, struct = struct_type, n = nrow(df_sub), k = ncol(df_sub), est = est_meth)
  })
  
  output$cfa_path_badges <- renderUI({
    info <- cfa_path_model_fit()
    req(info)
    struct_label <- switch(input$cfa_path_model_struct,
      "unidim" = "Unidimensional",
      "multi" = "Multidimensional",
      "second_order" = "Second-Order",
      "bifactor" = "Bi-Factor",
      "CFA Model"
    )
    style_label <- if (input$cfa_path_style == "kawe") "Academic Kawe Edition" else "Standar Klasik"
    tagList(
      span(class = "badge bg-primary me-1", struct_label),
      span(class = "badge bg-secondary me-1", paste(info$k, "Butir Aktif")),
      span(class = "badge bg-dark me-1", paste0("N = ", info$n)),
      span(class = ifelse(input$cfa_path_style == "kawe", "badge bg-danger", "badge bg-light text-dark border"), style_label)
    )
  })
  
  # Helper to build non-colliding layout according to APA 7th publication standards (Q1/Q2)
  build_cfa_custom_layout <- function(fit, struct_type, n_items, orientation) {
    if (orientation == "tree_h") {
      p <- semPaths(fit, layout = "tree", rotation = 2, thresholds = FALSE, intercepts = FALSE, DoNotPlot = TRUE)
      lay <- p$layout
      n_total_nodes <- nrow(lay)
      n_lat <- n_total_nodes - n_items
      
      lv_names <- fit@pta$vnames$lv[[1]]
      
      if (struct_type == "bifactor") {
        # Tata Letak Khusus Bi-Factor Standar APA 7th (Q1/Q2 Benchmark):
        # 3 Faktor Spesifik di KIRI (-1.0), Butir di TENGAH (0.0), General Factor di KANAN (1.0)
        y_items <- seq(-0.95, 0.95, length.out = n_items)
        for (i in 1:n_items) {
          lay[i, ] <- c(0.0, y_items[i])
        }
        
        gen_idx <- which(grepl("general|umum", lv_names, ignore.case = TRUE))
        if (length(gen_idx) == 0) gen_idx <- 1
        spec_indices <- setdiff(seq_along(lv_names), gen_idx)
        
        # General Factor di Kanan (x = 1.0)
        for (g in gen_idx) {
          lay[n_items + g, ] <- c(1.0, 0.0)
        }
        
        # Faktor-faktor Spesifik di Kiri (x = -1.0)
        y_spec <- seq(-0.75, 0.75, length.out = max(1, length(spec_indices)))
        for (s_i in seq_along(spec_indices)) {
          s_idx <- spec_indices[s_i]
          lay[n_items + s_idx, ] <- c(-1.0, y_spec[s_i])
        }
        
      } else if (struct_type == "second_order") {
        y_items <- seq(-0.95, 0.95, length.out = n_items)
        for (i in 1:n_items) {
          lay[i, ] <- c(1.0, y_items[i])
        }
        if (n_lat >= 2) {
          y_lat <- seq(-0.75, 0.75, length.out = n_lat - 1)
          for (i in 1:(n_lat - 1)) {
            lay[n_items + i, ] <- c(-0.25, y_lat[i])
          }
          lay[n_total_nodes, ] <- c(-1.1, 0.0)
        } else {
          lay[n_total_nodes, ] <- c(-1.0, 0.0)
        }
        
      } else if (struct_type == "multi") {
        # Standar Tata Letak Pohon Horizontal: Memanfaatkan kalkulasi kurva otomatis semPaths
        # sehingga kurva kovarians antar-laten non-berdampingan (F1 <-> F3) melengkung ke luar (curve = 2)
        # dan tidak tertutup oleh simpul laten di tengah (F2).
        return(list(layout = "tree", rotation = 2))
        
      } else { # unidim
        return(list(layout = "tree", rotation = 2))
      }
      
      return(list(layout = lay, rotation = 2))
    } else {
      rot <- 1
      lay_name <- if (orientation == "tree_v") "tree" else orientation
      return(list(layout = lay_name, rotation = rot))
    }
  }
  
  # --- CFA Path Diagram Rendering (Standar Publikasi APA 7th Q1/Q2) ---
  render_cfa_path_internal <- function() {
    model_obj <- cfa_path_model_fit()
    if (is.null(model_obj) || is.null(model_obj$fit)) {
      plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0, 0, "Model CFA belum dapat diestimasi atau belum konvergen.\nSilakan periksa konfigurasi data atau gunakan struktur model lain.", cex = 1.1, col = "#650915", font = 2)
      return()
    }
    fit <- model_obj$fit
    
    lay_pref <- input$cfa_path_layout
    style <- input$cfa_path_style
    st_type <- input$cfa_path_model_struct
    
    struct_title <- switch(st_type,
      "unidim" = "Model 1: Unidimensional (1-Faktor)",
      "multi" = "Model 2: Multidimensional (Faktor Berkorelasi)",
      "second_order" = "Model 3: Second-Order (Hierarkis)",
      "bifactor" = "Model 4: Bi-Factor Ortogonal",
      "Model Analisis Faktor Konfirmatori"
    )
    
    ov_names <- fit@pta$vnames$ov[[1]]
    lv_names <- fit@pta$vnames$lv[[1]]
    clean_lv <- gsub("_", "\n", lv_names)
    clean_labels <- c(ov_names, clean_lv)
    
    lay_info <- build_cfa_custom_layout(fit, st_type, length(ov_names), lay_pref)
    
    user_cex <- if (!is.null(input$cfa_path_cex)) input$cfa_path_cex else 0.85
    edge_w <- if (!is.null(input$cfa_path_edge_width)) input$cfa_path_edge_width else 1.2
    show_resid <- isTRUE(input$cfa_path_show_resid)
    
    # Margin dengan headroom atas aman dari tabrakan
    m_top <- 4.5
    m_bot <- 2.5
    m_side <- 4.0
    
    # Ukuran simpul node yang seimbang
    man_w <- if (lay_pref == "tree_h") (if (st_type == "bifactor") 6.4 else 7.0) else 6.0
    man_h <- if (lay_pref == "tree_h") 2.0 else 1.8
    lat_w <- if (lay_pref == "tree_h") 12.5 else 9.0
    lat_h <- if (lay_pref == "tree_h") 5.2 else 4.5
    
    if (style == "kawe") {
      # Academic Kawe Style (Maroon & Gold Edition - Standar APA 7th Q1/Q2)
      semPlot::semPaths(
        fit,
        what = ifelse(isTRUE(input$cfa_path_show_std), "std", "path"),
        whatLabels = ifelse(isTRUE(input$cfa_path_show_std), "std", "name"),
        layout = lay_info$layout,
        rotation = lay_info$rotation,
        weighted = FALSE,             # Standar Baku APA 7th: Garis langsing proporsional, bebas alay
        edge.width = edge_w,          # Garis ramping presisi tinggi
        asize = 1.2,                  # Mata panah tajam dan tegas
        residuals = show_resid,
        thresholds = FALSE,
        intercepts = FALSE,
        residScale = 14,
        nodeLabels = clean_labels,
        nCharNodes = 0,
        nCharEdges = 0,
        sizeMan = man_w,
        sizeMan2 = man_h,
        sizeLat = lat_w,
        sizeLat2 = lat_h,
        color = list(lat = "#FAF8F5", man = "#FFFFFF"),
        edge.color = "#650915",
        border.color = "#C5A059",
        border.width = 1.4,
        label.color = "#231618",
        edge.label.bg = "white",
        edge.label.cex = user_cex,
        mar = c(m_bot, m_side, m_top, m_side),
        curveAdjacent = TRUE
      )
      
      sub_txt <- if (st_type == "bifactor") {
        sprintf("Faktor Spesifik di Kiri | Butir di Tengah | General Factor di Kanan | Estimator: %s (N = %d, k = %d)", model_obj$est, model_obj$n, model_obj$k)
      } else {
        sprintf("Estimator: %s | Sampel N = %d | k = %d Butir Aktif | Standar Pelaporan APA 7th (Q1/Q2)", model_obj$est, model_obj$n, model_obj$k)
      }
      
      mtext(sprintf("%s - Academic Kawe Edition", struct_title), side = 3, line = 2.0, cex = 1.25, font = 2, col = "#650915")
      mtext(sub_txt, side = 3, line = 0.6, cex = 0.95, font = 3, col = "#5A4B4C")
      
    } else {
      # Path Diagram Biasa (Standar Klasik Publikasi Monokrom - APA 7th Q1/Q2)
      semPlot::semPaths(
        fit,
        what = ifelse(isTRUE(input$cfa_path_show_std), "std", "path"),
        whatLabels = ifelse(isTRUE(input$cfa_path_show_std), "std", "name"),
        layout = lay_info$layout,
        rotation = lay_info$rotation,
        weighted = FALSE,             # Standar Baku APA 7th: Garis langsing proporsional, bebas alay
        edge.width = edge_w,          # Garis ramping presisi tinggi
        asize = 1.1,                  # Mata panah klasik tajam
        residuals = show_resid,
        thresholds = FALSE,
        intercepts = FALSE,
        residScale = 14,
        nodeLabels = clean_labels,
        nCharNodes = 0,
        nCharEdges = 0,
        sizeMan = man_w,
        sizeMan2 = man_h,
        sizeLat = lat_w,
        sizeLat2 = lat_h,
        color = list(lat = "#F9F9F9", man = "#FFFFFF"),
        edge.color = "#222222",
        border.color = "#333333",
        border.width = 1.2,
        label.color = "#000000",
        edge.label.bg = "#FFFFFF",
        edge.label.cex = user_cex * 0.95,
        mar = c(m_bot, m_side, m_top, m_side),
        curveAdjacent = TRUE
      )
      
      sub_txt <- if (st_type == "bifactor") {
        sprintf("Specific Factors (Left) | Items (Center) | General Factor (Right) | Estimator: %s | N = %d | k = %d", model_obj$est, model_obj$n, model_obj$k)
      } else {
        sprintf("Estimator: %s | N = %d | k = %d | APA 7th Reporting Standard (Q1/Q2)", model_obj$est, model_obj$n, model_obj$k)
      }
      
      mtext(sprintf("CFA Path Diagram: %s (Standard Classic)", struct_title), side = 3, line = 2.0, cex = 1.20, font = 1, col = "#222222")
      mtext(sub_txt, side = 3, line = 0.6, cex = 0.92, font = 1, col = "#555555")
    }
  }

  output$cfa_path_plot_ui <- renderUI({
    h_px <- if (!is.null(input$cfa_path_height)) input$cfa_path_height else 820
    plotOutput("cfa_path_display", width = "100%", height = paste0(h_px, "px"))
  })
  
  output$cfa_path_display <- renderPlot({
    render_cfa_path_internal()
  }, res = 130)
  
  output$btn_dl_cfa_path_png <- downloadHandler(
    filename = function() {
      paste0("CFA_Path_Diagram_", input$cfa_path_model_struct, "_", input$cfa_path_style, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 3600, height = 2600, res = 300)
      render_cfa_path_internal()
      dev.off()
    }
  )
  output$btn_dl_cfa_path_pdf <- downloadHandler(
    filename = function() {
      paste0("CFA_Path_Diagram_", input$cfa_path_model_struct, "_", input$cfa_path_style, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 12, height = 8.5)
      render_cfa_path_internal()
      dev.off()
    }
  )

  # ============================================================================
  # CFA MODEL FIT COMPARISON SERVER LOGIC (UNI, MULTI, SECOND-ORDER, BI-FACTOR)
  # ============================================================================
  cfa_comp_results_store <- reactiveVal(NULL)
  
  observeEvent(input$btn_run_cfa_comp, {
    df <- cfa_data_store()
    items <- cfa_active_items()
    df_sub <- df[, intersect(items, colnames(df)), drop = FALSE]
    req(ncol(df_sub) >= 4, nrow(df_sub) >= 20)
    
    is_builtin_poly <- (input$cfa_source_type == "builtin") || all(paste0("Item.", 1:20) %in% items)
    est_meth <- isolate(input$cfa_estimator)
    
    withProgress(message = "Mengestimasi Komparasi 4 Model CFA...", detail = "Menyusun sintaks model & kalkulasi MLR...", value = 0.2, {
      
      # 1. Unidimensional Model
      syn_uni <- paste0("F_Uni =~ ", paste(colnames(df_sub), collapse = " + "))
      
      # 2. Multidimensional Model
      d1_it <- intersect(paste0("Item.", 1:5), colnames(df_sub))
      d2_it <- intersect(paste0("Item.", 6:10), colnames(df_sub))
      d3_it <- intersect(paste0("Item.", 11:15), colnames(df_sub))
      d4_it <- intersect(paste0("Item.", 16:20), colnames(df_sub))
      
      if (!is_builtin_poly || length(d1_it) < 2 || length(d2_it) < 2 || length(d3_it) < 2 || length(d4_it) < 2) {
        k_tot <- ncol(df_sub)
        p1 <- colnames(df_sub)[1:ceiling(k_tot/4)]
        p2 <- colnames(df_sub)[(ceiling(k_tot/4)+1):ceiling(2*k_tot/4)]
        p3 <- colnames(df_sub)[(ceiling(2*k_tot/4)+1):ceiling(3*k_tot/4)]
        p4 <- colnames(df_sub)[(ceiling(3*k_tot/4)+1):k_tot]
        d1_it <- p1; d2_it <- p2; d3_it <- p3; d4_it <- p4
      }
      
      syn_multi <- paste(
        paste0("Faktor_1 =~ ", paste(d1_it, collapse = " + ")),
        paste0("Faktor_2 =~ ", paste(d2_it, collapse = " + ")),
        paste0("Faktor_3 =~ ", paste(d3_it, collapse = " + ")),
        paste0("Faktor_4 =~ ", paste(d4_it, collapse = " + ")),
        sep = "\n"
      )
      
      # 3. Second-Order Model
      syn_second <- paste(
        syn_multi,
        "General =~ Faktor_1 + Faktor_2 + Faktor_3 + Faktor_4",
        sep = "\n"
      )
      
      # 4. Bi-Factor Model (General orthogonal to specific factors)
      syn_bi <- paste(
        paste0("General =~ ", paste(colnames(df_sub), collapse = " + ")),
        paste0("Faktor_1 =~ ", paste(d1_it, collapse = " + ")),
        paste0("Faktor_2 =~ ", paste(d2_it, collapse = " + ")),
        paste0("Faktor_3 =~ ", paste(d3_it, collapse = " + ")),
        paste0("Faktor_4 =~ ", paste(d4_it, collapse = " + ")),
        "General ~~ 0*Faktor_1",
        "General ~~ 0*Faktor_2",
        "General ~~ 0*Faktor_3",
        "General ~~ 0*Faktor_4",
        "Faktor_1 ~~ 0*Faktor_2",
        "Faktor_1 ~~ 0*Faktor_3",
        "Faktor_1 ~~ 0*Faktor_4",
        "Faktor_2 ~~ 0*Faktor_3",
        "Faktor_2 ~~ 0*Faktor_4",
        "Faktor_3 ~~ 0*Faktor_4",
        sep = "\n"
      )
      
      is_cat_est <- est_meth %in% c("WLSMV", "WLSM", "WLSMVS", "DWLS", "ULSMV", "ULSM", "ULSMVS", "PML")
      fit_comp_model <- function(syn, use_ortho = FALSE, use_std_lv = TRUE) {
        tryCatch({
          if (is_cat_est) {
            lavaan::cfa(syn, data = df_sub, estimator = est_meth, ordered = TRUE, orthogonal = use_ortho, std.lv = use_std_lv)
          } else {
            lavaan::cfa(syn, data = df_sub, estimator = est_meth, orthogonal = use_ortho, std.lv = use_std_lv)
          }
        }, error = function(e) {
          tryCatch({
            lavaan::cfa(syn, data = df_sub, estimator = est_meth, ordered = FALSE, orthogonal = use_ortho, std.lv = use_std_lv)
          }, error = function(e2) NULL)
        })
      }
      
      setProgress(value = 0.4, detail = "Menjalankan estimasi Unidimensi & Multidimensi...")
      fit_uni <- fit_comp_model(syn_uni, use_ortho = FALSE, use_std_lv = TRUE)
      fit_multi <- fit_comp_model(syn_multi, use_ortho = FALSE, use_std_lv = TRUE)
      
      setProgress(value = 0.7, detail = "Menjalankan estimasi Second-Order & Bi-Factor...")
      fit_second <- fit_comp_model(syn_second, use_ortho = FALSE, use_std_lv = FALSE)
      fit_bi <- fit_comp_model(syn_bi, use_ortho = TRUE, use_std_lv = TRUE)
      
      extract_cfa_row <- function(fit, label, desc) {
        if (is.null(fit) || !lavInspect(fit, "converged")) {
          return(data.frame(
            Model = label, Struktur = desc, Chisq = NA_real_, df = NA_integer_, p_value = NA_real_,
            CFI = NA_real_, TLI = NA_real_, RMSEA = NA_real_, SRMR = NA_real_, AIC = NA_real_, BIC = NA_real_,
            Evaluasi = "Tidak Konvergen", stringsAsFactors = FALSE
          ))
        }
        fm <- fitMeasures(fit)
        get_val <- function(k_sc, k_un) if (!is.null(fm[k_sc]) && !is.na(fm[k_sc])) unname(fm[k_sc]) else unname(fm[k_un])
        cs <- round(get_val("chisq.scaled", "chisq"), 3)
        df_v <- as.integer(get_val("df.scaled", "df"))
        pv <- round(get_val("pvalue.scaled", "pvalue"), 3)
        cfi_v <- round(get_val("cfi.scaled", "cfi"), 3)
        tli_v <- round(get_val("tli.scaled", "tli"), 3)
        rmsea_v <- round(get_val("rmsea.scaled", "rmsea"), 3)
        srmr_v <- round(unname(fm["srmr"]), 3)
        aic_v <- if (!is.na(fm["aic"])) round(unname(fm["aic"]), 3) else NA_real_
        bic_v <- if (!is.na(fm["bic"])) round(unname(fm["bic"]), 3) else NA_real_
        
        eval_txt <- if (!is.na(cfi_v) && !is.na(rmsea_v) && cfi_v >= 0.95 && rmsea_v <= 0.06) "Good Fit (Hu & Bentler, 1999)" else if (!is.na(cfi_v) && !is.na(rmsea_v) && cfi_v >= 0.90 && rmsea_v <= 0.08) "Acceptable Fit" else "Marginal / Poor Fit"
        
        data.frame(
          Model = label, Struktur = desc, Chisq = cs, df = df_v, p_value = pv,
          CFI = cfi_v, TLI = tli_v, RMSEA = rmsea_v, SRMR = srmr_v, AIC = aic_v, BIC = bic_v,
          Evaluasi = eval_txt, stringsAsFactors = FALSE
        )
      }
      
      row_uni <- extract_cfa_row(fit_uni, "1. Unidimensi", "1 Faktor Umum (General)")
      row_multi <- extract_cfa_row(fit_multi, "2. Multidimensi", "3 Faktor Orde-1 Berkorelasi")
      row_second <- extract_cfa_row(fit_second, "3. Second-Order", "Faktor Hierarkis Menjelaskan Orde-1")
      row_bi <- extract_cfa_row(fit_bi, "4. Bi-Factor", "Faktor Umum + Spesifik Ortogonal")
      
      comp_df <- rbind(row_uni, row_multi, row_second, row_bi)
      
      # LRT Satorra-Bentler (2001) berpasangan vs Unidimensional Baseline
      fit_list_comp <- list(fit_uni, fit_multi, fit_second, fit_bi)
      lrt_pair_res <- lapply(2:4, function(i) {
        f_alt <- fit_list_comp[[i]]
        if (is.null(fit_uni) || !lavInspect(fit_uni, "converged") || is.null(f_alt) || !lavInspect(f_alt, "converged")) {
          return(list(d_cs = NA_real_, d_df = NA_integer_, p_lrt = NA_real_, row_df = NULL))
        }
        lrt <- tryCatch(
          lavaan::lavTestLRT(fit_uni, f_alt, method = "satorra.bentler.2001"),
          error = function(e) tryCatch(lavaan::lavTestLRT(fit_uni, f_alt), error = function(e2) NULL)
        )
        if (!is.null(lrt) && nrow(lrt) >= 2 && !is.na(lrt$`Chisq diff`[2])) {
          d_cs <- round(as.numeric(lrt$`Chisq diff`[2]), 3)
          d_df <- as.integer(lrt$`Df diff`[2])
          p_lrt <- round(as.numeric(lrt$`Pr(>Chisq)`[2]), 3)
          lrt_df <- as.data.frame(lrt)
          lrt_df$Perbandingan <- c("Baseline: Unidimensi", paste0(comp_df$Model[i], " vs Unidimensi"))
          return(list(d_cs = d_cs, d_df = d_df, p_lrt = p_lrt, row_df = lrt_df[2, , drop = FALSE]))
        } else {
          d_cs <- if (!is.na(comp_df$Chisq[1]) && !is.na(comp_df$Chisq[i])) round(comp_df$Chisq[1] - comp_df$Chisq[i], 3) else NA_real_
          d_df <- if (!is.na(comp_df$df[1]) && !is.na(comp_df$df[i])) as.integer(comp_df$df[1] - comp_df$df[i]) else NA_integer_
          p_lrt <- if (!is.na(d_cs) && !is.na(d_df) && d_df > 0 && d_cs > 0) round(pchisq(d_cs, df = d_df, lower.tail = FALSE), 3) else NA_real_
          return(list(d_cs = d_cs, d_df = d_df, p_lrt = p_lrt, row_df = NULL))
        }
      })
      
      comp_df$`\\Delta_Chisq` <- c(NA_real_, sapply(lrt_pair_res, `[[`, "d_cs"))
      comp_df$`\\Delta_df` <- c(NA_integer_, sapply(lrt_pair_res, `[[`, "d_df"))
      comp_df$p_LRT <- c(NA_real_, sapply(lrt_pair_res, `[[`, "p_lrt"))
      
      # Nested Model LRT Table (Pairwise Satorra-Bentler vs Unidimensi)
      anova_tbl <- tryCatch({
        lrt_rows <- Filter(Negate(is.null), lapply(lrt_pair_res, `[[`, "row_df"))
        if (length(lrt_rows) > 0 && !is.null(fit_uni) && lavInspect(fit_uni, "converged")) {
          an_df <- do.call(rbind, lrt_rows)
          rownames(an_df) <- NULL
          cols_ord <- c("Perbandingan", setdiff(colnames(an_df), "Perbandingan"))
          an_df[, cols_ord, drop = FALSE]
        } else {
          NULL
        }
      }, error = function(e) NULL)
      
      best_aic_mod <- if (all(is.na(comp_df$AIC))) {
        ratio <- comp_df$Chisq / pmax(1, comp_df$df)
        if (all(is.na(ratio))) "-" else paste0(comp_df$Model[which.min(ratio)], " (Min \u03c7\u00b2/df)")
      } else {
        comp_df$Model[which.min(comp_df$AIC)]
      }
      best_bic_mod <- if (all(is.na(comp_df$BIC))) {
        if (all(is.na(comp_df$SRMR))) "-" else paste0(comp_df$Model[which.min(comp_df$SRMR)], " (Min SRMR)")
      } else {
        comp_df$Model[which.min(comp_df$BIC)]
      }
      best_cfi_mod <- if (all(is.na(comp_df$CFI))) "-" else comp_df$Model[which.max(comp_df$CFI)]
      
      setProgress(value = 1.0, detail = "Komparasi selesai!")
      
      res_cfa_comp <- list(
        summary = comp_df,
        anova = anova_tbl,
        best_aic = best_aic_mod,
        best_bic = best_bic_mod,
        best_cfi = best_cfi_mod,
        n_samp = nrow(df_sub),
        k_items = ncol(df_sub)
      )
      cfa_comp_results_store(res_cfa_comp)
      showNotification("Komparasi 4 Model CFA berhasil diestimasi! Menampilkan tabel fit.", type = "message")
    })
  })
  
  output$vb_cfa_comp_best_aic <- renderText({
    res <- cfa_comp_results_store()
    if (is.null(res)) return("-")
    res$best_aic
  })
  output$vb_cfa_comp_best_bic <- renderText({
    res <- cfa_comp_results_store()
    if (is.null(res)) return("-")
    res$best_bic
  })
  output$vb_cfa_comp_best_cfi <- renderText({
    res <- cfa_comp_results_store()
    if (is.null(res)) return("-")
    res$best_cfi
  })
  output$vb_cfa_comp_n_k <- renderText({
    res <- cfa_comp_results_store()
    if (is.null(res)) return(paste(length(cfa_active_items()), "Butir"))
    paste0("N = ", res$n_samp, " / k = ", res$k_items)
  })
  
  output$table_cfa_comp_summary <- renderDT({
    res <- cfa_comp_results_store()
    if (is.null(res)) {
      return(datatable(data.frame(Status = "Silakan klik tombol 'Bandingkan 4 Model CFA' untuk memulai komparasi."), rownames = FALSE))
    }
    datatable(res$summary, options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE, class = "display compact cell-border")
  })
  
  output$table_cfa_comp_lrt <- renderTable({
    res <- cfa_comp_results_store()
    if (is.null(res) || is.null(res$anova)) {
      return(data.frame(Status = "Uji Likelihood Ratio Test (LRT) akan ditampilkan setelah estimasi komparasi dijalankan."))
    }
    res$anova
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$cfa_comp_bar_plot <- renderPlot({
    res <- cfa_comp_results_store()
    req(res, !is.null(res$summary))
    d <- res$summary
    req(nrow(d) >= 2)
    
    par(mar = c(5, 5, 3, 2), bg = "#FFFFFF", mfrow = c(1, 2))
    
    # Plot 1: CFI & TLI
    cfi_tli_mat <- rbind(d$CFI, d$TLI)
    colnames(cfi_tli_mat) <- c("Uni", "Multi", "2nd-Order", "Bi-Factor")[1:ncol(cfi_tli_mat)]
    barplot(
      cfi_tli_mat,
      beside = TRUE,
      col = c("#650915", "#C5A059"),
      ylim = c(0.70, 1.0),
      xpd = FALSE,
      las = 1,
      ylab = "Nilai Indeks Fit",
      main = "Komparasi CFI & TLI (Ambang Good Fit >= 0.95)",
      col.main = "#650915"
    )
    abline(h = 0.95, col = "#27AE60", lty = 2, lwd = 1.8)
    abline(h = 0.90, col = "#E67E22", lty = 3, lwd = 1.5)
    legend("bottomright", legend = c("CFI", "TLI"), fill = c("#650915", "#C5A059"), bty = "n")
    
    # Plot 2: RMSEA & SRMR
    rmsea_srmr_mat <- rbind(d$RMSEA, d$SRMR)
    colnames(rmsea_srmr_mat) <- colnames(cfi_tli_mat)
    barplot(
      rmsea_srmr_mat,
      beside = TRUE,
      col = c("#800020", "#5A4B4C"),
      ylim = c(0, max(0.10, max(d$RMSEA, na.rm = TRUE) + 0.02)),
      las = 1,
      ylab = "Nilai Residual / Error",
      main = "Komparasi RMSEA & SRMR (Ambang Ideal <= 0.06)",
      col.main = "#650915"
    )
    abline(h = 0.06, col = "#27AE60", lty = 2, lwd = 1.8)
    abline(h = 0.08, col = "#E67E22", lty = 3, lwd = 1.5)
    legend("topright", legend = c("RMSEA", "SRMR"), fill = c("#800020", "#5A4B4C"), bty = "n")
  })
  
  output$btn_dl_cfa_comp_csv <- downloadHandler(
    filename = function() paste0("CFA_Model_Comparison_", Sys.Date(), ".csv"),
    content = function(file) {
      res <- cfa_comp_results_store()
      req(res)
      write.csv(res$summary, file, row.names = FALSE)
    }
  )

  # ============================================================================
  # 2. RASCH SERVER LOGIC (TAM ENGINE: ITEM/PERSON FIT & WRIGHT MAP)
  # ============================================================================
  
  # Fungsi Eksekusi Paralel Multi-Core (doParallel & parallel)
  run_parallel_rasch_tasks <- function(task_list, n_cores = 2) {
    if (n_cores <= 1 || length(task_list) <= 1 || isTRUE(getOption("webr.active", FALSE)) || exists(".webr_pkg_cache", envir = .GlobalEnv)) {
      return(lapply(task_list, function(f) tryCatch(f(), error = function(e) NULL)))
    }
    
    if (.Platform$OS.type == "unix") {
      parallel::mclapply(task_list, function(f) {
        tryCatch(f(), error = function(e) NULL)
      }, mc.cores = min(length(task_list), n_cores))
    } else {
      cl <- parallel::makeCluster(min(length(task_list), n_cores))
      on.exit(parallel::stopCluster(cl))
      parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages(library(TAM))
      })
      parallel::parLapply(cl, task_list, function(f) {
        tryCatch(f(), error = function(e) NULL)
      })
    }
  }
  
  # Reactive Data Storage for Rasch
  rasch_data_store <- reactiveVal(simulate_rasch_dataset(n_persons = 300, n_items = 10, n_dims = 1, item_type = "dichotomous"))
  rasch_data_label <- reactiveVal("Simulasi Default Rasch (Person N = 300, k = 10 Dikotomi)")
  rasch_fit_trigger <- reactiveVal(0)
  rasch_active_items <- reactiveVal(NULL)
  rasch_active_persons <- reactiveVal(NULL)
  rasch_reanalysis_results <- reactiveVal(NULL)
  
  # Handle Rasch Simulation Button
  observeEvent(input$btn_rasch_gen_sim, {
    req(input$rasch_sim_n, input$rasch_sim_k)
    rasch_active_items(NULL)
    rasch_active_persons(NULL)
    rasch_reanalysis_results(NULL)
    n_d <- switch(input$rasch_sim_dim_count, "uni" = 1, "multi_2" = 2, "multi_3" = 3)
    sim_df <- simulate_rasch_dataset(
      n_persons = input$rasch_sim_n,
      n_items = input$rasch_sim_k,
      n_dims = n_d,
      item_type = input$rasch_sim_format,
      dim_cor = input$rasch_sim_cor,
      n_categories = input$rasch_sim_n_cats
    )
    rasch_data_store(sim_df)
    rasch_data_label(sprintf("Simulasi Rasch Baru (N = %d, k = %d, %s, %s)",
                             input$rasch_sim_n, input$rasch_sim_k, input$rasch_sim_format, input$rasch_sim_dim_count))
    showNotification("Data simulasi Rasch TAM berhasil dibangkitkan!", type = "message", duration = 4)
  })
  
  # Handle Rasch Built-in Data
  observeEvent(input$rasch_source_type, {
    rasch_active_items(NULL)
    rasch_active_persons(NULL)
    rasch_reanalysis_results(NULL)
    if (input$rasch_source_type == "builtin_qq") {
      poly_path <- "02_Working_Data/Polytomous_Benchmark_Data.csv"
      if (file.exists(poly_path)) {
        df_poly <- read.csv(poly_path, check.names = FALSE)
      } else {
        set.seed(2026)
        df_poly <- as.data.frame(matrix(sample(1:7, 999 * 20, replace = TRUE), nrow = 999, ncol = 20))
        colnames(df_poly) <- paste0("Item.", 1:20)
      }
      rasch_data_store(df_poly)
      rasch_data_label(sprintf("Dataset Bawaan Politomi (PCM, N = %d, k = %d, 4 Dimensi, Skala 1-7)", nrow(df_poly), ncol(df_poly)))
      updateRadioButtons(session, "rasch_model_dim_spec", selected = "multi")
      updateSelectInput(session, "rasch_irt_model_type", selected = "PCM")
      updateSliderInput(session, "rasch_multi_q_dims", value = 4)
      showNotification("Dataset bawaan Politomi berhasil dimuat. Parameter disesuaikan ke Multidimensi (4 Dimensi, k = 20) & PCM.", type = "message")
    } else if (input$rasch_source_type == "builtin") {
      # Use benchmark items
      df_b <- as.data.frame(matrix(sample(0:1, 400 * 8, replace = TRUE), nrow = 400, ncol = 8))
      colnames(df_b) <- paste0("NeoUni.", 1:8)
      rasch_data_store(df_b)
      rasch_data_label("Dataset Bawaan Dikotomi Rasch (N = 400, k = 8)")
      showNotification("Dataset bawaan Rasch berhasil dimuat.", type = "message")
    }
  })
  
  # Handle Rasch File Upload
  observeEvent(input$rasch_upload_file, {
    req(input$rasch_upload_file)
    rasch_active_items(NULL)
    rasch_active_persons(NULL)
    rasch_reanalysis_results(NULL)
    fpath <- input$rasch_upload_file$datapath
    fname <- input$rasch_upload_file$name
    ext <- tools::file_ext(fname)
    tryCatch({
      if (ext == "csv") {
        df_up <- read.csv(fpath, check.names = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        df_up <- as.data.frame(readxl::read_excel(fpath))
      } else {
        stop("Format berkas harus .csv atau .xlsx")
      }
      num_cols <- names(df_up)[sapply(df_up, is.numeric)]
      if (length(num_cols) < 3) stop("Berkas harus memuat minimal 3 kolom butir numerik.")
      df_up <- df_up[, num_cols, drop = FALSE]
      rasch_data_store(df_up)
      rasch_data_label(sprintf("Unggahan Respons: %s (N = %d, k = %d)", fname, nrow(df_up), ncol(df_up)))
      showNotification(sprintf("Berhasil memuat %d butir dari berkas %s", ncol(df_up), fname), type = "message")
    }, error = function(e) {
      showNotification(paste("Gagal memuat berkas:", e$message), type = "error")
    })
  })

  # Dynamic UI for Interactive Q-Matrix Dimension Mapper
  output$rasch_qmatrix_interactive_ui <- renderUI({
    df <- rasch_data_store()
    items <- colnames(df)
    n_d <- input$rasch_multi_q_dims
    req(n_d, length(items) >= n_d)
    is_builtin_poly <- (input$rasch_source_type == "builtin_qq") || all(paste0("Item.", 1:20) %in% items)
    
    if (is_builtin_poly && n_d == 4) {
      poly_specs <- list(
        list(name = "Faktor_1", label = "Dimensi 1 (5 Butir: 1 s/d 5)", items = paste0("Item.", 1:5)),
        list(name = "Faktor_2", label = "Dimensi 2 (5 Butir: 6 s/d 10)", items = paste0("Item.", 6:10)),
        list(name = "Faktor_3", label = "Dimensi 3 (5 Butir: 11 s/d 15)", items = paste0("Item.", 11:15)),
        list(name = "Faktor_4", label = "Dimensi 4 (5 Butir: 16 s/d 20)", items = paste0("Item.", 16:20))
      )
      cards <- lapply(1:4, function(d) {
        div(
          class = "p-2 mb-2 bg-white rounded border",
          textInput(paste0("rasch_qdim_name_", d), sprintf("Nama Dimensi %d:", d), value = poly_specs[[d]]$name),
          span(class = "small text-muted d-block mb-1", poly_specs[[d]]$label),
          selectizeInput(
            paste0("rasch_qdim_items_", d),
            sprintf("Pilih Butir Dimensi %d:", d),
            choices = items,
            selected = intersect(poly_specs[[d]]$items, items),
            multiple = TRUE,
            options = list(plugins = list('remove_button'))
          )
        )
      })
      tagList(cards)
    } else {
      items_per_dim <- ceiling(length(items) / n_d)
      cards <- lapply(1:n_d, function(d) {
        start_i <- (d - 1) * items_per_dim + 1
        end_i <- min(d * items_per_dim, length(items))
        def_sel <- if (start_i <= length(items)) items[start_i:end_i] else items[1]
        div(
          class = "p-2 mb-2 bg-white rounded border",
          textInput(paste0("rasch_qdim_name_", d), sprintf("Nama Dimensi %d:", d), value = paste0("Dimensi_", d)),
          selectizeInput(
            paste0("rasch_qdim_items_", d),
            sprintf("Pilih Butir Dimensi %d:", d),
            choices = items,
            selected = def_sel,
            multiple = TRUE,
            options = list(plugins = list('remove_button'))
          )
        )
      })
      tagList(cards)
    }
  })

  # Handle Button to Load Standard Q-Matrix (4 Dimensi, k=20)
  observeEvent(input$btn_load_std_qmatrix, {
    df <- rasch_data_store()
    items <- colnames(df)
    updateSliderInput(session, "rasch_multi_q_dims", value = 4)
    updateTextInput(session, "rasch_qdim_name_1", value = "Faktor_1")
    updateSelectizeInput(session, "rasch_qdim_items_1", selected = intersect(paste0("Item.", 1:5), items))
    updateTextInput(session, "rasch_qdim_name_2", value = "Faktor_2")
    updateSelectizeInput(session, "rasch_qdim_items_2", selected = intersect(paste0("Item.", 6:10), items))
    updateTextInput(session, "rasch_qdim_name_3", value = "Faktor_3")
    updateSelectizeInput(session, "rasch_qdim_items_3", selected = intersect(paste0("Item.", 11:15), items))
    updateTextInput(session, "rasch_qdim_name_4", value = "Faktor_4")
    updateSelectizeInput(session, "rasch_qdim_items_4", selected = intersect(paste0("Item.", 16:20), items))
    showNotification("Desain Q-Matrix 4 Dimensi Standar (k = 20) berhasil dimuat.", type = "message")
  })

  # Reactive Q-Matrix Builder for TAM
  rasch_active_qmatrix <- reactive({
    df <- rasch_data_store()
    items <- colnames(df)
    k <- length(items)
    mode <- input$rasch_qmatrix_mode
    n_d <- input$rasch_multi_q_dims
    req(k >= 2, n_d >= 2)
    
    is_builtin_poly <- (input$rasch_source_type == "builtin_qq") || all(paste0("Item.", 1:20) %in% items)
    
    if (mode == "interactive") {
      Q <- matrix(0, nrow = k, ncol = n_d)
      rownames(Q) <- items
      c_names <- c()
      for (d in 1:n_d) {
        cn <- input[[paste0("rasch_qdim_name_", d)]]
        if (is.null(cn) || cn == "") cn <- paste0("Dimensi_", d)
        c_names <- c(c_names, cn)
        itms <- input[[paste0("rasch_qdim_items_", d)]]
        matched <- intersect(itms, items)
        if (length(matched) > 0) {
          Q[matched, d] <- 1
        }
      }
      colnames(Q) <- c_names
      if (all(Q == 0)) {
        if (is_builtin_poly && n_d == 4) {
          colnames(Q) <- c("Faktor_1", "Faktor_2", "Faktor_3", "Faktor_4")
          Q[intersect(paste0("Item.", 1:5), items), 1] <- 1
          Q[intersect(paste0("Item.", 6:10), items), 2] <- 1
          Q[intersect(paste0("Item.", 11:15), items), 3] <- 1
          Q[intersect(paste0("Item.", 16:20), items), 4] <- 1
        } else {
          for (i in 1:k) Q[i, min(n_d, ceiling(i / ceiling(k / n_d)))] <- 1
        }
      }
      Q
    } else if (mode == "syntax") {
      req(input$rasch_qmatrix_syntax_text)
      txt <- input$rasch_qmatrix_syntax_text
      q_mat <- tryCatch({
        env <- new.env()
        eval(parse(text = txt), envir = env)
        res_mat <- NULL
        for (obj in ls(env)) {
          if (is.matrix(env[[obj]])) {
            res_mat <- env[[obj]]
            break
          }
        }
        if (is.null(res_mat) && exists("Q", envir = env)) res_mat <- env$Q
        res_mat
      }, error = function(e) NULL)
      
      if (!is.null(q_mat) && nrow(q_mat) == k) {
        rownames(q_mat) <- items
        if (is.null(colnames(q_mat))) colnames(q_mat) <- paste0("Dimensi_", 1:ncol(q_mat))
        q_mat
      } else {
        Q <- matrix(0, nrow = k, ncol = n_d)
        rownames(Q) <- items
        colnames(Q) <- paste0("Dimensi_", 1:n_d)
        for (i in 1:k) Q[i, min(n_d, ceiling(i / ceiling(k / n_d)))] <- 1
        Q
      }
    } else {
      if (is_builtin_poly && n_d == 3) {
        Q <- matrix(0, nrow = 17, ncol = 3)
        rownames(Q) <- items
        colnames(Q) <- c("Faktor_1", "Faktor_2", "Faktor_3")
        Q[c("Item.1", "Item.4", "Item.8", "Item.11", "Item.13", "Item.16"), 1] <- 1
        Q[c("Item.2", "Item.6", "Item.7", "Item.9", "Item.12", "Item.15", "Item.17"), 2] <- 1
        Q[c("Item.3", "Item.5", "Item.10", "Item.14"), 3] <- 1
        Q
      } else {
        Q <- matrix(0, nrow = k, ncol = n_d)
        rownames(Q) <- items
        colnames(Q) <- paste0("Dimensi_", 1:n_d)
        for (i in 1:k) {
          d_assigned <- min(n_d, ceiling(i / ceiling(k / n_d)))
          Q[i, d_assigned] <- 1
        }
        Q
      }
    }
  })

  # Q-Matrix Live Table Display & Dimension Badge
  output$rasch_qmatrix_preview_table <- renderTable({
    Q <- rasch_active_qmatrix()
    req(Q)
    df_q <- as.data.frame(Q)
    cbind(Butir = rownames(df_q), df_q)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$rasch_qmatrix_badge_dims <- renderText({
    Q <- rasch_active_qmatrix()
    if (is.null(Q)) return("-")
    sprintf("%d Butir x %d Dimensi", nrow(Q), ncol(Q))
  })
  
  # Trigger estimation and show waiting status
  observeEvent(input$btn_run_rasch, {
    rasch_active_items(NULL)
    rasch_active_persons(NULL)
    rasch_reanalysis_results(NULL)
    showNotification(
      id = "tam_progress_notif",
      tagList(
        div(
          class = "d-flex align-items-center gap-2",
          tags$div(class = "spinner-border spinner-border-sm text-primary", role = "status"),
          tags$strong("Estimasi Model TAM Sedang Berjalan...")
        ),
        p(class = "small mb-0 mt-1", "Algoritma EM dan integrasi numerik sedang mengkalibrasi parameter butir dan kemampuan laten. Mohon tunggu beberapa saat.")
      ),
      type = "message",
      duration = NULL
    )
  })

  # Reactive TAM Fitted Model with Progress
  rasch_tam_fit_reactive <- reactive({
    input$btn_run_rasch
    rasch_fit_trigger()
    df_raw <- rasch_data_store()
    req(nrow(df_raw) >= 30, ncol(df_raw) >= 3)
    
    # Filter responden aktif (Person filter)
    act_p <- rasch_active_persons()
    df <- if (!is.null(act_p) && length(act_p) >= 30 && max(act_p) <= nrow(df_raw)) {
      df_raw[act_p, , drop = FALSE]
    } else {
      df_raw
    }
    
    # Filter butir aktif (Item filter)
    act_i <- rasch_active_items()
    if (!is.null(act_i) && length(act_i) >= 3) {
      common_i <- intersect(act_i, colnames(df))
      if (length(common_i) >= 3) {
        df <- df[, common_i, drop = FALSE]
      }
    }
    
    dim_spec <- isolate(input$rasch_model_dim_spec)
    irt_type <- isolate(input$rasch_irt_model_type)
    
    n_cores_alloc <- if (isTRUE(input$rasch_use_parallel)) {
      if (!is.null(input$rasch_n_cores)) as.integer(input$rasch_n_cores) else 2
    } else 1
    
    # Konfigurasi Threading Multi-Core doParallel & RcppParallel
    if (requireNamespace("RcppParallel", quietly = TRUE)) {
      tryCatch(RcppParallel::setThreadOptions(numThreads = n_cores_alloc), error = function(e) NULL)
    }
    
    withProgress(message = sprintf("Estimasi TAM (IRT) Sedang Berlangsung [%d Cores]...", n_cores_alloc), 
                 detail = "Menyiapkan matriks respons & estimasi EM (doParallel)...", value = 0.2, {
      
      df_mat <- as.matrix(df)
      min_val <- min(df_mat, na.rm = TRUE)
      if (min_val > 0) df_mat <- df_mat - min_val
      
      setProgress(value = 0.4, detail = "Menyusun matriks spesifikasi Q & parameter butir...")
      
      Q <- NULL
      if (dim_spec == "multi") {
        Q_raw <- isolate(rasch_active_qmatrix())
        if (!is.null(Q_raw)) {
          common_q <- intersect(rownames(Q_raw), colnames(df_mat))
          if (length(common_q) >= 2) {
            Q_cand <- Q_raw[common_q, , drop = FALSE]
            v_dims <- which(colSums(Q_cand) > 0)
            if (length(v_dims) >= 2) {
              Q <- Q_cand[, v_dims, drop = FALSE]
            }
          }
        }
      }
      
      setProgress(value = 0.6, detail = "Menjalankan iterasi EM pada engine TAM...")
      
      n_dims_cur <- if (!is.null(Q)) ncol(Q) else 1
      ctrl_opts <- list(maxiter = 35)
      if (n_dims_cur >= 3) {
        ctrl_opts$nodes <- seq(-3, 3, length.out = 8)
      } else if (n_dims_cur == 2) {
        ctrl_opts$nodes <- seq(-3.5, 3.5, length.out = 12)
      }
      
      mod <- tryCatch({
        if (irt_type %in% c("1PL", "PCM", "PCM2", "RSM")) {
          if (!is.null(Q)) {
            TAM::tam.mml(resp = df_mat, Q = Q, irtmodel = irt_type, verbose = FALSE, control = ctrl_opts)
          } else {
            TAM::tam.mml(resp = df_mat, irtmodel = irt_type, verbose = FALSE, control = ctrl_opts)
          }
        } else if (irt_type %in% c("2PL", "GPCM", "2PL.groups", "GPCM.groups")) {
          if (!is.null(Q)) {
            TAM::tam.mml.2pl(resp = df_mat, Q = Q, irtmodel = irt_type, verbose = FALSE, control = ctrl_opts)
          } else {
            TAM::tam.mml.2pl(resp = df_mat, irtmodel = irt_type, verbose = FALSE, control = ctrl_opts)
          }
        } else if (irt_type == "3PL") {
          TAM::tam.mml.3pl(resp = df_mat, Q = Q, est.guess = 1:ncol(df_mat), verbose = FALSE, control = list(maxiter = 30))
        } else if (irt_type == "JML") {
          TAM::tam.jml(resp = df_mat, verbose = FALSE, control = list(maxiter = 30))
        } else {
          TAM::tam.mml(resp = df_mat, irtmodel = "1PL", verbose = FALSE, control = ctrl_opts)
        }
      }, error = function(e) {
        NULL
      })
      
      setProgress(value = 1.0, detail = "Estimasi selesai!")
      removeNotification("tam_progress_notif")
      
      if (!is.null(mod)) {
        showNotification(
          sprintf("Estimasi TAM (%s) berhasil konvergen!", irt_type),
          type = "message",
          duration = 4
        )
      } else {
        showNotification(
          "Estimasi TAM tidak dapat konvergen. Periksa kembali spesifikasi model atau butir respons.",
          type = "error",
          duration = 5
        )
      }
      
      mod
    })
  })
  
  # Helper Komputasi Eksak Paket TAM (Mendukung 1D, Multidimensi Q-Matrix & JML)
  safe_tam_wle <- function(m) {
    if (inherits(m, "tam.jml")) {
      th <- if (!is.null(m$WLE)) as.vector(m$WLE) else as.vector(m$theta)
      err <- if (!is.null(m$errorWLE)) as.vector(m$errorWLE) else if (!is.null(m$errorP)) rep(NA_real_, length(th)) else rep(NA_real_, length(th))
      rel <- if (!is.null(m$WLEreliability)) m$WLEreliability[1] else NA_real_
      return(data.frame(theta = th, error = err, WLE.rel = rel))
    }
    tryCatch(TAM::tam.wle(m, progress = FALSE), error = function(e) NULL)
  }

  calc_rasch_item_se <- function(m) {
    nitems <- ncol(m$resp)
    item_names <- colnames(m$resp)
    if (!inherits(m, "tam.jml") && !is.null(m$xsi) && is.data.frame(m$xsi) && all(item_names %in% rownames(m$xsi))) {
      return(as.numeric(m$xsi[item_names, "se.xsi"]))
    }
    se_raw <- if (inherits(m, "tam.jml")) m$errorP else if (!is.null(m$xsi$se.xsi)) m$xsi$se.xsi else NULL
    if (is.null(se_raw) || is.null(m$A) || is.null(m$AXsi)) return(rep(NA_real_, nitems))
    maxK_vec <- apply(m$AXsi, 1, function(r) {
      v <- which(!is.na(r))
      if (length(v) > 0) max(v) else 1
    })
    sapply(seq_len(nitems), function(i) {
      k_idx <- maxK_vec[i]
      K_i <- max(1, k_idx - 1)
      c_vec <- m$A[i, k_idx, ] / K_i
      sqrt(sum((c_vec * se_raw)^2, na.rm = TRUE))
    })
  }

  calc_rasch_personfit <- function(m, wle_res = NULL, trim_val = 10) {
    if (inherits(m, "tam.jml")) {
      return(tryCatch(TAM::tam.fit(m)$fit.person, error = function(e) NULL))
    }
    nd <- if (!is.null(m$ndim)) m$ndim else 1
    if (nd <= 1) {
      p_1d <- tryCatch(TAM::tam.personfit(m), error = function(e) NULL)
      if (!is.null(p_1d)) return(p_1d)
    }
    # Estimasi Person Fit Multidimensi (Wright & Masters, 1982) berbasis WLE vektor theta_dim
    tryCatch({
      if (is.null(wle_res)) wle_res <- safe_tam_wle(m)
      if (is.null(wle_res)) return(NULL)
      theta_cols <- grep("^theta\\.Dim", colnames(wle_res))
      theta_mat <- if (length(theta_cols) > 0) as.matrix(wle_res[, theta_cols, drop = FALSE]) else matrix(wle_res$theta, ncol = 1)
      resp <- m$resp
      resp.ind <- m$resp.ind
      A <- m$A
      B <- m$B
      nstud <- nrow(resp)
      nitems <- ncol(resp)
      maxK <- dim(A)[2]
      xsi <- if (is.data.frame(m$xsi)) m$xsi$xsi else m$xsi
      AXsi <- matrix(0, nrow = nitems, ncol = maxK)
      for (k in 1:maxK) AXsi[, k] <- A[, k, ] %*% xsi
      B0 <- B
      B0[is.na(B0)] <- 0
      W_ik <- apply(B0, c(1, 2), sum)

      B_bari <- matrix(0, nrow = nstud, ncol = nitems)
      BB_bari <- matrix(0, nrow = nstud, ncol = nitems)
      C4 <- matrix(0, nrow = nstud, ncol = nitems)

      for (i in seq_len(nitems)) {
        eta <- matrix(0, nrow = nstud, ncol = maxK)
        for (k in 1:maxK) {
          b_vec <- B0[i, k, ]
          eta[, k] <- as.vector(theta_mat %*% b_vec) + AXsi[i, k]
        }
        valid_k <- !is.na(m$AXsi[i, ])
        eta[, !valid_k] <- -Inf
        max_eta <- apply(eta, 1, max)
        exp_eta <- exp(eta - max_eta)
        probs <- exp_eta / rowSums(exp_eta)
        w_i <- W_ik[i, ]
        e_i <- as.vector(probs %*% w_i)
        e2_i <- as.vector(probs %*% (w_i^2))
        B_bari[, i] <- e_i
        BB_bari[, i] <- e2_i
        for (k in which(valid_k)) {
          C4[, i] <- C4[, i] + ((w_i[k] - e_i)^4) * probs[, k]
        }
      }

      B_bari[resp.ind == 0] <- NA
      BB_bari[resp.ind == 0] <- NA
      C4[resp.ind == 0] <- NA
      resp_na <- resp
      resp_na[resp.ind == 0] <- NA

      B_var <- BB_bari - (B_bari^2)
      B_var[B_var <= 1e-8] <- NA
      z_sq <- ((resp_na - B_bari)^2) / B_var
      zsd <- stats::sd(as.vector(z_sq), na.rm = TRUE)
      if (!is.na(zsd) && zsd > 0) z_sq[z_sq > trim_val * zsd] <- trim_val * zsd

      outfitPerson <- rowMeans(z_sq, na.rm = TRUE)
      infitPerson <- rowSums((resp_na - B_bari)^2, na.rm = TRUE) / rowSums(B_var, na.rm = TRUE)
      n_valid <- pmax(1, rowSums(!is.na(z_sq)))
      var_outfit <- pmax(1e-6, rowSums(C4 / (B_var^2), na.rm = TRUE) / (n_valid^2) - 1 / n_valid)
      outfitPerson_t <- (outfitPerson^(1/3) - 1) * (3 / sqrt(var_outfit)) + sqrt(var_outfit) / 3
      var_infit <- pmax(1e-6, rowSums(C4 - B_var^2, na.rm = TRUE) / (pmax(1e-6, rowSums(B_var, na.rm = TRUE))^2))
      infitPerson_t <- (infitPerson^(1/3) - 1) * (3 / sqrt(var_infit)) + sqrt(var_infit) / 3

      data.frame(
        outfitPerson = outfitPerson,
        outfitPerson_t = outfitPerson_t,
        infitPerson = infitPerson,
        infitPerson_t = infitPerson_t
      )
    }, error = function(e) NULL)
  }

  # Rasch Model Summary & Value Boxes
  output$vb_rasch_dev <- renderText({
    m <- rasch_tam_fit_reactive()
    if (is.null(m) || is.null(m$deviance)) return("-")
    sprintf("%.3f", m$deviance)
  })
  output$vb_rasch_aic_bic <- renderText({
    m <- rasch_tam_fit_reactive()
    if (is.null(m)) return("-")
    aic_v <- if (!is.null(m$ic$AIC)) m$ic$AIC else if (!is.null(m$aic)) m$aic else NA
    bic_v <- if (!is.null(m$ic$BIC)) m$ic$BIC else if (!is.null(m$bic)) m$bic else NA
    if (is.na(aic_v) || is.na(bic_v)) return("-")
    sprintf("%.3f / %.3f", aic_v, bic_v)
  })
  output$vb_rasch_rel <- renderText({
    m <- rasch_tam_fit_reactive()
    if (is.null(m)) return("-")
    rel_val <- if (!is.null(m$EAP.rel)) m$EAP.rel[1] else if (!is.null(m$WLEreliability)) m$WLEreliability[1] else NA
    if (is.na(rel_val)) "-" else sprintf("%.3f", rel_val)
  })
  output$vb_rasch_dims <- renderText({
    m <- rasch_tam_fit_reactive()
    if (is.null(m)) return("-")
    nd <- if (!is.null(m$ndim)) m$ndim else 1
    paste(nd, "Dimensi Laten")
  })
  
  output$table_rasch_latent_cov <- renderTable({
    m <- rasch_tam_fit_reactive()
    req(m)
    if (is.null(m$variance)) return(data.frame(Keterangan = "Varians laten teridentifikasi tetap pada model (fixed)"))
    v <- as.data.frame(round(m$variance, 3))
    dim_names <- if (!is.null(colnames(m$variance))) colnames(m$variance) else paste0("Dimensi_", 1:ncol(v))
    colnames(v) <- dim_names
    rownames(v) <- dim_names
    cbind(Dimensi = rownames(v), v)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$table_rasch_latent_cor <- renderTable({
    m <- rasch_tam_fit_reactive()
    req(m)
    if (is.null(m$variance)) return(data.frame(Keterangan = "Varians laten teridentifikasi tetap pada model (fixed)"))
    cr <- tryCatch(cov2cor(m$variance), error = function(e) m$variance)
    v <- as.data.frame(round(cr, 3))
    dim_names <- if (!is.null(colnames(m$variance))) colnames(m$variance) else paste0("Dimensi_", 1:ncol(v))
    colnames(v) <- dim_names
    rownames(v) <- dim_names
    cbind(Dimensi = rownames(v), v)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$table_rasch_reliability_dims <- renderTable({
    m <- rasch_tam_fit_reactive()
    req(m)
    wle_res <- safe_tam_wle(m)
    nd <- if (!is.null(m$ndim)) m$ndim else 1
    dim_names <- if (!is.null(colnames(m$variance))) colnames(m$variance) else paste("Dimensi", 1:nd)
    
    eap_rels <- if (!is.null(m$EAP.rel)) round(m$EAP.rel, 3) else rep(NA, nd)
    wle_rels <- rep(NA, nd)
    if (!is.null(wle_res)) {
      for (i in 1:nd) {
        cands <- c(paste0("WLE.rel.Dim0", i), paste0("WLE.rel.Dim", i), paste0("WLE.rel.Dim0", i, "[1]"))
        matched <- intersect(cands, colnames(wle_res))
        if (length(matched) > 0) {
          wle_rels[i] <- round(wle_res[[matched[1]]][1], 3)
        } else if (i == 1 && !is.null(wle_res$WLE.rel)) {
          wle_rels[i] <- round(wle_res$WLE.rel[1], 3)
        }
      }
    }
    
    rel_check <- ifelse(!is.na(eap_rels), eap_rels, wle_rels)
    data.frame(
      Dimensi = dim_names,
      WLE_Reliability = ifelse(is.na(wle_rels), "-", sprintf("%.3f", wle_rels)),
      EAP_Reliability = ifelse(is.na(eap_rels), "-", sprintf("%.3f", eap_rels)),
      Evaluasi = ifelse(!is.na(rel_check) & rel_check >= 0.70, "Reliabilitas Konstruk Sangat Baik (>= 0.70)", "Dapat Diterima / Butuh Penguatan Butir"),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$table_rasch_thurstonian_thresholds <- renderTable({
    m <- rasch_tam_fit_reactive()
    req(m)
    thr <- tryCatch(TAM::tam.threshold(m), error = function(e) NULL)
    if (is.null(thr)) return(data.frame(Keterangan = "Thurstonian thresholds hanya dapat dihitung pada model politomi (PCM/RSM)."))
    thr_df <- as.data.frame(round(thr, 3))
    cbind(Butir = rownames(thr_df), thr_df)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$table_rasch_person_fit_dim_summary <- renderTable({
    m <- rasch_tam_fit_reactive()
    req(m)
    df <- rasch_data_store()
    min_val <- min(df, na.rm = TRUE)
    if (min_val > 0) df <- df - min_val
    
    is_multi_mod <- (!is.null(m$ndim) && m$ndim > 1)
    Q_mat <- if (is_multi_mod) tryCatch(isolate(rasch_active_qmatrix()), error = function(e) NULL) else NULL
    is_builtin_poly <- all(paste0("Item.", 1:20) %in% colnames(df))
    sub_irt <- if (max(df, na.rm = TRUE) <= 1) "1PL" else "PCM"
    
    eval_mnsq_vec <- function(inf_vec, out_vec) {
      inf_m <- mean(inf_vec, na.rm = TRUE)
      out_m <- mean(out_vec, na.rm = TRUE)
      if (!is.finite(inf_m)) return("-")
      if (inf_m >= 0.50 && inf_m <= 1.50 && out_m >= 0.50 && out_m <= 1.50) {
        "Produktif (0.5 - 1.5 MNSQ)"
      } else {
        "Perlu Evaluasi Respons Aberan"
      }
    }
    
    n_cores_alloc <- if (isTRUE(input$rasch_use_parallel)) {
      if (!is.null(input$rasch_n_cores)) as.integer(input$rasch_n_cores) else 2
    } else 1
    
    if (is_multi_mod && is_builtin_poly) {
      d1_items <- paste0("Item.", 1:5)
      d2_items <- paste0("Item.", 6:10)
      d3_items <- paste0("Item.", 11:15)
      d4_items <- paste0("Item.", 16:20)
      
      dim_tasks <- list(
        d1 = function() TAM::tam.personfit(TAM::tam.mml(resp = df[, d1_items], irtmodel = sub_irt, verbose = FALSE)),
        d2 = function() TAM::tam.personfit(TAM::tam.mml(resp = df[, d2_items], irtmodel = sub_irt, verbose = FALSE)),
        d3 = function() TAM::tam.personfit(TAM::tam.mml(resp = df[, d3_items], irtmodel = sub_irt, verbose = FALSE)),
        d4 = function() TAM::tam.personfit(TAM::tam.mml(resp = df[, d4_items], irtmodel = sub_irt, verbose = FALSE))
      )
      pfit_list <- run_parallel_rasch_tasks(dim_tasks, n_cores = n_cores_alloc)
      pfit_d1 <- pfit_list$d1
      pfit_d2 <- pfit_list$d2
      pfit_d3 <- pfit_list$d3
      pfit_d4 <- pfit_list$d4
      
      data.frame(
        Dimensi = c("D1 (Faktor 1)", "D2 (Faktor 2)", "D3 (Faktor 3)", "D4 (Faktor 4)"),
        Jumlah_Butir = c(length(d1_items), length(d2_items), length(d3_items), length(d4_items)),
        Mean_Infit_MNSQ = c(
          if (!is.null(pfit_d1)) sprintf("%.3f", mean(pfit_d1$infitPerson, na.rm = TRUE)) else "-",
          if (!is.null(pfit_d2)) sprintf("%.3f", mean(pfit_d2$infitPerson, na.rm = TRUE)) else "-",
          if (!is.null(pfit_d3)) sprintf("%.3f", mean(pfit_d3$infitPerson, na.rm = TRUE)) else "-",
          if (!is.null(pfit_d4)) sprintf("%.3f", mean(pfit_d4$infitPerson, na.rm = TRUE)) else "-"
        ),
        Mean_Outfit_MNSQ = c(
          if (!is.null(pfit_d1)) sprintf("%.3f", mean(pfit_d1$outfitPerson, na.rm = TRUE)) else "-",
          if (!is.null(pfit_d2)) sprintf("%.3f", mean(pfit_d2$outfitPerson, na.rm = TRUE)) else "-",
          if (!is.null(pfit_d3)) sprintf("%.3f", mean(pfit_d3$outfitPerson, na.rm = TRUE)) else "-",
          if (!is.null(pfit_d4)) sprintf("%.3f", mean(pfit_d4$outfitPerson, na.rm = TRUE)) else "-"
        ),
        Evaluasi = c(
          if (!is.null(pfit_d1)) eval_mnsq_vec(pfit_d1$infitPerson, pfit_d1$outfitPerson) else "-",
          if (!is.null(pfit_d2)) eval_mnsq_vec(pfit_d2$infitPerson, pfit_d2$outfitPerson) else "-",
          if (!is.null(pfit_d3)) eval_mnsq_vec(pfit_d3$infitPerson, pfit_d3$outfitPerson) else "-",
          if (!is.null(pfit_d4)) eval_mnsq_vec(pfit_d4$infitPerson, pfit_d4$outfitPerson) else "-"
        ),
        stringsAsFactors = FALSE
      )
    } else if (is_multi_mod && !is.null(Q_mat) && ncol(Q_mat) > 1) {
      dim_cols <- colnames(Q_mat)
      dim_tasks <- lapply(dim_cols, function(d_col) {
        force(d_col)
        function() {
          items_d <- rownames(Q_mat)[Q_mat[, d_col] == 1]
          items_in_df <- intersect(items_d, colnames(df))
          if (length(items_in_df) >= 2) {
            sub_mod <- tryCatch(TAM::tam.mml(resp = df[, items_in_df], irtmodel = sub_irt, verbose = FALSE), error = function(e) NULL)
            sub_pfit <- if (!is.null(sub_mod)) tryCatch(TAM::tam.personfit(sub_mod), error = function(e) NULL) else NULL
            data.frame(
              Dimensi = d_col,
              Jumlah_Butir = length(items_in_df),
              Mean_Infit_MNSQ = if (!is.null(sub_pfit)) sprintf("%.3f", mean(sub_pfit$infitPerson, na.rm = TRUE)) else "-",
              Mean_Outfit_MNSQ = if (!is.null(sub_pfit)) sprintf("%.3f", mean(sub_pfit$outfitPerson, na.rm = TRUE)) else "-",
              Evaluasi = if (!is.null(sub_pfit)) eval_mnsq_vec(sub_pfit$infitPerson, sub_pfit$outfitPerson) else "-",
              stringsAsFactors = FALSE
            )
          } else {
            data.frame(Dimensi = d_col, Jumlah_Butir = length(items_in_df), Mean_Infit_MNSQ = "-", Mean_Outfit_MNSQ = "-", Evaluasi = "Butir < 2", stringsAsFactors = FALSE)
          }
        }
      })
      dim_res_list <- run_parallel_rasch_tasks(dim_tasks, n_cores = n_cores_alloc)
      do.call(rbind, dim_res_list)
    } else {
      p_fit <- calc_rasch_personfit(m)
      data.frame(
        Dimensi = "Keseluruhan Instrumen (Unidimensi)",
        Jumlah_Butir = ncol(df),
        Mean_Infit_MNSQ = if (!is.null(p_fit)) sprintf("%.3f", mean(p_fit$infitPerson, na.rm = TRUE)) else "-",
        Mean_Outfit_MNSQ = if (!is.null(p_fit)) sprintf("%.3f", mean(p_fit$outfitPerson, na.rm = TRUE)) else "-",
        Evaluasi = if (!is.null(p_fit)) eval_mnsq_vec(p_fit$infitPerson, p_fit$outfitPerson) else "Tidak Dapat Diestimasi",
        stringsAsFactors = FALSE
      )
    }
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$rasch_tam_console_summary <- renderPrint({
    m <- rasch_tam_fit_reactive()
    req(m)
    summary(m)
  })
  
  output$rasch_data_badge_ui <- renderUI({
    act_i <- rasch_active_items()
    act_p <- rasch_active_persons()
    is_filtered <- (!is.null(act_i) && length(act_i) > 0) || (!is.null(act_p) && length(act_p) > 0)
    tagList(
      span(class = "badge bg-primary me-2", rasch_data_label()),
      if (is_filtered) {
        span(class = "badge bg-warning text-dark", 
             sprintf("Filter Pemurnian Aktif: %d Butir, %d Responden", 
                     ifelse(!is.null(act_i), length(act_i), ncol(rasch_data_store())), 
                     ifelse(!is.null(act_p), length(act_p), nrow(rasch_data_store()))))
      }
    )
  })
  output$table_rasch_resp_preview <- renderDT({
    df <- rasch_data_store()
    req(df)
    datatable(
      df,
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 20, 25, 50, 100),
        scrollX = TRUE,
        autoWidth = FALSE,
        dom = "lfrtip",
        columnDefs = list(
          list(className = "dt-center", targets = "_all")
        )
      ),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })
  
  # --- Rasch Item Fit Table & Plot ---
  rasch_itemfit_df_reactive <- reactive({
    m <- rasch_tam_fit_reactive()
    req(m)
    
    it_fit <- if (inherits(m, "tam.jml")) {
      tf <- tryCatch(TAM::tam.fit(m)$fit.item, error = function(e) NULL)
      if (!is.null(tf)) {
        data.frame(
          item = as.character(tf$item),
          Outfit = tf$outfitItem,
          Outfit_t = tf$outfitItem_t,
          Infit = tf$infitItem,
          Infit_t = tf$infitItem_t,
          stringsAsFactors = FALSE
        )
      } else NULL
    } else {
      tryCatch(TAM::msq.itemfit(m)$itemfit, error = function(e) NULL)
    }
    
    # Strictly prioritize Rasch Item Difficulty Delta (xsi.item) over empirical mean score (item$M)
    betas <- if (!is.null(m$item$xsi.item)) m$item$xsi.item else if (!is.null(m$xsi) && is.data.frame(m$xsi)) m$xsi$xsi else if (!is.null(m$item$M)) m$item$M else rep(0, ncol(m$resp))
    ses <- calc_rasch_item_se(m)
    
    is_multi_mod <- (!is.null(m$ndim) && m$ndim > 1)
    Q_mat <- if (is_multi_mod) tryCatch(isolate(rasch_active_qmatrix()), error = function(e) NULL) else NULL
    get_dim_label <- function(it_name) {
      if (!is_multi_mod) return("Unidimensi (Dimensi 1)")
      if (!is.null(Q_mat) && it_name %in% rownames(Q_mat)) {
        d_assigned <- colnames(Q_mat)[which(Q_mat[it_name, ] == 1)]
        if (length(d_assigned) > 0) return(paste(d_assigned, collapse = ", "))
      }
      if (it_name %in% c("Item.1", "Item.4", "Item.8", "Item.11", "Item.13", "Item.16")) return("D1 (Pembatasan Peran)")
      if (it_name %in% c("Item.2", "Item.6", "Item.7", "Item.9", "Item.12", "Item.15", "Item.17")) return("D2 (Pelepasan Keterikatan)")
      if (it_name %in% c("Item.3", "Item.5", "Item.10", "Item.14")) return("D3 (Ketidaksepadanan Usaha)")
      "Dimensi 1"
    }
    
    if (is.null(it_fit)) {
      data.frame(
        Butir = colnames(m$resp),
        Dimensi = sapply(colnames(m$resp), get_dim_label),
        `Item_Difficulty_\\delta` = round(betas[1:ncol(m$resp)], 3),
        SE = round(ses[1:ncol(m$resp)], 3),
        Infit_MNSQ = 1.0,
        Infit_t = 0.0,
        Outfit_MNSQ = 1.0,
        Outfit_t = 0.0,
        Evaluasi_Linacre = "Fit",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      eval_fit <- sapply(it_fit$Infit, function(x) {
        if (is.na(x)) "N/A"
        else if (x >= 0.50 && x <= 1.50) "Produktif / Fit (0.5 - 1.5)"
        else if (x > 1.50) "Underfit / Noise (> 1.5)"
        else "Overfit / Redundan (< 0.5)"
      })
      
      data.frame(
        Butir = as.character(it_fit$item),
        Dimensi = sapply(as.character(it_fit$item), get_dim_label),
        `Item_Difficulty_\\delta` = round(betas[1:nrow(it_fit)], 3),
        SE = round(ses[1:nrow(it_fit)], 3),
        Infit_MNSQ = round(it_fit$Infit, 3),
        Infit_t = round(it_fit$Infit_t, 3),
        Outfit_MNSQ = round(it_fit$Outfit, 3),
        Outfit_t = round(it_fit$Outfit_t, 3),
        Evaluasi_Linacre = eval_fit,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
  })
  
  output$table_rasch_itemfit_full <- renderDT({
    df <- rasch_itemfit_df_reactive()
    req(df)
    datatable(
      df,
      options = list(
        pageLength = 50,
        lengthMenu = list(c(15, 25, 50, 100, -1), c("15", "25", "50", "100", "Semua Butir")),
        scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = "_all")
        )
      ),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })
  
  output$vb_rasch_item_fit_count <- renderText({
    d <- rasch_itemfit_df_reactive()
    if (is.null(d)) return("-")
    n_fit <- sum(d$Infit_MNSQ >= 0.50 & d$Infit_MNSQ <= 1.50, na.rm = TRUE)
    sprintf("%d / %d Butir (%.1f%%)", n_fit, nrow(d), (n_fit / nrow(d)) * 100)
  })
  output$vb_rasch_item_underfit_count <- renderText({
    d <- rasch_itemfit_df_reactive()
    if (is.null(d)) return("-")
    paste(sum(d$Infit_MNSQ > 1.50, na.rm = TRUE), "Butir")
  })
  output$vb_rasch_item_overfit_count <- renderText({
    d <- rasch_itemfit_df_reactive()
    if (is.null(d)) return("-")
    paste(sum(d$Infit_MNSQ < 0.50, na.rm = TRUE), "Butir")
  })
  output$vb_rasch_item_mean_beta <- renderText({
    d <- rasch_itemfit_df_reactive()
    if (is.null(d)) return("-")
    b_col <- if ("Item_Difficulty_\\delta" %in% colnames(d)) d$`Item_Difficulty_\\delta` else if ("Item_Difficulty_Delta" %in% colnames(d)) d$Item_Difficulty_Delta else d[[3]]
    sprintf("%.3f Logit", mean(b_col, na.rm = TRUE))
  })
  
  output$btn_dl_rasch_itemfit_csv <- downloadHandler(
    filename = function() paste0("Rasch_Item_Fit_", Sys.Date(), ".csv"),
    content = function(file) write.csv(rasch_itemfit_df_reactive(), file, row.names = FALSE)
  )
  
  output$rasch_itemfit_scatterplot <- renderPlot({
    d <- rasch_itemfit_df_reactive()
    req(d)
    
    par(mar = c(5, 5, 3, 2), bg = "#FFFFFF")
    plot(
      d$Infit_MNSQ, d$Outfit_MNSQ,
      pch = 21,
      bg = ifelse(d$Infit_MNSQ >= 0.5 & d$Infit_MNSQ <= 1.5 & d$Outfit_MNSQ >= 0.5 & d$Outfit_MNSQ <= 1.5,
                  "#650915", "#800020"),
      col = "#C5A059",
      cex = 1.8,
      xlab = "Item Infit MNSQ",
      ylab = "Item Outfit MNSQ",
      main = "Item Fit Diagnostic Map (Zona Hijau / Emas: 0.5 s/d 1.5 MNSQ Produktif)",
      xlim = c(min(0.2, min(d$Infit_MNSQ) - 0.1), max(2.0, max(d$Infit_MNSQ) + 0.2)),
      ylim = c(min(0.2, min(d$Outfit_MNSQ) - 0.1), max(2.0, max(d$Outfit_MNSQ) + 0.2)),
      las = 1
    )
    # Shaded acceptable zone (0.5 to 1.5)
    rect(0.5, 0.5, 1.5, 1.5, col = rgb(197/255, 160/255, 89/255, 0.20), border = "#C5A059", lty = 2)
    abline(h = c(0.5, 1.5), v = c(0.5, 1.5), col = "#C5A059", lty = 3)
    abline(h = 1.0, v = 1.0, col = "gray50", lty = 2)
    text(d$Infit_MNSQ, d$Outfit_MNSQ, labels = d$Butir, pos = 3, cex = 0.85, col = "#231618")
  })
  
  # --- Rasch Person Fit & Parameters ---
  rasch_personfit_df_reactive <- reactive({
    m <- rasch_tam_fit_reactive()
    req(m)
    
    wle_res <- safe_tam_wle(m)
    req(wle_res)
    p_fit <- calc_rasch_personfit(m, wle_res)
    
    # Ekstraksi Pola Jawaban Butir per Responden
    pola_jawaban <- if (!is.null(m$resp)) {
      apply(m$resp, 1, function(r) paste(ifelse(is.na(r), "-", as.character(r)), collapse = ""))
    } else {
      rep("-", nrow(wle_res))
    }
    
    thetas <- if (!is.null(wle_res$theta)) {
      wle_res$theta
    } else {
      dim_cols <- grep("^theta\\.Dim", colnames(wle_res))
      if (length(dim_cols) > 0) rowMeans(as.matrix(wle_res[, dim_cols, drop = FALSE]), na.rm = TRUE) else wle_res[, 1]
    }
    ses <- if (!is.null(wle_res$error)) {
      wle_res$error
    } else {
      err_cols <- grep("^error\\.Dim", colnames(wle_res))
      if (length(err_cols) > 0) {
        sqrt(rowSums(as.matrix(wle_res[, err_cols, drop = FALSE])^2, na.rm = TRUE)) / length(err_cols)
      } else {
        rep(NA_real_, length(thetas))
      }
    }
    
    if (is.null(p_fit)) {
      data.frame(
        Person_ID = paste0("P_", sprintf("%03d", 1:length(thetas))),
        Pola_Jawaban = pola_jawaban,
        `Kemampuan_\\theta` = round(thetas, 3),
        `SE_\\theta` = round(ses, 3),
        Infit_MNSQ = NA_real_,
        Infit_t = NA_real_,
        Outfit_MNSQ = NA_real_,
        Outfit_t = NA_real_,
        Status = "Tidak Dapat Diestimasi",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      out_df <- data.frame(
        Person_ID = paste0("P_", sprintf("%03d", 1:nrow(p_fit))),
        Pola_Jawaban = pola_jawaban,
        `Kemampuan_\\theta` = round(thetas, 3),
        `SE_\\theta` = round(ses, 3),
        Infit_MNSQ = round(p_fit$infitPerson, 3),
        Infit_t = round(p_fit$infitPerson_t, 3),
        Outfit_MNSQ = round(p_fit$outfitPerson, 3),
        Outfit_t = round(p_fit$outfitPerson_t, 3),
        Status = ifelse(is.na(p_fit$infitPerson) | is.na(p_fit$outfitPerson), "N/A",
                        ifelse(p_fit$infitPerson > 1.5 | p_fit$outfitPerson > 1.5, "Aberan / Tebakan", "Normal / Konsisten")),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      out_df
    }
  })
  
  output$table_rasch_personfit_full <- renderDT({
    d <- rasch_personfit_df_reactive()
    req(d)
    datatable(
      d,
      options = list(
        pageLength = 50,
        lengthMenu = list(c(15, 25, 50, 100, -1), c("15", "25", "50", "100", "Semua Responden")),
        scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = "_all"),
          list(className = "dt-center font-monospace fw-bold text-primary", targets = 1)
        )
      ),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })
  
  output$vb_rasch_person_mean_theta <- renderText({
    d <- rasch_personfit_df_reactive()
    if (is.null(d)) return("-")
    th_vec <- if ("Kemampuan_\\theta" %in% colnames(d)) d[["Kemampuan_\\theta"]] else if ("Kemampuan_θ" %in% colnames(d)) d[["Kemampuan_θ"]] else d[[3]]
    sprintf("%.3f Logit", mean(th_vec, na.rm = TRUE))
  })
  output$vb_rasch_person_sd_theta <- renderText({
    d <- rasch_personfit_df_reactive()
    if (is.null(d)) return("-")
    th_vec <- if ("Kemampuan_\\theta" %in% colnames(d)) d[["Kemampuan_\\theta"]] else if ("Kemampuan_θ" %in% colnames(d)) d[["Kemampuan_θ"]] else d[[3]]
    sprintf("%.3f", sd(th_vec, na.rm = TRUE))
  })
  output$vb_rasch_person_fit_pct <- renderText({
    d <- rasch_personfit_df_reactive()
    if (is.null(d)) return("-")
    n_norm <- sum(d$Status == "Normal / Konsisten", na.rm = TRUE)
    sprintf("%.1f%% (%d/%d)", (n_norm / nrow(d)) * 100, n_norm, nrow(d))
  })
  output$vb_rasch_person_misfit_pct <- renderText({
    d <- rasch_personfit_df_reactive()
    if (is.null(d)) return("-")
    n_mis <- sum(d$Status == "Aberan / Tebakan", na.rm = TRUE)
    sprintf("%.1f%% (%d/%d)", (n_mis / nrow(d)) * 100, n_mis, nrow(d))
  })
  
  output$btn_dl_rasch_personfit_csv <- downloadHandler(
    filename = function() paste0("Rasch_Person_Fit_", Sys.Date(), ".csv"),
    content = function(file) write.csv(rasch_personfit_df_reactive(), file, row.names = FALSE)
  )
  
  output$rasch_person_theta_dist_plot <- renderPlot({
    d <- rasch_personfit_df_reactive()
    req(d)
    
    th_vec <- if ("Kemampuan_\\theta" %in% colnames(d)) d[["Kemampuan_\\theta"]] else if ("Kemampuan_θ" %in% colnames(d)) d[["Kemampuan_θ"]] else d[[3]]
    thetas <- th_vec[is.finite(th_vec)]
    req(length(thetas) >= 2)
    par(mar = c(4.5, 4.5, 2.5, 1.5), bg = "#FFFFFF")
    
    h <- hist(thetas, breaks = 20, col = rgb(101/255, 9/255, 21/255, 0.70), border = "#C5A059",
              xlab = "Skala Kemampuan Laten Responden (θ Logit)", ylab = "Kepadatan Probabilitas",
              main = "Distribusi Kemampuan Responden (Person Ability θ)", freq = FALSE, las = 1)
    
    # Overlay theoretical normal
    x_seq <- seq(min(thetas) - 0.5, max(thetas) + 0.5, length.out = 200)
    lines(x_seq, dnorm(x_seq, mean = mean(thetas), sd = sd(thetas)), col = "#C5A059", lwd = 2.5, lty = 2)
    grid(col = "#E5D3B3")
  })
  
  # Sinkronisasi otomatis pilihan Wright Map berdasarkan dimensi model yang aktif
  observeEvent(rasch_tam_fit_reactive(), {
    m <- rasch_tam_fit_reactive()
    req(m)
    if (!is.null(m$ndim) && m$ndim > 1) {
      updateRadioButtons(session, "rasch_wm_view_mode", selected = "multi")
    } else {
      updateRadioButtons(session, "rasch_wm_view_mode", selected = "uni")
    }
  }, ignoreInit = TRUE)

  # --- Wright Map Rendering Engine (Uni & Multi) ---
  render_wright_map_internal <- function(custom_mod = NULL, custom_view_mode = NULL) {
    # Pastikan reset seluruh layar split.screen agar kanvas bersih 100%
    tryCatch(close.screen(all.screens = TRUE), error = function(e) NULL)
    on.exit(tryCatch(close.screen(all.screens = TRUE), error = function(e) NULL), add = TRUE)
    
    m <- if (!is.null(custom_mod)) custom_mod else rasch_tam_fit_reactive()
    req(m)
    
    wle_res <- safe_tam_wle(m)
    req(wle_res)
    view_mode <- if (!is.null(custom_view_mode)) {
      custom_view_mode
    } else if (!is.null(input$rasch_wm_view_mode)) {
      input$rasch_wm_view_mode
    } else {
      "uni"
    }
    wm_style <- if (!is.null(input$rasch_wm_style)) input$rasch_wm_style else "kawe"
    cex_val <- if (!is.null(input$rasch_wm_cex_scale)) input$rasch_wm_cex_scale else 0.95
    
    irt_type <- isolate(input$rasch_irt_model_type)
    
    thresholds <- if (irt_type %in% c("PCM", "PCM2", "RSM", "GPCM", "GPCM.groups", "JML") && max(m$resp, na.rm = TRUE) > 1) {
      tryCatch(TAM::tam.threshold(m), error = function(e) if (!is.null(m$item$xsi.item)) m$item$xsi.item else if (is.data.frame(m$xsi)) m$xsi$xsi else m$xsi)
    } else {
      if (!is.null(m$item$xsi.item)) m$item$xsi.item else if (is.data.frame(m$xsi)) m$xsi$xsi else m$xsi
    }
    
    item_diff <- if (!is.null(m$item$xsi.item)) m$item$xsi.item else if (!is.null(m$xsi$xsi)) m$xsi$xsi else rep(0, ncol(m$resp))
    names(item_diff) <- colnames(m$resp)
    
    is_multi <- (view_mode == "multi" && !is.null(m$ndim) && m$ndim > 1)
    
    if (is_multi) {
      theta_cols <- grep("^theta\\.Dim", colnames(wle_res))
      thetas_active <- as.matrix(wle_res[, theta_cols, drop = FALSE])
      colnames(thetas_active) <- paste("Dimensi", 1:ncol(thetas_active))
    } else {
      if (!is.null(wle_res$theta)) {
        thetas_active <- wle_res$theta
      } else {
        theta_cols <- grep("^theta\\.Dim", colnames(wle_res))
        if (length(theta_cols) > 0) {
          thetas_active <- rowMeans(as.matrix(wle_res[, theta_cols, drop = FALSE]), na.rm = TRUE)
        } else {
          thetas_active <- wle_res[[1]]
        }
      }
    }
    
    n_dim_count <- if (is.matrix(thetas_active)) ncol(thetas_active) else 1
    calculated_item_prop <- if (is_multi) max(0.60, 1 - (n_dim_count * 0.10)) else 0.82
    
    if (wm_style == "kawe") {
      itemSideCustom <- function(thr, yRange = NULL, axis.items = "Butir Instrumen", 
                                 axis.logits = "Logit Trait (Theta)", show.axis.logits = "R", 
                                 axis.logits.cex = cex_val * 1.05, oma = c(0, 0, 0, 0), 
                                 label.items = NULL, label.items.cex = cex_val * 0.95, ...) {
        thr <- as.matrix(thr)
        nI <- nrow(thr)
        nL <- ncol(thr)
        item_diff_clean <- if (ncol(thr) == 1) as.vector(thr) else rowMeans(thr, na.rm = TRUE)
        if (is.null(yRange)) yRange <- range(c(thr, item_diff_clean), na.rm = TRUE)
        if (is.null(label.items)) label.items <- rownames(thr)
        
        # Kunci margin luar (oma) dan margin dalam (mar) identik dengan personSide
        par(oma = c(0, 0, 0, 0), mar = c(5.1, 0.2, 4.1, 0.1), mgp = c(2.2, 1, 0))
        plot(seq_len(nI), rep(0, nI), type = "n", axes = FALSE, xlab = "", 
             ylab = "", ylim = yRange, xlim = c(0.5, nI + 0.5), font.lab = 2)
        box(bty = "o", col = "gray30", lwd = 1.5)
        abline(h = 0, col = "gray70", lty = 2, lwd = 1.2)
        
        if (show.axis.logits %in% c("R", TRUE)) {
          axis(4, las = 1, cex.axis = axis.logits.cex, font.axis = 2, col = "gray30", lwd = 1.2)
          mtext(axis.logits, side = 4, line = 2.4, cex = axis.logits.cex, font = 2)
        }
        
        # Backbone line
        for (i in 1:nI) {
          y_vals <- c(thr[i, !is.na(thr[i, ])], item_diff[i])
          if (length(y_vals) > 0) lines(c(i, i), range(y_vals), col = "gray55", lwd = 2.4, lty = 1)
        }
        
        # Threshold horizontal colored lines
        cat_cols <- c("#C0392B", "#E67E22", "#27AE60", "#2980B9", "#8E44AD", "#16A085")
        for (i in 1:nI) {
          for (k in 1:min(nL, length(cat_cols))) {
            val <- thr[i, k]
            if (!is.na(val)) {
              lines(c(i - 0.32, i + 0.32), c(val, val), col = cat_cols[k], lwd = 3.2)
              text(i + 0.34, val, labels = paste0("T", k), col = cat_cols[k], cex = cex_val * 0.75, font = 2, adj = c(0, 0.5))
            }
          }
        }
        
        # Item Difficulty Delta points (Solid Orange with Dark Border)
        for (i in 1:nI) {
          diff_val <- if (!is.na(item_diff[i])) item_diff[i] else item_diff_clean[i]
          points(i, diff_val, pch = 21, bg = "#F39C12", col = "#1C2833", cex = cex_val * 2.2, lwd = 2.2)
        }
        
        # Bottom Item Labels
        text(seq_len(nI), y = par("usr")[3] - 0.16, labels = label.items, 
             srt = 45, adj = c(1, 1), xpd = TRUE, cex = label.items.cex, font = 2)
        axis(1, at = seq_len(nI), labels = FALSE, line = NA, tcl = -0.35, col = "gray30", lwd = 1.2)
        mtext(axis.items, side = 1, line = 3.6, cex = cex_val * 1.15, font = 2)
      }
      
      personSideCustom <- function(thetas_in, yRange = NULL, oma = c(0, 0, 0, 0), ...) {
        thetas_in <- as.matrix(thetas_in)
        nD <- ncol(thetas_in)
        if (is.null(yRange)) yRange <- range(thetas_in, na.rm = TRUE)
        
        # Netralkan injeksi par(oma = c(0,0,3,0)) dan margin agar tidak memicu figure margins too large
        par(oma = c(0, 0, 0, 0), mar = c(5.1, 0.2, 4.1, 0.1))
        
        old.screens <- split.screen()
        on.exit({
          curr.screens <- split.screen()
          new.screens <- curr.screens[!(curr.screens %in% old.screens)]
          if (length(new.screens) > 0) tryCatch(close.screen(new.screens), error = function(e) NULL)
        }, add = TRUE)
        split.screen(c(1, nD), erase = FALSE)
        
        bar_cols <- c("#FADBD8", "#D4E6F1", "#D5F5E3", "#FCF3CF", "#E8DAEF", "#D1F2EB")
        border_cols <- c("#C0392B", "#2980B9", "#27AE60", "#D68910", "#8E44AD", "#16A085")
        
        dim_names <- colnames(thetas_in)
        if (is.null(dim_names)) dim_names <- if (nD == 1) "" else paste("Dim", 1:nD)
        
        for (i in 1:nD) {
          # Margin dalam identik 100% dengan sisi butir (itemSideCustom)
          par(mar = c(5.1, 0.2, 4.1, 0.1))
          vals <- thetas_in[, i]
          b_min <- min(c(yRange[1], vals)) - 0.1
          b_max <- max(c(yRange[2], vals)) + 0.1
          brks <- seq(b_min, b_max, length.out = 24)
          h <- hist(vals, breaks = brks, plot = FALSE)
          
          plot(c(0, max(h$counts) * 1.15), yRange, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(max(h$counts) * 1.15, 0))
          box(bty = "c", col = "gray30", lwd = 1.5)
          abline(h = 0, col = "gray70", lty = 2, lwd = 1.2)
          
          for (j in 1:length(h$mids)) {
            ybot <- h$breaks[j]
            ytop <- h$breaks[j + 1]
            xright <- h$counts[j]
            rect(xright, ybot, 0, ytop, col = bar_cols[((i - 1) %% length(bar_cols)) + 1], 
                 border = border_cols[((i - 1) %% length(border_cols)) + 1], lwd = 1.4)
          }
          
          if (nzchar(dim_names[i])) mtext(dim_names[i], side = 1, line = 1.4, cex = cex_val * 0.85, font = 2)
          if (i == 1) mtext("Frekuensi", side = 2, line = 1.6, cex = cex_val * 0.95, font = 3)
          if (i < nD) screen(screen() + 1)
        }
        
        curr.screens <- split.screen()
        new.screens <- curr.screens[!(curr.screens %in% old.screens)]
        if (length(new.screens) > 0) tryCatch(close.screen(new.screens), error = function(e) NULL)
      }
      
      WrightMap::wrightMap(
        thetas = thetas_active,
        thresholds = thresholds,
        item.side = itemSideCustom,
        person.side = personSideCustom,
        item.prop = calculated_item_prop,
        main.title = ifelse(is_multi, "Multidimensional Wright Map (Gaya Akademis KAWE)", "Unidimensional Wright Map (Gaya Akademis KAWE)"),
        show.thr.lab = FALSE,
        show.thr.sym = FALSE
      )
    } else {
      # Standard WrightMap
      WrightMap::wrightMap(
        thetas = thetas_active,
        thresholds = thresholds,
        item.prop = calculated_item_prop,
        main.title = ifelse(is_multi, "Multidimensional Person-Item Wright Map (TAM)", "Unidimensional Person-Item Wright Map (TAM)"),
        cex.label = cex_val * 0.9,
        cex.axis = cex_val * 0.8
      )
    }
  }
  
  output$rasch_wright_map_display <- renderPlot({
    render_wright_map_internal()
  }, res = 120)
  
  output$btn_dl_wright_png <- downloadHandler(
    filename = function() paste0("Wright_Map_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 3200, height = 2400, res = 300)
      render_wright_map_internal()
      dev.off()
    }
  )
  output$btn_dl_wright_pdf <- downloadHandler(
    filename = function() paste0("Wright_Map_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 12, height = 9)
      render_wright_map_internal()
      dev.off()
    }
  )

  # ============================================================================
  # RASCH RE-ANALYSIS & PURIFICATION SERVER LOGIC (ITEM & PERSON MISFIT)
  # ============================================================================
  
  # Observasi Pilihan Butir Manual untuk Re-Analisis
  observe({
    df <- rasch_data_store()
    req(df)
    items <- colnames(df)
    updateSelectizeInput(session, "rasch_re_manual_drop_items", choices = items, server = TRUE)
  })

  # Navigasi Cepat dari Tab Item Fit & Person Fit
  observeEvent(input$btn_goto_rasch_reanalysis_from_item, {
    nav_select("main_navbar", "tab_rasch_reanalysis")
  })
  observeEvent(input$btn_goto_rasch_reanalysis_from_person, {
    nav_select("main_navbar", "tab_rasch_reanalysis")
  })

  # Eksekusi Analisis Ulang Rasch (Pemurnian Butir & Responden)
  observeEvent(input$btn_run_rasch_reanalysis, {
    m_base <- rasch_tam_fit_reactive()
    if (is.null(m_base)) {
      showNotification("Model dasar Rasch belum diestimasi. Silakan estimasi model pada Tab Model terlebih dahulu.", type = "warning", duration = 5)
      return()
    }
    
    df_raw <- rasch_data_store()
    all_items <- colnames(m_base$resp)
    n_total_persons <- nrow(m_base$resp)
    
    # 1. Identifikasi Butir Misfit
    dropped_items <- character(0)
    drop_item_reasons <- character(0)
    drop_item_dims <- character(0)
    drop_item_infits <- numeric(0)
    drop_item_outfits <- numeric(0)
    
    it_df <- tryCatch(rasch_itemfit_df_reactive(), error = function(e) NULL)
    
    if (isTRUE(input$rasch_re_filter_items) && !is.null(it_df) && nrow(it_df) > 0) {
      up_th <- input$rasch_re_item_upper_mnsq
      lo_th <- input$rasch_re_item_lower_mnsq
      crit <- input$rasch_re_item_criteria
      
      for (i in seq_len(nrow(it_df))) {
        it_name <- as.character(it_df$Butir[i])
        in_val <- it_df$Infit_MNSQ[i]
        out_val <- it_df$Outfit_MNSQ[i]
        dim_val <- if ("Dimensi" %in% colnames(it_df)) as.character(it_df$Dimensi[i]) else "-"
        
        is_under <- switch(crit,
          "infit" = (!is.na(in_val) && in_val > up_th),
          "outfit" = (!is.na(out_val) && out_val > up_th),
          "either" = ((!is.na(in_val) && in_val > up_th) || (!is.na(out_val) && out_val > up_th))
        )
        is_over <- switch(crit,
          "infit" = (!is.na(in_val) && in_val < lo_th),
          "outfit" = (!is.na(out_val) && out_val < lo_th),
          "either" = ((!is.na(in_val) && in_val < lo_th) || (!is.na(out_val) && out_val < lo_th))
        )
        
        if (isTRUE(is_under)) {
          dropped_items <- c(dropped_items, it_name)
          drop_item_dims <- c(drop_item_dims, dim_val)
          drop_item_infits <- c(drop_item_infits, in_val)
          drop_item_outfits <- c(drop_item_outfits, out_val)
          drop_item_reasons <- c(drop_item_reasons, sprintf("Underfit / Noise (Infit=%.2f, Outfit=%.2f > %.2f)", in_val, out_val, up_th))
        } else if (isTRUE(is_over)) {
          dropped_items <- c(dropped_items, it_name)
          drop_item_dims <- c(drop_item_dims, dim_val)
          drop_item_infits <- c(drop_item_infits, in_val)
          drop_item_outfits <- c(drop_item_outfits, out_val)
          drop_item_reasons <- c(drop_item_reasons, sprintf("Overfit / Redundan (Infit=%.2f, Outfit=%.2f < %.2f)", in_val, out_val, lo_th))
        }
      }
    }
    
    # Tambahan Butir Manual
    man_drops <- input$rasch_re_manual_drop_items
    if (length(man_drops) > 0) {
      for (md in man_drops) {
        if (!(md %in% dropped_items) && md %in% all_items) {
          dropped_items <- c(dropped_items, md)
          dim_val <- if (!is.null(it_df) && md %in% it_df$Butir) as.character(it_df$Dimensi[it_df$Butir == md][1]) else "-"
          in_v <- if (!is.null(it_df) && md %in% it_df$Butir) it_df$Infit_MNSQ[it_df$Butir == md][1] else NA_real_
          out_v <- if (!is.null(it_df) && md %in% it_df$Butir) it_df$Outfit_MNSQ[it_df$Butir == md][1] else NA_real_
          drop_item_dims <- c(drop_item_dims, dim_val)
          drop_item_infits <- c(drop_item_infits, in_v)
          drop_item_outfits <- c(drop_item_outfits, out_v)
          drop_item_reasons <- c(drop_item_reasons, "Pilihan Manual Peneliti")
        }
      }
    }
    
    retained_items <- setdiff(all_items, dropped_items)
    
    # 2. Identifikasi Responden Aberan
    dropped_persons <- integer(0)
    drop_p_ids <- character(0)
    drop_p_patterns <- character(0)
    drop_p_thetas <- numeric(0)
    drop_p_infits <- numeric(0)
    drop_p_outfits <- numeric(0)
    drop_p_reasons <- character(0)
    
    p_df <- tryCatch(rasch_personfit_df_reactive(), error = function(e) NULL)
    
    if (isTRUE(input$rasch_re_filter_persons) && !is.null(p_df) && nrow(p_df) > 0) {
      p_th <- input$rasch_re_person_mnsq_thresh
      p_crit <- input$rasch_re_person_criteria
      
      th_col <- if ("Kemampuan_\\theta" %in% colnames(p_df)) "Kemampuan_\\theta" else if ("Kemampuan_θ" %in% colnames(p_df)) "Kemampuan_θ" else 3
      
      for (p in seq_len(nrow(p_df))) {
        pin <- p_df$Infit_MNSQ[p]
        pout <- p_df$Outfit_MNSQ[p]
        if (is.na(pin) || is.na(pout)) next
        
        is_aberrant <- switch(p_crit,
          "infit" = (pin > p_th),
          "outfit" = (pout > p_th),
          "either" = (pin > p_th || pout > p_th)
        )
        
        if (isTRUE(is_aberrant)) {
          dropped_persons <- c(dropped_persons, p)
          drop_p_ids <- c(drop_p_ids, as.character(p_df$Person_ID[p]))
          drop_p_patterns <- c(drop_p_patterns, as.character(p_df$Pola_Jawaban[p]))
          drop_p_thetas <- c(drop_p_thetas, round(as.numeric(p_df[[th_col]][p]), 3))
          drop_p_infits <- c(drop_p_infits, round(pin, 3))
          drop_p_outfits <- c(drop_p_outfits, round(pout, 3))
          
          r_desc <- if (pin > p_th && pout > p_th) {
            sprintf("Infit (%.2f) & Outfit (%.2f) > %.2f [Aberansi Total / Satisficing]", pin, pout, p_th)
          } else if (pout > p_th) {
            sprintf("Outfit (%.2f) > %.2f [Tebakan Tak Terduga / Lucky Guessing]", pout, p_th)
          } else {
            sprintf("Infit (%.2f) > %.2f [Ketidakkonsistenan Pola Respons Terdekat]", pin, p_th)
          }
          drop_p_reasons <- c(drop_p_reasons, r_desc)
        }
      }
    }
    
    # Mapping person index to original dataframe row index
    act_p_curr <- rasch_active_persons()
    orig_person_indices <- if (!is.null(act_p_curr)) act_p_curr else seq_len(nrow(df_raw))
    
    dropped_orig_indices <- if (length(dropped_persons) > 0) orig_person_indices[dropped_persons] else integer(0)
    retained_orig_indices <- setdiff(orig_person_indices, dropped_orig_indices)
    
    # Validasi Ambang Batas Jumlah Butir & Responden
    if (length(retained_items) < 3) {
      showNotification("Gagal: Terlalu banyak butir yang dieliminasi. Minimal 3 butir harus dipertahankan untuk kalibrasi Rasch.", type = "error", duration = 6)
      return()
    }
    if (length(retained_orig_indices) < 30) {
      showNotification("Gagal: Terlalu banyak responden yang dieliminasi. Minimal 30 responden harus dipertahankan untuk kalibrasi Rasch.", type = "error", duration = 6)
      return()
    }
    
    # 3. Estimasi Model TAM Pemurnian
    withProgress(message = "Kalibrasi Ulang TAM Sedang Berjalan...", detail = "Mengestimasi model dengan data dan butir termurnikan...", value = 0.3, {
      sub_mat <- as.matrix(df_raw[retained_orig_indices, retained_items, drop = FALSE])
      min_val <- min(sub_mat, na.rm = TRUE)
      if (min_val > 0) sub_mat <- sub_mat - min_val
      
      dim_spec <- isolate(input$rasch_model_dim_spec)
      irt_type <- isolate(input$rasch_irt_model_type)
      
      Q_sub <- NULL
      if (dim_spec == "multi") {
        Q_raw <- isolate(rasch_active_qmatrix())
        if (!is.null(Q_raw)) {
          common_q <- intersect(rownames(Q_raw), retained_items)
          if (length(common_q) >= 2) {
            Q_cand <- Q_raw[common_q, , drop = FALSE]
            v_dims <- which(colSums(Q_cand) > 0)
            if (length(v_dims) >= 2) {
              Q_sub <- Q_cand[, v_dims, drop = FALSE]
            }
          }
        }
      }
      
      n_dims_cur <- if (!is.null(Q_sub)) ncol(Q_sub) else 1
      ctrl_opts <- list(maxiter = 35)
      if (n_dims_cur >= 3) {
        ctrl_opts$nodes <- seq(-3, 3, length.out = 8)
      } else if (n_dims_cur == 2) {
        ctrl_opts$nodes <- seq(-3.5, 3.5, length.out = 12)
      }
      
      setProgress(value = 0.6, detail = "Menjalankan iterasi EM pada model pemurnian...")
      
      mod_purified <- tryCatch({
        if (irt_type %in% c("1PL", "PCM", "PCM2", "RSM")) {
          if (!is.null(Q_sub)) TAM::tam.mml(resp = sub_mat, Q = Q_sub, irtmodel = irt_type, verbose = FALSE, control = ctrl_opts)
          else TAM::tam.mml(resp = sub_mat, irtmodel = irt_type, verbose = FALSE, control = ctrl_opts)
        } else if (irt_type %in% c("2PL", "GPCM", "2PL.groups", "GPCM.groups")) {
          if (!is.null(Q_sub)) TAM::tam.mml.2pl(resp = sub_mat, Q = Q_sub, irtmodel = irt_type, verbose = FALSE, control = ctrl_opts)
          else TAM::tam.mml.2pl(resp = sub_mat, irtmodel = irt_type, verbose = FALSE, control = ctrl_opts)
        } else if (irt_type == "3PL") {
          TAM::tam.mml.3pl(resp = sub_mat, Q = Q_sub, est.guess = 1:ncol(sub_mat), verbose = FALSE, control = list(maxiter = 30))
        } else if (irt_type == "JML") {
          TAM::tam.jml(resp = sub_mat, verbose = FALSE, control = list(maxiter = 30))
        } else {
          TAM::tam.mml(resp = sub_mat, irtmodel = "1PL", verbose = FALSE, control = ctrl_opts)
        }
      }, error = function(e) NULL)
      
      if (is.null(mod_purified)) {
        showNotification("Estimasi model pemurnian tidak konvergen. Periksa kembali subset butir dan responden.", type = "error", duration = 6)
        return()
      }
      
      setProgress(value = 0.85, detail = "Mengekstraksi metrik kesesuaian dan diagnostik fit...")
      
      # 4. Ekstraksi Metrik Model Awal vs Pemurnian
      # Model Awal
      k_orig <- ncol(m_base$resp)
      n_orig <- nrow(m_base$resp)
      dev_orig <- if (!is.null(m_base$deviance)) m_base$deviance else NA_real_
      aic_orig <- if (!is.null(m_base$aic)) m_base$aic else NA_real_
      bic_orig <- if (!is.null(m_base$bic)) m_base$bic else NA_real_
      eap_orig <- if (!is.null(m_base$EAP.rel)) mean(m_base$EAP.rel, na.rm = TRUE) else NA_real_
      wle_base <- safe_tam_wle(m_base)
      wle_rel_orig <- if (!is.null(wle_base$WLE.rel)) mean(wle_base$WLE.rel, na.rm = TRUE) else NA_real_
      
      mean_infit_it_orig <- if (!is.null(it_df)) mean(it_df$Infit_MNSQ, na.rm = TRUE) else NA_real_
      mean_outfit_it_orig <- if (!is.null(it_df)) mean(it_df$Outfit_MNSQ, na.rm = TRUE) else NA_real_
      pct_it_fit_orig <- if (!is.null(it_df)) mean(it_df$Infit_MNSQ >= 0.50 & it_df$Infit_MNSQ <= 1.50, na.rm = TRUE) * 100 else NA_real_
      pct_p_fit_orig <- if (!is.null(p_df)) mean(p_df$Infit_MNSQ <= 1.50 & p_df$Outfit_MNSQ <= 1.50, na.rm = TRUE) * 100 else NA_real_
      
      # Model Purified
      k_pur <- ncol(mod_purified$resp)
      n_pur <- nrow(mod_purified$resp)
      dev_pur <- if (!is.null(mod_purified$deviance)) mod_purified$deviance else NA_real_
      aic_pur <- if (!is.null(mod_purified$aic)) mod_purified$aic else NA_real_
      bic_pur <- if (!is.null(mod_purified$bic)) mod_purified$bic else NA_real_
      eap_pur <- if (!is.null(mod_purified$EAP.rel)) mean(mod_purified$EAP.rel, na.rm = TRUE) else NA_real_
      wle_pur <- safe_tam_wle(mod_purified)
      wle_rel_pur <- if (!is.null(wle_pur$WLE.rel)) mean(wle_pur$WLE.rel, na.rm = TRUE) else NA_real_
      
      it_pur_fit <- if (inherits(mod_purified, "tam.jml")) tryCatch(TAM::tam.fit(mod_purified)$fit.item, error = function(e) NULL) else tryCatch(TAM::msq.itemfit(mod_purified)$itemfit, error = function(e) NULL)
      in_it_pur <- if (!is.null(it_pur_fit$infitItem)) it_pur_fit$infitItem else if (!is.null(it_pur_fit$Infit)) it_pur_fit$Infit else rep(1, k_pur)
      out_it_pur <- if (!is.null(it_pur_fit$outfitItem)) it_pur_fit$outfitItem else if (!is.null(it_pur_fit$Outfit)) it_pur_fit$Outfit else rep(1, k_pur)
      
      mean_infit_it_pur <- mean(in_it_pur, na.rm = TRUE)
      mean_outfit_it_pur <- mean(out_it_pur, na.rm = TRUE)
      pct_it_fit_pur <- mean(in_it_pur >= 0.50 & in_it_pur <= 1.50, na.rm = TRUE) * 100
      
      p_pur_fit <- calc_rasch_personfit(mod_purified, wle_pur)
      pct_p_fit_pur <- if (!is.null(p_pur_fit)) mean(p_pur_fit$infitPerson <= 1.50 & p_pur_fit$outfitPerson <= 1.50, na.rm = TRUE) * 100 else NA_real_
      
      # 5. Tabel Komparasi
      comp_df <- data.frame(
        Indeks_Fit = c(
          "Jumlah Butir Aktif (k)",
          "Jumlah Responden (N)",
          "Deviance (-2LL)",
          "AIC (Akaike Information Criterion)",
          "BIC (Bayesian Information Criterion)",
          "Reliabilitas Laten (EAP)",
          "Reliabilitas Laten (WLE)",
          "Rata-rata Item Infit MNSQ",
          "Rata-rata Item Outfit MNSQ",
          "Persentase Butir Produktif / Fit (%)",
          "Persentase Responden Fit Normal (%)"
        ),
        Model_Awal = c(
          sprintf("%d Butir", k_orig),
          sprintf("%d Responden", n_orig),
          sprintf("%.2f", dev_orig),
          sprintf("%.2f", aic_orig),
          sprintf("%.2f", bic_orig),
          sprintf("%.3f", eap_orig),
          sprintf("%.3f", wle_rel_orig),
          sprintf("%.3f", mean_infit_it_orig),
          sprintf("%.3f", mean_outfit_it_orig),
          sprintf("%.1f%%", pct_it_fit_orig),
          sprintf("%.1f%%", pct_p_fit_orig)
        ),
        Model_Setelah_Pemurnian = c(
          sprintf("%d Butir", k_pur),
          sprintf("%d Responden", n_pur),
          sprintf("%.2f", dev_pur),
          sprintf("%.2f", aic_pur),
          sprintf("%.2f", bic_pur),
          sprintf("%.3f", eap_pur),
          sprintf("%.3f", wle_rel_pur),
          sprintf("%.3f", mean_infit_it_pur),
          sprintf("%.3f", mean_outfit_it_pur),
          sprintf("%.1f%%", pct_it_fit_pur),
          sprintf("%.1f%%", pct_p_fit_pur)
        ),
        Perubahan_Delta = c(
          sprintf("%+d Butir", k_pur - k_orig),
          sprintf("%+d Responden", n_pur - n_orig),
          sprintf("%+.2f", dev_pur - dev_orig),
          sprintf("%+.2f", aic_pur - aic_orig),
          sprintf("%+.2f", bic_pur - bic_orig),
          sprintf("%+.3f", eap_pur - eap_orig),
          sprintf("%+.3f", wle_rel_pur - wle_rel_orig),
          sprintf("%+.3f", mean_infit_it_pur - mean_infit_it_orig),
          sprintf("%+.3f", mean_outfit_it_pur - mean_outfit_it_orig),
          sprintf("%+.1f%%", pct_it_fit_pur - pct_it_fit_orig),
          sprintf("%+.1f%%", pct_p_fit_pur - pct_p_fit_orig)
        ),
        Evaluasi_Psikometris = c(
          "Penyaringan butir non-produktif (noise / redundansi)",
          "Pembersihan responden aberan memitigasi distorsi parameter",
          "Penurunan deviasi mencerminkan konsistensi kalibrasi",
          "Penurunan AIC mengindikasikan efisiensi parsimoni model",
          "Penurunan BIC mengonfirmasi kecocokan di bawah penalti kompleksitas",
          "Presisi konsistensi estimasi kemampuan laten (EAP)",
          "Presisi estimasi kemampuan tertimbang tanpa bias (WLE)",
          "Mendekati nilai ekspektasi teoretis 1.00 (model fit optimal)",
          "Pereduksian outliers dan respons tak terduga pada butir",
          "Peningkatan proporsi butir produktif dalam batas Linacre (2002)",
          "Peningkatan integritas kelompok partisipan bebas tebakan acak"
        ),
        stringsAsFactors = FALSE
      )
      colnames(comp_df)[4] <- "Perubahan (\\Delta)"
      
      # 6. Tabel Elemen Tereliminasi
      dropped_items_df <- if (length(dropped_items) > 0) {
        data.frame(
          No = seq_along(dropped_items),
          Butir = dropped_items,
          Dimensi = drop_item_dims,
          Infit_MNSQ = round(drop_item_infits, 3),
          Outfit_MNSQ = round(drop_item_outfits, 3),
          Alasan_Eliminasi = drop_item_reasons,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          Status = "Tidak ada butir yang dieliminasi. Seluruh butir memenuhi kriteria fit produktif (0.50 s/d 1.50 MNSQ).",
          stringsAsFactors = FALSE
        )
      }
      
      dropped_persons_df <- if (length(dropped_persons) > 0) {
        data.frame(
          No = seq_along(dropped_persons),
          Person_ID = drop_p_ids,
          Pola_Jawaban = drop_p_patterns,
          `Kemampuan_\\theta` = drop_p_thetas,
          Infit_MNSQ = drop_p_infits,
          Outfit_MNSQ = drop_p_outfits,
          Alasan_Eliminasi = drop_p_reasons,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          Status = "Tidak ada responden yang dieliminasi. Seluruh responden memenuhi batas toleransi kesesuaian fit (MNSQ <= ambang batas).",
          stringsAsFactors = FALSE
        )
      }
      
      # 7. Tabel Pergeseran Parameter Butir (Item Difficulty & Infit/Outfit Shifts)
      item_diff_base <- if (!is.null(m_base$item$xsi.item)) m_base$item$xsi.item else if (!is.null(m_base$xsi) && is.data.frame(m_base$xsi)) m_base$xsi$xsi else rep(0, ncol(m_base$resp))
      names(item_diff_base) <- colnames(m_base$resp)
      
      item_diff_pur <- if (!is.null(mod_purified$item$xsi.item)) mod_purified$item$xsi.item else if (!is.null(mod_purified$xsi) && is.data.frame(mod_purified$xsi)) mod_purified$xsi$xsi else rep(0, ncol(mod_purified$resp))
      names(item_diff_pur) <- colnames(mod_purified$resp)
      
      infit_pur_named <- in_it_pur
      names(infit_pur_named) <- if (!is.null(it_pur_fit$item)) as.character(it_pur_fit$item) else colnames(mod_purified$resp)
      outfit_pur_named <- out_it_pur
      names(outfit_pur_named) <- if (!is.null(it_pur_fit$item)) as.character(it_pur_fit$item) else colnames(mod_purified$resp)
      
      it_df_indexed <- it_df
      if (!is.null(it_df_indexed) && "Butir" %in% colnames(it_df_indexed)) {
        rownames(it_df_indexed) <- as.character(it_df_indexed$Butir)
      }
      
      items_comp_list <- lapply(retained_items, function(it) {
        d_orig <- if (it %in% names(item_diff_base)) item_diff_base[it] else NA_real_
        d_pur  <- if (it %in% names(item_diff_pur)) item_diff_pur[it] else NA_real_
        
        in_orig <- if (!is.null(it_df_indexed) && it %in% rownames(it_df_indexed)) it_df_indexed[it, "Infit_MNSQ"] else NA_real_
        in_pur_val <- if (it %in% names(infit_pur_named)) infit_pur_named[it] else NA_real_
        
        out_orig <- if (!is.null(it_df_indexed) && it %in% rownames(it_df_indexed)) it_df_indexed[it, "Outfit_MNSQ"] else NA_real_
        out_pur_val <- if (it %in% names(outfit_pur_named)) outfit_pur_named[it] else NA_real_
        
        dim_val <- if (!is.null(it_df_indexed) && it %in% rownames(it_df_indexed) && "Dimensi" %in% colnames(it_df_indexed)) as.character(it_df_indexed[it, "Dimensi"]) else "-"
        
        eval_str <- if (!is.na(in_pur_val) && !is.na(out_pur_val)) {
          if (in_pur_val >= 0.70 && in_pur_val <= 1.30 && out_pur_val >= 0.70 && out_pur_val <= 1.30) {
            "Optimal Fit (0.70 - 1.30 MNSQ)"
          } else if (in_pur_val >= 0.50 && in_pur_val <= 1.50 && out_pur_val >= 0.50 && out_pur_val <= 1.50) {
            "Produktif Fit (0.50 - 1.50 MNSQ)"
          } else if (in_pur_val > 1.50 || out_pur_val > 1.50) {
            "Underfit Residual"
          } else {
            "Overfit Residual"
          }
        } else "-"
        
        data.frame(
          Butir = it,
          Dimensi = dim_val,
          `Kesulitan_Awal_\\delta` = round(d_orig, 3),
          `Kesulitan_Baru_\\delta` = round(d_pur, 3),
          `\\Delta_Kesulitan` = round(d_pur - d_orig, 3),
          Infit_Awal = round(in_orig, 3),
          Infit_Baru = round(in_pur_val, 3),
          `\\Delta_Infit` = round(in_pur_val - in_orig, 3),
          Outfit_Awal = round(out_orig, 3),
          Outfit_Baru = round(out_pur_val, 3),
          `\\Delta_Outfit` = round(out_pur_val - out_orig, 3),
          Kategori_Fit_Baru = eval_str,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      })
      items_comp_df <- if (length(items_comp_list) > 0) do.call(rbind, items_comp_list) else data.frame()
      
      res <- list(
        mod_base = m_base,
        mod_purified = mod_purified,
        comp_table = comp_df,
        dropped_items_df = dropped_items_df,
        dropped_persons_df = dropped_persons_df,
        items_comp_df = items_comp_df,
        retained_items = retained_items,
        retained_persons = retained_orig_indices,
        k_orig = k_orig,
        k_pur = k_pur,
        n_orig = n_orig,
        n_pur = n_pur,
        eap_pur = eap_pur,
        pct_p_fit_pur = pct_p_fit_pur,
        delta_k = k_pur - k_orig,
        delta_n = n_pur - n_orig,
        delta_eap = eap_pur - eap_orig,
        delta_p_fit = pct_p_fit_pur - pct_p_fit_orig
      )
      
      rasch_reanalysis_results(res)
      setProgress(value = 1.0, detail = "Selesai!")
      showNotification("Analisis ulang Rasch berhasil diselesaikan! Periksa hasil komparasi di bawah.", type = "message", duration = 5)
    })
  })

  # Terapkan Model Pemurnian ke Seluruh Tab
  observeEvent(input$btn_apply_rasch_reanalysis, {
    res <- rasch_reanalysis_results()
    if (is.null(res)) {
      showNotification("Silakan jalankan analisis ulang terlebih dahulu sebelum menerapkan model.", type = "warning", duration = 4)
      return()
    }
    rasch_active_items(res$retained_items)
    rasch_active_persons(res$retained_persons)
    rasch_fit_trigger(rasch_fit_trigger() + 1)
    showNotification(
      sprintf("Model pemurnian berhasil diterapkan (%d butir, %d responden)! Seluruh tab Item Fit, Person Fit, dan Wright Map kini menampilkan model hasil pemurnian.",
              length(res$retained_items), length(res$retained_persons)),
      type = "message",
      duration = 6
    )
  })

  # Reset Model ke Data Lengkap Awal
  observeEvent(input$btn_reset_rasch_reanalysis, {
    rasch_active_items(NULL)
    rasch_active_persons(NULL)
    rasch_reanalysis_results(NULL)
    rasch_fit_trigger(rasch_fit_trigger() + 1)
    showNotification("Seluruh filter telah di-reset. Model Rasch kembali menggunakan seluruh butir dan responden lengkap awal.", type = "message", duration = 5)
  })

  # Value Boxes Output Re-Analisis
  output$vb_rasch_re_items <- renderText({
    res <- rasch_reanalysis_results()
    if (is.null(res)) {
      m <- rasch_tam_fit_reactive()
      if (is.null(m)) return("-")
      return(sprintf("%d Butir", ncol(m$resp)))
    }
    sprintf("%d Butir (%+d)", res$k_pur, res$delta_k)
  })

  output$vb_rasch_re_persons <- renderText({
    res <- rasch_reanalysis_results()
    if (is.null(res)) {
      m <- rasch_tam_fit_reactive()
      if (is.null(m)) return("-")
      return(sprintf("%d Responden", nrow(m$resp)))
    }
    sprintf("%d Responden (%+d)", res$n_pur, res$delta_n)
  })

  output$vb_rasch_re_eap <- renderText({
    res <- rasch_reanalysis_results()
    if (is.null(res)) {
      m <- rasch_tam_fit_reactive()
      if (is.null(m) || is.null(m$EAP.rel)) return("-")
      return(sprintf("%.3f", mean(m$EAP.rel, na.rm = TRUE)))
    }
    sprintf("%.3f (%+.3f)", res$eap_pur, res$delta_eap)
  })

  output$vb_rasch_re_fit_pct <- renderText({
    res <- rasch_reanalysis_results()
    if (is.null(res)) {
      d <- rasch_personfit_df_reactive()
      if (is.null(d)) return("-")
      pct <- mean(d$Infit_MNSQ <= 1.50 & d$Outfit_MNSQ <= 1.50, na.rm = TRUE) * 100
      return(sprintf("%.1f%%", pct))
    }
    sprintf("%.1f%% (%+.1f%%)", res$pct_p_fit_pur, res$delta_p_fit)
  })

  # Tabel Output Re-Analisis
  output$table_rasch_reanalysis_comp <- renderTable({
    res <- rasch_reanalysis_results()
    if (is.null(res) || is.null(res$comp_table)) {
      return(data.frame(Status = "Silakan atur parameter eliminasi butir dan responden di panel kiri, lalu klik 'Jalankan Analisis Ulang Rasch'."))
    }
    res$comp_table
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$table_rasch_re_items_comp <- renderDT({
    res <- rasch_reanalysis_results()
    if (is.null(res) || is.null(res$items_comp_df) || nrow(res$items_comp_df) == 0) {
      return(datatable(data.frame(Status = "Belum ada analisis ulang yang dijalankan. Silakan klik 'Jalankan Analisis Ulang Rasch'.")))
    }
    datatable(
      res$items_comp_df,
      options = list(pageLength = 15, scrollX = TRUE, autoWidth = FALSE),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })

  output$table_rasch_re_dropped_items <- renderDT({
    res <- rasch_reanalysis_results()
    if (is.null(res) || is.null(res$dropped_items_df)) {
      return(datatable(data.frame(Status = "Belum ada analisis ulang yang dijalankan.")))
    }
    datatable(
      res$dropped_items_df,
      options = list(pageLength = 10, scrollX = TRUE, autoWidth = FALSE),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })

  output$table_rasch_re_dropped_persons <- renderDT({
    res <- rasch_reanalysis_results()
    if (is.null(res) || is.null(res$dropped_persons_df)) {
      return(datatable(data.frame(Status = "Belum ada analisis ulang yang dijalankan.")))
    }
    datatable(
      res$dropped_persons_df,
      options = list(pageLength = 10, scrollX = TRUE, autoWidth = FALSE),
      rownames = FALSE,
      class = "display compact cell-border stripe hover"
    )
  })

  # Plot Wright Map Model Pemurnian
  output$plot_rasch_reanalysis_wright_map <- renderPlot({
    res <- rasch_reanalysis_results()
    req(res, res$mod_purified)
    v_mode <- if (!is.null(res$mod_purified$ndim) && res$mod_purified$ndim > 1) "multi" else "uni"
    render_wright_map_internal(custom_mod = res$mod_purified, custom_view_mode = v_mode)
  }, res = 120)

  output$btn_dl_rasch_re_wright_png <- downloadHandler(
    filename = function() paste0("Wright_Map_Pemurnian_", Sys.Date(), ".png"),
    content = function(file) {
      res <- rasch_reanalysis_results()
      req(res, res$mod_purified)
      png(file, width = 3200, height = 2400, res = 300)
      v_mode <- if (!is.null(res$mod_purified$ndim) && res$mod_purified$ndim > 1) "multi" else "uni"
      render_wright_map_internal(custom_mod = res$mod_purified, custom_view_mode = v_mode)
      dev.off()
    }
  )

  output$btn_dl_rasch_reanalysis_comp_csv <- downloadHandler(
    filename = function() paste0("Komparasi_Reanalisis_Rasch_", Sys.Date(), ".csv"),
    content = function(file) {
      res <- rasch_reanalysis_results()
      req(res, res$comp_table)
      write.csv(res$comp_table, file, row.names = FALSE)
    }
  )
  
  # ============================================================================
  # RASCH MODEL FIT COMPARISON SERVER LOGIC (RSM, PCM, MRSM, MPCM)
  # ============================================================================
  rasch_comp_results_store <- reactiveVal(NULL)
  
  observeEvent(input$btn_run_rasch_comp, {
    df <- rasch_data_store()
    req(nrow(df) >= 30, ncol(df) >= 4)
    
    withProgress(message = "Estimasi Komparasi 4 Model Rasch (TAM)...", detail = "Menyiapkan matriks respons politomi...", value = 0.15, {
      
      df_mat <- as.matrix(df)
      min_val <- min(df_mat, na.rm = TRUE)
      if (min_val > 0) df_mat <- df_mat - min_val
      
      Q <- isolate(rasch_active_qmatrix())
      k_items <- ncol(df_mat)
      is_builtin_poly <- all(paste0("Item.", 1:20) %in% colnames(df_mat))
      
      if (is.null(Q) || nrow(Q) != k_items) {
        if (is_builtin_poly) {
          Q <- matrix(0, nrow = 20, ncol = 4)
          rownames(Q) <- paste0("Item.", 1:20)
          colnames(Q) <- c("Faktor_1", "Faktor_2", "Faktor_3", "Faktor_4")
          Q[paste0("Item.", 1:5), "Faktor_1"] <- 1
          Q[paste0("Item.", 6:10), "Faktor_2"] <- 1
          Q[paste0("Item.", 11:15), "Faktor_3"] <- 1
          Q[paste0("Item.", 16:20), "Faktor_4"] <- 1
        } else {
          Q <- matrix(0, nrow = k_items, ncol = 2)
          rownames(Q) <- colnames(df_mat)
          colnames(Q) <- c("Dimensi.1", "Dimensi.2")
          h <- ceiling(k_items / 2)
          Q[1:h, "Dimensi.1"] <- 1
          Q[(h+1):k_items, "Dimensi.2"] <- 1
        }
      }
      
      ctrl_nodes <- list(nodes = seq(-3.5, 3.5, len = 11), maxiter = 45)
      
      n_cores_comp <- if (isTRUE(input$rasch_use_parallel)) {
        if (!is.null(input$rasch_n_cores)) as.integer(input$rasch_n_cores) else 2
      } else 1
      
      setProgress(value = 0.35, detail = sprintf("Mengeksekusi 4 Model (RSM, PCM, MRSM, MPCM) Secara Paralel (%d Cores doParallel)...", n_cores_comp))
      
      model_tasks <- list(
        rsm  = function() TAM::tam.mml(resp = df_mat, irtmodel = "RSM", verbose = FALSE),
        pcm  = function() TAM::tam.mml(resp = df_mat, irtmodel = "PCM", verbose = FALSE),
        mrsm = function() TAM::tam.mml(resp = df_mat, Q = Q, irtmodel = "RSM", control = ctrl_nodes, verbose = FALSE),
        mpcm = function() TAM::tam.mml(resp = df_mat, Q = Q, irtmodel = "PCM", control = ctrl_nodes, verbose = FALSE)
      )
      
      fits <- run_parallel_rasch_tasks(model_tasks, n_cores = n_cores_comp)
      mod_rsm  <- fits$rsm
      mod_pcm  <- fits$pcm
      mod_mrsm <- fits$mrsm
      mod_mpcm <- fits$mpcm
      
      setProgress(value = 0.85, detail = "Menyusun matriks evaluasi kriteria informasi & LRT...")
      
      extract_tam_row <- function(m, name, dim_desc, thresh_desc) {
        if (is.null(m)) {
          return(data.frame(
            Model = name, Dimensi = dim_desc, Threshold = thresh_desc,
            Deviance = NA_real_, Npars = NA_integer_, AIC = NA_real_, BIC = NA_real_, aBIC = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
        ic <- m$ic
        dev <- round(ic$deviance, 3)
        np <- ic$Npars
        aic_val <- round(ic$AIC, 3)
        bic_val <- round(ic$BIC, 3)
        abic_val <- round(ic$aBIC, 3)
        
        data.frame(
          Model = name, Dimensi = dim_desc, Threshold = thresh_desc,
          Deviance = dev, Npars = np, AIC = aic_val, BIC = bic_val, aBIC = abic_val,
          stringsAsFactors = FALSE
        )
      }
      
      comp_df <- rbind(
        extract_tam_row(mod_rsm, "1. RSM (Rating Scale Model)", "Unidimensi (D=1)", "Invarian (Tetap per Butir)"),
        extract_tam_row(mod_pcm, "2. PCM (Partial Credit Model)", "Unidimensi (D=1)", "Bebas (Spesifik per Butir)"),
        extract_tam_row(mod_mrsm, "3. MRSM (Multidimensional RSM)", paste0("Multidimensi (D=", ncol(Q), ")"), "Invarian (Tetap per Butir)"),
        extract_tam_row(mod_mpcm, "4. MPCM (Multidimensional PCM)", paste0("Multidimensi (D=", ncol(Q), ")"), "Bebas (Spesifik per Butir)")
      )
      
      # Delta Deviance vs RSM Baseline (positif = reduksi deviance / statistik chi-square LRT)
      base_dev <- comp_df$Deviance[1]
      base_df <- comp_df$Npars[1]
      
      comp_df$`\\Delta_Deviance` <- sapply(1:nrow(comp_df), function(i) {
        if (i == 1 || is.na(comp_df$Deviance[i]) || is.na(base_dev)) return(NA_real_)
        round(base_dev - comp_df$Deviance[i], 3)
      })
      comp_df$`\\Delta_df` <- sapply(1:nrow(comp_df), function(i) {
        if (i == 1 || is.na(comp_df$Npars[i]) || is.na(base_df)) return(NA_integer_)
        as.integer(comp_df$Npars[i] - base_df)
      })
      comp_df$p_LRT <- sapply(1:nrow(comp_df), function(i) {
        if (i == 1 || is.na(comp_df$`\\Delta_Deviance`[i])) return(NA_real_)
        dd <- comp_df$`\\Delta_Deviance`[i]
        d_df <- comp_df$`\\Delta_df`[i]
        if (d_df > 0 && dd > 0) round(pchisq(dd, df = d_df, lower.tail = FALSE), 3) else NA_real_
      })
      
      min_aic <- min(comp_df$AIC, na.rm = TRUE)
      min_bic <- min(comp_df$BIC, na.rm = TRUE)
      min_dev <- min(comp_df$Deviance, na.rm = TRUE)
      
      comp_df$Keputusan_Parsimoni <- sapply(1:nrow(comp_df), function(i) {
        if (is.na(comp_df$AIC[i])) "Tidak Konvergen"
        else if (comp_df$AIC[i] == min_aic && comp_df$BIC[i] == min_bic) "Model Terbaik (AIC & BIC Minimum)"
        else if (comp_df$AIC[i] == min_aic) "Model Terbaik (AIC Minimum)"
        else if (comp_df$BIC[i] == min_bic) "Model Terbaik (BIC Minimum - Lebih Hemat Parameter)"
        else if (comp_df$Deviance[i] == min_dev) "Fit Data Tertinggi (-2LL Terendah)"
        else "Model Pembanding"
      })
      
      # ANOVA Comparisons: RSM vs PCM and MRSM vs MPCM
      anova_tbl <- tryCatch({
        a_uni <- if (!is.null(mod_rsm) && !is.null(mod_pcm)) anova(mod_rsm, mod_pcm) else NULL
        a_multi <- if (!is.null(mod_mrsm) && !is.null(mod_mpcm)) anova(mod_mrsm, mod_mpcm) else NULL
        
        rows <- list()
        if (!is.null(a_uni)) {
          a_u_df <- as.data.frame(a_uni)
          a_u_df$Perbandingan <- c("RSM (Invarian)", "PCM (Bebas) [vs RSM]")
          rows[[1]] <- a_u_df
        }
        if (!is.null(a_multi)) {
          a_m_df <- as.data.frame(a_multi)
          a_m_df$Perbandingan <- c("MRSM (Invarian)", "MPCM (Bebas) [vs MRSM]")
          rows[[2]] <- a_m_df
        }
        if (length(rows) > 0) do.call(rbind, rows) else NULL
      }, error = function(e) NULL)
      
      best_aic_name <- comp_df$Model[which.min(comp_df$AIC)]
      best_bic_name <- comp_df$Model[which.min(comp_df$BIC)]
      best_dev_name <- comp_df$Model[which.min(comp_df$Deviance)]
      
      setProgress(value = 1.0, detail = "Komparasi 4 Model Selesai!")
      
      res_tam_comp <- list(
        summary = comp_df,
        anova = anova_tbl,
        best_aic = best_aic_name,
        best_bic = best_bic_name,
        best_dev = best_dev_name
      )
      rasch_comp_results_store(res_tam_comp)
      showNotification("Komparasi RSM, PCM, MRSM, dan MPCM selesai!", type = "message")
    })
  })
  
  output$vb_rasch_comp_best_aic <- renderText({
    res <- rasch_comp_results_store()
    if (is.null(res)) return("-")
    res$best_aic
  })
  output$vb_rasch_comp_best_bic <- renderText({
    res <- rasch_comp_results_store()
    if (is.null(res)) return("-")
    res$best_bic
  })
  output$vb_rasch_comp_best_dev <- renderText({
    res <- rasch_comp_results_store()
    if (is.null(res)) return("-")
    res$best_dev
  })
  output$vb_rasch_comp_parsimony <- renderText({
    res <- rasch_comp_results_store()
    if (is.null(res)) return("-")
    "BIC Lebih Efektif Menghukum Kompleksitas Parameter"
  })
  
  output$table_rasch_comp_summary <- renderDT({
    res <- rasch_comp_results_store()
    if (is.null(res)) {
      return(datatable(data.frame(Status = "Silakan klik tombol 'Bandingkan 4 Model Rasch' untuk memulai estimasi komparasi."), rownames = FALSE))
    }
    datatable(res$summary, options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE, class = "display compact cell-border")
  })
  
  output$table_rasch_comp_anova <- renderTable({
    res <- rasch_comp_results_store()
    if (is.null(res) || is.null(res$anova)) {
      return(data.frame(Status = "Uji Likelihood Ratio Test (ANOVA) akan ditampilkan setelah estimasi komparasi dijalankan."))
    }
    res$anova
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$rasch_comp_bar_plot <- renderPlot({
    res <- rasch_comp_results_store()
    req(res, !is.null(res$summary))
    d <- res$summary
    req(nrow(d) >= 2)
    
    par(mar = c(5, 5, 3, 2), bg = "#FFFFFF", mfrow = c(1, 2))
    
    # Plot 1: Deviance (-2LL)
    devs <- d$Deviance
    names(devs) <- c("RSM", "PCM", "MRSM", "MPCM")
    barplot(
      devs,
      col = c("#FAF8F5", "#C5A059", "#800020", "#650915"),
      border = "#231618",
      las = 1,
      ylim = c(min(devs, na.rm = TRUE) * 0.95, max(devs, na.rm = TRUE) * 1.02),
      xpd = FALSE,
      ylab = "Deviance (-2LL)",
      main = "Perbandingan Deviance (-2 Log-Likelihood)",
      col.main = "#650915"
    )
    
    # Plot 2: AIC vs BIC
    aic_bic_mat <- rbind(d$AIC, d$BIC)
    colnames(aic_bic_mat) <- names(devs)
    barplot(
      aic_bic_mat,
      beside = TRUE,
      col = c("#650915", "#C5A059"),
      ylim = c(min(aic_bic_mat, na.rm = TRUE) * 0.95, max(aic_bic_mat, na.rm = TRUE) * 1.02),
      xpd = FALSE,
      las = 1,
      ylab = "Nilai Indeks Informasi",
      main = "Kriteria Informasi AIC vs BIC",
      col.main = "#650915"
    )
    legend("bottomright", legend = c("AIC", "BIC"), fill = c("#650915", "#C5A059"), bty = "n")
  })
  
  output$btn_dl_rasch_comp_csv <- downloadHandler(
    filename = function() paste0("Rasch_Model_Comparison_", Sys.Date(), ".csv"),
    content = function(file) {
      res <- rasch_comp_results_store()
      req(res)
      write.csv(res$summary, file, row.names = FALSE)
    }
  )

  
  # Reactive Data Storage
  benchmark_storage <- reactiveVal(create_benchmark_data())
  calc_compiled_data <- reactiveVal(data.frame()) # starts clean, users can load sample or compute
  active_data <- reactiveVal(create_benchmark_data())
  active_source_label <- reactiveVal("Benchmark RCT (k = 14)")
  calc_result <- reactiveVal(NULL)
  
  # Filtered Active Dataset (supports Quick Sensitivity Filter in Tab 4)
  filtered_active_data <- reactive({
    df <- active_data()
    if (nrow(df) == 0) return(df)
    flt <- input$quick_sensitivity_filter
    if (!is.null(flt) && flt != "all") {
      if (flt == "low_rob" && "Risk_of_Bias" %in% colnames(df)) {
        df <- df[df$Risk_of_Bias == "Low Risk", , drop = FALSE]
      } else if (flt == "no_waitlist" && "Control_Type" %in% colnames(df)) {
        df <- df[df$Control_Type != "Waitlist", , drop = FALSE]
      } else if (flt == "clinical_only" && "Population" %in% colnames(df)) {
        df <- df[grepl("^Clinical", df$Population, ignore.case = TRUE), , drop = FALSE]
      }
    }
    df
  })
  
  # Hero & Workflow Navigation Buttons (CFA, Rasch, and Meta-Analysis)
  observeEvent(input$hero_btn_cfa_model, { nav_select("main_navbar", "tab_cfa_model") })
  observeEvent(input$hero_btn_cfa_elim, { nav_select("main_navbar", "tab_cfa_elim") })
  observeEvent(input$hero_btn_cfa_fit, { nav_select("main_navbar", "tab_cfa_fit") })
  observeEvent(input$hero_btn_cfa_comp, { nav_select("main_navbar", "tab_cfa_comp") })
  observeEvent(input$hero_btn_cfa_path, { nav_select("main_navbar", "tab_cfa_path") })
  
  observeEvent(input$hero_btn_rasch_model, { nav_select("main_navbar", "tab_rasch_model") })
  observeEvent(input$hero_btn_rasch_item, { nav_select("main_navbar", "tab_rasch_itemfit") })
  observeEvent(input$hero_btn_rasch_person, { nav_select("main_navbar", "tab_rasch_personfit") })
  observeEvent(input$hero_btn_rasch_reanalysis, { nav_select("main_navbar", "tab_rasch_reanalysis") })
  observeEvent(input$hero_btn_rasch_wright, { nav_select("main_navbar", "tab_rasch_wright") })
  observeEvent(input$hero_btn_rasch_comp, { nav_select("main_navbar", "tab_rasch_comp") })
  
  observeEvent(input$hero_btn_calc, { nav_select("main_navbar", "tab_calc") })
  observeEvent(input$hero_btn_data, { nav_select("main_navbar", "tab_model") })
  observeEvent(input$hero_btn_subgroup, { nav_select("main_navbar", "tab_forest") })
  observeEvent(input$hero_btn_reg, { nav_select("main_navbar", "tab_reg") })
  observeEvent(input$hero_btn_pub, { nav_select("main_navbar", "tab_bias") })
  observeEvent(input$btn_goto_compiled, {
    nav_select("main_navbar", "tab_compiled")
  })
  observeEvent(input$btn_data_to_model, {
    nav_select("main_navbar", "tab_model")
  })
  observeEvent(input$btn_model_to_subgroup, {
    nav_select("main_navbar", "tab_subgroup")
  })
  observeEvent(input$btn_model_to_forest, {
    nav_select("main_navbar", "tab_forest")
  })
  observeEvent(input$btn_model_to_bias, {
    nav_select("main_navbar", "tab_bias")
  })
  observeEvent(input$btn_calc_direct_model, {
    df_c <- calc_compiled_data()
    if (nrow(df_c) >= 2 && !isTRUE(input$calc_auto_push_active)) {
      active_data(df_c)
      active_source_label(paste0("Kompilasi Kalkulator (k = ", nrow(df_c), ")"))
      updateRadioButtons(session, "dataset_source_choice", selected = "calculator")
    }
    nav_select("main_navbar", "tab_model")
  })
  observeEvent(input$btn_calc_to_model, {
    df <- calc_compiled_data()
    if (nrow(df) < 2) {
      showNotification("Minimal dibutuhkan 2 studi dalam matriks kalkulator untuk menjalankan meta-analisis.", type = "warning")
      return()
    }
    active_data(df)
    active_source_label(paste0("Kompilasi Kalkulator (k = ", nrow(df), ")"))
    updateRadioButtons(session, "dataset_source_choice", selected = "calculator")
    nav_select("main_navbar", "tab_model")
    showNotification("Matriks kalkulator diaktifkan. Menampilkan hasil pemodelan meta-analisis.", type = "message")
  })
  
  # Tab 1: Normal Distribution Engine Plot
  output$kawe_normal_demo_plot <- renderPlot({
    mu_val <- input$demo_mu
    sigma_val <- input$demo_sigma
    show_zones <- input$demo_show_zones
    show_scales <- input$demo_show_scales
    
    x <- seq(mu_val - 4 * sigma_val, mu_val + 4 * sigma_val, length.out = 500)
    y <- dnorm(x, mean = mu_val, sd = sigma_val)
    
    par(mar = c(4.8, 4.5, 2.5, 1.5), bg = "#FFFFFF")
    plot(
      x, y, type = "n",
      xlab = "Psychometric Theta Scale / Latent Metric (x)",
      ylab = "Probability Density f(x)",
      main = "Normal Distribution Metric — KAWE Psychometric Engine",
      bty = "l",
      las = 1,
      col.axis = "#231618",
      col.lab = "#650915",
      col.main = "#650915",
      font.main = 2
    )
    
    grid(col = "#E5D3B3", lty = 3)
    
    if (isTRUE(show_zones)) {
      idx2 <- x >= (mu_val - 2 * sigma_val) & x <= (mu_val + 2 * sigma_val)
      polygon(c(x[idx2], rev(x[idx2])), c(y[idx2], rep(0, sum(idx2))), col = rgb(197/255, 160/255, 89/255, 0.20), border = NA)
      idx1 <- x >= (mu_val - 1 * sigma_val) & x <= (mu_val + 1 * sigma_val)
      polygon(c(x[idx1], rev(x[idx1])), c(y[idx1], rep(0, sum(idx1))), col = rgb(197/255, 160/255, 89/255, 0.45), border = NA)
    }
    
    lines(x, y, col = "#650915", lwd = 3.2)
    abline(v = mu_val, col = "#C5A059", lty = 2, lwd = 1.8)
    
    usr <- par("usr")
    text(mu_val, usr[4] * 0.92, sprintf("mu = %.2f, sigma = %.2f", mu_val, sigma_val), col = "#650915", font = 2, pos = 4)
    
    if (isTRUE(show_scales)) {
      text(usr[1] + 0.2, usr[4] * 0.85, "Standard Psychometric Equivalences:\n• z-score = (x - mu) / sigma\n• T-score = 50 + 10z\n• Wechsler IQ = 100 + 15z", adj = c(0, 1), cex = 0.82, col = "#5A4B4C", font = 3)
    }
  })
  
  # Core Function to Compute Effect Size from Inputs
  compute_current_es <- function() {
    itype <- input$input_type
    metric <- input$es_metric
    
    apply_hedges_var <- function(r_obj) {
      if (metric == "g" && !is.null(r_obj) && !is.null(r_obj$totaln) && is.finite(r_obj$totaln) && r_obj$totaln > 3 && is.finite(r_obj$var)) {
        j_val <- 1 - (3 / (4 * r_obj$totaln - 9))
        r_obj$var <- r_obj$var * (j_val^2)
        r_obj$se <- sqrt(r_obj$var)
        r_obj$ci.lo <- r_obj$es - qnorm(0.975) * r_obj$se
        r_obj$ci.hi <- r_obj$es + qnorm(0.975) * r_obj$se
        r_obj$w <- 1 / r_obj$var
      }
      r_obj
    }
    
    if (itype == "mean_gain") {
      gsub <- input$gain_sub_type
      gdir <- input$gain_outcome_dir
      r_val <- input$gain_r
      
      if (gsub == "pre_post_full") {
        p1_m <- input$pre1_m
        p1_sd <- input$pre1_sd
        po1_m <- input$post1_m
        po1_sd <- input$post1_sd
        n1 <- input$pre1_n
        
        p2_m <- input$pre2_m
        p2_sd <- input$pre2_sd
        po2_m <- input$post2_m
        po2_sd <- input$post2_sd
        n2 <- input$pre2_n
        
        if (gdir == "decrease_better") {
          g1_m <- p1_m - po1_m
          g2_m <- p2_m - po2_m
        } else {
          g1_m <- po1_m - p1_m
          g2_m <- po2_m - p2_m
        }
        
        s_c1 <- sqrt(pmax(1e-8, p1_sd^2 + po1_sd^2 - 2 * r_val * p1_sd * po1_sd))
        s_c2 <- sqrt(pmax(1e-8, p2_sd^2 + po2_sd^2 - 2 * r_val * p2_sd * po2_sd))
      } else {
        g1_m <- input$gain1_m
        s_c1 <- input$gain1_sd
        n1 <- input$gain1_n
        
        g2_m <- input$gain2_m
        s_c2 <- input$gain2_sd
        n2 <- input$gain2_n
        
        if (gdir == "decrease_better" && g1_m < 0) {
          g1_m <- -g1_m
          g2_m <- -g2_m
        }
      }
      
      sp_change <- sqrt(((n1 - 1) * s_c1^2 + (n2 - 1) * s_c2^2) / (n1 + n2 - 2))
      d_val <- (g1_m - g2_m) / sp_change
      j_val <- 1 - (3 / (4 * (n1 + n2) - 9))
      vd_val <- (2 * (1 - r_val)) / n1 + (2 * (1 - r_val)) / n2 + (d_val^2) / (2 * (n1 + n2))
      es_out <- if (metric == "g") (d_val * j_val) else d_val
      var_out <- if (metric == "g") (j_val^2 * vd_val) else vd_val
      se_out <- sqrt(var_out)
      z_crit <- qnorm(0.975)
      
      res <- list(
        es = es_out,
        var = var_out,
        se = se_out,
        ci.lo = es_out - z_crit * se_out,
        ci.hi = es_out + z_crit * se_out,
        w = 1 / var_out,
        totaln = n1 + n2,
        measure = metric
      )
    } else if (itype == "mean_sd") {
      res <- esc_mean_sd(
        grp1m = input$grp1m_sd, grp1sd = input$grp1sd_sd, grp1n = input$grp1n_sd,
        grp2m = input$grp2m_sd, grp2sd = input$grp2sd_sd, grp2n = input$grp2n_sd,
        es.type = metric
      )
      if (!is.null(input$sd_outcome_dir) && input$sd_outcome_dir == "decrease_better") {
        res$es <- -res$es
        temp_lo <- -res$ci.hi
        res$ci.hi <- -res$ci.lo
        res$ci.lo <- temp_lo
      }
      res <- apply_hedges_var(res)
    } else if (itype == "mean_se") {
      s1_est <- input$grp1se_se * sqrt(input$grp1n_se)
      s2_est <- input$grp2se_se * sqrt(input$grp2n_se)
      res <- esc_mean_sd(
        grp1m = input$grp1m_se, grp1sd = s1_est, grp1n = input$grp1n_se,
        grp2m = input$grp2m_se, grp2sd = s2_est, grp2n = input$grp2n_se,
        es.type = metric
      )
      if (!is.null(input$se_outcome_dir) && input$se_outcome_dir == "decrease_better") {
        res$es <- -res$es
        temp_lo <- -res$ci.hi
        res$ci.hi <- -res$ci.lo
        res$ci.lo <- temp_lo
      }
      res <- apply_hedges_var(res)
    } else if (itype == "t_test") {
      res <- apply_hedges_var(esc_t(t = input$t_val, grp1n = input$t_n1, grp2n = input$t_n2, es.type = metric))
    } else if (itype == "p_val") {
      df_deg <- input$p_n1 + input$p_n2 - 2
      t_equiv <- qt(1 - (input$p_value / 2), df_deg)
      if (input$p_direction == "neg") t_equiv <- -t_equiv
      res <- apply_hedges_var(esc_t(t = t_equiv, grp1n = input$p_n1, grp2n = input$p_n2, es.type = metric))
    } else if (itype == "f_test") {
      res <- esc_f(f = input$f_val, grp1n = input$f_n1, grp2n = input$f_n2, es.type = metric)
      if (input$f_direction == "neg") {
        res$es <- -abs(res$es)
        temp_lo <- -res$ci.hi
        res$ci.hi <- -res$ci.lo
        res$ci.lo <- temp_lo
      }
      res <- apply_hedges_var(res)
    } else if (itype == "table_2x2") {
      res <- apply_hedges_var(esc_2x2(
        grp1yes = input$tab_e1, grp1no = input$tab_n1 - input$tab_e1,
        grp2yes = input$tab_e2, grp2no = input$tab_n2 - input$tab_e2,
        es.type = metric
      ))
    } else if (itype == "correlation") {
      res <- apply_hedges_var(esc_rpb(r = input$corr_r, totaln = input$corr_n, es.type = metric))
    }
    res
  }
  
  # Live Reactive Auto-Calculation for Calculator Inputs
  observe({
    res <- tryCatch(compute_current_es(), error = function(e) NULL)
    if (!is.null(res) && length(res$es) == 1 && is.finite(res$es) && is.finite(res$var) && res$var > 0) {
      calc_result(res)
    }
  })
  
  # Observe Calculate Button
  observeEvent(input$btn_calculate, {
    tryCatch({
      res <- compute_current_es()
      calc_result(res)
      showNotification("Estimasi besaran efek berhasil dikalkulasi ulang.", type = "message")
    }, error = function(e) {
      showNotification(paste("Kesalahan perhitungan:", e$message), type = "error")
    })
  })
  
  # Observe Append Study Button (Smart auto-compute + add)
  observeEvent(input$btn_append_study, {
    res <- calc_result()
    if (is.null(res)) {
      res <- tryCatch(compute_current_es(), error = function(e) NULL)
      if (is.null(res)) {
        showNotification("Gagal menghitung besaran efek. Periksa kembali input numerik.", type = "error")
        return()
      }
      calc_result(res)
    }
    
    itype <- input$input_type
    n1_val <- 50
    n2_val <- 50
    if (itype == "mean_gain") {
      if (input$gain_sub_type == "pre_post_full") {
        n1_val <- input$pre1_n
        n2_val <- input$pre2_n
      } else {
        n1_val <- input$gain1_n
        n2_val <- input$gain2_n
      }
    } else if (itype == "mean_sd") {
      n1_val <- input$grp1n_sd
      n2_val <- input$grp2n_sd
    } else if (itype == "mean_se") {
      n1_val <- input$grp1n_se
      n2_val <- input$grp2n_se
    } else if (itype == "t_test") {
      n1_val <- input$t_n1
      n2_val <- input$t_n2
    } else if (itype == "p_val") {
      n1_val <- input$p_n1
      n2_val <- input$p_n2
    } else if (itype == "f_test") {
      n1_val <- input$f_n1
      n2_val <- input$f_n2
    } else if (itype == "table_2x2") {
      n1_val <- input$tab_n1
      n2_val <- input$tab_n2
    } else if (itype == "correlation") {
      n1_val <- round(input$corr_n / 2)
      n2_val <- round(input$corr_n / 2)
    }
    
    active_comp_val <- if (input$study_ctrl == "Active Control") {
      input$study_active_comparator
    } else if (input$study_ctrl == "Waitlist") {
      "None (Waitlist)"
    } else if (input$study_ctrl == "TAU") {
      "Treatment as Usual"
    } else if (input$study_ctrl == "Placebo") {
      "Attention Placebo"
    } else {
      "-"
    }
    
    new_row <- data.frame(
      Study = input$study_name,
      Year = input$study_year,
      Intervention = input$study_interv,
      Format = input$study_format,
      Control_Type = input$study_ctrl,
      Active_Comparator = active_comp_val,
      Population = input$study_pop,
      Age_Group = input$study_age,
      Therapist_Level = input$study_therapist,
      Risk_of_Bias = input$study_rob,
      Duration_Weeks = input$study_dur,
      Total_Sessions = input$study_sessions,
      N_Intervention = n1_val,
      N_Control = n2_val,
      yi = round(res$es, 4),
      vi = round(res$var, 5),
      sei = round(res$se, 4),
      stringsAsFactors = FALSE
    )
    
    # Append to compiled calculator batch
    comp_df <- calc_compiled_data()
    calc_compiled_data(rbind(comp_df, new_row))
    
    # If auto-push is checked, also append to active analysis dataset
    if (isTRUE(input$calc_auto_push_active)) {
      curr_act <- active_data()
      active_data(rbind(curr_act, new_row))
    }
    
    showNotification(paste0("Studi '", input$study_name, "' berhasil disimpan ke matriks gabungan kalkulator."), type = "message")
  })
  
  # Reset Calculator Form Button
  observeEvent(input$btn_reset_calc_form, {
    updateTextInput(session, "study_name", value = paste("Studi Baru", nrow(calc_compiled_data()) + 1))
    calc_result(NULL)
    showNotification("Formulir telah disiapkan untuk entri studi berikutnya.", type = "default")
  })
  
  # Calculator Outputs
  output$vb_calc_es <- renderText({
    res <- calc_result()
    if (is.null(res)) return("-")
    sprintf("%.3f", res$es)
  })
  output$vb_calc_var <- renderText({
    res <- calc_result()
    if (is.null(res)) return("-")
    sprintf("%.3f", res$var)
  })
  output$vb_calc_se <- renderText({
    res <- calc_result()
    if (is.null(res)) return("-")
    sprintf("%.3f", res$se)
  })
  
  output$calc_visual_preview_plot <- renderPlot({
    res <- calc_result()
    req(!is.null(res), is.finite(res$es), is.finite(res$ci.lo), is.finite(res$ci.hi))
    par(mar = c(3.2, 1.5, 1.2, 1.5), bg = "#FFFFFF")
    xlim_range <- c(min(-0.6, res$ci.lo - 0.25), max(1.4, res$ci.hi + 0.25))
    plot(
      0, 0, type = "n",
      xlim = xlim_range, ylim = c(0.4, 1.6),
      yaxt = "n", ylab = "",
      xlab = "Besaran Efek Terstandardisasi (g > 0 = Memihak Intervensi)",
      bty = "n", cex.lab = 0.85, col.lab = "#650915"
    )
    # Shaded magnitude zones
    rect(0.2, 0.45, 0.5, 1.55, col = rgb(197/255, 160/255, 89/255, 0.10), border = NA)
    rect(0.5, 0.45, 0.8, 1.55, col = rgb(197/255, 160/255, 89/255, 0.20), border = NA)
    rect(0.8, 0.45, xlim_range[2], 1.55, col = rgb(101/255, 9/255, 21/255, 0.08), border = NA)
    abline(v = 0, lty = 2, lwd = 1.5, col = "#5A4B4C")
    abline(v = c(0.2, 0.5, 0.8), lty = 3, col = "#C5A059")
    segments(res$ci.lo, 1.0, res$ci.hi, 1.0, col = "#650915", lwd = 3.0)
    arrows(res$ci.lo, 1.0, res$ci.hi, 1.0, code = 3, angle = 90, length = 0.06, col = "#650915", lwd = 2.2)
    points(res$es, 1.0, pch = 18, cex = 2.2, col = "#650915")
    text(
      res$es, 1.32,
      sprintf("g = %.3f [%.3f, %.3f]", res$es, res$ci.lo, res$ci.hi),
      col = "#650915", font = 2, cex = 0.88
    )
  })
  
  output$calc_summary_table <- renderTable({
    res <- calc_result()
    if (is.null(res)) {
      return(data.frame(Status = "Masukkan parameter numerik dan klik 'Hitung Parameter Efek'."))
    }
    w_val <- ifelse(is.finite(res$var) && res$var > 0, sprintf("%.3f (1/vi)", 1 / res$var), "-")
    
    data.frame(
      Parameter = c(
        "Metrik & Arah Besaran Efek",
        "Titik Estimasi Efek (Point Estimate)",
        "Varians Sampling (vi)",
        "Standar Error (sei)",
        "Bobot Studi (Weight wi = 1/vi)",
        "Interval Kepercayaan 95% CI",
        "Kategori Cohen (1988)"
      ),
      Nilai = c(
        paste0(ifelse(input$es_metric == "g", "Hedges' g (Unbiased)", "Cohen's d"), " (", ifelse(res$es >= 0, "Memihak Intervensi", "Memihak Kontrol"), ")"),
        sprintf("%.3f", res$es),
        sprintf("%.3f", res$var),
        sprintf("%.3f", res$se),
        w_val,
        sprintf("[%.3f, %.3f]", res$ci.lo, res$ci.hi),
        ifelse(abs(res$es) < 0.2, "Dapat Diabaikan (< 0.20)",
               ifelse(abs(res$es) < 0.5, "Kecil (0.20 - 0.49)",
                      ifelse(abs(res$es) < 0.8, "Sedang (0.50 - 0.79)", "Besar (>= 0.80)")))
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Active Method Badge in Card Header
  output$calc_active_formula_badge <- renderUI({
    itype <- input$input_type
    metric_lbl <- ifelse(input$es_metric == "g", "Hedges' g (Unbiased)", "Cohen's d")
    type_lbl <- switch(itype,
      "mean_gain" = "Pre-Post Mean Gain",
      "mean_sd" = "Post-Test Means & SDs",
      "mean_se" = "Means & SE",
      "t_test" = "Independent t-test",
      "p_val" = "p-value & Sample Sizes",
      "f_test" = "One-Way ANOVA F-test",
      "table_2x2" = "2x2 Contingency Table (OR)",
      "correlation" = "Correlation Coefficient (r)",
      itype
    )
    span(class = "badge bg-primary", paste("Metode:", type_lbl, "|", metric_lbl))
  })
  
  # Quick Formula Callout in Output Card
  output$calc_quick_formula_callout <- renderUI({
    itype <- input$input_type
    metric_lbl <- ifelse(input$es_metric == "g", "Hedges' g", "Cohen's d")
    eq_str <- switch(itype,
      "mean_sd" = "$$d = \\frac{\\bar{X}_1 - \\bar{X}_2}{s_{\\text{pooled}}}, \\quad s_{\\text{pooled}} = \\sqrt{\\frac{(n_1-1)s_1^2 + (n_2-1)s_2^2}{n_1+n_2-2}}, \\quad g = J \\cdot d$$",
      "mean_se" = "$$s_i = SE_i \\sqrt{n_i}, \\quad d = \\frac{\\bar{X}_1 - \\bar{X}_2}{s_{\\text{pooled}}}, \\quad g = J \\cdot d$$",
      "mean_gain" = "$$\\text{Gain}_i = \\Delta\\bar{X}_i, \\quad s_{\\text{change}} = \\sqrt{s_{\\text{pre}}^2 + s_{\\text{post}}^2 - 2 r s_{\\text{pre}}s_{\\text{post}}}, \\quad d = \\frac{\\Delta\\text{Gain}}{s_{\\text{pooled}}}, \\quad g = J \\cdot d$$",
      "t_test" = "$$d = t \\sqrt{\\frac{n_1 + n_2}{n_1 n_2}}, \\quad J = 1 - \\frac{3}{4(n_1+n_2)-9}, \\quad g = J \\cdot d$$",
      "f_test" = "$$|t| = \\sqrt{F}, \\quad d = \\pm \\sqrt{F\\left(\\frac{n_1 + n_2}{n_1 n_2}\\right)}, \\quad g = J \\cdot d$$",
      "p_val" = "$$t = F_{t, df}^{-1}(1 - p/2), \\quad d = \\pm t \\sqrt{\\frac{n_1 + n_2}{n_1 n_2}}, \\quad g = J \\cdot d$$",
      "table_2x2" = "$$OR = \\frac{a \\cdot d}{b \\cdot c}, \\quad d = \\ln(OR) \\times \\frac{\\sqrt{3}}{\\pi} \\approx 0.5513 \\ln(OR), \\quad g = J \\cdot d$$",
      "correlation" = "$$d = \\frac{2r}{\\sqrt{1 - r^2}}, \\quad z = \\frac{1}{2}\\ln\\left(\\frac{1+r}{1-r}\\right), \\quad g = J \\cdot d$$",
      "$$d = \\dots, \\quad g = J \\cdot d$$"
    )
    div(
      class = "kawe-formula-step mt-3",
      div(
        class = "d-flex justify-content-between align-items-center mb-1",
        strong(class = "text-primary", icon("square-root-alt"), " Rumus:"),
        span(class = "badge bg-secondary", metric_lbl)
      ),
      div(class = "kawe-eq-box my-1", HTML(eq_str)),
      div(class = "small text-muted", "Formulasi lengkap untuk seluruh variasi metodologis luaran klinis dapat dipelajari pada tab Kumpulan Rumus dan Panduan Perhitungan.")
    )
  })

  # Navigation Observers inside calc_results_navset
  observeEvent(input$btn_view_calc_compendium, {
    nav_select("calc_results_navset", "panel_formula_catalog")
  })
  observeEvent(input$btn_back_to_calc_results2, {
    nav_select("calc_results_navset", "panel_calc_results")
  })

  # Dynamic Active Formula Display
  output$calc_active_formula_display <- renderUI({
    res <- calc_result()
    build_active_formula_ui(input$input_type, input$es_metric, input, res)
  })
  
  output$calc_saved_count_badge <- renderUI({
    df <- calc_compiled_data()
    span(class = "badge bg-primary", paste(nrow(df), "Studi Tersimpan"))
  })
  
  output$table_calc_recent_preview <- renderTable({
    df <- calc_compiled_data()
    if (nrow(df) == 0) {
      return(data.frame(Status = "Belum ada studi yang disimpan pada sesi kalkulator ini."))
    }
    tail(df[, c("Study", "Year", "Intervention", "Control_Type", "yi", "sei")], 5)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Tab 3: Gabungan Hasil Kalkulator Server Handlers
  output$vb_calc_comp_k <- renderText({
    df <- calc_compiled_data()
    if (nrow(df) == 0) return("0")
    paste(nrow(df), "Studi")
  })
  output$vb_calc_comp_mean_g <- renderText({
    df <- calc_compiled_data()
    if (nrow(df) == 0) return("-")
    sprintf("%.3f", mean(df$yi, na.rm = TRUE))
  })
  output$vb_calc_comp_range_g <- renderText({
    df <- calc_compiled_data()
    if (nrow(df) == 0) return("-")
    sprintf("%.3f s/d %.3f", min(df$yi, na.rm = TRUE), max(df$yi, na.rm = TRUE))
  })
  output$vb_calc_comp_n <- renderText({
    df <- calc_compiled_data()
    if (nrow(df) == 0) return("0")
    total_n <- sum(df$N_Intervention + df$N_Control, na.rm = TRUE)
    sprintf("%s Partisipan", format(total_n, big.mark = "."))
  })
  
  output$calc_compiled_status_alert <- renderUI({
    df <- calc_compiled_data()
    act <- active_data()
    is_active <- nrow(df) > 0 && identical(df$Study, act$Study)
    if (is_active) {
      div(
        class = "alert alert-success border my-2 py-2 d-flex align-items-center",
        icon("check-circle", class = "me-2 fs-5"),
        span("Matriks kalkulator ini saat ini ", strong("SEDANG AKTIF DIGUNAKAN"), " sebagai data analisis utama pada modul Meta-Analisis, Subgroup, dan Forest Plot.")
      )
    } else {
      div(
        class = "alert alert-light border my-2 py-2 d-flex align-items-center justify-content-between",
        div(
          icon("info-circle", class = "text-primary me-2 fs-5"),
          span("Studi hasil kalkulator terkumpul di sini. Klik ", strong("Jadikan Dataset Analisis Utama"), " untuk menganalisis batch ini secara langsung.")
        )
      )
    }
  })
  
  observeEvent(input$btn_use_calc_as_active, {
    df <- calc_compiled_data()
    if (nrow(df) == 0) {
      showNotification("Matriks kalkulator masih kosong. Silakan hitung studi terlebih dahulu.", type = "warning")
      return()
    }
    active_data(df)
    active_source_label(paste0("Kompilasi Kalkulator (k = ", nrow(df), ")"))
    updateRadioButtons(session, "dataset_source_choice", selected = "calculator")
    showNotification("Matriks kalkulator berhasil ditetapkan sebagai dataset analisis utama.", type = "message")
  })
  
  observeEvent(input$btn_append_calc_to_active, {
    df <- calc_compiled_data()
    if (nrow(df) == 0) {
      showNotification("Matriks kalkulator masih kosong.", type = "warning")
      return()
    }
    curr <- active_data()
    # avoid exact study name duplicates
    new_studies <- df[!(df$Study %in% curr$Study), ]
    if (nrow(new_studies) == 0) {
      showNotification("Semua studi dalam matriks kalkulator sudah ada dalam dataset aktif.", type = "warning")
      return()
    }
    combined <- rbind(curr, new_studies)
    active_data(combined)
    active_source_label(paste0("Dataset Gabungan (k = ", nrow(combined), ")"))
    showNotification(paste("Berhasil menambahkan", nrow(new_studies), "studi baru ke dataset aktif."), type = "message")
  })
  
  observeEvent(input$btn_load_calc_sample, {
    calc_compiled_data(create_calc_sample_batch())
    showNotification("Contoh batch studi kalkulator (k = 5) berhasil dimuat.", type = "message")
  })
  
  observeEvent(input$btn_delete_calc_row, {
    sel <- input$table_calc_compiled_rows_selected
    df <- calc_compiled_data()
    if (is.null(sel) || length(sel) == 0 || nrow(df) == 0) {
      showNotification("Pilih satu atau beberapa baris pada tabel terlebih dahulu untuk dihapus.", type = "warning")
      return()
    }
    calc_compiled_data(df[-sel, , drop = FALSE])
    showNotification(paste("Berhasil menghapus", length(sel), "studi dari matriks kompilasi."), type = "message")
  })
  
  observeEvent(input$btn_clear_calc_compiled, {
    calc_compiled_data(data.frame())
    showNotification("Matriks kompilasi kalkulator telah dikosongkan.", type = "warning")
  })
  
  output$btn_download_calc_csv <- downloadHandler(
    filename = function() paste0("Compiled_Calculator_Effects_", Sys.Date(), ".csv"),
    content = function(file) write.csv(calc_compiled_data(), file, row.names = FALSE)
  )
  
  output$table_calc_compiled <- renderDT({
    df <- calc_compiled_data()
    if (nrow(df) == 0) return(NULL)
    inv_var <- 1 / df$vi
    df$Bobot_Persen <- round((inv_var / sum(inv_var)) * 100, 2)
    datatable(
      df,
      selection = "multiple",
      options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE),
      rownames = FALSE,
      class = "display compact cell-border"
    ) %>%
      formatRound(columns = c("yi", "vi", "sei"), digits = 4) %>%
      formatRound(columns = "Bobot_Persen", digits = 2)
  })
  
  # Tab 4: Dataset Sintesis Bukti Server Handlers
  observeEvent(input$dataset_source_choice, {
    src <- input$dataset_source_choice
    if (src == "benchmark") {
      active_data(benchmark_storage())
      active_source_label("Benchmark RCT (k = 14)")
    } else if (src == "calculator") {
      df_calc <- calc_compiled_data()
      if (nrow(df_calc) > 0) {
        active_data(df_calc)
        active_source_label(paste0("Kompilasi Kalkulator (k = ", nrow(df_calc), ")"))
      } else {
        showNotification("Matriks kompilasi kalkulator masih kosong. Menggunakan data aktif sebelumnya.", type = "warning")
      }
    }
  })
  
  observeEvent(input$btn_load_benchmark, {
    active_data(create_benchmark_data())
    benchmark_storage(create_benchmark_data())
    active_source_label("Benchmark RCT (k = 14)")
    updateRadioButtons(session, "dataset_source_choice", selected = "benchmark")
    updateSelectInput(session, "quick_sensitivity_filter", selected = "all")
    showNotification("Data benchmark intervensi klinis (k = 14) berhasil dimuat.", type = "message")
  })
  
  observeEvent(input$btn_delete_active_row, {
    sel <- input$table_active_dataset_rows_selected
    df_flt <- filtered_active_data()
    df_all <- active_data()
    if (is.null(sel) || length(sel) == 0 || nrow(df_flt) == 0) {
      showNotification("Pilih satu atau beberapa baris studi pada tabel untuk dihapus.", type = "warning")
      return()
    }
    studies_to_drop <- df_flt$Study[sel]
    df_new <- df_all[!(df_all$Study %in% studies_to_drop), , drop = FALSE]
    active_data(df_new)
    active_source_label(paste0("Dataset Kustom (k = ", nrow(df_new), ")"))
    showNotification(paste("Berhasil menghapus", length(studies_to_drop), "studi dari dataset aktif."), type = "message")
  })
  
  observeEvent(input$table_active_dataset_cell_edit, {
    info <- input$table_active_dataset_cell_edit
    df <- active_data()
    row_idx <- info$row
    col_idx <- info$col + 1 # rownames = FALSE in DT
    if (row_idx >= 1 && row_idx <= nrow(df) && col_idx >= 1 && col_idx <= ncol(df)) {
      val <- info$value
      if (is.numeric(df[[col_idx]])) {
        num_val <- suppressWarnings(as.numeric(val))
        if (!is.na(num_val)) {
          df[row_idx, col_idx] <- num_val
          if (colnames(df)[col_idx] == "vi" && num_val > 0) {
            df$sei[row_idx] <- round(sqrt(num_val), 4)
          } else if (colnames(df)[col_idx] == "sei" && num_val > 0) {
            df$vi[row_idx] <- round(num_val^2, 5)
          }
          active_data(df)
        }
      } else {
        df[row_idx, col_idx] <- as.character(val)
        active_data(df)
      }
    }
  })
  
  observeEvent(input$btn_clear_data, {
    empty_df <- data.frame(
      Study = character(0), Year = numeric(0), Intervention = character(0), Format = character(0),
      Control_Type = character(0), Active_Comparator = character(0), Population = character(0),
      Age_Group = character(0), Therapist_Level = character(0), Risk_of_Bias = character(0),
      Duration_Weeks = numeric(0), Total_Sessions = numeric(0), N_Intervention = numeric(0),
      N_Control = numeric(0), yi = numeric(0), vi = numeric(0), sei = numeric(0), stringsAsFactors = FALSE
    )
    active_data(empty_df)
    active_source_label("Kosong (k = 0)")
    showNotification("Dataset analisis aktif telah direset.", type = "warning")
  })
  
  observeEvent(input$upload_file, {
    req(input$upload_file)
    file_path <- input$upload_file$datapath
    ext <- tools::file_ext(file_path)
    
    tryCatch({
      if (ext == "csv") {
        df_uploaded <- read.csv(file_path, stringsAsFactors = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        if (!requireNamespace("readxl", quietly = TRUE)) stop("Paket 'readxl' diperlukan untuk membaca Excel.")
        df_uploaded <- as.data.frame(readxl::read_excel(file_path))
      } else {
        stop("Format berkas tidak didukung. Harap unggah berkas .csv atau .xlsx.")
      }
      
      colnames_lower <- tolower(colnames(df_uploaded))
      yi_col <- which(colnames_lower %in% c("yi", "es", "hedges_g", "cohen_d", "effectsize", "g", "d"))[1]
      if (is.na(yi_col)) stop("Kolom besaran efek (yi / es / g) tidak ditemukan.")
      df_uploaded$yi <- as.numeric(df_uploaded[[yi_col]])
      
      vi_col <- which(colnames_lower %in% c("vi", "var", "variance", "samp_var"))[1]
      sei_col <- which(colnames_lower %in% c("sei", "se", "standarderror", "std_err"))[1]
      
      if (!is.na(vi_col)) {
        df_uploaded$vi <- as.numeric(df_uploaded[[vi_col]])
        df_uploaded$sei <- sqrt(df_uploaded$vi)
      } else if (!is.na(sei_col)) {
        df_uploaded$sei <- as.numeric(df_uploaded[[sei_col]])
        df_uploaded$vi <- (df_uploaded$sei)^2
      } else {
        stop("Kolom varians (vi) atau standard error (sei) wajib tersedia.")
      }
      
      study_col <- which(colnames_lower %in% c("study", "author", "author_year", "peneliti", "label"))[1]
      if (!is.na(study_col)) {
        df_uploaded$Study <- as.character(df_uploaded[[study_col]])
      } else {
        df_uploaded$Study <- paste("Study", seq_len(nrow(df_uploaded)))
      }
      
      df_valid <- df_uploaded[!is.na(df_uploaded$yi) & !is.na(df_uploaded$vi) & df_uploaded$vi > 0, ]
      if (nrow(df_valid) == 0) stop("Tidak ditemukan baris data yang valid (yi numerik dan vi > 0).")
      
      active_data(df_valid)
      active_source_label(paste0("Berkas Unggahan (k = ", nrow(df_valid), ")"))
      updateSelectInput(session, "quick_sensitivity_filter", selected = "all")
      showNotification(paste("Berhasil memuat", nrow(df_valid), "studi dari berkas kustom."), type = "message")
    }, error = function(e) {
      showNotification(paste("Gagal memuat berkas:", e$message), type = "error")
    })
  })
  
  output$table_active_dataset <- renderDT({
    df <- filtered_active_data()
    if (nrow(df) == 0) return(NULL)
    inv_var <- 1 / df$vi
    df$Bobot_Persen <- round((inv_var / sum(inv_var)) * 100, 2)
    datatable(
      df,
      editable = TRUE,
      selection = "multiple",
      options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE),
      rownames = FALSE,
      class = "display compact cell-border"
    ) %>%
      formatRound(columns = c("yi", "vi", "sei"), digits = 4) %>%
      formatRound(columns = "Bobot_Persen", digits = 2)
  })
  
  output$btn_download_csv <- downloadHandler(
    filename = function() paste0("Intervention_Meta_Dataset_", Sys.Date(), ".csv"),
    content = function(file) write.csv(filtered_active_data(), file, row.names = FALSE)
  )
  
  output$btn_download_template_csv <- downloadHandler(
    filename = function() "Template_Meta_Analysis_KAWE.csv",
    content = function(file) write.csv(head(create_benchmark_data(), 4), file, row.names = FALSE)
  )
  
  # Dynamic Covariates & Subgroup Dropdowns
  observe({
    df <- filtered_active_data()
    if (nrow(df) > 0) {
      num_cols <- names(df)[sapply(df, is.numeric)]
      valid_covars <- setdiff(num_cols, c("yi", "vi", "sei", "Year", "N_Intervention", "N_Control"))
      if (length(valid_covars) > 0) {
        choices <- c("Tanpa Kovariat Kontinu" = "none", valid_covars)
      } else {
        choices <- c("Tanpa Kovariat Kontinu" = "none")
      }
      updateSelectInput(session, "meta_continuous_covar", choices = choices)
      
      # Categorical columns for Subgroup & Multilevel Cluster
      cat_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
      cat_cols_sub <- setdiff(cat_cols, c("Study", "Active_Comparator"))
      if (length(cat_cols_sub) > 0) {
        label_map <- c(
          "Control_Type" = "Tipe Kelompok Kontrol (Control_Type)",
          "Format" = "Format Penghantaran (Format)",
          "Intervention" = "Pendekatan Terapi (Intervention)",
          "Risk_of_Bias" = "Risiko Bias (Risk_of_Bias)",
          "Population" = "Populasi Sasaran (Population)",
          "Therapist_Level" = "Kualifikasi Terapis (Therapist_Level)",
          "Age_Group" = "Kelompok Usia (Age_Group)"
        )
        sub_choices <- setNames(cat_cols_sub, ifelse(cat_cols_sub %in% names(label_map), label_map[cat_cols_sub], cat_cols_sub))
        curr_sub <- isolate(input$subgroup_var_select)
        sel_sub <- if (!is.null(curr_sub) && curr_sub %in% cat_cols_sub) curr_sub else cat_cols_sub[1]
        updateSelectInput(session, "subgroup_var_select", choices = sub_choices, selected = sel_sub)
      }
    }
  })
  
  # Reactive Meta-Analysis Model (Supports Standard z and Hartung-Knapp knha)
  meta_model <- reactive({
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    method <- input$meta_estimator
    test_type <- ifelse(is.null(input$meta_test_type) || method == "FE", "z", input$meta_test_type)
    conf_level <- as.numeric(input$conf_level)
    covar <- input$meta_continuous_covar
    
    if (!is.null(covar) && covar != "none" && covar %in% colnames(df)) {
      form <- as.formula(paste("~", covar))
      res <- rma(yi = yi, vi = vi, mods = form, data = df, method = method, test = test_type, level = conf_level * 100)
    } else {
      res <- rma(yi = yi, vi = vi, data = df, method = method, test = test_type, level = conf_level * 100)
    }
    res
  })
  
  # Unconditional Pooled Model (for Clinical Translation, Sensitivity, and Forest/Funnel)
  pooled_only_model <- reactive({
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    method <- input$meta_estimator
    test_type <- ifelse(is.null(input$meta_test_type) || method == "FE", "z", input$meta_test_type)
    conf_level <- as.numeric(input$conf_level)
    rma(yi = yi, vi = vi, data = df, method = method, test = test_type, level = conf_level * 100)
  })
  
  # Tab 5: Value Boxes
  output$vb_pooled_es <- renderText({
    res <- tryCatch(pooled_only_model(), error = function(e) NULL)
    if (is.null(res)) return("-")
    sprintf("%.3f [%.3f, %.3f]", res$b[1], res$ci.lb[1], res$ci.ub[1])
  })
  output$vb_hetero_i2 <- renderText({
    res <- tryCatch(pooled_only_model(), error = function(e) NULL)
    if (is.null(res)) return("-")
    if (input$meta_estimator == "FE") return("N/A (FE)")
    sprintf("%.3f%% (Q = %.3f, p = %.3f)", res$I2, res$QE, res$QEp)
  })
  output$vb_hetero_tau2 <- renderText({
    res <- tryCatch(pooled_only_model(), error = function(e) NULL)
    if (is.null(res)) return("-")
    if (input$meta_estimator == "FE") return("0 (FE)")
    sprintf("%.3f (\\tau = %.3f)", res$tau2, sqrt(res$tau2))
  })
  output$vb_pred_interval <- renderText({
    res <- tryCatch(pooled_only_model(), error = function(e) NULL)
    if (is.null(res) || input$meta_estimator == "FE") return("N/A (FE Model)")
    pred <- predict(res)
    sprintf("[%.3f, %.3f]", pred$pi.lb, pred$pi.ub)
  })
  
  output$meta_model_summary <- renderPrint({
    res <- tryCatch(meta_model(), error = function(e) e)
    if (inherits(res, "error")) {
      cat("Gagal mengestimasi model meta-analisis:\n", res$message)
    } else {
      summary(res)
    }
  })
  
  output$meta_continuous_covar_ui <- renderUI({
    covar <- input$meta_continuous_covar
    df <- filtered_active_data()
    if (is.null(covar) || covar == "none" || !covar %in% colnames(df)) {
      return(
        div(
          class = "alert alert-light border my-2",
          h6(class = "fw-bold text-primary", icon("info-circle"), " Meta-Regresi Kontinu"),
          p("Pilih variabel kontinu (misal: ", tags$strong("Duration_Weeks"), " atau ", tags$strong("Total_Sessions"), ") pada bilah samping untuk memodelkan hubungan linear kovariat dengan besaran efek intervensi.")
        )
      )
    }
    
    tagList(
      plotOutput("meta_regplot_display", height = "420px"),
      hr(),
      verbatimTextOutput("meta_covar_print")
    )
  })
  
  output$meta_regplot_display <- renderPlot({
    df <- filtered_active_data()
    covar <- input$meta_continuous_covar
    req(covar != "none", covar %in% colnames(df))
    res <- tryCatch(meta_model(), error = function(e) NULL)
    req(!is.null(res))
    
    x_var <- df[[covar]]
    weights <- 1 / res$vi
    cex_scaled <- 1.0 + 3.0 * (weights - min(weights)) / (max(weights) - min(weights) + 1e-5)
    x_seq <- seq(min(x_var, na.rm = TRUE), max(x_var, na.rm = TRUE), length.out = 100)
    preds <- predict(res, newmods = x_seq)
    
    par(mar = c(5, 5, 3, 2), bg = "#FFFFFF")
    plot(
      x_var, res$yi,
      pch = 21,
      bg = rgb(101/255, 9/255, 21/255, 0.65),
      col = "#650915",
      cex = cex_scaled,
      xlab = paste("Kovariat Kontinu:", covar),
      ylab = "Besaran Efek Intervensi (Hedges' g)",
      main = paste("Meta-Regression Bubble Plot (Ukuran Titik = Presisi Studi 1/vi)"),
      las = 1,
      bty = "l",
      ylim = c(min(c(res$yi, preds$ci.lb)) - 0.1, max(c(res$yi, preds$ci.ub)) + 0.1)
    )
    polygon(c(x_seq, rev(x_seq)), c(preds$ci.lb, rev(preds$ci.ub)), col = rgb(197/255, 160/255, 89/255, 0.25), border = NA)
    lines(x_seq, preds$pred, col = "#650915", lwd = 2.5)
    grid(nx = NULL, ny = NULL, lty = 2, col = "#E5D3B3")
    points(x_var, res$yi, pch = 21, bg = rgb(101/255, 9/255, 21/255, 0.65), col = "#650915", cex = cex_scaled)
    abline(h = 0, lty = 3, col = "gray40")
  })
  
  output$meta_covar_print <- renderPrint({
    res <- tryCatch(meta_model(), error = function(e) e)
    if (inherits(res, "error")) cat("Gagal menjalankan meta-regresi:\n", res$message) else print(res)
  })
  
  # Baujat Plot & Leave-One-Out Sensitivity Analysis
  output$meta_baujat_plot <- renderPlot({
    df <- filtered_active_data()
    req(nrow(df) >= 3)
    res <- tryCatch(rma(yi = yi, vi = vi, slab = df$Study, data = df, method = input$meta_estimator), error = function(e) NULL)
    req(!is.null(res))
    par(mar = c(4.5, 4.5, 2.5, 1.5), bg = "#FFFFFF")
    baujat(
      res,
      xlab = "Kontribusi terhadap Heterogenitas Cochran Q",
      ylab = "Pengaruh terhadap Estimasi Gabungan (Squared Standardized Diff)",
      col = "#650915",
      bg = rgb(197/255, 160/255, 89/255, 0.45),
      bty = "l"
    )
    title(main = "Baujat Diagnostic Plot", col.main = "#650915", cex.main = 0.95)
  })
  
  output$table_leave1out <- renderTable({
    df <- filtered_active_data()
    req(nrow(df) >= 3)
    res <- tryCatch(rma(yi = yi, vi = vi, slab = df$Study, data = df, method = input$meta_estimator), error = function(e) NULL)
    req(!is.null(res))
    l1o <- tryCatch(leave1out(res), error = function(e) NULL)
    req(!is.null(l1o))
    data.frame(
      `Studi Dikecualikan` = l1o$slab,
      `Pooled g` = sprintf("%.3f", l1o$estimate),
      `95% CI` = sprintf("[%.3f, %.3f]", l1o$ci.lb, l1o$ci.ub),
      `p-value` = sprintf("%.4f", l1o$pval),
      `I^2 (%)` = sprintf("%.1f%%", l1o$I2),
      `tau^2` = sprintf("%.4f", l1o$tau2),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Multilevel Three-Level Meta-Analysis (rma.mv)
  multilevel_fit <- reactive({
    df <- filtered_active_data()
    req(nrow(df) >= 3)
    clust_col <- input$ml_cluster_var
    if (is.null(clust_col) || !(clust_col %in% colnames(df))) clust_col <- "Study"
    df$Cluster_ID <- as.factor(df[[clust_col]])
    df$Effect_ID <- seq_len(nrow(df))
    method_mv <- ifelse(input$meta_estimator %in% c("REML", "ML"), input$meta_estimator, "REML")
    test_mv <- ifelse(!is.null(input$meta_test_type) && input$meta_test_type == "knha", "t", "z")
    res_mv <- tryCatch(
      rma.mv(yi = yi, V = vi, random = ~ 1 | Cluster_ID / Effect_ID, data = df, method = method_mv, test = test_mv),
      error = function(e) NULL
    )
    if (is.null(res_mv)) return(NULL)
    # Cheung (2014) typical sampling variance
    wi <- 1 / df$vi
    k <- nrow(df)
    v_typ <- ((k - 1) * sum(wi)) / (sum(wi)^2 - sum(wi^2))
    s2_between <- res_mv$sigma2[1] # Level 3 (Between-Cluster)
    s2_within <- res_mv$sigma2[2]  # Level 2 (Within-Cluster)
    total_var <- v_typ + s2_within + s2_between
    i2_level1 <- (v_typ / total_var) * 100
    i2_level2 <- (s2_within / total_var) * 100
    i2_level3 <- (s2_between / total_var) * 100
    list(
      model = res_mv,
      clust_col = clust_col,
      n_clusters = length(unique(df$Cluster_ID)),
      v_typ = v_typ,
      s2_within = s2_within,
      s2_between = s2_between,
      i2_level1 = i2_level1,
      i2_level2 = i2_level2,
      i2_level3 = i2_level3
    )
  })
  
  output$table_multilevel_decomp <- renderTable({
    ml <- multilevel_fit()
    if (is.null(ml)) return(data.frame(Status = "Model multilevel tiga aras belum dapat diestimasi."))
    data.frame(
      `Aras Varians (Multilevel Level)` = c(
        "Level 1: Kesalahan Sampling (Sampling Error)",
        paste0("Level 2: Varians Dalam-Klaster (Within-", ml$clust_col, ")"),
        paste0("Level 3: Varians Antarklaster (Between-", ml$clust_col, ")"),
        "Total Heterogenitas Sejati (I^2 Level 2 + Level 3)"
      ),
      `Komponen Varians (\\sigma^2)` = c(
        sprintf("%.3f (Typical v)", ml$v_typ),
        sprintf("%.3f", ml$s2_within),
        sprintf("%.3f", ml$s2_between),
        sprintf("%.3f", ml$s2_within + ml$s2_between)
      ),
      `Proporsi Varians (%)` = c(
        sprintf("%.3f%%", ml$i2_level1),
        sprintf("%.3f%%", ml$i2_level2),
        sprintf("%.3f%%", ml$i2_level3),
        sprintf("%.3f%%", ml$i2_level2 + ml$i2_level3)
      ),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$multilevel_model_summary <- renderPrint({
    ml <- multilevel_fit()
    if (is.null(ml)) cat("Gagal mengestimasi rma.mv.") else summary(ml$model)
  })
  
  output$multilevel_variance_plot <- renderPlot({
    ml <- multilevel_fit()
    req(!is.null(ml))
    vals <- c(ml$i2_level1, ml$i2_level2, ml$i2_level3)
    labs <- c(
      sprintf("Level 1 (Sampling)\n%.1f%%", ml$i2_level1),
      sprintf("Level 2 (Within)\n%.1f%%", ml$i2_level2),
      sprintf("Level 3 (Between)\n%.1f%%", ml$i2_level3)
    )
    par(mar = c(4.5, 4.5, 2.5, 1.5), bg = "#FFFFFF")
    bp <- barplot(
      vals,
      names.arg = labs,
      ylim = c(0, 105),
      col = c("#E5D3B3", "#C5A059", "#650915"),
      border = "#650915",
      ylab = "Proporsi Varians Total (%)",
      main = paste0("Dekomposisi Varians 3-Level (Klaster: ", ml$clust_col, ")"),
      col.main = "#650915",
      las = 1
    )
    grid(nx = NA, ny = NULL, lty = 2, col = "#E5D3B3")
    text(bp, vals + 4, sprintf("%.1f%%", vals), font = 2, col = "#231618", cex = 0.9)
  })
  
  # Tab 6: Analisis Subkelompok (Subgroup Analysis)
  subgroup_models <- reactive({
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    grp_var <- input$subgroup_var_select
    req(grp_var %in% colnames(df))
    
    df_clean <- df[!is.na(df[[grp_var]]) & df[[grp_var]] != "", ]
    req(nrow(df_clean) >= 2)
    
    method <- input$meta_estimator
    test_type <- ifelse(is.null(input$meta_test_type) || method == "FE", "z", input$meta_test_type)
    conf_level <- if (!is.null(input$conf_level)) as.numeric(input$conf_level) else 0.95
    res_all <- tryCatch(rma(yi = yi, vi = vi, data = df_clean, method = method, test = test_type, level = conf_level * 100), error = function(e) NULL)
    
    form <- as.formula(paste("~ factor(", grp_var, ")"))
    res_mod <- tryCatch(rma(yi = yi, vi = vi, mods = form, data = df_clean, method = method, test = test_type, level = conf_level * 100), error = function(e) NULL)
    
    list(data = df_clean, grp_var = grp_var, res_all = res_all, res_mod = res_mod, method = method, test_type = test_type, conf_level = conf_level)
  })
  
  output$vb_sub_qb <- renderText({
    m <- tryCatch(subgroup_models(), error = function(e) NULL)
    if (is.null(m) || is.null(m$res_mod)) return("-")
    sprintf("%.3f (df = %d, p = %.3f)", m$res_mod$QM, m$res_mod$QMdf[1], m$res_mod$QMp)
  })
  output$vb_sub_qe <- renderText({
    m <- tryCatch(subgroup_models(), error = function(e) NULL)
    if (is.null(m) || is.null(m$res_mod)) return("-")
    sprintf("%.3f (df = %d, p = %.3f)", m$res_mod$QE, m$res_mod$QEdf, m$res_mod$QEp)
  })
  output$vb_sub_r2 <- renderText({
    m <- tryCatch(subgroup_models(), error = function(e) NULL)
    if (is.null(m) || is.null(m$res_mod)) return("-")
    if (input$meta_estimator == "FE") return("N/A (FE)")
    if (!is.null(m$res_mod$R2) && !is.na(m$res_mod$R2)) sprintf("%.1f%%", m$res_mod$R2) else "0.0%"
  })
  output$vb_sub_k_cats <- renderText({
    m <- tryCatch(subgroup_models(), error = function(e) NULL)
    if (is.null(m)) return("-")
    cats <- unique(m$data[[m$grp_var]])
    sprintf("%d Kategori (k = %d)", length(cats), nrow(m$data))
  })
  
  # Helper: Publication-Grade meta::forest Style Renderer (Author, Year, g, SE | Standardised Mean Difference | SMD, 95%-CI, Weight)
  draw_meta_publication_forest <- function(df, res, show_pi = TRUE, show_weights = TRUE, cex_val = 0.95, method = "REML", test_type = "z") {
    k <- nrow(df)
    sei_vec <- if ("sei" %in% colnames(df)) df$sei else sqrt(df$vi)
    yr_vec <- if ("Year" %in% colnames(df)) as.character(df$Year) else rep("-", k)
    
    # In metafor, res$level is alpha (e.g. 0.05 for 95% CI)
    alpha_ci <- ifelse(!is.null(res$level) && res$level < 0.5, res$level, 0.05)
    ci_pct <- round((1 - alpha_ci) * 100)
    z_crit <- qnorm(1 - alpha_ci / 2)
    ci_lb <- df$yi - z_crit * sei_vec
    ci_ub <- df$yi + z_crit * sei_vec
    w_raw <- tryCatch(weights(res), error = function(e) (1 / df$vi) / sum(1 / df$vi) * 100)
    
    # Prediction Interval matching meta::metagen & metafor::predict.rma
    if (k >= 3 && method != "FE") {
      t_pi <- qt(1 - alpha_ci / 2, df = k - 2)
      se_pi <- sqrt(res$tau2 + res$se[1]^2)
      pi_lb <- res$b[1] - t_pi * se_pi
      pi_ub <- res$b[1] + t_pi * se_pi
    } else {
      pred <- tryCatch(predict(res), error = function(e) list(pi.lb = NA_real_, pi.ub = NA_real_))
      pi_lb <- pred$pi.lb
      pi_ub <- pred$pi.ub
    }
    
    has_pi <- isTRUE(show_pi) && method != "FE" && !is.na(pi_lb) && !is.na(pi_ub)
    
    core_x <- c(df$yi, res$ci.lb[1], res$ci.ub[1], 0)
    if (has_pi) core_x <- c(core_x, pi_lb, pi_ub)
    max_core <- max(1.0, max(abs(core_x), na.rm = TRUE))
    ticks <- pretty(c(-max_core, max_core), n = 6)
    x_lim <- range(ticks)
    
    map_x <- function(val) {
      0.460 + (val - x_lim[1]) / (x_lim[2] - x_lim[1]) * (0.695 - 0.460)
    }
    
    fmt_ci <- function(lb, ub) {
      s_lb <- ifelse(lb >= 0, sprintf(" %.2f", lb), sprintf("%.2f", lb))
      s_ub <- ifelse(ub >= 0, sprintf(" %.2f", ub), sprintf("%.2f", ub))
      paste0("[", s_lb, "; ", s_ub, "]")
    }
    
    # Row y-coordinates (from top to bottom)
    y_header_top <- k + 5.25
    y_header <- k + 4.45
    y_studies <- rev(seq_len(k)) + 3.0
    y_re <- 2.75
    y_pi <- if (has_pi) 1.85 else NA_real_
    y_het <- if (has_pi) 0.95 else 1.85
    y_axis <- if (has_pi) 0.95 else 1.65
    y_max <- k + 6.1
    
    par(mar = c(1.2, 0.4, 0.6, 0.4), bg = "#FFFFFF", family = "sans")
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(-0.4, y_max), xaxs = "i", yaxs = "i")
    
    # Column X positions
    x_auth <- 0.015
    x_year <- 0.255
    x_g    <- 0.340
    x_se   <- 0.420
    x_plot_mid <- (0.455 + 0.705) / 2
    x_smd  <- 0.775
    x_ci   <- if (isTRUE(show_weights)) 0.905 else 0.975
    x_wt   <- 0.985
    
    # Headers
    text(x_auth, y_header, "Author", adj = c(0, 0.5), font = 2, cex = cex_val)
    text(x_year, y_header, "Year", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_g, y_header, "g", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_se, y_header, "SE", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_plot_mid, y_header_top, "Standardised Mean", adj = c(0.5, 0.5), font = 2, cex = cex_val)
    text(x_plot_mid, y_header, "Difference", adj = c(0.5, 0.5), font = 2, cex = cex_val)
    text(x_smd, y_header, "SMD", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_ci, y_header, paste0(ci_pct, "%-CI"), adj = c(1, 0.5), font = 2, cex = cex_val)
    if (isTRUE(show_weights)) {
      text(x_wt, y_header, "Weight", adj = c(1, 0.5), font = 2, cex = cex_val)
    }
    
    # Reference vertical line at 0 and dotted pooled line
    x0_px <- map_x(0)
    xp_px <- map_x(res$b[1])
    segments(x0_px, y_axis, x0_px, max(y_studies) + 0.65, col = "black", lwd = 1.2)
    segments(xp_px, y_re + 0.32, xp_px, max(y_studies) + 0.65, col = "black", lty = 3, lwd = 1.1)
    
    # Study Rows
    w_sqrt <- sqrt(pmax(w_raw, 1))
    w_norm <- 0.14 + 0.22 * (w_sqrt - min(w_sqrt)) / (max(w_sqrt) - min(w_sqrt) + 1e-6)
    
    for (i in seq_len(k)) {
      yi_r <- y_studies[i]
      text(x_auth, yi_r, df$Study[i], adj = c(0, 0.5), cex = cex_val)
      text(x_year, yi_r, yr_vec[i], adj = c(1, 0.5), cex = cex_val)
      text(x_g, yi_r, sprintf("%.4f", df$yi[i]), adj = c(1, 0.5), cex = cex_val)
      text(x_se, yi_r, sprintf("%.4f", sei_vec[i]), adj = c(1, 0.5), cex = cex_val)
      
      # Square & CI in center panel (square first, then CI line & center tick)
      xc <- map_x(df$yi[i])
      xl <- map_x(ci_lb[i])
      xr <- map_x(ci_ub[i])
      sq_half_y <- w_norm[i]
      sq_half_x <- sq_half_y * 0.032
      
      rect(xc - sq_half_x, yi_r - sq_half_y, xc + sq_half_x, yi_r + sq_half_y, col = "#C4C4C4", border = "#B0B0B0")
      segments(xl, yi_r, xr, yi_r, col = "black", lwd = 1.25)
      segments(xc, yi_r - 0.09, xc, yi_r + 0.09, col = "black", lwd = 1.4)
      
      # Right columns
      text(x_smd, yi_r, sprintf("%.2f", df$yi[i]), adj = c(1, 0.5), cex = cex_val)
      text(x_ci, yi_r, fmt_ci(ci_lb[i], ci_ub[i]), adj = c(1, 0.5), cex = cex_val)
      if (isTRUE(show_weights)) {
        text(x_wt, yi_r, sprintf("%.1f%%", w_raw[i]), adj = c(1, 0.5), cex = cex_val)
      }
    }
    
    # Random / Fixed effects summary row
    re_label <- if (method == "FE") {
      "Common effect model"
    } else if (test_type == "knha") {
      "Random effects model (HK)"
    } else {
      paste0("Random effects model (", method, ")")
    }
    text(x_auth, y_re, re_label, adj = c(0, 0.5), font = 2, cex = cex_val)
    
    # Summary Diamond
    d_xc <- map_x(res$b[1])
    d_xl <- map_x(res$ci.lb[1])
    d_xr <- map_x(res$ci.ub[1])
    polygon(
      x = c(d_xl, d_xc, d_xr, d_xc),
      y = c(y_re, y_re + 0.30, y_re, y_re - 0.30),
      col = "#BEBEBE",
      border = "black",
      lwd = 1.3
    )
    
    text(x_smd, y_re, sprintf("%.2f", res$b[1]), adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_ci, y_re, fmt_ci(res$ci.lb[1], res$ci.ub[1]), adj = c(1, 0.5), font = 2, cex = cex_val)
    if (isTRUE(show_weights)) {
      text(x_wt, y_re, "100.0%", adj = c(1, 0.5), font = 2, cex = cex_val)
    }
    
    # Prediction Interval row (Red bar)
    if (has_pi) {
      text(x_auth, y_pi, "Prediction interval", adj = c(0, 0.5), font = 2, cex = cex_val)
      pi_xl <- map_x(pi_lb)
      pi_xr <- map_x(pi_ub)
      rect(pi_xl, y_pi - 0.11, pi_xr, y_pi + 0.11, col = "#EE0000", border = "black", lwd = 1.15)
      text(x_ci, y_pi, fmt_ci(pi_lb, pi_ub), adj = c(1, 0.5), font = 2, cex = cex_val)
    }
    
    # Heterogeneity footer (mengikuti I^2 estimator aktif dengan fallback Cochran Q)
    q_val <- ifelse(is.null(res$QE) || is.na(res$QE), 0, res$QE)
    i2_q <- if (!is.null(res$I2) && !is.na(res$I2)) res$I2 else if (q_val > (k - 1) && q_val > 0) ((q_val - (k - 1)) / q_val) * 100 else 0
    tau2_val <- ifelse(is.null(res$tau2) || is.na(res$tau2), 0, res$tau2)
    p_het <- ifelse(is.null(res$QEp) || is.na(res$QEp), 1, res$QEp)
    het_expr <- bquote(
      paste(
        "Heterogeneity: ",
        italic(I)^2, " = ", .(sprintf("%.1f%%", i2_q)), ", ",
        tau^2, " = ", .(sprintf("%.4f", tau2_val)), ", ",
        italic(p), " = ", .(sprintf("%.4f", p_het))
      )
    )
    text(x_auth, y_het, het_expr, adj = c(0, 0.5), cex = cex_val * 0.92)
    
    # Bottom Axis under Center Plot
    ax_left <- map_x(x_lim[1])
    ax_right <- map_x(x_lim[2])
    segments(ax_left, y_axis, ax_right, y_axis, col = "black", lwd = 1.4)
    for (tk in ticks) {
      xtk <- map_x(tk)
      segments(xtk, y_axis, xtk, y_axis - 0.26, col = "black", lwd = 1.3)
      text(xtk, y_axis - 0.65, as.character(tk), adj = c(0.5, 0.5), cex = cex_val * 0.95)
    }
  }
  
  # Helper: Subgroup Forest Plot in meta::forest Publication Style
  draw_meta_subgroup_forest <- function(m, cex_val = 0.90) {
    df <- m$data
    grp_var <- m$grp_var
    method <- m$method
    test_type <- m$test_type
    res_all <- m$res_all
    res_mod <- m$res_mod
    
    ord_idx <- order(df[[grp_var]], df$yi, decreasing = FALSE)
    df <- df[ord_idx, ]
    cats <- unique(df[[grp_var]])
    
    conf_level <- if (!is.null(m$conf_level)) m$conf_level else 0.95
    alpha_ci <- 1 - conf_level
    sei_vec <- if ("sei" %in% colnames(df)) df$sei else sqrt(df$vi)
    yr_vec <- if ("Year" %in% colnames(df)) as.character(df$Year) else rep("-", nrow(df))
    z_crit <- qnorm(1 - alpha_ci / 2)
    ci_lb <- df$yi - z_crit * sei_vec
    ci_ub <- df$yi + z_crit * sei_vec
    w_raw <- tryCatch(weights(res_all), error = function(e) (1 / df$vi) / sum(1 / df$vi) * 100)
    k_all <- nrow(df)
    if (k_all >= 3 && method != "FE") {
      t_pi_all <- qt(1 - alpha_ci / 2, df = k_all - 2)
      se_pi_all <- sqrt(res_all$tau2 + res_all$se[1]^2)
      pi_lb_all <- res_all$b[1] - t_pi_all * se_pi_all
      pi_ub_all <- res_all$b[1] + t_pi_all * se_pi_all
    } else {
      pred_all <- tryCatch(predict(res_all), error = function(e) list(pi.lb = res_all$ci.lb[1], pi.ub = res_all$ci.ub[1]))
      pi_lb_all <- pred_all$pi.lb
      pi_ub_all <- pred_all$pi.ub
    }
    
    core_x <- c(df$yi, res_all$ci.lb[1], res_all$ci.ub[1], pi_lb_all, pi_ub_all, 0)
    max_abs <- max(1.0, min(4.0, max(abs(core_x), na.rm = TRUE)))
    ticks <- pretty(c(-max_abs, max_abs), n = 6)
    x_lim <- range(ticks)
    
    map_x <- function(val) {
      0.460 + (val - x_lim[1]) / (x_lim[2] - x_lim[1]) * (0.695 - 0.460)
    }
    fmt_ci <- function(lb, ub) {
      s_lb <- ifelse(lb >= 0, sprintf(" %.2f", lb), sprintf("%.2f", lb))
      s_ub <- ifelse(ub >= 0, sprintf(" %.2f", ub), sprintf("%.2f", ub))
      paste0("[", s_lb, "; ", s_ub, "]")
    }
    
    # Build row positions bottom-up
    row_cursor <- 4.2
    cat_het_rows <- list()
    cat_poly_rows <- list()
    cat_header_rows <- list()
    study_rows <- numeric(nrow(df))
    
    for (g_name in rev(cats)) {
      idx <- which(df[[grp_var]] == g_name)
      cat_het_rows[[g_name]] <- row_cursor
      row_cursor <- row_cursor + 0.85
      cat_poly_rows[[g_name]] <- row_cursor
      row_cursor <- row_cursor + 0.95
      for (i in rev(seq_along(idx))) {
        study_rows[idx[i]] <- row_cursor
        row_cursor <- row_cursor + 0.90
      }
      cat_header_rows[[g_name]] <- row_cursor
      row_cursor <- row_cursor + 1.35
    }
    
    y_header <- row_cursor + 0.4
    y_header_top <- y_header + 0.8
    y_re <- 2.85
    y_pi <- 2.00
    y_het <- 1.15
    y_qb <- 0.35
    y_axis <- 1.15
    
    par(mar = c(1.2, 0.4, 0.6, 0.4), bg = "#FFFFFF", family = "sans")
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(-0.8, y_header_top + 0.8), xaxs = "i", yaxs = "i")
    
    x_auth <- 0.015
    x_year <- 0.255
    x_g    <- 0.340
    x_se   <- 0.420
    x_plot_mid <- (0.455 + 0.705) / 2
    x_smd  <- 0.775
    x_ci   <- 0.905
    x_wt   <- 0.985
    
    text(x_auth, y_header, "Author", adj = c(0, 0.5), font = 2, cex = cex_val)
    text(x_year, y_header, "Year", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_g, y_header, "g", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_se, y_header, "SE", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_plot_mid, y_header_top, "Standardised Mean", adj = c(0.5, 0.5), font = 2, cex = cex_val)
    text(x_plot_mid, y_header, "Difference", adj = c(0.5, 0.5), font = 2, cex = cex_val)
    text(x_smd, y_header, "SMD", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_ci, y_header, "95%-CI", adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_wt, y_header, "Weight", adj = c(1, 0.5), font = 2, cex = cex_val)
    
    x0_px <- map_x(0)
    xp_px <- map_x(res_all$b[1])
    segments(x0_px, y_axis, x0_px, y_header - 0.5, col = "black", lwd = 1.2)
    segments(xp_px, y_re + 0.28, xp_px, y_header - 0.5, col = "black", lty = 3, lwd = 1.1)
    
    w_sqrt <- sqrt(pmax(w_raw, 1))
    w_norm <- 0.12 + 0.18 * (w_sqrt - min(w_sqrt)) / (max(w_sqrt) - min(w_sqrt) + 1e-6)
    
    for (i in seq_len(nrow(df))) {
      yi_r <- study_rows[i]
      text(x_auth, yi_r, df$Study[i], adj = c(0, 0.5), cex = cex_val)
      text(x_year, yi_r, yr_vec[i], adj = c(1, 0.5), cex = cex_val)
      text(x_g, yi_r, sprintf("%.4f", df$yi[i]), adj = c(1, 0.5), cex = cex_val)
      text(x_se, yi_r, sprintf("%.4f", sei_vec[i]), adj = c(1, 0.5), cex = cex_val)
      
      xc <- map_x(df$yi[i])
      xl <- map_x(ci_lb[i])
      xr <- map_x(ci_ub[i])
      sq_half_y <- w_norm[i]
      sq_half_x <- sq_half_y * 0.030
      rect(xc - sq_half_x, yi_r - sq_half_y, xc + sq_half_x, yi_r + sq_half_y, col = "#C4C4C4", border = "#B0B0B0")
      segments(xl, yi_r, xr, yi_r, col = "black", lwd = 1.2)
      segments(xc, yi_r - 0.08, xc, yi_r + 0.08, col = "black", lwd = 1.3)
      
      text(x_smd, yi_r, sprintf("%.2f", df$yi[i]), adj = c(1, 0.5), cex = cex_val)
      text(x_ci, yi_r, fmt_ci(ci_lb[i], ci_ub[i]), adj = c(1, 0.5), cex = cex_val)
      text(x_wt, yi_r, sprintf("%.1f%%", w_raw[i]), adj = c(1, 0.5), cex = cex_val)
    }
    
    for (g_name in cats) {
      idx <- which(df[[grp_var]] == g_name)
      text(x_auth, cat_header_rows[[g_name]], paste0(grp_var, " = ", g_name), adj = c(0, 0.5), font = 2, cex = cex_val)
      sub_test <- ifelse(length(idx) >= 2, test_type, "z")
      sub_res <- tryCatch(rma(df$yi[idx], df$vi[idx], method = method, test = sub_test, level = conf_level * 100), error = function(e) NULL)
      if (!is.null(sub_res)) {
        yp <- cat_poly_rows[[g_name]]
        yh <- cat_het_rows[[g_name]]
        text(x_auth, yp, "Random effects model", adj = c(0, 0.5), font = 2, cex = cex_val)
        d_xc <- map_x(sub_res$b[1])
        d_xl <- map_x(sub_res$ci.lb[1])
        d_xr <- map_x(sub_res$ci.ub[1])
        polygon(c(d_xl, d_xc, d_xr, d_xc), c(yp, yp + 0.25, yp, yp - 0.25), col = "#BEBEBE", border = "black", lwd = 1.2)
        text(x_smd, yp, sprintf("%.2f", sub_res$b[1]), adj = c(1, 0.5), font = 2, cex = cex_val)
        text(x_ci, yp, fmt_ci(sub_res$ci.lb[1], sub_res$ci.ub[1]), adj = c(1, 0.5), font = 2, cex = cex_val)
        text(x_wt, yp, sprintf("%.1f%%", sum(w_raw[idx])), adj = c(1, 0.5), font = 2, cex = cex_val)
        
        sub_het <- bquote(
          paste(
            "Heterogeneity: ",
            italic(I)^2, " = ", .(sprintf("%.1f%%", ifelse(is.na(sub_res$I2), 0, sub_res$I2))), ", ",
            tau^2, " = ", .(sprintf("%.4f", ifelse(is.na(sub_res$tau2), 0, sub_res$tau2))), ", ",
            italic(p), " = ", .(sprintf("%.4f", ifelse(is.na(sub_res$QEp), 1, sub_res$QEp)))
          )
        )
        text(x_auth, yh, sub_het, adj = c(0, 0.5), cex = cex_val * 0.88)
      }
    }
    
    re_label <- ifelse(test_type == "knha", "Random effects model (HK)", paste0("Random effects model (", method, ")"))
    text(x_auth, y_re, re_label, adj = c(0, 0.5), font = 2, cex = cex_val)
    d_xc <- map_x(res_all$b[1])
    d_xl <- map_x(res_all$ci.lb[1])
    d_xr <- map_x(res_all$ci.ub[1])
    polygon(c(d_xl, d_xc, d_xr, d_xc), c(y_re, y_re + 0.28, y_re, y_re - 0.28), col = "#BEBEBE", border = "black", lwd = 1.3)
    text(x_smd, y_re, sprintf("%.2f", res_all$b[1]), adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_ci, y_re, fmt_ci(res_all$ci.lb[1], res_all$ci.ub[1]), adj = c(1, 0.5), font = 2, cex = cex_val)
    text(x_wt, y_re, "100.0%", adj = c(1, 0.5), font = 2, cex = cex_val)
    
    if (method != "FE" && !is.na(pi_lb_all)) {
      text(x_auth, y_pi, "Prediction interval", adj = c(0, 0.5), font = 2, cex = cex_val)
      rect(map_x(pi_lb_all), y_pi - 0.10, map_x(pi_ub_all), y_pi + 0.10, col = "#EE0000", border = "black", lwd = 1.1)
      text(x_ci, y_pi, fmt_ci(pi_lb_all, pi_ub_all), adj = c(1, 0.5), font = 2, cex = cex_val)
    }
    
    q_all <- ifelse(is.null(res_all$QE) || is.na(res_all$QE), 0, res_all$QE)
    i2_q_all <- if (!is.null(res_all$I2) && !is.na(res_all$I2)) res_all$I2 else if (q_all > (k_all - 1) && q_all > 0) ((q_all - (k_all - 1)) / q_all) * 100 else 0
    het_all <- bquote(
      paste(
        "Heterogeneity: ",
        italic(I)^2, " = ", .(sprintf("%.1f%%", i2_q_all)), ", ",
        tau^2, " = ", .(sprintf("%.4f", ifelse(is.na(res_all$tau2), 0, res_all$tau2))), ", ",
        italic(p), " = ", .(sprintf("%.4f", ifelse(is.na(res_all$QEp), 1, res_all$QEp)))
      )
    )
    text(x_auth, y_het, het_all, adj = c(0, 0.5), cex = cex_val * 0.88)
    
    if (!is.null(res_mod)) {
      qb_expr <- bquote(
        paste(
          "Test for subgroup differences: ",
          chi[.(res_mod$QMdf[1])]^2, " = ", .(sprintf("%.2f", res_mod$QM)), ", ",
          "df = ", .(res_mod$QMdf[1]), " (", italic(p), " = ", .(sprintf("%.4f", res_mod$QMp)), ")"
        )
      )
      text(x_auth, y_qb, qb_expr, adj = c(0, 0.5), cex = cex_val * 0.88)
    }
    
    ax_left <- map_x(x_lim[1])
    ax_right <- map_x(x_lim[2])
    segments(ax_left, y_axis, ax_right, y_axis, col = "black", lwd = 1.4)
    for (tk in ticks) {
      xtk <- map_x(tk)
      segments(xtk, y_axis, xtk, y_axis - 0.24, col = "black", lwd = 1.2)
      text(xtk, y_axis - 0.58, as.character(tk), adj = c(0.5, 0.5), cex = cex_val * 0.92)
    }
  }
  
  # Subgroup Forest Plot Internal Renderer
  render_subgroup_forest_internal <- function() {
    m <- subgroup_models()
    req(!is.null(m), !is.null(m$res_all))
    
    if (is.null(input$subgroup_forest_style) || input$subgroup_forest_style == "meta_style") {
      draw_meta_subgroup_forest(m, cex_val = input$subgroup_forest_cex)
      return(invisible(NULL))
    }
    
    df <- m$data
    grp_var <- m$grp_var
    method <- m$method
    
    ord_idx <- order(df[[grp_var]], df$yi, decreasing = FALSE)
    df <- df[ord_idx, ]
    cats <- unique(df[[grp_var]])
    
    row_cursor <- 3
    cat_poly_rows <- list()
    cat_header_rows <- list()
    study_rows <- numeric(nrow(df))
    
    for (g_name in rev(cats)) {
      idx <- which(df[[grp_var]] == g_name)
      cat_poly_rows[[g_name]] <- row_cursor
      row_cursor <- row_cursor + 1
      
      for (i in seq_along(idx)) {
        study_rows[idx[i]] <- row_cursor
        row_cursor <- row_cursor + 1
      }
      
      cat_header_rows[[g_name]] <- row_cursor
      row_cursor <- row_cursor + 2
    }
    
    top_row <- row_cursor
    
    if (grp_var == "Control_Type" && "Active_Comparator" %in% colnames(df)) {
      slabs <- paste0(df$Study, " (", df$Year, ") [", ifelse(df$Control_Type == "Active Control" & !is.na(df$Active_Comparator), df$Active_Comparator, df$Control_Type), "]")
    } else if ("Year" %in% colnames(df)) {
      slabs <- paste0(df$Study, " (", df$Year, ")")
    } else {
      slabs <- df$Study
    }
    
    cex_val <- input$subgroup_forest_cex
    
    par(mar = c(5, 5, 3, 2), bg = "#FFFFFF")
    forest(
      df$yi,
      df$vi,
      slab = slabs,
      rows = study_rows,
      ylim = c(-3.0, top_row + 1.5),
      xlab = "Besaran Efek Intervensi (Hedges' g)",
      header = c("Studi Primer & Subkelompok", "Hedges' g [95% CI]"),
      cex = cex_val,
      col = "#650915",
      refline = 0,
      alim = c(min(-0.5, min(df$yi - 1.96 * sqrt(df$vi)) - 0.2), max(1.5, max(df$yi + 1.96 * sqrt(df$vi)) + 0.2))
    )
    
    addpoly(
      m$res_all,
      row = 1,
      mlab = paste0("Total Pooled (", method, ")"),
      cex = cex_val,
      col = "#650915",
      border = "#650915"
    )
    
    for (g_name in cats) {
      idx <- which(df[[grp_var]] == g_name)
      k_cat <- length(idx)
      sub_test <- ifelse(k_cat >= 2, m$test_type, "z")
      conf_lvl <- if (!is.null(m$conf_level)) m$conf_level * 100 else 95
      sub_res <- tryCatch(rma(df$yi[idx], df$vi[idx], method = method, test = sub_test, level = conf_lvl), error = function(e) NULL)
      if (!is.null(sub_res)) {
        addpoly(
          sub_res,
          row = cat_poly_rows[[g_name]],
          mlab = paste0("Subtotal: ", g_name, " (k = ", k_cat, ")"),
          cex = cex_val,
          col = "#C5A059",
          border = "#650915"
        )
      }
      
      text(
        par("usr")[1],
        cat_header_rows[[g_name]],
        paste0("Subkelompok: ", g_name, " (k = ", k_cat, ")"),
        font = 2,
        pos = 4,
        cex = cex_val,
        col = "#650915"
      )
    }
    
    if (!is.null(m$res_mod)) {
      text(
        par("usr")[1],
        -1.8,
        sprintf("Uji Beda Subkelompok (QB): QM = %.2f, df = %d, p = %.4f  |  Heterogenitas Residual: QE = %.2f, p = %.4f",
                m$res_mod$QM, m$res_mod$QMdf[1], m$res_mod$QMp, m$res_mod$QE, m$res_mod$QEp),
        pos = 4,
        cex = cex_val * 0.85,
        font = 3,
        col = "#5A4B4C"
      )
    }
  }
  
  output$subgroup_forest_plot_display <- renderPlot({
    render_subgroup_forest_internal()
  })
  
  output$btn_dl_subgroup_png <- downloadHandler(
    filename = function() paste0("Subgroup_Forest_Plot_", input$subgroup_var_select, "_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 3000, height = 2400, res = 300)
      render_subgroup_forest_internal()
      dev.off()
    }
  )
  
  output$btn_dl_subgroup_pdf <- downloadHandler(
    filename = function() paste0("Subgroup_Forest_Plot_", input$subgroup_var_select, "_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 11, height = 8.5)
      render_subgroup_forest_internal()
      dev.off()
    }
  )
  
  subgroup_summary_df <- reactive({
    m <- tryCatch(subgroup_models(), error = function(e) NULL)
    if (is.null(m)) return(data.frame(Status = "Data tidak mencukupi untuk analisis subkelompok."))
    
    df <- m$data
    grp_var <- m$grp_var
    method <- m$method
    cats <- sort(unique(df[[grp_var]]))
    
    res_list <- lapply(cats, function(cat) {
      sub_df <- df[df[[grp_var]] == cat, ]
      k_sub <- nrow(sub_df)
      total_n <- if ("N_Intervention" %in% colnames(sub_df) && "N_Control" %in% colnames(sub_df)) {
        as.character(sum(sub_df$N_Intervention + sub_df$N_Control, na.rm = TRUE))
      } else {
        "-"
      }
      sub_test <- ifelse(k_sub >= 2, m$test_type, "z")
      conf_lvl <- if (!is.null(m$conf_level)) m$conf_level * 100 else 95
      sub_res <- tryCatch(rma(yi = yi, vi = vi, data = sub_df, method = method, test = sub_test, level = conf_lvl), error = function(e) NULL)
      if (is.null(sub_res)) return(NULL)
      pred <- predict(sub_res)
      
      data.frame(
        Subkelompok = as.character(cat),
        `k (Studi)` = as.character(k_sub),
        `Total N` = total_n,
        `Pooled Hedges' g` = sprintf("%.3f", sub_res$b[1]),
        SE = sprintf("%.3f", sub_res$se[1]),
        `95% CI` = sprintf("[%.3f, %.3f]", sub_res$ci.lb[1], sub_res$ci.ub[1]),
        `95% PI` = if (method != "FE" && k_sub >= 3) sprintf("[%.3f, %.3f]", pred$pi.lb, pred$pi.ub) else "-",
        `p-value` = sprintf("%.4f", sub_res$pval[1]),
        `I^2 (%)` = sprintf("%.1f%%", sub_res$I2),
        `tau^2` = sprintf("%.4f", sub_res$tau2),
        `Heterogenitas Qw` = sprintf("%.2f (df = %d, p = %.4f)", sub_res$QE, sub_res$QEdf, sub_res$QEp),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, res_list)
  })
  
  output$table_subgroup_comprehensive <- renderTable({
    subgroup_summary_df()
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$btn_dl_subgroup_csv <- downloadHandler(
    filename = function() paste0("Subgroup_Comparison_", input$subgroup_var_select, "_", Sys.Date(), ".csv"),
    content = function(file) write.csv(subgroup_summary_df(), file, row.names = FALSE)
  )
  
  observeEvent(input$btn_copy_sub_apa, {
    m <- tryCatch(subgroup_models(), error = function(e) NULL)
    if (is.null(m) || is.null(m$res_mod)) return()
    res_mod <- m$res_mod
    txt <- sprintf(
      "Analisis moderator subkelompok pada variabel %s menunjukkan statistik uji beda antarsubkelompok QB(df = %d) = %.2f, p = %.4f (R^2 = %.1f%%), dengan heterogenitas residual QE(df = %d) = %.2f, p = %.4f.",
      m$grp_var, res_mod$QMdf[1], res_mod$QM, res_mod$QMp,
      ifelse(!is.null(res_mod$R2) && !is.na(res_mod$R2), res_mod$R2, 0),
      res_mod$QEdf, res_mod$QE, res_mod$QEp
    )
    session$sendCustomMessage("copy_to_clipboard", list(text = txt))
    showNotification("Naratif subkelompok APA 7th berhasil disalin ke Clipboard.", type = "message")
  })
  
  output$subgroup_interpretation_narrative <- renderUI({
    m <- tryCatch(subgroup_models(), error = function(e) NULL)
    if (is.null(m) || is.null(m$res_mod)) return(p("Model subkelompok belum tersedia."))
    
    df <- m$data
    grp_var <- m$grp_var
    res_mod <- m$res_mod
    cats <- sort(unique(df[[grp_var]]))
    qb_sig <- res_mod$QMp < 0.05
    
    qb_text <- if (qb_sig) {
      paste0(
        "Hasil uji omnibus moderator membuktikan adanya variasi besaran efek antarsubkelompok yang bermakna secara statistik (",
        "<em>QB</em> = ", sprintf("%.2f, df = %d, p = %.4f", res_mod$QM, res_mod$QMdf[1], res_mod$QMp), "). ",
        "Temuan ini mengindikasikan bahwa karakteristik <strong>", grp_var,
        "</strong> berperan sebagai moderator substantif yang menjelaskan heterogenitas efektivitas intervensi antarstudi klinis, ",
        "dengan proporsi varians sejati antarstudi yang berhasil dijelaskan mencapai <strong>",
        sprintf("%.1f%%", ifelse(!is.null(res_mod$R2) && !is.na(res_mod$R2), res_mod$R2, 0)), "</strong>."
      )
    } else {
      paste0(
        "Hasil uji omnibus moderator mengindikasikan bahwa perbedaan besaran efek antarsubkelompok belum mencapai signifikansi statistik konvensional (",
        "<em>QB</em> = ", sprintf("%.2f, df = %d, p = %.4f", res_mod$QM, res_mod$QMdf[1], res_mod$QMp), "). ",
        "Meskipun demikian, variasi parameter deskriptif antarkategori tetap memberikan gambaran diferensial yang penting secara klinis."
      )
    }
    
    sub_details <- lapply(cats, function(cat) {
      sub_df <- df[df[[grp_var]] == cat, ]
      k_sub <- nrow(sub_df)
      sub_test <- ifelse(k_sub >= 2, m$test_type, "z")
      conf_lvl <- if (!is.null(m$conf_level)) m$conf_level * 100 else 95
      sub_res <- tryCatch(rma(yi = yi, vi = vi, data = sub_df, method = m$method, test = sub_test, level = conf_lvl), error = function(e) NULL)
      if (is.null(sub_res)) return(NULL)
      mag <- ifelse(abs(sub_res$b[1]) < 0.2, "dapat diabaikan", ifelse(abs(sub_res$b[1]) < 0.5, "kecil", ifelse(abs(sub_res$b[1]) < 0.8, "sedang", "besar")))
      tags$li(
        strong(paste0(cat, " (k = ", k_sub, "): ")),
        sprintf("Pooled Hedges' g = %.3f [%.3f, %.3f], p = %.4f (kategori efek %s). ", sub_res$b[1], sub_res$ci.lb[1], sub_res$ci.ub[1], sub_res$pval[1], mag),
        sprintf("Tingkat heterogenitas internal I^2 = %.1f%% (tau^2 = %.4f).", sub_res$I2, sub_res$tau2)
      )
    })
    
    psy_context <- if (grp_var == "Control_Type") {
      p(
        "Secara khusus pada telaah tipe kelompok kontrol, variasi besaran efek mencerminkan fenomena atenuasi akibat komparator aktif: ",
        "uji coba dengan kelompok pembanding ", tags$em("Waitlist"), " secara konsisten menghasilkan estimasi besaran efek intervensi yang paling tinggi. ",
        "Sebaliknya, pengujian terhadap ", tags$em("Active Control"), " (misal: ", tags$em("Relaxation Training"), " atau ", tags$em("Supportive Counseling"), ") ",
        "menghasilkan estimasi besaran efek yang lebih teratenuasi secara substansial. ",
        "Hal tersebut menegaskan bahwa kelompok kontrol aktif berhasil mengendalikan faktor-faktor terapeutik non-spesifik seperti ekspektasi pemulihan dan kontak relasional terapis (Cuijpers et al., 2023)."
      )
    } else if (grp_var == "Format") {
      p(
        "Secara khusus pada telaah format penghantaran intervensi, format ", tags$em("Individual"), " dan ", tags$em("Group"), " tatap muka ",
        "memperlihatkan besaran efek yang relatif lebih tinggi dibandingkan format ", tags$em("Guided Digital"), ". ",
        "Meskipun intervensi digital terpandu menawarkan skalabilitas implementasi yang luas, perbedaan intensitas relasi terapeutik (", tags$em("therapeutic alliance"), ") ",
        "memberikan kontribusi nyata terhadap variabilitas respons klinis partisipan."
      )
    } else {
      NULL
    }
    
    tagList(
      h5("Sintesis Komparatif Subkelompok & Implikasi Klinis", class = "text-primary fw-bold"),
      p(HTML(qb_text)),
      h6("Rincian Estimasi per Subkelompok:", class = "fw-bold mt-3"),
      tags$ul(sub_details),
      psy_context,
      p(
        "Evaluasi terhadap heterogenitas residual menghasilkan ",
        strong(sprintf("QE = %.2f (df = %d, p = %.4f)", res_mod$QE, res_mod$QEdf, res_mod$QEp)), ". ",
        ifelse(res_mod$QEp < 0.05,
               "Heterogenitas residual yang masih signifikan mengindikasikan adanya faktor perancu atau karakteristik metodologis lain yang belum sepenuhnya terakomodasi dalam pengelompokan ini.",
               "Nilai p residual yang tidak signifikan mengindikasikan bahwa sebagian besar heterogenitas antarstudi telah berhasil dieliminasi melalui stratifikasi subkelompok ini.")
      )
    )
  })
  
  # Tab 7: Forest Plot Internal Renderer
  render_forest_internal <- function() {
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    
    if (input$forest_sort == "yi_asc") {
      df <- df[order(df$yi), ]
    } else if (input$forest_sort == "yi_desc") {
      df <- df[order(df$yi, decreasing = TRUE), ]
    } else if (input$forest_sort == "year" && "Year" %in% colnames(df)) {
      df <- df[order(df$Year), ]
    } else if (input$forest_sort == "control" && "Control_Type" %in% colnames(df)) {
      df <- df[order(df$Control_Type), ]
    } else if (input$forest_sort == "format" && "Format" %in% colnames(df)) {
      df <- df[order(df$Format), ]
    } else if (input$forest_sort == "interv" && "Intervention" %in% colnames(df)) {
      df <- df[order(df$Intervention), ]
    } else if (input$forest_sort == "rob" && "Risk_of_Bias" %in% colnames(df)) {
      df <- df[order(df$Risk_of_Bias), ]
    }
    
    test_type <- ifelse(is.null(input$meta_test_type) || input$meta_estimator == "FE", "z", input$meta_test_type)
    res <- rma(yi = yi, vi = vi, data = df, method = input$meta_estimator, test = test_type, level = as.numeric(input$conf_level) * 100)
    
    if (is.null(input$forest_style) || input$forest_style == "meta_style") {
      draw_meta_publication_forest(
        df = df,
        res = res,
        show_pi = input$forest_show_pi,
        show_weights = input$forest_show_weights,
        cex_val = input$forest_cex,
        method = input$meta_estimator,
        test_type = test_type
      )
      return(invisible(NULL))
    }
    
    if (input$forest_sort %in% c("format", "control", "interv", "rob")) {
      mod_tag <- switch(
        input$forest_sort,
        "format" = df$Format,
        "control" = ifelse(df$Control_Type == "Active Control" & "Active_Comparator" %in% colnames(df) & !is.na(df$Active_Comparator),
                           paste0("Active: ", df$Active_Comparator), df$Control_Type),
        "interv" = df$Intervention,
        "rob" = df$Risk_of_Bias
      )
      slabs <- paste0(df$Study, " (", df$Year, ") [", mod_tag, "]")
    } else if ("Year" %in% colnames(df)) {
      slabs <- paste0(df$Study, " (", df$Year, ")")
    } else {
      slabs <- df$Study
    }
    
    par(mar = c(5, 5, 3, 2), bg = "#FFFFFF")
    forest(
      res,
      slab = slabs,
      xlab = "Besaran Efek Intervensi (Hedges' g)",
      mlab = paste0("Model Efek Acak (", input$meta_estimator, ")"),
      header = c("Studi Primer & Tahun", "Hedges' g [95% CI]"),
      addpred = input$forest_show_pi && input$meta_estimator != "FE",
      showweights = input$forest_show_weights,
      cex = input$forest_cex,
      col = "#650915",
      ylim = c(-3.0, nrow(df) + 3.5),
      alim = c(min(-0.5, min(df$yi - 1.96 * sqrt(df$vi)) - 0.2), max(1.5, max(df$yi + 1.96 * sqrt(df$vi)) + 0.2)),
      refline = 0
    )
    
    # Safe Axis labels positioned inside plot bounds
    usr <- par("usr")
    text(usr[1] + 0.1, -2.0, input$forest_label_left, pos = 4, cex = input$forest_cex * 0.85, font = 3, col = "#650915")
    text(usr[2] - 0.1, -2.0, input$forest_label_right, pos = 2, cex = input$forest_cex * 0.85, font = 3, col = "#650915")
  }
  
  output$forest_plot_display <- renderPlot({
    render_forest_internal()
  })
  
  output$btn_dl_forest_png <- downloadHandler(
    filename = function() paste0("Forest_Plot_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 2800, height = 2000, res = 300)
      render_forest_internal()
      dev.off()
    }
  )
  
  output$btn_dl_forest_pdf <- downloadHandler(
    filename = function() paste0("Forest_Plot_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 7.5)
      render_forest_internal()
      dev.off()
    }
  )
  
  # Tab 8: Funnel Plot Internal Renderer
  render_funnel_internal <- function() {
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    res <- rma(yi = yi, vi = vi, data = df, method = input$meta_estimator)
    
    par(mar = c(5, 5, 3, 2), bg = "#FFFFFF")
    if (input$funnel_type == "contour") {
      funnel(
        res,
        level = c(90, 95, 99),
        shade = c("#FFFFFF", "#FAF8F5", "#F2EBE1"),
        refline = 0,
        legend = TRUE,
        xlab = "Besaran Efek Intervensi (Hedges' g)",
        ylab = "Galat Baku (SE)"
      )
    } else if (input$funnel_type == "trimfill") {
      tf <- tryCatch(trimfill(res), error = function(e) res)
      funnel(
        tf,
        legend = TRUE,
        xlab = "Besaran Efek Intervensi (Hedges' g) — Imputasi Trim & Fill",
        ylab = "Galat Baku (SE)"
      )
    } else {
      funnel(
        res,
        xlab = "Besaran Efek Intervensi (Hedges' g)",
        ylab = "Galat Baku (SE)"
      )
    }
  }
  
  output$funnel_plot_display <- renderPlot({
    render_funnel_internal()
  })
  
  output$btn_dl_funnel_png <- downloadHandler(
    filename = function() paste0("Funnel_Plot_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 2400, height = 1800, res = 300)
      render_funnel_internal()
      dev.off()
    }
  )
  
  output$btn_dl_funnel_pdf <- downloadHandler(
    filename = function() paste0("Funnel_Plot_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 8, height = 6)
      render_funnel_internal()
      dev.off()
    }
  )
  
  # Publication Bias Tables (Forced Computation at k.min = 5)
  output$table_egger_test <- renderTable({
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    res <- tryCatch(rma(yi = yi, vi = vi, data = df, method = input$meta_estimator), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    eg_lm <- tryCatch(regtest(res, model = "lm", predictor = "sei"), error = function(e) NULL)
    eg_rma <- tryCatch(regtest(res, model = "rma", predictor = "sei"), error = function(e) NULL)
    eg <- if (!is.null(eg_lm)) eg_lm else eg_rma
    if (is.null(eg)) return(data.frame(Status = "Uji Egger tidak dapat dikomputasi (minimal dibutuhkan 3 studi)."))
    
    stat_label <- if (!is.null(eg$dfs) && !is.na(eg$dfs)) {
      sprintf("t = %.3f (df = %d)", eg$zval, eg$dfs)
    } else {
      sprintf("z = %.3f", eg$zval)
    }
    
    data.frame(
      Metrik = c(
        "Statistik Regresi Linear Egger (k.min = 5)",
        "Nilai p (p-value)",
        "Limit Estimate (Intersep b0 saat SE -> 0)",
        "Batas 95% CI Limit Estimate",
        "Kesimpulan Bias Publikasi"
      ),
      Nilai = c(
        stat_label,
        sprintf("%.3f", eg$pval),
        sprintf("%.3f", eg$est),
        sprintf("[%.3f, %.3f]", eg$ci.lb, eg$ci.ub),
        ifelse(eg$pval < 0.05, "Terindikasi Asimetri Bermakna (Potensi Bias Publikasi)", "Tidak Terindikasi Asimetri Bermakna (Simetris)")
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$table_begg_test <- renderTable({
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    res <- tryCatch(rma(yi = yi, vi = vi, data = df, method = input$meta_estimator), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    bg <- tryCatch(ranktest(res), error = function(e) NULL)
    if (is.null(bg)) return(data.frame(Status = "Uji Begg tidak dapat dikomputasi."))
    
    data.frame(
      Metrik = c("Statistik Kendall's Tau", "Nilai p (p-value)", "Kesimpulan"),
      Nilai = c(
        sprintf("%.3f", bg$tau),
        sprintf("%.3f", bg$pval),
        ifelse(bg$pval < 0.05, "Terindikasi Korelasi Asimetri Signifikan", "Korelasi Peringkat Tidak Signifikan")
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$table_trimfill_summary <- renderTable({
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    res <- tryCatch(rma(yi = yi, vi = vi, data = df, method = input$meta_estimator), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    tf <- tryCatch(trimfill(res), error = function(e) NULL)
    if (is.null(tf)) return(data.frame(Status = "Analisis Trim-and-Fill tidak dapat dijalankan."))
    
    data.frame(
      Parameter = c(
        "Estimasi Studi Hilang (k0)",
        "Sisi Asimetri Imputasi",
        "Pooled Effect Asli (Sebelum Imputasi)",
        "Adjusted Pooled Effect (Setelah Imputasi)",
        "Pergeseran Besaran Efek"
      ),
      Nilai = c(
        as.character(tf$k0),
        tf$side,
        sprintf("%.3f [%.3f, %.3f]", res$b[1], res$ci.lb[1], res$ci.ub[1]),
        sprintf("%.3f [%.3f, %.3f]", tf$b[1], tf$ci.lb[1], tf$ci.ub[1]),
        sprintf("%.3f (Atenuasi Efek)", tf$b[1] - res$b[1])
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$table_failsafe_n <- renderTable({
    df <- filtered_active_data()
    req(nrow(df) >= 2)
    k_min_val <- ifelse(is.null(input$bias_k_min), 5, as.numeric(input$bias_k_min))
    target_orw <- ifelse(is.null(input$fsn_orwin_target), 0.20, as.numeric(input$fsn_orwin_target))
    
    # Compatible with both metafor >= 5.0 (x = ...) and older metafor (yi = ...)
    run_fsn <- function(type_name, target_val = NULL) {
      tryCatch({
        if (is.null(target_val)) {
          fsn(x = df$yi, vi = df$vi, type = type_name, alpha = 0.05)
        } else {
          fsn(x = df$yi, vi = df$vi, type = type_name, alpha = 0.05, target = target_val)
        }
      }, error = function(e) {
        tryCatch({
          if (is.null(target_val)) {
            fsn(yi = df$yi, vi = df$vi, type = type_name, alpha = 0.05)
          } else {
            fsn(yi = df$yi, vi = df$vi, type = type_name, alpha = 0.05, target = target_val)
          }
        }, error = function(e2) NULL)
      })
    }
    
    fs_ros <- run_fsn("Rosenthal")
    fs_orw <- run_fsn("Orwin", target_val = target_orw)
    fs_rsb <- run_fsn("Rosenberg")
    
    if (is.null(fs_ros)) return(data.frame(Status = "Fail-Safe N tidak dapat dikomputasi."))
    threshold_ros <- 5 * nrow(df) + 10
    
    data.frame(
      `Metode & Parameter (k.min = 5)` = c(
        sprintf("Status Komputasi Paksa (k = %d, k.min = %d)", nrow(df), k_min_val),
        "Rosenthal Fail-Safe N (Target alpha = 0.05)",
        "Signifikansi Observasi Rosenthal (Observed p)",
        "Ambang Ketahanan Rosenthal (5k + 10)",
        sprintf("Orwin Fail-Safe N (Target g = %.3f)", target_orw),
        "Rata-Rata Efek Observasi Orwin (Average ES)",
        "Rosenberg Fail-Safe N (Weighted, alpha = 0.05)",
        "Signifikansi Observasi Rosenberg (Observed p)",
        "Kesimpulan Ketahanan Bukti (File-Drawer)"
      ),
      `Nilai Estimasi` = c(
        ifelse(nrow(df) >= k_min_val, sprintf("Aktif & Terpenuhi (k = %d >= %d)", nrow(df), k_min_val), sprintf("Dipaksa Aktif (k = %d < %d)", nrow(df), k_min_val)),
        paste(fs_ros$fsnum, "Studi"),
        sprintf("p = %.3f (Target \\alpha = %.3f)", fs_ros$pval, fs_ros$alpha),
        paste(threshold_ros, "Studi"),
        if (!is.null(fs_orw)) paste(fs_orw$fsnum, "Studi") else "-",
        if (!is.null(fs_orw)) sprintf("g = %.3f (Target = %.3f)", fs_orw$meanes, fs_orw$target) else "-",
        if (!is.null(fs_rsb)) paste(fs_rsb$fsnum, "Studi") else "-",
        if (!is.null(fs_rsb)) sprintf("p = %.3f (Average g = %.3f)", fs_rsb$pval, fs_rsb$meanes) else "-",
        ifelse(
          fs_ros$fsnum > threshold_ros,
          "Kuat / Tahan Bias File-Drawer (N_fs > 5k + 10)",
          ifelse(fs_ros$fsnum == 0, "N_fs = 0 (Efek gabungan belum melampaui ambang signifikansi target)", "Rentan terhadap Studi Tidak Terpublikasi (N_fs <= 5k + 10)")
        )
      ),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
}

# Run Application
shinyApp(ui = ui, server = server)
