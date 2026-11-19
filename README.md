# Kawe: Quantitative and Psychometrics

[English](#english) | [Bahasa Indonesia](#bahasa-indonesia)

---

<a name="english"></a>
## English

A unified interactive platform and toolkit built with **R Shiny** for psychometric modeling, item response theory, and evidence synthesis. It integrates three core quantitative methodologies:

1. **Confirmatory Factor Analysis (CFA & Polychoric Correlations)**
2. **Rasch / Item Response Theory (TAM Engine & Instrument Purification)**
3. **Clinical Evidence Synthesis (Multilevel Meta-Analysis)**

The web application runs entirely in modern browsers via WebAssembly (**Shinylive / WebR**) without requiring local R installation.

---

### Web Application Access
- **Latest Release**: [https://adityakusumawijaya.github.io/quantitative-and-psychometrics/](https://adityakusumawijaya.github.io/quantitative-and-psychometrics/)
- **Release v0.11**: [https://adityakusumawijaya.github.io/quantitative-and-psychometrics/v0.11/](https://adityakusumawijaya.github.io/quantitative-and-psychometrics/v0.11/)

---

### Key Features

#### 1. Confirmatory Factor Analysis (CFA) & Polychoric Correlations
- **Polychoric Matrix Estimation**: Two-step Olsson estimation for ordinal Likert indicators.
- **Attenuation Analysis ($\Delta$)**: Quantification of artificial attenuation between Pearson and Polychoric correlations ($\Delta = r_{\text{poly}} - r_{\text{pearson}}$).
- **Robust Model Estimation (`lavaan`)**: WLSMV (DWLS with robust standard errors and mean-and-variance adjusted test statistics) and MLR.
- **Scale Reliability Decomposition**: McDonald's Omega ($\omega_t$, $\omega_h$) and Cronbach's Alpha.
- **Path Diagrams (`semPlot`)**: Interactive structural visualization with standardized factor loadings.

#### 2. Rasch & Item Response Theory (`TAM`)
- **Polytomous & Dichotomous Models**: Partial Credit Model (PCM) and Rating Scale Model (RSM) via Marginal Maximum Likelihood (MML / EM algorithm).
- **Wright Map Visualizations**: Clear item-person targeting maps for unidimensional and multidimensional constructs.
- **Scale Purification**: Automated outlier/misfit detection (Outfit MNSQ $> 1.50$ or $t > 2.0$) with re-estimation on purified samples to evaluate parameter stability.
- **Q-Matrix Configuration**: Interactive matrix assignment supporting multidimensional IRT structures.

#### 3. Clinical Evidence Synthesis & Meta-Analysis
- **Effect Size Calculator (`esc`)**: Converts 8 reporting formats (pre-post gains, means & SDs, $t$-tests, $p$-values, $F$-tests, $2 \times 2$ contingency tables, correlation coefficients) into standardized Hedges' $g$ and sampling variance $v_i$.
- **Random-Effects Modeling (`metafor`)**: Pooled effect estimation using Restricted Maximum Likelihood (REML), with comprehensive heterogeneity diagnostics ($Q$, $\tau^2$, $I^2$).
- **Publication Bias & Diagnostics**: Funnel plot, Egger's regression test, and Duval & Tweedie's *Trim and Fill* analysis.
- **Publication-Ready Figures**: Forest plots, Subgroup analyses, and Meta-regression bubble plots (`ggplot2`).

---

### Running Locally

#### Prerequisites
Make sure R (>= 4.2.0) is installed. Install the required R packages:

```r
install.packages(c(
  "shiny",
  "bslib",
  "lavaan",
  "TAM",
  "WrightMap",
  "metafor",
  "esc",
  "DT",
  "ggplot2",
  "misty",
  "psych",
  "semTools",
  "readxl",
  "semPlot",
  "Matrix",
  "MASS"
))
```

#### Launch the App
Clone the repository and run:

```r
# Launch from the repository root:
shiny::runApp(".", launch.browser = TRUE)

# Or run via script launcher:
source("01_Analysis/01_Run_Shiny_App.R")
```

---

<a name="bahasa-indonesia"></a>
## Bahasa Indonesia

Aplikasi interaktif dan perangkat analisis terpadu berbasis **R Shiny** untuk pemodelan psikometri, teori respons butir (*item response theory* / IRT), dan sintesis bukti penelitian. Aplikasi ini mengintegrasikan tiga metodologi kuantitatif utama:

1. **Analisis Faktor Konfirmatori (CFA & Korelasi Polikorik)**
2. **Pemodelan Rasch / IRT (TAM Engine & Pemurnian Skala)**
3. **Sintesis Bukti Klinis (Meta-Analisis Multilevel)**

Aplikasi web ini berjalan langsung di peramban modern melalui teknologi WebAssembly (**Shinylive / WebR**) tanpa memerlukan instalasi R lokal.

---

### Akses Aplikasi Web
- **Versi Terbaru**: [https://adityakusumawijaya.github.io/quantitative-and-psychometrics/](https://adityakusumawijaya.github.io/quantitative-and-psychometrics/)
- **Rilis v0.11**: [https://adityakusumawijaya.github.io/quantitative-and-psychometrics/v0.11/](https://adityakusumawijaya.github.io/quantitative-and-psychometrics/v0.11/)

---

### Fitur Utama

#### 1. Analisis Faktor Konfirmatori (CFA) & Korelasi Polikorik
- **Estimasi Matriks Polikorik**: Estimasi dua-tahap Olsson untuk indikator ordinal berskala Likert.
- **Analisis Atenuasi ($\Delta$)**: Kuantifikasi atenuasi artifisial antara korelasi Pearson dan Polikorik ($\Delta = r_{\text{poly}} - r_{\text{pearson}}$).
- **Estimasi Model Robust (`lavaan`)**: WLSMV (DWLS dengan standar error robust dan statistik uji yang disesuaikan rata-rata dan varians) serta MLR.
- **Dekomposisi Reliabilitas Pengukuran**: Koefisien Omega McDonald ($\omega_t$, $\omega_h$) dan Alpha Cronbach.
- **Diagram Jalur Interaktif (`semPlot`)**: Visualisasi struktural model pengukuran dengan muatan faktor terstandarisasi.

#### 2. Pemodelan Rasch & Teori Respons Butir (`TAM`)
- **Model Politomi & Dikotomi**: *Partial Credit Model* (PCM) dan *Rating Scale Model* (RSM) berbasis Marginal Maximum Likelihood (algoritma MML / EM).
- **Visualisasi Wright Map**: Pemetaan keselarasan tingkat kesulitan butir dan kemampuan responden (*person-item map*) untuk model unidimensi maupun multidimensi.
- **Pemurnian Skala**: Deteksi otomatis data aberan/misfit (Outfit MNSQ $> 1,50$ atau $t > 2,0$) serta re-estimasi pada sampel murni untuk menguji stabilitas parameter.
- **Konfigurasi Matriks Q**: Pengaturan interaktif matriks spesifikasi butir-faktor untuk pemodelan IRT multidimensi.

#### 3. Sintesis Bukti Klinis & Meta-Analisis
- **Kalkulator Ukuran Efek (`esc`)**: Mengonversi 8 format laporan studi primer (*pre-post gains*, rerata & SD, uji-$t$, nilai-$p$, uji-$F$, tabel kontingensi $2 \times 2$, koefisien korelasi) menjadi standar *Hedges' g* dan varians sampling $v_i$.
- **Model Efek Acak (`metafor`)**: Estimasi efek gabungan menggunakan Restricted Maximum Likelihood (REML) lengkap dengan diagnostik heterogenitas ($Q$, $\tau^2$, $I^2$).
- **Deteksi Bias Publikasi**: Funnel plot, uji regresi Egger, dan analisis *Trim and Fill* Duval & Tweedie.
- **Visualisasi Standar Publikasi**: Forest plot, analisis subkelompok, dan diagram sebar meta-regresi (*bubble plot* via `ggplot2`).

---

### Menjalankan Secara Lokal

#### Prasyarat
Pastikan R (>= 4.2.0) telah terpasang. Pasang paket-paket dependensi berikut:

```r
install.packages(c(
  "shiny",
  "bslib",
  "lavaan",
  "TAM",
  "WrightMap",
  "metafor",
  "esc",
  "DT",
  "ggplot2",
  "misty",
  "psych",
  "semTools",
  "readxl",
  "semPlot",
  "Matrix",
  "MASS"
))
```

#### Menjalankan Aplikasi
Klon repositori ini, lalu jalankan di konsol R atau RStudio:

```r
# Jalankan langsung dari direktori repositori:
shiny::runApp(".", launch.browser = TRUE)

# Atau melalui skrip launcher:
source("01_Analysis/01_Run_Shiny_App.R")
```

---

## References / Daftar Pustaka

### Methodological References / Rujukan Metodologis
- Borenstein, M., Higgins, J. P., Hedges, L. V., & Rothstein, H. R. (2017). Basics of meta-analysis: $I^2$ is not an absolute measure of heterogeneity. *Research Synthesis Methods*, 8(1), 5–18. https://doi.org/10.1002/jrsm.1230
- Cuijpers, P., Karyotaki, E., Eckshtain, D., Ng, M. Y., Cortes-Garcia, L., Lan, G., ... & Weisz, J. R. (2020). The effects of psychotherapies for depression on response, remission, reliable change, and deterioration: A meta-analysis. *Acta Psychiatrica Scandinavica*, 141(4), 312–322. https://doi.org/10.1111/acps.13155
- Duval, S., & Tweedie, R. (2000). Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. *Biometrics*, 56(2), 455–463. https://doi.org/10.1111/j.0006-341X.2000.00455.x
- Egger, M., Davey Smith, G., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. *BMJ*, 315(7109), 629–634. https://doi.org/10.1136/bmj.315.7109.629
- Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. *Journal of Educational Statistics*, 6(2), 107–128. https://doi.org/10.3102/10769986006002107

### R Packages / Sitasi Paket R
- **shiny**: Chang, W., Cheng, J., Allaire, J., Sievert, C., Schloerke, B., Aden-Buie, G., Xie, Y., Allen, J., McPherson, J., Dipert, A., & Borges, B. (2026). *shiny: Web Application Framework for R*. R package version 1.14.0. https://doi.org/10.32614/CRAN.package.shiny
- **bslib**: Sievert, C., Cheng, J., & Aden-Buie, G. (2026). *bslib: Custom 'Bootstrap' 'Sass' Themes for 'shiny' and 'rmarkdown'*. R package version 0.12.0. https://doi.org/10.32614/CRAN.package.bslib
- **lavaan**: Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. *Journal of Statistical Software*, 48(2), 1–36. https://doi.org/10.18637/jss.v048.i02
- **TAM**: Robitzsch, A., Kiefer, T., & Wu, M. (2025). *TAM: Test Analysis Modules*. R package version 4.3-25. https://doi.org/10.32614/CRAN.package.TAM
- **WrightMap**: Torres Irribarra, D., & Freund, R. (2026). *WrightMap: IRT Item-Person Map with 'ConQuest' Integration*. R package version 1.5. https://doi.org/10.32614/CRAN.package.WrightMap
- **metafor**: Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. *Journal of Statistical Software*, 36(3), 1–48. https://doi.org/10.18637/jss.v036.i03
- **esc**: Lüdecke, D. (2019). *esc: Effect Size Computation for Meta Analysis* (Version 0.5.1). https://doi.org/10.5281/zenodo.1249218
- **DT**: Xie, Y., Cheng, J., Tan, X., & Aden-Buie, G. (2025). *DT: A Wrapper of the JavaScript Library 'DataTables'*. R package version 0.34.0. https://doi.org/10.32614/CRAN.package.DT
- **ggplot2**: Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. https://ggplot2.tidyverse.org
- **misty**: Yanagida, T. (2024). *misty: Miscellaneous Functions 'T. Yanagida'*. R package version 0.6.6. https://CRAN.R-project.org/package=misty
- **psych**: Revelle, W. (2026). *psych: Procedures for Psychological, Psychometric, and Personality Research*. Northwestern University, Evanston, Illinois. R package version 2.6.9. https://CRAN.R-project.org/package=psych
- **semTools**: Jorgensen, T. D., Pornprasertmanit, S., Schoemann, A. M., & Rosseel, Y. (2026). *semTools: Useful tools for structural equation modeling*. R package version 0.5-10. https://CRAN.R-project.org/package=semTools
- **readxl**: Wickham, H., & Bryan, J. (2026). *readxl: Read Excel Files*. R package version 1.5.0.1. https://doi.org/10.32614/CRAN.package.readxl
- **semPlot**: Epskamp, S. (2026). *semPlot: Path Diagrams and Visual Analysis of Various SEM Packages' Output*. R package version 1.2.0. https://doi.org/10.32614/CRAN.package.semPlot
- **Matrix**: Bates, D., Maechler, M., & Jagan, M. (2026). *Matrix: Sparse and Dense Matrix Classes and Methods*. R package version 1.7-6. https://doi.org/10.32614/CRAN.package.Matrix
- **MASS**: Venables, W. N., & Ripley, B. D. (2002). *Modern Applied Statistics with S* (4th ed.). Springer, New York. https://doi.org/10.1007/978-0-387-21706-2
