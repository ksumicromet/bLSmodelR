############### TEST SUBSET - bLS WEST SIDE ###################################

library(bLSmodelR)
library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

# ============================================================
# DIRECTORIES
# ============================================================

dirPath  <- "C:/Users/esantos/OneDrive - Kansas State University/AgroComb shared/SEK feedyard experiment/bls_R/Bls_R_test"
CatPath  <- file.path(dirPath, "Catalog_test")
outPath  <- file.path(dirPath, "resultados_bls_test")

dir.create(CatPath,  showWarnings = FALSE, recursive = TRUE)
dir.create(outPath,  showWarnings = FALSE, recursive = TRUE)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

ajustar_timestamp <- function(x) {
  x <- trimws(as.character(x))
  parse_date_time(x,
    orders = c("mdy HMS","mdy HM","ymd HMS","ymd HM","dmy HMS","dmy HM"),
    tz = "America/Chicago")
}

converter_numero <- function(x) as.numeric(gsub(",", ".", trimws(as.character(x))))

# ============================================================
# READ METEOROLOGICAL FILES
# ============================================================

met_1 <- fread(file.path(dirPath, "windtrax_westpens_sync_2024-03-14_1330.csv"),
               sep = ",", dec = ".", na.strings = c("", "NA", "NaN"))

met_2 <- fread(file.path(dirPath, "westpens_sync_2025-01-15_10-18.csv"),
               sep = ",", dec = ".", na.strings = c("", "NA", "NaN"))

cat("met_1 columns:", paste(names(met_1), collapse=", "), "\n")
cat("met_2 columns:", paste(names(met_2), collapse=", "), "\n")

# Align columns
if (!setequal(names(met_1), names(met_2))) {
  stop("Meteorological files have different columns.")
}
met_2 <- met_2[, names(met_1), with = FALSE]

# Parse timestamps
met_1 <- met_1 %>%
  mutate(TimeStamp = ajustar_timestamp(TimeStamp)) %>%
  filter(!is.na(TimeStamp)) %>% arrange(TimeStamp) %>% distinct(TimeStamp, .keep_all = TRUE)

met_2 <- met_2 %>%
  mutate(TimeStamp = ajustar_timestamp(TimeStamp)) %>%
  filter(!is.na(TimeStamp)) %>% arrange(TimeStamp) %>% distinct(TimeStamp, .keep_all = TRUE)

dados <- bind_rows(met_1, anti_join(met_2, met_1 %>% select(TimeStamp), by = "TimeStamp")) %>%
  arrange(TimeStamp) %>% distinct(TimeStamp, .keep_all = TRUE)

cat("Rows in complete met file:", nrow(dados), "\n")

# ============================================================
# READ SEK FILE
# ============================================================

SEK_raw <- fread(file.path(dirPath, "SEK_west_sync_all_2025-03-14.csv"),
                 sep = ",", dec = ".", na.strings = c("", "NA", "NaN"),
                 check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

SEK <- as.data.frame(SEK_raw[48:1185, ])

# Rename columns
rename_if <- function(df, old, new) { if (old %in% names(df)) names(df)[names(df)==old] <- new; df }
SEK <- SEK %>%
  rename_if("V2",  "Ustar") %>% rename_if("V3", "L") %>% rename_if("V4", "Zo") %>%
  rename_if("V5",  "WD")   %>% rename_if("V6", "Q_WindTrax_g_m2_day") %>%
  rename_if("V12", "CH4_LineSensor1_ppmv") %>% rename_if("V13", "CH4_LineSensor2_ppmv")

# Parse timestamp
ts_col <- intersect(c("WindTrax Simulation Results","TimeStamp"), names(SEK))[1]
if (is.na(ts_col)) stop("No timestamp column found in SEK file.")
SEK <- SEK %>% rename(TimeStamp_raw = !!ts_col) %>%
  mutate(
    TimeStamp_raw = trimws(as.character(TimeStamp_raw)),
    TimeStamp_raw = if_else(grepl("^Missing data", TimeStamp_raw), NA_character_, TimeStamp_raw),
    TimeStamp = ajustar_timestamp(TimeStamp_raw)
  ) %>%
  filter(!is.na(TimeStamp)) %>%
  select(-TimeStamp_raw)

# Convert numeric columns
SEK <- SEK %>% mutate(across(c(Ustar, L, Zo, WD, Q_WindTrax_g_m2_day,
                                CH4_LineSensor1_ppmv, CH4_LineSensor2_ppmv), converter_numero))

# ============================================================
# MERGE SEK + MET
# ============================================================

met_extra <- dados %>%
  select(TimeStamp, u_var, v_var, w_var, air_temperature, air_pressure) %>%
  mutate(across(c(u_var, v_var, w_var, air_temperature, air_pressure), converter_numero)) %>%
  filter(!is.na(TimeStamp)) %>% arrange(TimeStamp) %>% distinct(TimeStamp, .keep_all = TRUE)

dados_modelo <- SEK %>%
  select(TimeStamp, Ustar, L, Zo, WD, Q_WindTrax_g_m2_day,
         CH4_LineSensor1_ppmv, CH4_LineSensor2_ppmv) %>%
  inner_join(met_extra, by = "TimeStamp") %>%
  mutate(d = 0.067, z_sWu = 2.10, N0 = 3E4, MaxFetch = 500,
         sigma_u = sqrt(u_var), sigma_v = sqrt(v_var), sigma_w = sqrt(w_var),
         sUu = sigma_u / Ustar, sVu = sigma_v / Ustar, sWu = sigma_w / Ustar)

# ============================================================
# FILTER VALID ROWS
# ============================================================

dados_bls <- dados_modelo %>%
  filter(!is.na(TimeStamp), !is.na(Ustar), !is.na(L), !is.na(Zo), !is.na(WD),
         !is.na(u_var), !is.na(v_var), !is.na(w_var),
         !is.na(sUu), !is.na(sVu), !is.na(sWu),
         Ustar > 0, WD >= 0, WD <= 360,
         is.finite(L), is.finite(sUu), is.finite(sVu), is.finite(sWu))

cat("Valid rows for bLS:", nrow(dados_bls), "\n")

# ============================================================
# *** SUBSET: first 5 intervals only ***
# ============================================================

dados_bls <- head(dados_bls, 5)
cat("Subset rows used for test:", nrow(dados_bls), "\n")
print(dados_bls[, c("TimeStamp","Ustar","L","Zo","WD","sUu","sVu","sWu")])

# ============================================================
# SENSORS & SOURCES
# ============================================================

Sensors <- genSensors(
  LineSensor1 = list(x = c(56.554, 5.523), y = c(33.988, 46.169), z = c(1.50, 2.36), d = 1.81),
  LineSensor2 = list(x = c(54.936, 0.0),   y = c(4.784,  0.0),   z = c(1.33, 2.06), d = 1.90)
)

Sources <- genSources(
  Feedlot = list(Class = "p",
                 x = c(4.57, 35.15, 36.42, 6.12),
                 y = c(4.08, 4.08,  34.38, 34.67161695))
)

png(file.path(outPath, "test_sensors_sources.png"), width = 800, height = 600)
plot(Sources, Sensors, main = "Source and Sensors")
dev.off()
cat("Plot saved to:", file.path(outPath, "test_sensors_sources.png"), "\n")

SourceArea_m2 <- as.numeric(getArea(Sources)["Feedlot"])
cat("Source area (m2):", SourceArea_m2, "\n")

# ============================================================
# BUILD INTERVALS
# ============================================================

IntervalData <- dados_bls %>%
  mutate(SensorNames = "LineSensor1,LineSensor2", SourceNames = "Feedlot") %>%
  select(TimeStamp, Ustar, L, Zo, sUu, sVu, sWu, z_sWu, WD, d, N0, MaxFetch,
         SensorNames, SourceNames)

Intervals <- genInterval(Data = IntervalData)
cat("Intervals created:", length(Intervals), "\n")

# ============================================================
# MODEL CONFIG  (reduced N0 & cores for speed)
# ============================================================

Model <- genModel(
  TDwrite     = TRUE,
  overwriteTD = TRUE,
  TDread      = TRUE,
  TDonly      = FALSE,
  ncores      = 4       # reduced for test
)

Tolerances <- genTolerances()

InputList <- genInputList(Intervals, Sensors, Sources, Model, Tolerances)

# ============================================================
# RUN bLS
# ============================================================

cat("\nRunning bLS on", length(Intervals), "intervals...\n")
Results_bls <- runbLS(InputList, Cat.Path = CatPath, asDT = FALSE)
Results_bls <- switchNames(Results_bls)

cat("bLS result columns:", paste(names(Results_bls), collapse=", "), "\n")
print(head(Results_bls))

# ============================================================
# COMPUTE EMISSION
# ============================================================

Results_CE <- Results_bls %>%
  as.data.frame() %>%
  select(TimeStamp, Sensor, Source, CE, CE_se, CE_lo, CE_hi, everything()) %>%
  filter(!is.na(CE), CE > 0)

cat("Rows with valid CE:", nrow(Results_CE), "\n")

conc_df <- dados_modelo %>%
  filter(TimeStamp %in% dados_bls$TimeStamp,
         !is.na(CH4_LineSensor1_ppmv), !is.na(CH4_LineSensor2_ppmv)) %>%
  mutate(
    Cb_CH4_ppmv  = pmin(CH4_LineSensor1_ppmv, CH4_LineSensor2_ppmv, na.rm = TRUE),
    Sensor_downwind = case_when(
      CH4_LineSensor1_ppmv > CH4_LineSensor2_ppmv ~ "LineSensor1",
      CH4_LineSensor2_ppmv > CH4_LineSensor1_ppmv ~ "LineSensor2",
      TRUE ~ NA_character_)
  ) %>%
  pivot_longer(cols = c(CH4_LineSensor1_ppmv, CH4_LineSensor2_ppmv),
               names_to = "Sensor_col", values_to = "CH4_ppmv") %>%
  mutate(
    Sensor = case_when(Sensor_col == "CH4_LineSensor1_ppmv" ~ "LineSensor1",
                       Sensor_col == "CH4_LineSensor2_ppmv" ~ "LineSensor2"),
    CH4_delta_ppmv = CH4_ppmv - Cb_CH4_ppmv
  ) %>%
  filter(Sensor == Sensor_downwind) %>%
  select(TimeStamp, Sensor, CH4_ppmv, Cb_CH4_ppmv, CH4_delta_ppmv)

met_conv <- dados_modelo %>%
  filter(TimeStamp %in% dados_bls$TimeStamp) %>%
  select(TimeStamp, air_temperature, air_pressure) %>%
  mutate(T_K = converter_numero(air_temperature), P_Pa = converter_numero(air_pressure)) %>%
  select(TimeStamp, T_K, P_Pa)

MW_CH4 <- 16.04
R_gas  <- 8.314462618

Results_emission <- Results_CE %>%
  inner_join(conc_df, by = c("TimeStamp","Sensor")) %>%
  inner_join(met_conv, by = "TimeStamp") %>%
  mutate(
    CH4_delta_ppmv  = ifelse(CH4_delta_ppmv <= 0, NA, CH4_delta_ppmv),
    CH4_delta_ug_m3 = CH4_delta_ppmv * P_Pa * MW_CH4 / (R_gas * T_K),
    Q_ug_m2_s       = CH4_delta_ug_m3 / CE,
    Q_g_m2_day      = Q_ug_m2_s * 86400 / 1E6,
    SourceArea_m2   = SourceArea_m2,
    Q_total_kg_day  = Q_ug_m2_s * SourceArea_m2 * 86400 / 1E9
  )

cat("\n===== SUBSET TEST RESULTS =====\n")
print(Results_emission[, c("TimeStamp","Sensor","CE","CH4_delta_ppmv","Q_g_m2_day","Q_total_kg_day")])

cat("\nSummary Q_g_m2_day:\n")
print(summary(Results_emission$Q_g_m2_day))

cat("\nSummary Q_total_kg_day:\n")
print(summary(Results_emission$Q_total_kg_day))

# Save test results
fwrite(Results_emission %>% mutate(TimeStamp = format(TimeStamp, "%Y-%m-%d %H:%M:%S")),
       file.path(outPath, "test_subset_results.csv"), sep = ",", dec = ".")

cat("\nTest completed successfully. Results saved to:", outPath, "\n")
