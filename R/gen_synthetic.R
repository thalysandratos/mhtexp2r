################################################################################
# gen_synthetic.R
#
# Generates a synthetic dataset that mirrors the structure of the Karlan-List
# (2007) data used in mhtexp2.do, for step-by-step debugging of R vs Stata.
#
# Design choices:
#   - n = 5000 rows before missingness
#   - ~3% of rows get at least one covariate set to NA -> dropped -> ~4850
#   - All variables match the types and roles in mhtexp2.do
#   - Seed is fixed so the dataset is fully reproducible
#   - Exported to CSV for import into Stata
################################################################################

set.seed(20260219)

n_raw <- 5000

# ---------------------------------------------------------------------------
# 1. Treatment variables
# ---------------------------------------------------------------------------

# Binary treatment (0/1), roughly 50/50
treatment <- sample(0:1, n_raw, replace = TRUE, prob = c(0.5, 0.5))

# Match ratio: 0 = control, 1 = 1:1, 2 = 2:1, 3 = 3:1
# Roughly 1/4 each (like Karlan-List: control + 3 match arms)
ratio <- sample(0:3, n_raw, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))

# ---------------------------------------------------------------------------
# 2. Geographic indicators (used to construct groupid)
# ---------------------------------------------------------------------------

# redcty: county voted Republican in last election (0/1)
redcty <- rbinom(n_raw, 1, 0.5)

# red0: state voted Republican in last election (0/1)
red0   <- rbinom(n_raw, 1, 0.5)

# groupid: exactly as in mhtexp2.do
#   1 = Red County, Red State
#   2 = Blue County, Red State
#   3 = Blue County, Blue State
#   4 = Red County, Blue State
#   0 (impossible by construction, but set to NA as in Stata)
groupid <- (redcty == 1 & red0 == 1) * 1 +
           (redcty == 0 & red0 == 1) * 2 +
           (redcty == 0 & red0 == 0) * 3 +
           (redcty == 1 & red0 == 0) * 4
groupid[groupid == 0] <- NA   # mirrors Stata's replace groupid = . if groupid == 0

# ---------------------------------------------------------------------------
# 3. Outcome variables
# ---------------------------------------------------------------------------

# gave: donated (binary). True effect of treatment ~ +8pp
p_gave <- 0.20 + 0.08 * treatment + 0.02 * (ratio > 0)
p_gave <- pmin(pmax(p_gave, 0), 1)
gave   <- rbinom(n_raw, 1, p_gave)

# amount: dollars given (continuous, 0 for non-donors)
# Among donors, lognormal with mild treatment effect
amount_raw <- rlnorm(n_raw, meanlog = 3.5 + 0.15 * treatment, sdlog = 1.2)
amount     <- amount_raw * gave   # zero for non-donors

# amountmat: amount * (1 + ratio) — exactly as in mhtexp2.do
amountmat <- amount * (1 + ratio)

# amountchange: change in donation amount (can be negative)
amountchange <- rnorm(n_raw, mean = 5 * treatment + 2 * ratio, sd = 40)

# ---------------------------------------------------------------------------
# 4. Baseline covariates (10 variables, matching mhtexp2.do controls)
# ---------------------------------------------------------------------------

# female: fraction female (0/1, slight imbalance)
female     <- rbinom(n_raw, 1, 0.55)

# pwhite: proportion white in zip code (0-1)
pwhite     <- rbeta(n_raw, 8, 2)   # right-skewed toward 1

# pblack: proportion black in zip code (0-1)
pblack     <- rbeta(n_raw, 2, 8)   # right-skewed toward 0

# page18_39: proportion aged 18-39 in zip code
page18_39  <- rbeta(n_raw, 4, 6)

# ave_hh_sz: average household size (roughly 2-3.5)
ave_hh_sz  <- rnorm(n_raw, mean = 2.6, sd = 0.4)

# years: years on mailing list (1-30)
years      <- round(runif(n_raw, min = 1, max = 30), 1)

# couple: household is couple (0/1)
couple     <- rbinom(n_raw, 1, 0.35)

# dormant: lapsed donor (0/1)
dormant    <- rbinom(n_raw, 1, 0.30)

# nonlit: non-litigation donor (0/1)
nonlit     <- rbinom(n_raw, 1, 0.60)

# cases: number of cases supported (count, 0-20)
cases      <- rpois(n_raw, lambda = 4)

# ---------------------------------------------------------------------------
# 5. Assemble raw data frame
# ---------------------------------------------------------------------------

df_raw <- data.frame(
  treatment    = treatment,
  ratio        = ratio,
  redcty       = redcty,
  red0         = red0,
  groupid      = groupid,
  gave         = gave,
  amount       = amount,
  amountmat    = amountmat,
  amountchange = amountchange,
  female       = female,
  pwhite       = pwhite,
  pblack       = pblack,
  page18_39    = page18_39,
  ave_hh_sz    = ave_hh_sz,
  years        = years,
  couple       = couple,
  dormant      = dormant,
  nonlit       = nonlit,
  cases        = cases
)

# ---------------------------------------------------------------------------
# 6. Introduce ~3% missingness on covariate/groupid columns
#    (mirrors Stata's "drop if v == ." loop)
#
#    We scatter NAs across the 11 columns (10 covariates + groupid).
#    Each row independently gets a NA on one of these columns with prob 3%.
#    A row may get multiple NAs but that still just results in one drop.
# ---------------------------------------------------------------------------

vars_with_missing <- c("female", "pwhite", "pblack", "page18_39",
                       "ave_hh_sz", "years", "couple", "dormant",
                       "nonlit", "cases", "groupid")

set.seed(20260219 + 1)   # separate seed for missingness pattern

for (v in vars_with_missing) {
  # Each variable independently gets ~0.5% NA -> together ~3-5% of rows
  # have at least one NA
  miss_idx <- which(runif(n_raw) < 0.005)
  df_raw[miss_idx, v] <- NA
}

cat("Raw dataset: ", nrow(df_raw), " rows\n")
cat("Missing counts per variable:\n")
print(colSums(is.na(df_raw[, vars_with_missing])))

# ---------------------------------------------------------------------------
# 7. Drop rows with any NA in the relevant columns
#    (mirrors Stata's foreach loop over allvars)
# ---------------------------------------------------------------------------

complete_mask <- complete.cases(df_raw[, vars_with_missing])
df <- df_raw[complete_mask, ]

# Reset row index
rownames(df) <- NULL

cat("\nAfter dropping incomplete rows: ", nrow(df), " rows\n")
cat("Fraction dropped: ", round(1 - nrow(df) / n_raw, 4), "\n\n")

# ---------------------------------------------------------------------------
# 8. Quick sanity checks
# ---------------------------------------------------------------------------

cat("=== Treatment proportions ===\n")
print(prop.table(table(df$treatment)))

cat("\n=== Ratio proportions ===\n")
print(prop.table(table(df$ratio)))

cat("\n=== GroupID proportions ===\n")
print(prop.table(table(df$groupid)))

cat("\n=== Gave rate by treatment ===\n")
print(tapply(df$gave, df$treatment, mean))

cat("\n=== Amount mean by treatment (donors only) ===\n")
print(tapply(df$amount[df$gave == 1], df$treatment[df$gave == 1], mean))

cat("\n=== Variable summary ===\n")
print(summary(df))

# ---------------------------------------------------------------------------
# 9. Export to CSV for Stata
# ---------------------------------------------------------------------------

write.csv(df, "synthetic_data.csv", row.names = FALSE, na = "")

cat("\nExported to synthetic_data.csv\n")
cat("Rows:", nrow(df), "\n")
cat("Columns:", ncol(df), "\n")
cat("Column names:", paste(names(df), collapse = ", "), "\n")
