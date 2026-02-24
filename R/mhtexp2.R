#' Multiple Hypothesis Testing with Covariate Adjustment
#'
#' Implements the stepwise multiple testing procedure described in
#' List, Shaikh, and Vayalinkal (2023). The procedure controls the familywise
#' error rate (FWER) while exploiting the joint dependence structure of test
#' statistics and observed baseline covariates for improved power.
#'
#' @param outcomes A data frame or matrix of outcome variables (n x K).
#'   Each column is a separate outcome.
#' @param treatment A vector or single-column matrix of treatment assignments.
#'   Should contain integer-coded group labels (0 = control, 1, 2, ... for
#'   treatment groups).
#' @param subgroup An optional vector of subgroup identifiers. If `NULL`
#'   (default), all units are treated as belonging to a single group.
#' @param combo A string specifying the type of treatment comparisons:
#'   `"treatmentcontrol"` (default) compares each treatment group to the
#'   control; `"pairwise"` performs all pairwise comparisons.
#' @param exclude An optional matrix with 3 columns (outcome, subgroup,
#'   comparison) specifying hypotheses to exclude. Default `NULL`.
#' @param only An optional matrix with 3 columns (outcome, subgroup,
#'   comparison) specifying which hypotheses to include. If provided, only
#'   these hypotheses are tested. Default `NULL`.
#' @param bootstrap Number of bootstrap iterations (default: 3000).
#' @param controls An optional data frame or matrix of baseline covariates
#'   for covariate-adjusted inference. Default `NULL`.
#' @param treatnames An optional character vector of names for the result rows.
#'   Default `NULL`.
#' @param studentized Logical; whether to use studentized test statistics
#'   (default: `TRUE`). See Remark 3.5 in the paper.
#' @param idbootmat An optional pre-specified bootstrap sampling matrix
#'   (n x B integer matrix of row indices). Useful for reproducibility or
#'   for sharing bootstrap draws across analyses. Default `NULL`.
#' @param transitivitycheck Logical; whether to apply the transitivity-based
#'   refinement from Remark 3.8 (default: `TRUE`). Only relevant for pairwise
#'   comparisons with 3+ groups.
#'
#' @return A numeric matrix with one row per tested hypothesis and columns:
#' \describe{
#'   \item{outcome}{Index of the outcome variable}
#'   \item{subgroup}{Index of the subgroup}
#'   \item{treatment1}{First treatment group in the comparison (typically control = 0)}
#'   \item{treatment2}{Second treatment group in the comparison}
#'   \item{coefficient}{Estimated treatment effect (covariate-adjusted if controls provided)}
#'   \item{Remark3_2}{Single-hypothesis bootstrap p-value (Remark 3.2)}
#'   \item{Thm3_1}{Stepdown-adjusted p-value (Theorem 3.1 / Algorithm 1)}
#'   \item{Remark3_8}{Transitivity-refined p-value (Remark 3.8; equals Thm3_1 if transitivitycheck = FALSE)}
#'   \item{Bonf}{Bonferroni-adjusted p-value}
#'   \item{Holm}{Holm-adjusted p-value}
#' }
#'
#' @details
#' The procedure tests null hypotheses of the form
#' \eqn{H_s: E[Y_k(d) - Y_k(d') | Z = z] = 0}, where \eqn{d} and \eqn{d'}
#' are treatment conditions, \eqn{k} indexes outcomes, and \eqn{z} indexes
#' subgroups.
#'
#' The five p-value columns represent, in order of increasing power:
#' Bonferroni (most conservative), Holm, Theorem 3.1 (stepdown using bootstrap
#' joint distribution), Remark 3.8 (transitivity refinement), and Remark 3.2
#' (single-hypothesis, no multiplicity adjustment).
#'
#' When covariates are supplied via `controls`, the treatment effect estimator
#' follows Ye et al. (2022), which guarantees asymptotic efficiency gains
#' over the simple difference in means.
#'
#' @references
#' List, J. A., Shaikh, A. M., & Vayalinkal, A. (2023). Multiple testing with
#' covariate adjustment in experimental economics. *Journal of Applied
#' Econometrics*, 38(6), 920--939. \doi{10.1002/jae.2985}
#'
#' Romano, J. P., & Wolf, M. (2010). Balanced control of generalized error rates.
#' *Annals of Statistics*, 38(1), 598--633.
#'
#' Ye, T., Shao, J., Yi, Y., & Zhao, Q. (2022). Toward better practice of
#' covariate adjustment in analyzing randomized clinical trials. *Journal of
#' the American Statistical Association*, 118(544), 2370--2382.
#'
#' @examples
#' # Simple two-group comparison with one outcome
#' set.seed(42)
#' n <- 200
#' D <- sample(0:1, n, replace = TRUE)
#' Y <- matrix(rnorm(n) + 1.5 * D, ncol = 1)
#' result <- mhtexp2(Y, D, bootstrap = 200)
#' print(result)
#'
#' # Multiple outcomes
#' Y2 <- cbind(rnorm(n) + 2 * D, rnorm(n))
#' result2 <- mhtexp2(Y2, D, bootstrap = 200)
#' print(result2)
#'
#' @export
mhtexp2 <- function(outcomes, treatment, subgroup = NULL,
                    combo = "treatmentcontrol", exclude = NULL, only = NULL,
                    bootstrap = 3000, controls = NULL, treatnames = NULL,
                    studentized = TRUE, idbootmat = NULL, transitivitycheck = TRUE) {

  if (!is.null(combo) && !combo %in% c("pairwise", "treatmentcontrol")) {
    stop("INVALID combo: choose either pairwise or treatmentcontrol")
  }

  Y <- as.matrix(outcomes)
  D <- as.matrix(treatment)

  if (!is.null(controls)) {
    X <- as.matrix(controls)
    DX <- cbind(D, X)
  } else {
    X <- 0
    DX <- D
  }

  if (is.null(exclude)) {
    excludemat <- matrix(c(NA, NA, NA), nrow = 1)
  } else {
    excludemat <- exclude
  }

  if (is.null(only)) {
    onlymat <- matrix(c(NA, NA, NA), nrow = 1)
  } else {
    onlymat <- only
  }

  if (is.null(subgroup)) {
    sub <- matrix(1, nrow = nrow(D), ncol = 1)
  } else {
    sub <- as.matrix(subgroup)
  }

  numoc <- ncol(Y)
  numsub <- length(unique(sub[!is.na(sub)]))
  numg <- length(unique(D[!is.na(D)])) - 1

  if (combo == "pairwise") {
    combo_matrix <- t(utils::combn(0:numg, 2))
  } else {
    combo_matrix <- cbind(rep(0, numg), 1:numg)
  }

  numpc <- nrow(combo_matrix)

  select <- .build_select(onlymat, excludemat, numoc, numsub, numpc)

  if (is.null(idbootmat)) {
    set.seed(0)
    idbootmat <- matrix(sample(1:nrow(Y), nrow(Y) * bootstrap, replace = TRUE),
                        nrow = nrow(Y), ncol = bootstrap)
  }

  results <- .seidelxu2(Y, sub, D, combo_matrix, select, bootstrap, DX, X,
                        studentized, idbootmat, transitivitycheck)

  if (!is.null(treatnames) && length(treatnames) == nrow(results)) {
    rownames(results) <- treatnames
  } else if (!is.null(treatnames) && length(treatnames) != nrow(results) && length(treatnames) > 0) {
    message("treatnames length does not match number of treatments!")
  }

  return(results)
}
