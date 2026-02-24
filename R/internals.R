# -----------------------------------------

# ---- helpers to match Stata/Mata edge behavior ----
.safe_var_vec <- function(x) {
  x <- as.vector(x)
  if (length(x) <= 1 || all(is.na(x))) return(0)
  stats::var(x, na.rm = TRUE)
}

.safe_var_mat <- function(M) {
  # stats::var() on a 1-row matrix yields NA in R; Stata/Mata behavior is effectively 0
  if (is.null(dim(M)) || nrow(M) <= 1) {
    return(matrix(0, nrow = ncol(M), ncol = ncol(M)))
  }
  stats::var(M)
}

.solve_ols <- function(A, b) {
  # Prefer QR solve; fall back to generalized inverse if singular
  out <- tryCatch(as.vector(qr.solve(A, b)), error = function(e) NULL)
  if (!is.null(out) && !anyNA(out)) return(out)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Singular system in .run_reg and MASS not available for ginv().")
  }
  as.vector(MASS::ginv(A) %*% b)
}

# Build selection matrix
# select[k, i, j] = 1 if hypothesis (outcome=i, subgroup=j, comparison=k) is selected
.build_select <- function(only, exclude, numoc, numsub, numpc) {
  select <- array(0, dim = c(numpc, numoc, numsub))

  if (!all(is.na(only[1, ]))) {
    for (r in 1:nrow(only)) {
      select[only[r, 3], only[r, 1], only[r, 2]] <- 1
    }
  } else {
    select[] <- 1
  }

  if (!all(is.na(exclude[1, ]))) {
    for (r in 1:nrow(exclude)) {
      select[exclude[r, 3], exclude[r, 1], exclude[r, 2]] <- 0
    }
  }

  return(select)
}

# Main analysis engine implementing Algorithm 1
.seidelxu2 <- function(Y, sub, D, combo, select, bootstrap, DX, X,
                       studentized, idbootmat, transitivitycheck) {

  n <- nrow(Y)
  B <- bootstrap
  numoc <- ncol(Y)
  numsub <- length(unique(sub[!is.na(sub)]))
  numg <- length(unique(D[!is.na(D)])) - 1
  numpc <- nrow(combo)
  stud <- as.numeric(studentized)

  # --- Actual data statistics ---
  regact <- array(NA, dim = c(numoc, numsub, numpc))
  abregact <- array(NA, dim = c(numoc, numsub, numpc))

  for (i in 1:numoc) {
    for (j in 1:numsub) {
      if (ncol(DX) > 1) {
        sg <- (sub == j)
        cursgX <- X[sg, , drop = FALSE]
        barXz <- colMeans(cursgX)
        varXz <- .safe_var_mat(cursgX)
      } else {
        sg <- (sub == j)
        barXz <- 0
        varXz <- 0
      }

      for (l in 1:numpc) {
        w <- (sub == j & (D == combo[l, 1] | D == combo[l, 2]))
        pi_2 <- sum(sub == j & D == combo[l, 2]) / n
        pi_1 <- sum(sub == j & D == combo[l, 1]) / n
        pi_z <- sum(sub == j) / n

        if (sum(w) > 0) {
          curDX <- DX[w, , drop = FALSE]
          curD <- D[w, , drop = FALSE]
          curY <- Y[w, i, drop = FALSE]
          curDX[, 1] <- as.numeric(curD == combo[l, 2])

          regres <- .run_reg(curDX, curY, barXz, varXz, pi_2, pi_1, pi_z, sum(sg))
          regact[i, j, l] <- regres$coef
          abregact[i, j, l] <- abs(regres$coef) / (stud * regres$se + (1 - stud))
        }
      }
    }
  }

  # --- Bootstrap ---
  cat("Bootstrap iteration\n")
  abregboot <- array(0, dim = c(B, numoc, numsub, numpc))

  for (i in 1:B) {
    if (i %% 10 == 0) cat("-")
    if (i %% 1000 == 0) cat(sprintf("%g\n", i))

    Yboot <- Y[idbootmat[, i], , drop = FALSE]
    subboot <- sub[idbootmat[, i], , drop = FALSE]
    Dboot <- D[idbootmat[, i], , drop = FALSE]
    DXboot <- DX[idbootmat[, i], , drop = FALSE]

    if (ncol(DXboot) > 1) {
      Xboot <- X[idbootmat[, i], , drop = FALSE]
    }

    for (j in 1:numoc) {
      for (k in 1:numsub) {
        if (ncol(DXboot) > 1) {
          sg <- (subboot == k)
          cursgX <- Xboot[sg, , drop = FALSE]
          barXz <- colMeans(cursgX)
          varXz <- .safe_var_mat(cursgX)
        } else {
          sg <- (sub == k)
          barXz <- 0
          varXz <- 0
        }

        for (l in 1:numpc) {
          w <- (subboot == k & (Dboot == combo[l, 1] | Dboot == combo[l, 2]))
          pi_2 <- sum(sub == k & Dboot == combo[l, 2]) / n
          pi_1 <- sum(sub == k & Dboot == combo[l, 1]) / n
          pi_z <- sum(sub == k) / n

          if (sum(w) > 0) {
            curDX <- DXboot[w, , drop = FALSE]
            curD <- Dboot[w, , drop = FALSE]
            curY <- Yboot[w, j, drop = FALSE]
            curDX[, 1] <- as.numeric(curD == combo[l, 2])

            regres <- .run_reg(curDX, curY, barXz, varXz, pi_2, pi_1, pi_z, sum(sg))
            abregboot[i, j, k, l] <- abs(regres$coef - regact[j, k, l]) / (stud * regres$se + (1 - stud))
          }
        }
      }
    }
  }

  # --- P-values ---
  pact <- array(0, dim = c(numoc, numsub, numpc))
  pboot <- array(0, dim = c(B, numoc, numsub, numpc))

  for (i in 1:numoc) {
    for (j in 1:numsub) {
      for (k in 1:numpc) {
        pact[i, j, k] <- 1 - sum(abregboot[, i, j, k] >= abregact[i, j, k]) / B
        for (l in 1:B) {
          pboot[l, i, j, k] <- 1 - sum(abregboot[, i, j, k] >= abregboot[l, i, j, k]) / B
        }
      }
    }
  }

  # --- Single hypothesis testing (Remark 3.2) ---
  alphasin <- array(0, dim = c(numoc, numsub, numpc))

  for (i in 1:numoc) {
    for (j in 1:numsub) {
      for (k in 1:numpc) {
        sortp <- sort(pboot[, i, j, k], decreasing = TRUE)
        v <- (pact[i, j, k] >= sortp)
        indx <- which(v)[1]
        alphasin[i, j, k] <- if (is.na(indx)) 1 else indx / B
      }
    }
  }

  psin <- alphasin

  # --- Multiple hypothesis testing ---
  nh <- 0
  for (k in 1:numpc) {
    nh <- nh + sum(select[k, , ])
  }

  statsall <- matrix(0, nrow = nh, ncol = 8 + B)

  counter <- 1
  for (i in 1:numoc) {
    for (j in 1:numsub) {
      for (k in 1:numpc) {
        if (select[k, i, j] == 1) {
          statsall[counter, ] <- c(counter, i, j, combo[k, ],
                                   regact[i, j, k], psin[i, j, k],
                                   pact[i, j, k], pboot[, i, j, k])
          counter <- counter + 1
        }
      }
    }
  }

  statsrank <- statsall[order(statsall[, 7]), , drop = FALSE]

  alphamul <- rep(0, nh)
  alphamulm <- rep(0, nh)

  # --- Theorem 3.1 ---
  for (i in 1:nh) {
    if (i < nh) {
      maxstats <- apply(statsrank[i:nh, 9:(8 + B), drop = FALSE], 2, max)
    } else {
      maxstats <- statsrank[nh, 9:(8 + B)]
    }
    sortmaxstats <- sort(maxstats, decreasing = TRUE)
    indx <- which(statsrank[i, 8] >= sortmaxstats)[1]
    alphamul[i] <- if (is.na(indx)) 1 else indx / B

    if (i == 1 || !transitivitycheck) {
      alphamulm[i] <- alphamul[i]
    } else {
      # --- Remark 3.8: transitivity check ---
      sortmaxstatsm <- rep(0, B)

      for (j in (nh - i + 1):1) {
        remaining_indices <- statsrank[i:nh, 1]
        subsets <- utils::combn(remaining_indices, j)

        sumcont <- 0

        for (k_sub in 1:ncol(subsets)) {
          subset_indices <- subsets[, k_sub]
          cont <- 0

          for (l in 1:(i - 1)) {
            tempA <- statsall[subset_indices, 2:3, drop = FALSE]
            tempB <- matrix(rep(statsrank[l, 2:3], each = nrow(tempA)), ncol = 2)

            same_mask <- (tempA[, 1] == tempB[, 1]) & (tempA[, 2] == tempB[, 2])
            sameocsub <- subset_indices[same_mask]

            if (length(sameocsub) >= 1) {
              tran_list <- list()
              for (m in seq_along(sameocsub)) {
                tran_list[[m]] <- statsall[sameocsub[m], 4:5]
              }
              trantemp_list <- tran_list
            }

            if (length(sameocsub) <= 1) {
              cont <- 0
              maxstatsm <- apply(statsall[subset_indices, 9:(8 + B), drop = FALSE], 2, max)
              sortmaxstatsm <- pmax(sortmaxstatsm, sort(maxstatsm, decreasing = TRUE))
              break
            } else {
              # Transitivity merging
              counter_merge <- 1
              repeat {
                tran_list <- trantemp_list
                if (length(tran_list) <= 1) break
                trantemp_list <- list(tran_list[[1]])
                counter_merge <- counter_merge + 1
                merged_any <- FALSE

                for (m in 2:length(tran_list)) {
                  belong <- 0

                  for (nn in seq_along(trantemp_list)) {
                    trantempn <- trantemp_list[[nn]]
                    tranm <- tran_list[[m]]

                    unq <- sort(unique(c(trantempn, tranm)))
                    if (all(unq < (length(trantempn) + length(tranm)))) {
                      trantemp_list[[nn]] <- unq
                      belong <- belong + 1
                      merged_any <- TRUE
                    }
                  }

                  if (belong == 0) {
                    trantemp_list[[length(trantemp_list) + 1]] <- tranm
                  }
                }

                if (!merged_any || length(trantemp_list) >= length(tran_list)) {
                  break
                }
              }

              # Check if rejected hypothesis contradicts subset via transitivity
              for (p in seq_along(trantemp_list)) {
                rejected_pair <- statsrank[l, 4:5]
                if (sum(rejected_pair %in% trantemp_list[[p]]) == 2) {
                  cont <- 1
                  break
                }
              }
            }

            if (cont == 1) break
          }

          sumcont <- sumcont + cont

          if (cont == 0) {
            maxstatsm <- apply(statsall[subset_indices, 9:(8 + B), drop = FALSE], 2, max)
            sortmaxstatsm <- pmax(sortmaxstatsm, sort(maxstatsm, decreasing = TRUE))
          }
        }

        if (sumcont == 0) {
          break
        }
      }

      indx <- which(statsrank[i, 8] >= sortmaxstatsm)[1]
      alphamulm[i] <- if (is.na(indx)) 1 else indx / B
    }
  }

  # Bonferroni and Holm
  bon <- pmin(statsrank[, 7] * nh, 1)
  holm <- pmin(statsrank[, 7] * (nh:1), 1)
  for (i in 2:nh) holm[i] <- max(holm[i], holm[i-1])


  # Sort back to original order
  output <- cbind(statsrank[, 2:7, drop = FALSE], alphamul, alphamulm, bon, holm)
  output <- output[order(statsrank[, 1]), , drop = FALSE]

  colnames(output) <- c("outcome", "subgroup", "treatment1", "treatment2",
                         "coefficient", "Remark3_2", "Thm3_1", "Remark3_8",
                         "Bonf", "Holm")

  return(output)
}

# Regression for covariate-adjusted treatment effect estimation.
# Implements the estimator from Ye et al. (2022).
.run_reg <- function(X, y, Xbar, Xvar, pi_2, pi_1, pi_z, full_n) {
  y <- as.vector(y)

  if (is.matrix(X) && ncol(X) > 1) {
    D <- X[, 1]
  } else {
    D <- as.vector(X)
  }

  Ddummy <- (D == 1)
  Cdummy <- (D == 0)

  y1 <- y[Ddummy]
  y0 <- y[Cdummy]

  if (is.matrix(X) && ncol(X) > 1) {
    X_covs <- X[, -1, drop = FALSE]
    X1 <- X_covs[Ddummy, , drop = FALSE]
    X0 <- X_covs[Cdummy, , drop = FALSE]

    DX1 <- cbind(1, X1)
    DX0 <- cbind(1, X0)

    b1 <- tryCatch(
      as.vector(solve(crossprod(DX1)) %*% crossprod(DX1, y1)),
      error = function(e) rep(0, ncol(DX1))
    )
    b0 <- tryCatch(
      as.vector(solve(crossprod(DX0)) %*% crossprod(DX0, y0)),
      error = function(e) rep(0, ncol(DX0))
    )

    bX1 <- b1[-1]
    bX0 <- b0[-1]

    e1 <- y1 - X1 %*% bX1
    s1 <- .safe_var_vec(e1)

    e0 <- y0 - X0 %*% bX0
    s0 <- .safe_var_vec(e0)

    est_ATE <- (b1[1] - b0[1]) + sum(Xbar * (bX1 - bX0))
    est_VAR <- (1 / pi_2) * s1 + (1 / pi_1) * s0 +
      (1 / pi_z) * as.numeric(t(bX1 - bX0) %*% Xvar %*% (bX1 - bX0))

  } else {
    s1 <- if (length(y1) > 1) stats::var(y1) else 0
    s0 <- if (length(y0) > 1) stats::var(y0) else 0

    est_ATE <- mean(y1) - mean(y0)
    est_VAR <- (1 / pi_2) * s1 + (1 / pi_1) * s0
  }

  est_SE <- sqrt(abs(est_VAR))

  return(list(coef = est_ATE, se = est_SE))
}
