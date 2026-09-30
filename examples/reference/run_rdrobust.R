#==============================================================================
# run_rdrobust.R -- reproduce the _test_vs_wls.do comparison against the
# R package rdrobust.
#
# The estimator of didc is a local polynomial regression of the differenced
# outcome on the (centered) running variable.  That is exactly what
# rdrobust() does in R when it is handed dy and z, so the two numbers must
# agree to machine precision once the bandwidth is held fixed.
#
# The in-repository test (examples/_test_vs_wls.do) performs the same check
# without R, using a hand-rolled weighted least squares fit.  Use this script
# when R is available, for a check that involves no Stata code at all.
#
# Usage
#   1. In Stata, from the examples directory:
#          use "../data/didc_sim1.dta", clear
#          import ... (see below)
#      or simply run the export block in examples/_test_vs_rdrobust.do
#   2. Rscript run_rdrobust.R
#
# Expected input : didc_sim1.csv with columns id, t, z, y
# Output         : rdrobust_reference.csv with the estimates and bandwidths
#==============================================================================

if (!requireNamespace("rdrobust", quietly = TRUE)) {
  install.packages("rdrobust")
}
library(rdrobust)

d <- read.csv("didc_sim1.csv")

# form Delta Y within unit, independently of anything Stata did
wide <- reshape(d, idvar = "id", timevar = "t", direction = "wide")
# reshape() names the columns y.0 / y.1 by the values of t
dy <- wide[["y.1"]] - wide[["y.0"]]
z  <- wide[["z.1"]]

fit <- rdrobust(y = dy, x = z, c = 0, p = 1, q = 2)

ref <- data.frame(
  tau_cl  = fit$Estimate[1, 1],
  tau_bc  = fit$Estimate[1, 2],
  se_rb   = fit$se[3],
  ci_l_rb = fit$ci[3, 1],
  ci_r_rb = fit$ci[3, 2],
  h_l     = fit$bws["h", "left"],
  h_r     = fit$bws["h", "right"],
  b_l     = fit$bws["b", "left"],
  b_r     = fit$bws["b", "right"]
)

write.csv(ref, "rdrobust_reference.csv", row.names = FALSE)
print(ref)

cat("\nCompare tau_cl, se_rb and the robust interval with didc's e(tau_didc),\n")
cat("e(se_didc) and e(ci_didc_*).  With the same bandwidth they must agree.\n")
