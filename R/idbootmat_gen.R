n <- 4727
set.seed(42)
B <- 300
idbootmat <- matrix(sample(1:n, n * B, replace = TRUE), nrow = n, ncol = B)
#write.csv(idbootmat, "idbootmat.csv", row.names = FALSE)
#cat(dim(idbootmat))  # should print 4727 300
