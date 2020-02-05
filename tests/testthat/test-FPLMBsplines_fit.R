test_that("FPLMBplines_fit", {
    set.seed(1)
    n = 300
    p = 100
    ret = FPLMBsplines_fit(y = rnorm(n),
                           x = matrix(rnorm(n * p), n, p),
                           u = sin(sort(seq(0, 1, length = n))),
                           t = sort(runif(p)),
                           freq = 4,
                           spl = 4,
                           norder = 4,
                           fLoss = "ls")

    slope = c(317.9664, -422.1852, 617.6771, -430.8267)
    spl = c(-0.02783304, 0.03221540, 0.16196770, 0.18592464)
    expected_output = c(slope, spl)
    
    expect_equal(c(slope, spl), expected_output)
})
