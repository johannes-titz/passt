weights_sedlmeier14_2 <- matrix(c(.36, .03, .08, .29,
                                  .27, .53, .29, .15), nrow = 2)

sedlmeier14_2 <- updt_winner_weights(c(1, 0, 0, 1),
                                     weights_sedlmeier14_2,
                                     lrate = .1)

test_that("updating of winner weights works via Sedlmeier1999 fig.
          14.2 (p. 169)", {
            expect_equal(sedlmeier14_2,
                         matrix(c(.374, .072, .243, .311,
                                  .03, .29, .53, .15), byrow = T,
                                nrow = 2))
          })

set.seed(20191016)
test_that("updt_winner_weights work if there are equal output
          activations", {
            expect_equal(updt_winner_weights(c(1, 0),
                                             matrix(rep(0.5, 4),
                                                    nrow = 2),
                                             0.5),
                         matrix(c(0.5, 0.75, 0.50, 0.25), nrow = 2))
          })

test_that("calculating output sum for Sedlmeier fig 14.2 works
          (p. 169)",{
  expect_equal(calc_output_sum(c(1, 0, 0, 1),
                               weights_sedlmeier14_2), .83)
})

test_that("get_attention works", {
  expect_equal(get_attention(10, 1, 2, .5)[[1]], c(1, rep(.5, 9)))
  expect_equal(get_attention(10, 1, 10, .5)[[1]], c(rep(1, 9), .5))
  expect_equal(get_attention(c(10, 5), 1, 6, .2),
               list(c(rep(1, 5), rep(.2, 5)), c(rep(1, 5))))
})


presentation_matrix <- create_pres_matrix(diag(10), 1:10, 10:1,
                                          1, 2, 0.5,
                                          pulses_per_second = 1)[[1]]

test_that("create_pres_matrix works for orthogonal stimuli", {
  expect_equal(colSums(presentation_matrix), c(1:10 * 10:1))
})
