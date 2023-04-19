test_that("compute_ci_confidence works", {

    ref = NEMOcode::refDist
    obs = ref$adm_mu
    control = ref$adm_mu

    tumor_mat = cbind(obs, obs)
    control_mat = cbind(control, control)

    rownames(tumor_mat) <- rownames(control_mat) <- ref$reg_id

    state = data.frame(reg_id = rownames(tumor_mat), meth_state = sample(x = c(-1,1), size = nrow(control_mat), replace = T))

    res = compute_ci_confidence(tumor_mat, control_mat, reg_df = state, nsub = 10, quant_prob = 0.05, frac_sub = 0.5)

    expect_true(typeof(res) == "list")
    expect_false(any(is.na(res)))

})
