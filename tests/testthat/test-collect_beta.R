test_that("collect beta works", {
    fpath <- system.file("extdata", "PM155_P-run_pilotCells_unif.RData", package="NEMOcode")
    exfile = local(get(load(fpath)))
    pdes = NEMOcode::panel_design

    colldes = collect_beta_panel(cov_file = exfile, panel_design = pdes)
    expect_true(typeof(colldes) == "list")

    statsdes = collapse_stats(colldes, ssid = "PM155_P")
    expect_true(typeof(statsdes) == "list")

})


