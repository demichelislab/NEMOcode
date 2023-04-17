#' Compute the fraction of non-WBC derived cell-free DNA in a sample
#'
#' @param tumor_mat Matrix of tumor samples in which estimate IC fraction, regions on rows, samples on columns
#' @param control_mat Matrix of PBMC profiles, regions on rows, samples on columns
#' @param reg_df Dataframe of informative regions with 2 columns reg_id and meth state (-1/+1) for hypo and hyper respectively
#'
#' @import tidyr dplyr
#' @return A vector with TC estimation for input samples
#' @export
#'
compute_ic_meth <- function(tumor_mat, control_mat, reg_df) {

    rownames(reg_df) <- reg_df$reg_id

    ## Checks
    assertthat::assert_that((max(control_mat, na.rm = T) - min(control_mat, na.rm = T) <= 1))
    assertthat::assert_that((max(tumor_mat, na.rm = T) - min(tumor_mat, na.rm = T) <= 1))
    assertthat::assert_that(all(reg_df$reg_id %in% rownames(tumor_mat)))
    assertthat::assert_that(all(reg_df$reg_id == rownames(control_mat)))

    ## Based on the selection strategy,
    ## We assume that their discrete methylation state is opposite in PBMCs and CRPCs # nolint

    ## Separate Hyper and Hypo
    reg_hyper <- reg_df %>% filter(meth_state == +1) %>% pull(reg_id)
    reg_hypo <- reg_df %>% filter(meth_state == -1) %>% pull(reg_id)


    ## Extract rows and homogenize beta vals - (1-Beta) for hypo sites
    cancer_mat = tumor_mat
    cancer_mat[reg_hypo, ] = 1 - cancer_mat[reg_hypo, ]

    ctrl_mat = control_mat
    ctrl_mat[reg_hypo, ] = 1 - ctrl_mat[reg_hypo, ]

    tc_est_df = cancer_mat - apply(ctrl_mat, 1, median, na.rm = T)


    ## Offsets:
    ## using PBMCs/cfDNA we compute the expected lower and upper values (technical offsets).

    if (length(reg_hyper) >= 1) {
        offset_hyper <-
            median(apply(ctrl_mat[reg_hyper,], 1, mean, na.rm = T))
    } else {
        offset_hyper = 0
    }

    if (length(reg_hypo) >= 1) {
        offset_hypo <-
            median(apply(ctrl_mat[reg_hypo,], 1, mean, na.rm = T))
    } else {
        offset_hypo = 0
    }

    tc_est = apply(tc_est_df, 2, median, na.rm = T)

    ## Stretch the result based on offsets, computed on PBMCs
    tc_est <- tc_est/ ((1 - offset_hypo) - offset_hyper)

    ## If some values exceed the boundaries, correct
    ## This can happen due to noise and as the correction is
    ## based on the average value

    tc_est[which(tc_est < 0)] <- 0
    tc_est[which(tc_est > 1)] <- 1

    tc_res <- list(tc_est = tc_est, local_est = tc_est_df)

    return(tc_res)

}





#' Compute the fraction of unexpected epithelial cell-free DNA in a sample with subsampling confidence interval
#'
#' @param tumor_mat Matrix of tumor samples in which estimate IC fraction, regions on rows, samples on columns
#' @param control_mat Matrix of PBMC profiles, regions on rows, samples on columns
#' @param reg_df Dataframe of informative regions with 2 columns reg_id and meth state (-1/+1) for hypo and hyper respectively
#' @param nsub number of subsampling iterations to estimate stability (default = 100)
#' @param quant_prob the quantile of probabilities to produce stability interval via subsampling (0.05 = 95% CI)
#' @param frac_sub the fraction of informative regions to use for subsampling iterations
#'
#' @import dplyr tidyr
#' @return A dataframe with tumor content estimation
#' @export
compute_ci_confidence = function(tumor_mat,
                                 control_mat,
                                 reg_df,
                                 nsub = 100,
                                 quant_prob = .05,
                                 frac_sub = .5) {
    requireNamespace(tidyverse)

    ll_bt = list()

    for (i in seq(1:nsub)) {

        message(i)

        keep_ = sample(1:nrow(tumor_mat),
                       size = round(nrow(reg_df) * frac_sub),
                       replace = F)

        tumor_bt = tumor_mat[keep_,]
        control_bt = control_mat[keep_,]
        reg_bt = reg_df[keep_,]

        tc_bt = compute_ic_meth(tumor_bt, control_bt, reg_bt)

        tc_bt = data.frame(
            SampleName = names(tc_bt$tc_est),
            meth_est = tc_bt$tc_est,
            bt = i)

        rownames(tc_bt) = NULL
        ll_bt[[paste0("bt_", i)]] = tc_bt

    }

    df_bt = bind_rows(ll_bt)

    tc_est_ci = df_bt %>%
        group_by(SampleName) %>%
        summarise(
            est_mu = mean(meth_est, na.rm = T),
            est_sd = sd(meth_est, na.rm = T),
            est_min = min(meth_est, na.rm = T),
            est_max = max(meth_est, na.rm = T),
            q025_tc = quantile(
                meth_est,
                probs = quant_prob,
                na.rm = T,
                names = F
            ),
            q975_tc = quantile(
                meth_est,
                probs = 1 - quant_prob,
                na.rm = T,
                names = F
            )
        ) %>%
        mutate(
            est_mu = .clamp_val(est_mu, 0, 1),
            ci_lower = .clamp_val(q025_tc, 0, 1),
            ci_upper = .clamp_val(q975_tc, 0, 1),
            est_min = .clamp_val(est_min, 0, 1),
            est_max = .clamp_val(est_max, 0, 1))

    return(tc_est_ci)

}







#' Decompose tumor component into NE and Adeno modeling observed signal as a linear combination of reference profiles,
#' using a strong prior based on informed tumor content estimation
#'
#' @param obs Observed DNA methylation (Beta) for each region, Vector
#' @param cf_mu Mean DNA methylation for each region in cfDNA
#' @param cf_sd SD DNA methylation for each region in cfDNA
#' @param ad_mu Mean DNA methylation for each region in CRPC-Adeno
#' @param ad_sd SD DNA methylation for each region in CRPC-Adeno
#' @param ne_mu Mean DNA methylation for each region in CRPC-NE
#' @param ne_sd SD DNA methylation for each region in CRPC-NE
#' @param tc Estimation of tumor content for the current sample
#' @param tc_lw lower condifence of TC (not currently used)
#' @param tc_up upper confidence of TC (not currently used)
#'
#' @import brms tidyr dplyr
#' @return The estimated components for the observed sample
#' @export
compute_evidence_brms = function(obs,
                                 cf_mu,
                                 cf_sd,
                                 ad_mu,
                                 ad_sd,
                                 ne_mu,
                                 ne_sd,
                                 tc,
                                 tc_lw,
                                 tc_up) {
    if (is.na(tc) | is.nan(tc) | tc == 0) {

        res = data.frame(
            immune = 1 - tc,
            adeno = NA,
            ne = NA,
            rel_error = NA,
            var_score = NA
        )

    } else {

        refm = as.matrix(data.frame(cf_mu, ad_mu, ne_mu))
        rownames(refm) = paste0("idx", seq(1, nrow(refm)))

        wgh = 1 / (ad_sd + ne_sd)
        wgh = wgh*(1/max(wgh))

        dfin = data.frame(obs, refm, wgh)

        suppressPackageStartupMessages(base::library(brms))

        priors <- c(
            prior(normal(1 - tc, 0.01), coef = "cf_mu"), # Strongly informative, as TC is known
            prior(normal(tc / 2, .5), coef = "ad_mu"), # Uninformative
            prior(normal(tc / 2, .5), coef = "ne_mu") # Uninformative
        )

        idata = na.omit(dfin)

        if (TRUE) {

            ## Linear or quadratic weighting based on variance in the reference
            brm.1 <- brm(
                obs | weights(wgh) ~ cf_mu + ad_mu + ne_mu - 1,
                family = brmsfamily("gaussian"),
                data = idata,
                chains = 1,
                cores = getOption("mc.cores", 1),
                iter = 2000,
                warmup = 500,
                thin = 5,
                prior = priors,
                stanvars = stanvar(tc, name = "tc"),
                silent = 2
            )
        } else {
            brm.1 <- brm(
                obs ~ cf_mu + ad_mu + ne_mu - 1,
                family = brmsfamily("gaussian"),
                data = idata,
                chains = 1,
                cores = getOption("mc.cores", 1),
                iter = 2000,
                warmup = 500,
                thin = 5,
                prior = priors,
                stanvars = stanvar(tc, name = "tc"),
                silent = 2
            )


        }

        outp = fixef(brm.1, robust = T)[, 1]

        ## Negative weights should go to zero, weights of more than 1 are normalized to sum to 1

        p_ne = .clamp_val(outp["ne_mu"], 0, 2)
        p_adeno = .clamp_val(outp["ad_mu"], 0 , 2)
        cc = .clamp_val(outp["cf_mu"], 0, 2)

        nrm_out = norm_one(c(p_ne, p_adeno, cc))

        p_ne = nrm_out[1]
        p_adeno = nrm_out[2]
        cc = nrm_out[3]

        ## The relative error is computed as the maximum deviation across the three parameters,
        ## divided by the total tumor content

        var_score =  2*max(fixef(brm.1, robust = T)[, 2])
        # rel_error =  var_score / (p_ne + p_adeno)

        ## Compute the reconstruction loss as a goodness of fit metric
        pred = (idata[,4]*p_ne)+(idata[,3]*p_adeno)+(idata[,2]*cc)
        rmse_recons = Metrics::rmse(idata[,1], pred)
        rel_error = rmse_recons

        res = data.frame(
            immune = cc,
            adeno = p_adeno,
            ne = p_ne,
            rel_error,
            var_score
        )

    }

    return(res)

}




#' Compute PE score wrapper
#'
#' @param dd A matrix of DNA methylation with informative regions on rows and samples on columns
#' @param atlas_tc A data.frame with previously computed tumor content estimation for each sample
#' @param ref_mat A matrix with mean and sd DNA methylation for each cell-type of interest and each region
#' @param atlas_samples List of samples for which the TC should be considered 1 regardless of the estimation
#' @param cfdna_ids List of samples for which healthy cfDNA should be used as background, instead of PBMCs
#' @param mode Deconvolution strategy used, default to brms Bayesian regression (Stan backend)
#' @param sequential execute the computation sequentially (for debugging purpose)
#' @param nclust Number of cluster for parallel computation
#'
#' @import parallel doParallel
#' @return A data.frame with PE score
#' @export
compute_all = function(dd, atlas_tc, ref_mat, atlas_samples, cfdna_ids, mode = "brms", sequential = FALSE, nclust = 4) {


    if (!sequential) {
        cl <- makeCluster(nclust)
        registerDoParallel(cl)
    } else {
        registerDoSEQ()
    }

    est_all = foreach(
        i = 1:ncol(dd),
        .combine = rbind,
        .packages = c("dplyr"), .export = c("deconv_fun"),
    ) %dopar% {

        # source("/shares/CIBIO-Storage/CO/SPICE/gianmarco/meth_collection/script/F_estimates.R")
        `%!in%` <- Negate(`%in%`)

        nn = colnames(dd)[i]
        message(nn)

        ## If the sample is in Atlas, consider 100% purity
        if (nn %in% atlas_samples) {
            tcs = 1
            tc_up = 1
            tc_lw = 1
        } else {
            tcs = atlas_tc %>% filter(SampleName == nn) %>% pull(est_mu)
            tc_up = atlas_tc %>% filter(SampleName == nn) %>% pull(ci_upper)
            tc_lw = atlas_tc %>% filter(SampleName == nn) %>% pull(ci_lower)
        }


        if (mode == "brms"){
            deconv_fun = compute_evidence_brms
        } else {
            message("Current deconvolution options include only bayesian regression (brms). Future options might be implemented.")
            stop()
        }


        if (!length(tcs) == 0) {

            if (nn %!in% cfdna_ids) {
                tern = deconv_fun(
                    obs = dd[, i],
                    cf_mu = ref_mat$adm_mu,
                    cf_sd = ref_mat$adm_sd,
                    ad_mu = ref_mat$adeno_mu,
                    ad_sd = ref_mat$adeno_sd,
                    ne_mu = ref_mat$ne_mu,
                    ne_sd = ref_mat$ne_sd,
                    tc = tcs,
                    tc_lw = tc_lw,
                    tc_up = tc_up
                )

                tern$SampleName = nn
                return(tern)


            } else {

                tern = deconv_fun(
                    obs = dd[, i],
                    cf_mu = ref_mat$cf_mu,
                    cf_sd = ref_mat$cf_sd,
                    ad_mu = ref_mat$adeno_mu,
                    ad_sd = ref_mat$adeno_sd,
                    ne_mu = ref_mat$ne_mu,
                    ne_sd = ref_mat$ne_sd,
                    tc = tcs,
                    tc_lw = tc_lw,
                    tc_up = tc_up
                )

                tern$SampleName = nn
                return(tern)

            }

        } else{
            tern = data.frame(
                immune = NA,
                adeno = NA,
                ne = NA,
                rel_error = NA,
                var_score = NA,
                SampleName = nn
            )

            return(tern)
        }
    }

    stopCluster(cl)
    return(est_all)

}







#' Compute final scores from deconvolved fractions
#'
#' @param score_df input data.frame
#' @param limit lower limit of detection for tumor content
#'
#' @import dplyr tidyr
#' @return A data.frame with complete estimations
#' @export
#'
format_score = function(score_df, limit = 0.05){

    # Compute final evidence

    score_ss = score_df %>%
        mutate(
            pes = ne/(ne+adeno),
            pes_lw = ne/(ne+adeno),
            pes_up = ne/(ne+adeno),
            quality_flag = est_mu > limit) %>%
        select(
            SampleName,
            file_id,
            Dataset,
            Class,
            Type,
            tc_est = est_mu,
            tc_lw = ci_lower,
            tc_up = ci_upper,
            evidence_ne = ne,
            evidence_adeno = adeno,
            pes,
            pes_lw,
            pes_up,
            quality_flag,
            TumorContent,
            var_score,
            rel_error
        )

    ## If score is below the limit of detection, set to NA

    score_ss$pes = with(score_ss, ifelse(tc_est < limit, NA, pes))
    score_ss$pes_lw = with(score_ss, ifelse(tc_est < limit, NA, pes_lw))
    score_ss$pes_up = with(score_ss, ifelse(tc_est < limit, NA, pes_up))

    return(score_ss)

}
