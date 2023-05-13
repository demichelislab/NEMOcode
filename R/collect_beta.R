#' Collect Beta Panel
#'
#' @param cov_file Uniformed output file
#' @param panel_design Data frame with CONCATENATED design and region_set and reg_id columns
#' @param sets_analyze The panel submodules that should be processed
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @import dplyr tidyr
#' @return A data.frame with collapsed beta values by region
#' @export
collect_beta_panel = function(cov_file,
                              panel_design,
                              sets_analyze = c(
                                  "mCRPC_PBMC",
                                  "NE_Adeno",
                                  "rockermeth_NEvsAdeno",
                                  "DMS_genesOfInterest",
                                  "Promoter_DEgenes",
                                  "CpG_site_DNAmeth_classifierJCI_2020"
                              )) {

    `%!in%` <- Negate(`%in%`)

    assertthat::are_equal(ncol(cov_file), 6)
    colnames(cov_file) = c("chr", "start", "end", "beta", "meth", "unmeth")

    cov_file$cov = with(cov_file, meth + unmeth)
    cov_gr = makeGRangesFromDataFrame(cov_file, keep.extra.columns = T)

    assertthat::are_equal(all(c("region_set", "reg_id") %in% colnames(panel_design)), TRUE)

    coll_rr = list()


    for (ss in unique(panel_design$region_set)) {
        if (ss %in% sets_analyze) {
            message(ss)

            ss_gr = panel_design %>%
                filter(region_set == ss)

            reg_reg = data.frame(ss_gr) %>% select(reg_id, region_set)

            ## Enlarge of 1 bp for cases with length = 1,
            ## as HOMER annotation caused start and end to invert

            if (ss == "CpG_site_DNAmeth_classifierJCI_2020") {
                ss_gr$end = ss_gr$end + 1
                ss_gr$start = ss_gr$start - 1
            }

            ss_gr = ss_gr %>%
                makeGRangesFromDataFrame(keep.extra.columns = T)

            ovv = suppressWarnings(findOverlaps(cov_gr, ss_gr))

            if (length(ovv) > 0) {
                rr = data.frame(cov_gr[queryHits(ovv),])
                rr$reg_id = ss_gr[subjectHits(ovv),]$reg_id
                rr$region_set = ss

                ## If there is no intersection, fill all values and return them
            } else {
                rr = reg_reg %>%
                    mutate(
                        seqnames = NA,
                        start = NA,
                        end = NA,
                        width = NA,
                        strand = NA,
                        beta = NA,
                        meth = NA,
                        unmeth = NA,
                        cov = NA
                    ) %>%
                    separate(
                        reg_id,
                        into = c("seqnames", "pos"),
                        sep = c(":"),
                        remove = F
                    ) %>%
                    separate(pos,
                             into = c("start", "end"),
                             sep = c("-")) %>%
                    mutate(start = as.numeric(start),
                           end = as.numeric(end)) %>%
                    select(
                        seqnames,
                        start,
                        end,
                        width,
                        strand,
                        beta,
                        meth,
                        unmeth,
                        cov,
                        reg_id,
                        region_set
                    )

                coll_rr[[ss]] = rr

            }

            ## Add regions that are not covered as NA values
            add_reg = reg_reg %>%
                mutate(
                    seqnames = NA,
                    start = NA,
                    end = NA,
                    width = NA,
                    strand = NA,
                    beta = NA,
                    meth = NA,
                    unmeth = NA,
                    cov = NA
                ) %>%
                filter(reg_id %!in% rr$reg_id) %>%
                separate(
                    reg_id,
                    into = c("seqnames", "pos"),
                    sep = c(":"),
                    remove = F
                ) %>%
                separate(pos,
                         into = c("start", "end"),
                         sep = c("-")) %>%
                mutate(start = as.numeric(start), end = as.numeric(end)) %>%
                select(
                    seqnames,
                    start,
                    end,
                    width,
                    strand,
                    beta,
                    meth,
                    unmeth,
                    cov,
                    reg_id,
                    region_set
                )

            rr = rbind(rr, add_reg) %>% arrange(reg_id)
            coll_rr[[ss]] = rr
        }

    }

    res = bind_rows(coll_rr)
    return(res)

}



#' Collapse statistics of regions
#'
#' @param res_collected the result of collect_beta_panel function
#'
#' @import dplyr tidyr
#' @return A data.frame of collapsed regions
#' @export
collapse_stats = function(res_collected, ssid) {

    assertthat::are_equal(
        colnames(res_collected),
        c(
            "seqnames",
            "start",
            "end",
            "width",
            "strand",
            "beta",
            "meth",
            "unmeth",
            "cov",
            "reg_id",
            "region_set"
        )
    )

    stats_collected = res_collected %>%
        mutate(beta = meth / cov) %>%
        group_by(region_set, reg_id) %>%
        summarise(
            mean_meth = mean(beta, na.rm = T),
            sd_meth = sd(beta, na.rm = T),
            n_cpgs = n(),
            tot_cov = sum(cov),
            tot_meth = sum(meth)
        ) %>%
        mutate(sampleID = ssid)

    return(stats_collected)

}
