#' Uniform input
#'
#' @param file_in path to the input file with DNA methylation counts
#' @param format the format in which files are expected
#' @param all_cpg the position of all CpGs in the human genome
#' @param min_cov the minimum coverage required to retain a CpG site
#' @param ch_file chain file for hg38 to hg19 conversion - used for datasets with hg38 data
#'
#' @importFrom data.table fread
#' @import dplyr tidyr
#' @return A uniformed version of the count from the target file
#' @export
uniform_input = function(file_in,
                         format,
                         all_cpg,
                         min_cov = 10,
                         ch_file) {

    format = match.arg(
        format,
        choices = c(
            "bismark",
            "epicore",
            "cfdna",
            "epicore_old",
            "feng",
            "blueprint",
            "cpgea",
            "atlas"
        )
    )

    ## The uniform out is chr, start, end, meth_count, coverage
    ## The single position is 1 based and indicated the C site in the CpG dinucleotide,
    ## considering the TOP strand

    ## The format of epicore facility (Beltran2016 and organoids) and FH shared WGBS
    if (format == "epicore") {
        df_in = fread(
            file_in,
            sep = "\t",
            quote = "",
            header = T,
            data.table = F,
            stringsAsFactors = FALSE,
            nrows = read_mode
        )

        assertthat::are_equal(
            colnames(df_in),
            c(
                "chrBase",
                "chr",
                "base",
                "strand",
                "coverage",
                "freqC",
                "freqT"
            )
        )

        df_in = df_in %>%
            mutate(
                start = base,
                end = base,
                beta = as.numeric(freqC / 100),
                cov = coverage,
                strand = ifelse(strand == "F", "+", "-")
            ) %>%
            mutate(meth = round(beta * cov)) %>%
            select(chr, start, end, strand, meth, cov) %>%
            arrange(chr, start)

        uniformed = df_in %>%
            filter(end == start) %>%
            select(chr, start, end, meth, cov) %>%
            filter(cov >= min_cov)


    }

    ## The format of ctDNA cohorts (JCI2020 and unpub2019)
    if (format == "cfdna") {
        df_in = fread(
            file_in,
            sep = "\t",
            quote = "",
            header = F,
            data.table = F,
            stringsAsFactors = FALSE,
            nrows = read_mode
        )

        ## Single position, 1-based

        assertthat::are_equal(colnames(df_in), paste0("V", seq(1, 7)))

        colnames(df_in) = c("pos_id",
                            "chr",
                            "start",
                            "strand",
                            "cov",
                            "beta",
                            "minusBeta")
        df_in$end = df_in$start

        df_in$strand = ifelse(df_in$strand == "F", "+", "-")

        uniformed = df_in %>%
            mutate(meth = cov * beta / 100) %>%
            filter(end == start) %>%
            select(chr, start, end, meth, cov) %>%
            filter(cov >= min_cov)


    }

    ## The format of old tissue data (Lin2013)
    if (format == "epicore_old") {
        df_in = fread(
            file_in,
            sep = "\t",
            quote = "",
            header = T,
            data.table = F,
            stringsAsFactors = FALSE,
            nrows = read_mode
        )

        assertthat::are_equal(colnames(df_in),
                              c("chromosome", "position", "numCs", "numTs"))
        df_in = df_in %>% arrange(chromosome, position)

        ## C top strand

        uniformed = df_in %>%
            mutate(meth = numCs, cov = numCs + numTs) %>%
            select(chr = chromosome,
                   start = position,
                   end = position,
                   meth,
                   cov) %>%
            filter(cov >= min_cov)

    }

    ## The format of Felix cohort (NatGen2020)
    if (format == "feng") {
        df_in = fread(
            file_in,
            sep = "\t",
            quote = "",
            header = F,
            data.table = F,
            stringsAsFactors = FALSE,
            nrows = read_mode
        )


        assertthat::are_equal(colnames(df_in), paste0("V", seq(1, 8)))
        colnames(df_in) = c("mp",
                            "start",
                            "end",
                            "strand",
                            "meth",
                            "unmeth",
                            "cov",
                            "beta")


        uniformed = df_in %>%
            select(mp, start, end, strand, meth, cov)


        chrmap = NEMOcode::contigs_hg38

        chrmap = chrmap %>% select(mp = V4, chr = V1) %>%
            filter(chr != "unplaced")


        uniformed = uniformed %>%
            inner_join(chrmap, by = "mp") %>%
            select(chr, start, end, strand, meth, cov) %>%
            filter(cov >= min_cov)



        # Lifting

        uniformed = data.frame(unlist(
            rtracklayer::liftOver(
                x = makeGRangesFromDataFrame(uniformed, keep.extra.columns = T),
                chain = ch_file
            )
        ))


        ## Minus strand gets shifted by 1 and collapsed in the forward

        uniformed = uniformed %>%
            filter(end == start) %>%
            select(chr = seqnames, start, end, meth, cov)

        ## Shift by one as the position reported in the C in the bot strand



    }

    ## The Blueprint format
    if (format == "blueprint") {
        df_in = fread(
            file_in,
            sep = "\t",
            quote = "",
            header = F,
            data.table = F,
            stringsAsFactors = FALSE,
            nrows = read_mode
        )


        assertthat::are_equal(colnames(df_in), paste0("V", seq(1, 6)))
        colnames(df_in) = c("chr", "start", "end", "meth_str", "meth_mil", "strand")


        uniformed = df_in %>%
            mutate(meth_str = gsub("'", "", meth_str)) %>%
            separate(meth_str,
                     into = c("meth", "cov"),
                     sep = "/") %>%
            mutate(meth = as.integer(meth), cov = as.integer(cov)) %>%
            select(chr, start, end, strand, meth, cov) %>%
            filter(cov >= min_cov)

        # Lifting
        uniformed = rtracklayer::liftOver(x = makeGRangesFromDataFrame(uniformed, keep.extra.columns = T),
                                          chain = ch_file)

        ## The liftover might casue shifts in coordinates in some rare cases. Drop if end - start != 0.

        ## 0 based, forward strand. Convert to 1 based.

        uniformed = data.frame(unlist(uniformed)) %>%
            select(chr = seqnames, start, strand, end, meth, cov) %>%
            mutate(start = start + 1, end = end - 1) %>%
            filter(end == start) %>%
            group_by(chr, start, end) %>%
            mutate(meth = sum(meth), cov = sum(cov)) %>%
            select(chr, start, end, meth, cov)

    }

    ## The format of CPGEA cohort (Nat2020)
    if (format == "cpgea") {
        df_in = fread(
            file_in,
            sep = "\t",
            quote = "",
            header = F,
            data.table = F,
            stringsAsFactors = FALSE,
            nrows = read_mode
        )

        assertthat::are_equal(colnames(df_in), paste0("V", seq(1, 7)))

        colnames(df_in) = c("chr", "start", "end", "meth", "unmeth", "cov", "beta")
        df_in$chr = paste0("chr", df_in$chr)

        uniformed = df_in %>%
            select(chr, start, end, meth, cov) %>%
            filter(cov >= min_cov)

        uniformed = rtracklayer::liftOver(x = makeGRangesFromDataFrame(uniformed, keep.extra.columns = T),
                                          chain = ch_file)


        ## The liftover might casue shifts in coordinates in some rare cases. Drop if end - start != 0.
        ## These files do not have strand information, so we need to overlap with known CpGs and then collapse
        ## them to Top strand Cs

        uniformed = data.frame(unlist(uniformed)) %>%
            filter(end == start) %>%
            select(chr = seqnames, start, end, meth, cov) %>%
            makeGRangesFromDataFrame(keep.extra.columns = T)


        ovv = findOverlaps(uniformed, all_cpg)

        uniformed = cbind(data.frame(uniformed[queryHits(ovv)]),
                          pos = all_cpg[subjectHits(ovv)]$pos)


        uniformed = uniformed %>%
            group_by(pos) %>%
            summarise(meth = sum(meth), cov = sum(cov)) %>%
            separate(pos,
                     into = c("chr", "pos"),
                     sep = c("\\.")) %>%
            mutate(pos = as.numeric(pos)) %>%
            arrange(chr, pos) %>%
            select(chr, start = pos, end = pos, meth, cov)


    }

    ## The format of the DNA methylation Atlas
    if (format == "atlas") {
        df_in = fread(
            file_in,
            sep = "\t",
            quote = "",
            header = F,
            data.table = F,
            stringsAsFactors = FALSE,
            nrows = read_mode
        )

        assertthat::are_equal(colnames(df_in), paste0("V", seq(1, 5)))

        ## Coordinates indicate the CpG dinucleotide, centering on the C in the TOP strand.
        ## Map back to one based centered on C

        colnames(df_in) = c("chr", "start", "end", "meth", "cov")

        df_in$start = df_in$start + 1
        df_in$end = df_in$end - 1

        uniformed = df_in %>% filter(cov >= min_cov)

    }

    ## The format of our runs, obtained with bismark aligner
    if (format == "bismark") {
        df_in = fread(
            file_in,
            sep = "\t",
            quote = "",
            header = F,
            data.table = F,
            stringsAsFactors = FALSE,
            nrows = read_mode
        )

        ## This output is 1 based, centered on the C in the REVERSE strand.
        ## Subtract 1 to go in forward strand

        colnames(df_in) = c("chr", "start", "end", "beta", "meth", "unmeth")

        uniformed = df_in %>%
            mutate(cov = meth + unmeth) %>%
            filter(cov >= min_cov) %>%
            select(chr, start, end, meth, cov)

    }


    ## Overlap with CpG intervals

    gr_unif = makeGRangesFromDataFrame(uniformed, keep.extra.columns = T)
    ovv = findOverlaps(gr_unif, all_cpgs)

    uniformed = data.frame(gr_unif[queryHits(ovv)], pos = all_cpgs[subjectHits(ovv)]$pos)

    uniformed = uniformed %>%
        group_by(pos) %>%
        summarise(meth = sum(meth), cov = sum(cov)) %>%
        separate(
            pos,
            into = c("chr", "start"),
            sep = c("\\."),
            remove = T
        ) %>%
        mutate(start = as.numeric(start)) %>%
        mutate(end = start) %>%
        select(chr, start, end, meth, cov) %>%
        arrange(chr, start) %>%
        data.frame()

    ## Check that only 1 strand is present
    assertthat::assert_that(assertthat::are_equal(!any(
        paste0(uniformed$chr, ".", uniformed$start) %in%
            paste0(uniformed$chr, ".", uniformed$start + 1)
    ), TRUE))

    return(uniformed)

}

