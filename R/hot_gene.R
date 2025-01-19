#' Get hot insert genes by Chi-square test
#'
#' @param virus_info The virus information data frame contains at least three
#' columns.
#' The first column is the gene name, the second column is the start site,
#' and the third column is the stop site.
#' @param insert_info The virus insertion information data frame contains at
#' least four columns.
#' The first column is the chromosome,
#' the second column is the host insertion site,
#' the third column is the virus break site, and the fourth column is the
#' number of reads.
#' @param tssRegion Region Range of TSS
#' @param TxDb TxDb or EnsDb annotation object
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom ChIPseeker annotatePeak
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom stats chisq.test
#' @importFrom IRanges findOverlaps
#' @importFrom stats fisher.test
#' @importFrom stats na.omit
#' @importFrom stats aggregate
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#'
#' @return list
#' @export
#'
#' @examples
#' data(insert_info)
#' virus_info <- data.frame(
#'       gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
#'       start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
#'       end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))
#' hot_gene <- get_hot_gene(virus_info, insert_info)
get_hot_gene <- function(virus_info, insert_info, tssRegion = c(-3000, 3000), TxDb = NULL) {
    if (ncol(insert_info) == 4) {
        insert_info$sample <- "sample1"
    }
    if (ncol(insert_info) < 4) {
        stop("insert_info should has at leat 4 columns: chr, host_loc, virus_loc, reads")
    }
    insert_info <- insert_info[, 1:5]
    sample_length <- unique(insert_info[, 5]) |> length()


    # annotate host gene
    ranges <- GRanges(
        seqnames = insert_info$chr,
        ranges = IRanges(start = insert_info$host_loc, end = insert_info$host_loc)
    )
    if (is.null(TxDb)) {
        TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    
    peakAnno1 <- annotatePeak(ranges, level = "gene", # annotate_multiple_region = TRUE,
                             TxDb=TxDb, annoDb="org.Hs.eg.db", tssRegion = tssRegion,
                             verbose = FALSE) |>
                 suppressMessages()
    peakAnno1 <- as.data.frame(peakAnno1)
    peakAnno1$id <- paste(peakAnno1[, 1], peakAnno1[, 2], sep = "_")
    insert_info$id <- paste(insert_info[, 1], insert_info[, 2], sep = "_")
    insert_info$SYMBOL <- peakAnno1[match(insert_info$id, peakAnno1$id), "SYMBOL"]

    # annotate virus gene
    circlerange <- GRanges(
       seqnames = "virus",
       ranges = IRanges(virus_info$start, end = virus_info$end))

    insert_hpv_range <- GRanges(
       seqnames = "virus",
       ranges = IRanges(insert_info[, 3], end = insert_info[, 3]))
    aa <- findOverlaps(insert_hpv_range,circlerange, ignore.strand = T)
    bb <- as.data.frame(aa)
    insert_info2 <- matrix(0, nrow = nrow(bb), ncol = 8)
    for (i in 1:nrow(bb)) {
        insert_info2[i, ] <- c(insert_info[bb[i, 1], ], virus_info[bb[i, 2], "gene"]) |> unlist()
    }
    colnames(insert_info2) <- c(colnames(insert_info), "virus_gene")

    # fisher.test test for host genes
    host_line <- as.data.frame(insert_info2)[, c("host_loc", "SYMBOL")] |> unique()
    gene_num_host <- table(host_line$SYMBOL) |> as.data.frame()
    gene_num_host[, 1] <- as.character(gene_num_host[, 1])
    gene_cov$length <- gene_cov$length + sum(abs(tssRegion))
    exon_length <- sum(gene_cov[, 1])
    gene_num_host$length <- gene_cov[match(gene_num_host[, 1], gene_cov$symbol), "length"]
    insert_sum_host <- sum(gene_num_host[, 2])
    gene_num_host$expect <- insert_sum_host * gene_num_host$length / exon_length
    gene_num_host <- gene_num_host[order(gene_num_host[, 2], decreasing = TRUE), ]
    gene_num_host <- na.omit(gene_num_host)
    colnames(gene_num_host) <- c("gene", "insert", "length", "expect")
    # Contingency table:
    # (inserted sites in gene A,  non-inserted sites in gene)
    # (inserted sites on non-gene A, non-inserted sites on non-gene A)
    # Number of bits should be gene length * number of samples
    # But the number of bits is too large, to divide bins. Default to the inserted bins first.
    gene_num_host$pvalue <- 1
    gene_num_host$noinsert <- (gene_num_host$length * sample_length - gene_num_host$insert) / 1000
    gene_num_host$insert_noA <- insert_sum_host - gene_num_host$insert
    gene_num_host$insert_noA_noinsert <- ((exon_length - gene_num_host$length) * sample_length - gene_num_host$insert_noA) / 1000
    for (i in 1:nrow(gene_num_host)) {
        df <- matrix(c(gene_num_host[i, "insert"], gene_num_host[i, "noinsert"],gene_num_host[i, "insert_noA"],gene_num_host[i, "insert_noA_noinsert"] ), nrow = 2)
        gene_num_host[i, "pvalue"] <- fisher.test(df)$p.value
    }
    result_host <- gene_num_host[order(gene_num_host$pvalue), c("gene", "insert", "expect", "pvalue")]

    # Chi-square test for virus genes
    hpv_line <- as.data.frame(insert_info2)[, c(3, 8)] |> unique()
    colnames(hpv_line)[2] <- "gene"
    gene_num_hpv <- table(hpv_line$gene) |> as.data.frame()
    gene_num_hpv[, 1] <- as.character(gene_num_hpv[, 1])

    ## gene length
    virus_info$length <- virus_info$end - virus_info$start
    result <- aggregate(. ~ gene, data = virus_info[, c("gene", "length")], FUN = sum)
    virus_info$length <- result[match(virus_info$gene, result[, 1]), 2]
    exon_length <- sum(result[, 2])
    gene_num_hpv$length <- result[match(gene_num_hpv[, 1], result[, 1]), 2]
    insert_sum_hpv <- sum(gene_num_hpv[, 2])
    gene_num_hpv$expect <- insert_sum_hpv * gene_num_hpv$length / exon_length
    gene_num_hpv <- gene_num_hpv[order(gene_num_hpv[, 2], decreasing = TRUE), ]
    colnames(gene_num_hpv) <- c("gene", "insert", "length", "expect")
    # Contingency table:
    # (inserted sites in gene A,  non-inserted sites in gene)
    # (inserted sites on non-gene A, non-inserted sites on non-gene A)
    gene_num_hpv$pvalue <- 1
    gene_num_hpv$noinsert <- gene_num_hpv$length * sample_length - gene_num_hpv$insert
    gene_num_hpv$insert_noA <- insert_sum_hpv - gene_num_hpv$insert
    gene_num_hpv$insert_noA_noinsert <- (exon_length - gene_num_hpv$length) * sample_length - gene_num_hpv$insert_noA
    for (i in 1:nrow(gene_num_hpv)) {
        df <- matrix(c(gene_num_hpv[i, "insert"], gene_num_hpv[i, "noinsert"],gene_num_hpv[i, "insert_noA"],gene_num_hpv[i, "insert_noA_noinsert"] ), nrow = 2)
        gene_num_hpv[i, "pvalue"] <- chisq.test(df)$p.value
    }
    result_hpv <- gene_num_hpv[order(gene_num_hpv$pvalue), c("gene", "insert", "expect", "pvalue")]
    return(list(host = result_host, virus = result_hpv))
}

#' Title
#'
#' @param hot_gene_result result of get_hot_gene()
#' @param hot_gene_host hot insert genes of host to plot, can be a number of a vector of genes
#' @param hot_gene_virus hot insert genes of virus to plot, can be a number of a vector of genes
#' @param break_y If ture (the default), use ggbreak to set an y axis break.
#' @param observed_color color of observed value.
#' @param expected_color color of expected value.
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggbreak scale_y_break
#'
#' @return gg object
#' @export
#'
#' @examples
#' data(insert_info)
#' virus_info <- data.frame(
#'       gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
#'       start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
#'       end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))
#' hot_gene <- get_hot_gene(virus_info, insert_info)
#' insert_plot <- hot_gene_plot(hot_gene)
hot_gene_plot <- function(hot_gene_result, hot_gene_host = 5, hot_gene_virus = 8,
    break_y = TRUE, observed_color = "grey30", expected_color = "grey60") {
    gene <- group <- NULL
    if (is(hot_gene_host, "numeric")) {
        hot_gene_host <- min(hot_gene_host, nrow(hot_gene_result$host))
        result_host <- hot_gene_result$host[seq_len(hot_gene_host), ]
    } else {
        if (is(hot_gene_host, "character")) {
            result_host <- hot_gene_result$host[hot_gene_result$host[, 1] %in% hot_gene_host, ]
        } else {
            stop("hot_gene_host should be a number or a vector of genes")
        }
    }

    if (is(hot_gene_virus, "numeric")) {
        hot_gene_virus <- min(hot_gene_virus, nrow(hot_gene_result$virus))
        result_virus <- hot_gene_result$virus[seq_len(hot_gene_virus), ]
    } else {
        if (is(hot_gene_virus, "character")) {
            result_virus <- hot_gene_result$virus[hot_gene_result$virus[, 1] %in% hot_gene_virus, ]
        } else {
            stop("hot_gene_virus should be a number or a vector of genes")
        }
    }

    if (nrow(result_host) == 0) {
        stop("No hot insert gene in host")
    }

    if (nrow(result_virus) == 0) {
        stop("No hot insert gene in virus")
    }


    long <- result_host[, 1:3] %>%
      pivot_longer(cols = !gene,
                   names_to = "group",
                   values_to = "insert")
    long$group <- gsub("insert", "Observed", long$group)
    long$group <- gsub("expect", "Expected", long$group)


    rmsk_pvalue2 <- result_host
    rmsk_pvalue2$y <- 0
    for (i in 1:nrow(rmsk_pvalue2)) {
        rmsk_pvalue2$y[i] <- max(result_host[i, "insert"], result_host[i, "expect"]) + 2.8
    }
    # rmsk_pvalue2$y2 <- rmsk_pvalue2$y + 1.3

    # rmsk_pvalue2$xstart <- 1:nrow(rmsk_pvalue2) - 0.2
    # rmsk_pvalue2$xend <- 1:nrow(rmsk_pvalue2) + 0.2


    long$group <- factor(long$group, levels = c("Observed", "Expected"))
    p <- ggplot(long, aes_string(x = "gene", y = "insert"))+
      geom_bar(mapping = aes(fill = group), stat = 'identity',position = position_dodge(0.9)) +#使柱子并排放置
      scale_fill_manual(values=c(Observed = observed_color, Expected = expected_color)) +
      theme(text=element_text(family="Songti SC",size=12,face = "bold"),
            axis.text.x = element_text(size=10)) + # 设置X轴文字大小
      # annotate("text", x = rmsk_pvalue2$element, y = rmsk_pvalue2$y2, label = paste("P =", round(rmsk_pvalue2$pvalue, 9))) +
      annotate("text", x = rmsk_pvalue2$gene, y = rmsk_pvalue2$y, label = paste("P =", format(rmsk_pvalue2$pvalue, scientific = TRUE, digits = 3))) +
      theme_classic() +
      # theme_dose() +
      # geom_segment(mapping = aes(x = xstart, y = y, xend = xend, yend = y), data = rmsk_pvalue2) +
      xlab("Human") +
      ylab("Number of breakpoints")
    # p + scale_y_break(c(10, 15, 70, 120, 1300, 1400), scales="free") + ylim(0, 5250) +
        # theme(legend.text=element_text(size=15)) +
        # theme(legend.title=element_text(size=15))
    if (break_y) {
        b1 <- min(long$insert[long$group == "Expected"])
        b2 <- max(long$insert[long$group == "Expected"])
        p <- p + scale_y_break(c(b1, b2), scales="free")
    }
    p_host <- p + theme(legend.text=element_text(size=15)) +
        theme(legend.title=element_text(size=15))

    long <- result_virus[, 1:3] %>%
    pivot_longer(cols = !gene,
                 names_to = "group",
                 values_to = "insert")
    long$group <- gsub("insert", "Observed", long$group)
    long$group <- gsub("expect", "Expected", long$group)


    rmsk_pvalue2 <- result_virus
    rmsk_pvalue2$y <- 0
    for (i in 1:nrow(rmsk_pvalue2)) {
        rmsk_pvalue2$y[i] <- max(result_virus[i, "insert"], result_virus[i, "expect"]) + 10
    }
    # rmsk_pvalue2$y2 <- rmsk_pvalue2$y + 1.3

    # rmsk_pvalue2$xstart <- 1:nrow(rmsk_pvalue2) - 0.2
    # rmsk_pvalue2$xend <- 1:nrow(rmsk_pvalue2) + 0.2

    long$group <- factor(long$group, levels = c("Observed", "Expected"))

    p <- ggplot(long, aes_string(x = "gene", y = "insert"))+
      geom_bar(mapping = aes(fill = group), stat = 'identity',position = position_dodge(0.9)) +#使柱子并排放置
      scale_fill_manual(values=c(Observed = observed_color, Expected = expected_color)) +
      theme(text=element_text(family="Songti SC",size=12,face = "bold"),
            axis.text.x = element_text(size=10)) + # 设置X轴文字大小
      # annotate("text", x = rmsk_pvalue2$element, y = rmsk_pvalue2$y2, label = paste("P =", round(rmsk_pvalue2$pvalue, 9))) +
      annotate("text", x = rmsk_pvalue2$gene, y = rmsk_pvalue2$y, label = paste("P =", format(rmsk_pvalue2$pvalue, scientific = TRUE, digits = 3))) +
      theme_classic() +
      # theme_dose() +
      # geom_segment(mapping = aes(x = xstart, y = y, xend = xend, yend = y), data = rmsk_pvalue2) +
      # xlab("HPV16") +
      xlab("") +
      ylab("Number of breakpoints")
    p_virus <- p + theme(legend.text=element_text(size=15)) +
    theme(legend.title=element_text(size=15))
    return(list(p_host = p_host, p_virus = p_virus))

}
