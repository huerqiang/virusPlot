#' Strudel plot of virus insert
#'
#' @param virus_info The virus information data frame contains three columns.
#' The first column is the gene name, the second column is the start site,
#' and the third column is the stop site.
#' @param insert_info The virus insertion information data frame contains four columns.
#' The first column is the chromosome, the second column is the host insertion site,
#' the third column is the virus break site, and the fourth column is the number of reads.
#' @param virus_color color of virus rect.
#' @param host_color color of host rect.
#' @param label_virus label of virus.
#' @param label_host label of host.
#' @param hot_gene number of insert hot genes.
#' @param size_gene size of gene labels.
#' @param size_label size of label_virus and label_host.
#' @import ggplot2
#' @importFrom utils data
#' @importFrom methods is
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom ChIPseeker annotatePeak
#' @importFrom ChIPseeker annotatePeak
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom ggrepel geom_text_repel
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @return a gg object
#' @export
#'
#' @examples
#' data(insert_info)
#' virus_info <- data.frame(
#'       gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
#'       start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
#'       end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))
#' p <- strudel_plot(virus_info, insert_info)
strudel_plot <- function(virus_info, insert_info, virus_color = "#EAFEFF",
    host_color = "#EAFEFF", label_virus = "HPV16", label_host = "Host", hot_gene = 5,
    size_gene = 6.5, size_label = 6.5) {
    start <- num <- ymin <- ymax <- xmin <- xmax <- fill <- x <- gene <- SYMBOL <- NULL
    label <- hpv_loc <- y <- log10reads <- host_loc <- log10reads2 <- NULL
    virus_info <- as.data.frame(virus_info[, seq_len(3)])
    colnames(virus_info) <- c("gene", "start", "end")
    virus_info$start <- as.numeric(virus_info$start)
    virus_info$end <- as.numeric(virus_info$end)

    insert_info <- as.data.frame(insert_info[, seq_len(4)])
    colnames(insert_info) <- c("chr", "host_loc", "hpv_loc", "reads")
    insert_info$log10reads <- log10(insert_info$reads)
    # level
    virus_info <- virus_info[, seq_len(3)] |> as.data.frame()
    colnames(virus_info) <- c("gene", "start", "end")
    rownames(virus_info) <- virus_info$gene
    virus_info$start <- as.numeric(virus_info$start)
    virus_info$end <- as.numeric(virus_info$end)
    virus_info <- virus_info[order(virus_info$start), ]
    virus_info$level <- "level1"
    level2 <- 1 + which(virus_info$start[-1] -
        virus_info$end[seq_len(nrow(virus_info) - 1)] < 0)
    if (length(level2) > 0) {
        virus_info$level[level2] <- "level2"
        virus_info_level2 <- virus_info[level2, ]
        level3 <- 1 + which(virus_info_level2$start[-1] -
            virus_info_level2$end[seq_len(nrow(virus_info_level2) - 1)] < 0)

        if (length(level3) > 0) {
            virus_info_level2$level[level3] <- "level3"
        }
        virus_info[virus_info_level2$gene, "level"] <- virus_info_level2$level
    }
    max_reads <- ceiling(max(insert_info$log10reads))
    ymin <- max_reads + c(1, 0.6, 0.2)
    names(ymin) <- c("level1", "level2", "level3")

    rect_df <- data.frame(ymin = ymin[virus_info$level], xmin = virus_info$start, xmax = virus_info$end)
    rect_df$ymax <- rect_df$ymin + 0.38
    text_df <- data.frame(label = virus_info$gene,
        x = (virus_info$start + virus_info$end) / 2,
        y = (rect_df$ymin + rect_df$ymax) / 2)

    rect2_df <- data.frame(ymin = 0, ymax = max_reads, xmin = virus_info$start, xmax = virus_info$end)
    rect2_df <- rect2_df[seq(1, nrow(rect2_df), 2), ]
    p <- ggplot() + geom_rect(data = rect_df, mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), fill = "#FFBDBD") +
        geom_text(data = text_df, mapping = aes(x = x, y = y, label = label), color = "black", size = size_gene) +
        geom_rect(data = rect2_df, mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), fill = virus_color) +
        geom_point(data = insert_info, aes(x = hpv_loc, y = log10reads))  +
        annotate("text", x = max(virus_info$end) * 1.1, y = max_reads/2, label = label_virus, color = "black", size = size_label)
    # data(UCSC.HG38.Human.CytoBandIdeogram, package = "RCircos")
    # cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
    chr_length <- split(cyto.info$chromEnd, cyto.info$Chromosome) |> lapply(max) |> unlist()
    chr_length <- chr_length[paste("chr", c(1:22, "X"), sep = "")]
    len_sum <- sum(chr_length)

    insert_info_host_list <- split(insert_info[, ], insert_info$chr)

    insert_info_host_list <- insert_info_host_list[paste("chr", c(1:22, "X"), sep = "")]

    chr_length <- chr_length / (len_sum / max(virus_info$end))
    insert_info_host_list <- lapply(insert_info_host_list, function(x) {
        x$host_loc_old <- x$host_loc
        x$host_loc <- x$host_loc / (len_sum / max(virus_info$end))
        x
    })

    csum <- cumsum(chr_length)
    csum2 <- c(0, csum[1:22])
    names(csum2) <- names(csum)

    for (i in names(insert_info_host_list)) {
        insert_info_host_list[[i]][, "host_loc"] <- insert_info_host_list[[i]][, "host_loc"] + csum2[i]
    }
    insert_info_host <- do.call(rbind, insert_info_host_list) |> as.data.frame()


    insert_info_host$log10reads2 <- insert_info_host$log10reads - max_reads * 2.5

    hostfile <- data.frame(chr = names(csum), start = csum2, end = csum)
    colnames(hostfile) <- c("chr", "start", "end")
    hostfile$ymax <- -max_reads*1.5
    hostfile$ymin <- -max_reads*2.5
    hostfile$xmin <- hostfile$start
    hostfile$xmax <- hostfile$end

    hostfile2 <- hostfile[seq(1, nrow(hostfile), 2), ]
    p <- p + geom_rect(data = hostfile2, mapping = aes(ymin=ymin, ymax=ymax, xmin=xmin,  xmax = xmax), fill = host_color) +
        geom_point(data = insert_info_host, aes(x = host_loc, y = log10reads2))  +
        annotate("text", x = max(virus_info$end) * 1.1, y = -max_reads*2, label = label_host, color = "black", size = 6.5) +
        annotate("text", x = 0, y = -max_reads*2.65, label = "Chr.", color = "black", size = 6.5) +
        theme_classic() +
        geom_segment(data = insert_info_host,
            aes(x = hpv_loc, xend = host_loc),
            y = 0, yend = -max_reads*1.5,
            lty = 1, colour = "black", size = .2, alpha = 0.8)

    text_df <- data.frame(label = c(1:22, "x"),
        x = (hostfile$start + hostfile$end)/2,
        y = -max_reads*2.55)
    text_df[c(17, 19, 21), "y"] <- -max_reads*2.65
    y_df <- data.frame(label = as.character(c(rev(seq_len(max_reads)), 0, seq_len(max_reads), 0)),
        x = -max_reads,
        y = c(rev(seq_len(max_reads)), 0, -((1.5*max_reads) : (2.5*max_reads))))
    p <- p + geom_text(data = text_df, aes(label = label, x = x, y = y),
                  size = 5.2, color = "black") +
               geom_text(data = y_df, aes(label = label, x = x, y = y)) +
               ylab("log10(supporting reads)") +
               xlab("") +
               theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.line = element_blank(),
                     axis.title.y = element_text(size = 20))
    # hot gene anatation
    ranges <- GRanges(
        seqnames = insert_info_host$chr,
        ranges = IRanges(start = insert_info_host$host_loc_old, end = insert_info_host$host_loc_old)
    )
    txdb_38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno1 <- annotatePeak(ranges, level = "gene", annotate_multiple_region = TRUE,
                             TxDb=txdb_38, annoDb="org.Hs.eg.db", verbose = FALSE) |> 
                 suppressMessages()      						 
    peakAnno1 <- as.data.frame(peakAnno1)
    peakAnno1$id <- paste(peakAnno1[, 1], peakAnno1[, 2], sep = "_") 
    
    
    
    insert_info_host$id <- paste(insert_info_host$chr, insert_info_host$host_loc_old, sep = "_")
    insert_info_host <- insert_info_host[!duplicated(insert_info_host$id), ]
    rownames(insert_info_host) <- insert_info_host$id
    
    peakAnno1$log10reads <- insert_info_host[peakAnno1$id, "log10reads"] 
    peakAnno1$host_loc <- insert_info_host[peakAnno1$id, "host_loc"] 
    
    
    peakAnno1 <- peakAnno1[, c("seqnames", "log10reads", "host_loc", "SYMBOL")]
    peakAnno1 <- peakAnno1[order(peakAnno1$log10reads, decreasing = TRUE), ]
    
    genes <- NULL
    if (is(hot_gene, "numeric")) {
        genes <- unique(peakAnno1[, 4])[seq_len(hot_gene)]
    }
    if (is(hot_gene, "character")) {
        genes <- intersect(hot_gene, peakAnno1[, 4])
    }
    if (!is.null(genes) && length(genes) > 0) {
        peakAnno1 <- peakAnno1[peakAnno1$SYMBOL %in% genes, ]
        peakAnno1 <- peakAnno1[!duplicated(peakAnno1$SYMBOL), ] 
        # 怎么分层？同个染色体就分层？
        # 不行，相邻染色体还是会覆盖，那么相邻染色体再弄个分层？
        peakAnno1[, 1] <- as.character(peakAnno1[, 1])
        peakAnno1_list <- split(peakAnno1, peakAnno1[, 1])
        peakAnno1_list <- lapply(peakAnno1_list, function(x) {
            x[, 5] <- -10.2 - 0.8 * seq_len(nrow(x)) 
            return(x)
        })
        peakAnno1 <- do.call(rbind, peakAnno1_list)
        colnames(peakAnno1)[5] <- "y"
    }
    p + ggrepel::geom_text_repel(data = peakAnno1, aes(label = SYMBOL, x = host_loc, y = y),
              size = 5.2, color = "black")
}
