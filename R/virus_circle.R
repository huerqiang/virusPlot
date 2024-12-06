#' Circle plot of virus insert
#'
#' @param virus_info The virus information data frame contains three columns.
#' The first column is the gene name, the second column is the start site,
#' and the third column is the stop site.
#' @param insert_num The data frame of the number of breakpoints in each segment
#' contains three columns. The first column is the starting site of the segment,
#' the second column is the ending site of the segment, and the third column is
#' the number of insertion sites.
#' @param bar_color color of bar plot.
#' @param label a character of virus name.
#'
#' @return a gg object
#' @export
#'
#' @examples
#' virus_info <- data.frame(
#'          gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
#'          start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
#'          end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))
#'
#' insert_num <- data.frame(start = seq(1, 7801, 100),
#'     end = seq(101, 7901, 100),
#'     num = sample(1:28, 79, replace = TRUE))
#' p <- circle_virus(virus_info, insert_num)
circle_virus <- function(virus_info, insert_num, label = "HPV16",
                         bar_color = "#366CB4") {
    start <- num <- ymin <- ymax <- xmin <- xmax <- fill <- x <- gene <- NULL
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

    virus_info$x <- (virus_info$start + virus_info$end) / 2
    y <- c(-7.5, -13.5, -19.5)
    names(y) <- c("level1", "level2", "level3")
    virus_info$y <- y[virus_info$level]

    colnames(insert_num) <- c("start", "end", "num")

    col <- ggsci::pal_npg()(nrow(virus_info))
    ymax<- c(-5, -11, -17)
    names(ymax) <- c("level1", "level2", "level3")
    rect_df <- data.frame(ymax  = ymax[virus_info$level],
        xmin = virus_info[, 2], xmax = virus_info[, 3], fill = col)
    rect_df$ymin <- rect_df$ymax - 5
    text_df <- data.frame(label = as.character(seq(0, max(insert_num$end), 500)),
        x = seq(0, max(insert_num$end), 500), y = 38)
    # node为每隔10个，要超过最长的。
    n <- ceiling(max(insert_num$num) / 10)
    node_df <- data.frame(x=0, y = c(0,seq_len(n) * 10))
    y_df <- data.frame(label = as.character(node_df$y), x = 0, y = node_df$y)


    p <- ggplot() +
        geom_col(data = insert_num, mapping = aes(x = start, y = num),
                 fill = bar_color) +
        geom_rect(data = rect_df, mapping = aes(ymin=ymin, ymax=ymax,
                                                xmin = xmin, xmax = xmax,
                                                fill = I(fill))) +
        shadowtext::geom_shadowtext(data = virus_info, aes(x = x, y = y,
                                                        label = gene),
            colour = "black",bg.color = "white", bg.r = 0.1) +
        annotate("text", x = max(insert_num$end) / 2, y = -60, label = label,
            size = 8) +
        geom_rect(aes(ymin=max(insert_num$num), ymax=max(insert_num$num) + 1,
                      xmin=0, xmax = max(insert_num$end)),
            fill = "#DBDBDB")  +
        geom_rect(aes(ymin=max(insert_num$num) + 5,
                      ymax=max(insert_num$num) + 10,
                      xmin=0, xmax = max(insert_num$end)),
            fill = "#DBDBDB")  +
        geom_text(data = text_df, aes(label = label, x = x, y = y),
                  size = 4, color = "black") +
        geom_segment(x = 0, y = 0, xend = 0, yend = n * 10, size = 1,
                     color = "grey") +
        geom_point(data = node_df, aes(x = x, y = y), size = 1) +
        geom_text(data = y_df, aes(label = label, x = x, y = y),
                  size = 4, color = "black") +
        ylim(-60, NA) +
        xlim(0, max(virus_info[,3])) +
        theme_void() +
        coord_polar(theta='x', start = pi / 2)  +
        theme(legend.position = "none")
    p
}
