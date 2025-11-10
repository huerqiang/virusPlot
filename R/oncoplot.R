#' oncoplot_fill
#'
#' @param breaks breaks of scale_fill_manual
#' @param varis_color namesd vector, color of varis
#' @param name name of fill legend
#' @param na.value color of NA value.
#' @importFrom utils getFromNamespace
#'
#' @return a gg object
#' @export
oncoplot_fill <- function(breaks=NULL, varis_color = NULL,
                          name = NULL, na.value = "#bdbdbd") {
    # vc_col <- get_vcColors(websafe = FALSE)  # maftools color setting
    # if (is.null(values)) {
        # values <- vc_col
    # }
    if (is.null(breaks)) {
        vc_lev <- names(varis_color)
        breaks <- factor(vc_lev, levels = rev(vc_lev))
    }
    scale_fill_manual(name = NULL,
            breaks = breaks,
            values = varis_color,
            na.value = na.value)
}

oncoprint_row_order = function(count_matrix, n_mut) {
	order(rowSums(count_matrix), n_mut, decreasing = TRUE)
}

#' main plot
#'
#' @param mat matrix, each row is a gene,
#' and each column is a sample
#' @param varis_color namesd vector, color of varis
#' @param na.value color of na value.
#'
#' @return a gg object
#' @export
oncoplot_main <- function(mat, varis_color, na.value = "#F3F5F7") {
    oncoplot_scale <- getFromNamespace("oncoplot_scale", "aplotExtra")
    mat_01 <- matrix(0, nrow(mat), ncol(mat))
    colnames(mat_01) <- colnames(mat)
    rownames(mat_01) <- rownames(mat)
    mat_01[mat != " "] <- 1

    oncoprint_column_order = function(count_matrix, row_order) {
	    scoreCol = function(x) {
	    	score = 0
	    	for(i in 1:length(x)) {
	    		if(x[i]) {
	    			score = score + 2^(length(x)-i*1/x[i])
	    		}
	    	}
	    	return(score)
	    }
	    scores = apply(count_matrix[row_order, ,drop = FALSE], 2, scoreCol)
	    order(scores, decreasing=TRUE)
    }

    n_mut <- rowSums(mat_01)

    row_order <- oncoprint_row_order(count_matrix = mat_01, n_mut = n_mut)
    column_order <- oncoprint_column_order(count_matrix = mat_01, row_order = row_order)

    sample_level <- colnames(mat)[column_order]
    gene_level <- rownames(mat)[row_order]

    mat$Gene <- rownames(mat)
    d <- tidyr::pivot_longer(mat,cols =-c('Gene'), names_to = "Sample",values_to = "Type")
    d$Type[d$Type == " "] <- NA #"Non_Mut"
    d$Gene <- factor(d$Gene, levels = rev(gene_level))
    d$Sample <- factor(d$Sample, levels = sample_level)
    p <- ggplot(d, aes(x=.data$Sample, y=.data$Gene, fill=.data$Type)) +
        # geom_tile(colour="white", linewidth=.05) +
        geom_raster() +
        list(theme_minimal(), ggfun::theme_noxaxis(), theme(legend.position = "none", panel.grid.major = element_blank()),
          oncoplot_scale(continuous = FALSE, scale='y'),
          oncoplot_fill(varis_color = varis_color, na.value = na.value))+
        theme(legend.position = "bottom",
            axis.text.y.left=element_text(face='italic'))
    leg <- ggfun::get_legend(p)
    p + theme(legend.position = "none",
            plot.margin = margin(b = 60, r = 5)) +
        annotation_custom(leg, ymin=-3, ymax=-1.5) +
        coord_cartesian(clip='off') +
    ylab("") +
    xlab("")
}

oncoprint_row_order = function(count_matrix, n_mut) {
    order(rowSums(count_matrix), n_mut, decreasing = TRUE)
}
oncoprint_column_order = function(count_matrix, row_order) {
    scoreCol = function(x) {
    	score = 0
    	for(i in 1:length(x)) {
    		if(x[i]) {
    			score = score + 2^(length(x)-i*1/x[i])
    		}
    	}
    	return(score)
    }
    scores = apply(count_matrix[row_order, ,drop = FALSE], 2, scoreCol)
    order(scores, decreasing=TRUE)
}

#' oncoplot_sample
#'
#' @param mat matrix, each row is a gene,
#' and each column is a sample
#' @param varis_color namesd vector, color of varis
#'
#' @return a gg object
#' @export
oncoplot_sample <- function(mat, varis_color) {
    oncoplot_scale <- getFromNamespace("oncoplot_scale", "aplotExtra")
    mat_01 <- matrix(0, nrow(mat), ncol(mat))
    colnames(mat_01) <- colnames(mat)
    rownames(mat_01) <- rownames(mat)
    mat_01[mat != " "] <- 1

    n_mut <- rowSums(mat_01)

    row_order <- oncoprint_row_order(count_matrix = mat_01, n_mut = n_mut)
    column_order <- oncoprint_column_order(count_matrix = mat_01, row_order = row_order)

    sample_level <- colnames(mat)[column_order]

    vars <- mat |> as.matrix() |> as.vector() |> unique()
    dMat <- matrix(0, nrow = length(vars), ncol = ncol(mat))
    rownames(dMat) <- vars
    colnames(dMat) <- colnames(mat)
    for (i in colnames(mat)) {
        tableDf <- mat[, i] |> table() |> as.data.frame()
        tableDf[, 1] <- as.character(tableDf[, 1])
        dMat[tableDf[, 1], i] <- dMat[tableDf[, 1], i] + tableDf[, 2]
    }
    dMat <- as.data.frame(dMat)
    dMat$Type <- rownames(dMat)
    d <- tidyr::pivot_longer(dMat,cols =-c('Type'), names_to = "Sample",values_to = "Freq")
    d <- d[d$Type != " ", ] 


    d$Sample <- factor(d$Sample, levels = sample_level)


    p <- ggplot(d, aes(x=.data$Sample,y=.data$Freq,fill=.data$Type)) +
        geom_col(position="stack") +
        list(
            theme_minimal(),
            ggfun::theme_noxaxis(),
            oncoplot_fill(varis_color = varis_color),
            theme(legend.position = "none", panel.grid.major = element_blank()),
            oncoplot_scale(continuous = TRUE, scale = 'y'),
            xlab(NULL),
            ylab(NULL)
        )

    p  + annotation_custom(grob = textGrob(
            # label = "TMB",
            label = NULL,
            rot = 90,
            gp = gpar(fontsize=11),
            x = -30, default.units = "pt"
        )) +
        coord_cartesian(clip="off") +
        theme(axis.title.y = element_blank(),
              plot.margin = margin(l=30),
              panel.grid=element_blank())
}


#' oncoplot_gene
#'
#' @param mat matrix, each row is a gene,
#' and each column is a sample
#' @param ylab ylab
#' @param varis_color namesd vector, color of varis
#' @import grid
#'
#' @return a gg object
#' @export
oncoplot_gene <- function(mat, ylab = 'percentage', varis_color) {
    oncoplot_scale <- getFromNamespace("oncoplot_scale", "aplotExtra")
    mat_01 <- matrix(0, nrow(mat), ncol(mat))
    colnames(mat_01) <- colnames(mat)
    rownames(mat_01) <- rownames(mat)
    mat_01[mat != " "] <- 1
    n_mut <- rowSums(mat_01)

    row_order <- oncoprint_row_order(count_matrix = mat_01, n_mut = n_mut)
    column_order <- oncoprint_column_order(count_matrix = mat_01, row_order = row_order)

    sample_level <- colnames(mat)[column_order]
    gene_level <- rownames(mat)[row_order]

    mat$Gene <- rownames(mat)
    d <- tidyr::pivot_longer(mat,cols =-c('Gene'), names_to = "Sample",values_to = "Type")
    d$Type[d$Type == " "] <- NA #"Non_Mut"
    d$Gene <- factor(d$Gene, levels = rev(gene_level))
    d$Sample <- factor(d$Sample, levels = sample_level)
    d <- d[!is.na(d$Type), ]
    ylab <- match.arg(ylab, c("gene", "percentage"))


    p <- ggplot(d, aes(y = .data$Gene, fill = .data$Type)) +
        geom_bar(position='stack', orientation = 'y') +
        list(
            theme_minimal(),
            theme(legend.position = "none", panel.grid.major = element_blank()),
            oncoplot_scale(continuous = TRUE, scale='none'),
            oncoplot_fill(varis_color = varis_color),
            # oncoplot_fill(),
            xlab(NULL),
            ylab(NULL)
        )
        # oncoplot_setting(noxaxis = FALSE, scale='none') #+


    if (ylab == 'percentage') {
        # numMat <- get_oncoplot_numericMatrix(maf, genes)
        numMat <- mat_01
        # totSamps = as.numeric(maf@summary[3, 'summary'])
        totSamps = ncol(mat_01)
        percent <- apply(numMat, 1, function(x) sum(x >0))/totSamps
        percent_alt <- paste0(round(percent * 100), '%')
        p <- p + scale_y_discrete(breaks = rownames(numMat),
                                labels = percent_alt,
                                expand = c(0, 0))
    }

    p <- p + annotation_custom(grob = textGrob(
                label = NULL,
                gp = gpar(fontsize=11)
        ), ymin=-2, ymax=-1) +
    coord_cartesian(clip="off") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(b=40),
          panel.grid=element_blank())

    return(p)
}



#' oncoplot_clinical
#'
#' @param clinical_data dataframe of clinical data,
#' the first column is sample ID.
#' @param mapping mapping.
#' @param clinical_color namesd vector, color of clinical varis.
#' @param p_main p_main.
#' @importFrom rlang quo_name
#'
#' @return a gg object
#' @export
oncoplot_clinical <- function(clinical_data, mapping, clinical_color = NULL, p_main = p_main) {
    level <- levels(p_main$data$Sample)
    rownames(clinical_data) <- clinical_data[, 1]
    clinical_data <- clinical_data[level, ]
    clinical_data$ID <- factor(clinical_data$ID, levels = level)
    p <- ggplot(data = clinical_data) + geom_tile(mapping = mapping)+
        scale_y_discrete(position = "right", expand = c(0,0)) +
        theme(panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              plot.background = element_blank())
    if (!is.null(clinical_color)) {
        fill_var <- quo_name(mapping$fill)
        colors <- clinical_data[, fill_var] |> unique()
        p <- p + scale_fill_manual(values = clinical_color[colors])
    }
    p
}


insert_top <- function(.data, plot, height = 1) {
  dynamic_call("aplot", "insert_tb", 
               .data = .data, plot = plot, 
               height = height, side = "top")
}


insert_bottom <- function(.data, plot, height = 1) {
  dynamic_call("aplot", "insert_tb", 
               .data = .data, plot = plot, 
               height = height, side = "bottom")
}


insert_right <- function(.data, plot, width = 1) {
  dynamic_call("aplot", "insert_lr", 
               .data = .data, plot = plot, 
               width = width, side = "right")
}


insert_left <- function(.data, plot, width = 1) {
  dynamic_call("aplot", "insert_lr", 
               .data = .data, plot = plot, 
               width = width, side = "left")
}



#' main function of oncoplot
#'
#' @param mat matrix, each row is a gene,
#' and each column is a sample
#' @param varis_color namesd vector, color of varis
#' @param clinical dataframe of clinical data,
#' the first column is sample ID.
#' @param clinical_color namesd vector, color of clinical varis.
#' @param na.value color of na value.
# @import aplot
#'
#' @return an aplot object
#' @export
oncoplot <- function(mat, varis_color = NULL, clinical, clinical_color = NULL, na.value) {
    # library(aplot)
    p_main <- oncoplot_main(mat, varis_color = varis_color, na.value = na.value)
    p_top <- oncoplot_sample(mat, varis_color = varis_color)
    p_right <- oncoplot_gene(mat, ylab = 'percentage', varis_color = varis_color)
    p_spacer <- ggplot() + ggfun::theme_transparent()
    clinicals <- colnames(clinical)[-1]
    ID <- colnames(clinical)[1]
    clinical2 <- clinical
    for (i in clinicals) {
        clinical2[[paste0("y_", i)]] <- i
    }

    p_cli_list <- lapply(clinicals, function(k) {
        oncoplot_clinical(clinical_data = clinical2,
                          mapping = aes_string(x = ID, y = paste0("y_", k), fill = k),
                          clinical_color = clinical_color, p_main = p_main)
    })
    p <- p_main
    for (i in 1:length(p_cli_list)) {
        p <- p |> insert_top(p_cli_list[[i]], height=.02)
    }

    p |> insert_top(p_top, height=.2) |>
        insert_right(p_spacer, width=.01) |>
        insert_right(p_right, width=.2) |>
        insert_bottom(p_spacer, height=.1)
}


