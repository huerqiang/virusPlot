#' Get virus annotation from NCBI
#'
#' @param accession_number accession_number, such as NC_001526.2.
#' @param email email address.
#' @param output Output file name, if not NULL, will be exported to the file.
#' @importFrom httr status_code
#' @importFrom httr GET
#' @importFrom httr content
#'
#' @return character
#' @export
#'
#' @examples
#' gene_features <- get_virus_annotation(accession_number = "NC_001526.2",
#'     email = "13766876214@163.com")
get_virus_annotation <- function(accession_number, email, output = NULL) {
    db <- "nuccore"
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=", db, "&id=", accession_number, "&rettype=gp&retmode=text&email=", email)
    response <- GET(url)
    gene_features <- NULL
    if (status_code(response) == 200) {
        gene_features <- content(x = response, as = "text", encoding = "UTF-8")
    } 
    if (!is.null(output)) {
        writeLines(gene_features, output)
    }
    return(gene_features)
}


#' Get virus genome
#'
#' @param accession_number, such as NC_001526.2.
#' @param email email address.
#' @param output Output file name, if not NULL, will be exported to the file.
#' @importFrom rentrez entrez_fetch
#'
#' @return character
#' @export
#'
#' @examples
#' genome <- get_virus_genom(accession_number = "NC_001526.2",
#'     email = "13766876214@163.com")
get_virus_genom <- function(accession_number, email, output = NULL) {
    db <- "nuccore"
    handle <- entrez_fetch(db, id = accession_number, rettype = "fasta", retmode = "text", email = email)
    if (!is.null(output)) {
        writeLines(handle, output)
    }
    handle <- read.table(text = handle, sep = "^", header = FALSE)
    return(handle)
}

#' Deal virus annotation
#'
#' @param gene_features result of get_virus_annotation()
#' @importFrom utils read.table
#'
#' @return data.frame
#' @export
#'
#' @examples
#' gene_features <- get_virus_annotation(accession_number = "NC_001526.2",
#'     email = "13766876214@163.com")
#' virus_info <- deal_virus_annotation(gene_features)
deal_virus_annotation <- function(gene_features) {
    aa <- read.table(text = gene_features, sep = "^", header = FALSE)[, 1]
    gene_lines <- grep("gene=", aa)
    gene_lines2 <- gene_lines - 1
    bb <- data.frame(aa[gene_lines], aa[gene_lines2])
    bb <- bb[grep("gene", bb[, 2]), ]
    bb[, 1] <- gsub(".*=", "", bb[, 1])
    bb[, 2] <- gsub("[]+gene[]+", "", bb[,2])
    bb[, 2] <- gsub(" ", "", bb[,2])
    bb[, 2] <- gsub(".*\\((.*)\\).*", "\\1", bb[, 2])
    anno <- strsplit(bb[, 2], ",")
    names(anno) <- bb[, 1]
    anno <- stack(anno)
    se <- strsplit(anno[, 1], "\\.\\.") |> do.call(rbind, args = _)
    start <- gsub("[^0-9]", "", se[, 1]) |> as.numeric()
    end <- gsub("[^0-9]", "", se[, 2]) |> as.numeric()
    return(data.frame(gene = anno[, 2], start = start, end = end))
}

