#' Import a network text file
#'
#' This function reads in a TAB-delimited file with the global network data needed for NEA.
#' Checks the format of the given NETWORK file, processes unique genes from the given network file,  and counts the network links.
#' @param tbl Input network file.
#' @param Lowercase Render gene/protein IDs lower-case.
#' @param col.1 Number of the column in the input file that corresponds to node 1 (gene/protein ID).
#' @param col.2 Number of the column in the input file that corresponds to node 2 (gene/protein ID).
#' @param col.score Number of column where edge confidence score is found.
#' @param min.score Minimum confidence score for an edge to be included in the network.
#' @param n.top  Number of edges to be included in the network, top when ranked by confidence score (is alternative to min.score).
#' @param echo if execution progress should be reported.

#' @return \code{$Ntotal}, the total number of edges in the network and \code{$links}, an object with members (named by gene/protein IDs) that correspond to each node in the global network. Each entry of the list \code{$links} contains a vector of neighbours of the respective node.
#' @seealso \code{\link{mutations2ags}}, \code{\link{samples2ags}}
#' @examples
#' data(net.kegg)
#' net <- import.net(net.kegg)
#' summary(net$links)
#' print(paste(names(net$links)[100], net$links[[100]], sep=": "))
#'
#' @importFrom utils read.table
#' @export

import.net <- function(tbl, Lowercase = 1, col.1 = 1, col.2 = 2, col.score = 3, min.score=NA, n.top=NA, echo = 1) {
if (is.null(tbl)) {stop("No network file name given...");}
if (!is.na(min.score) & col.score < 1) {stop("Edge score/weight column col.score and/or cut-off are undefined...");}
if (!is.na(min.score) & !is.na(n.top)) {stop("Edge score/weight column col.score should be used either for filetering at min.score or for selectig n.top edges.\n\tDefining both min.score and n.top is not allowed...");}
if (is.data.frame(tbl)){
  net <- tbl
} else {
  net <- read.table(tbl, row.names = NULL, header = FALSE, sep = "\t", quote = "", dec = ".", na.strings = "", skip = 0, colClasses="character", comment.char = "");
}
if (is.null(net)) {stop(paste("Not a proper network file:", tbl, "...", sep=" "));}
if (nrow(net) < 10) {stop(paste("Not a proper network file:", tbl, "...", sep=" "));}
if (length(unique(net[,col.1])) < 2 & length(unique(net[,col.2])) < 2) {stop("Multiple node (gene) IDs are not found: check parameters 'col.1', 'col.2'... ");}

if (Lowercase > 0 ) {
for (i in c(col.1, col.2)) {
net[,i] <- tolower(net[,i]);
}}
net<-net[which(net[,col.1] != net[,col.2]),];
net <- unique(net); Net <- NULL; Net$links <- NULL;
if (col.score > 0) {
if (!is.na(min.score)) {
net <- net[which(net[,col.score] >= min.score),];
} else {
if (!is.na(n.top)) {
net <- net[order(net[,col.score], decreasing=T),][1:n.top,];
}}}
t1 <- c(
tapply(net[,col.1], factor(net[,col.2]), paste),
tapply(net[,col.2], factor(net[,col.1]), paste)
);
Net$links <- tapply(t1, factor(names(t1)), function (x) unique(c(x, recursive = T)));
for (i  in names(Net$links)) {
Net$links[[i]] <- unique(Net$links[[i]]);
}
Net$Ntotal <- sum(sapply(Net$links, length))/2;
if (echo>0) {
print(paste("Network of ", Net$Ntotal, " edges between ", length(Net$links), " nodes...", sep = ""));
}
return(Net)
}
