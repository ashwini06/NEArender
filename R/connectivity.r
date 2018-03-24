#' Connectivity
#'
#' Function for plotting node degree distribution in order to test if the network is scale-free

#' @param NET Input network file.
#' @param Lowercase If node IDs should be rendered lower-case (Default:1, i.e. 'yes').
#' @param col.1 Number of column where 1st node of each edge should be found (only needed when NET is a text file rather than a list).
#' @param col.2 Number of column where 2nd node of each edge should be found (only needed when NET is a text file rather than a list, i.e. passed down to import.net).
#' @param col.score Number of column  where edge confidence score is found (only needed when NET is a text file rather than a list, i.e. passed down to import.net).
#' @param min.score Minimum confidence score for an edge to be included in the network (is alternative to n.top and is only used when NET is a text file rather than a list, i.e. passed down to import.net).
#' @param n.top Number of edges to be included in the network, top when ranked by confidence score (is alternative to min.score and only used when NET is a text file rather than a list, i.e. passed down to import.net).
#' @param echo If messages about execution progress should appear.
#' @param main title name for the plot, default: "Connectivity plot"
#' @param hex If the node degree distribution should be presented as density plot using package hexbin.
#' @seealso \code{\link{benchmark}} and \code{\link{import.net}}
#' @examples
#' file <-  system.file("extdata", "Connectivity.FC1_full", package = "NEArender")
#' \donttest{
#' connect <- connectivity(file,hex=TRUE)
#' }
#' @import hexbin
#' @importFrom graphics axis plot.new

#' @export


connectivity <- function (NET, Lowercase = 1, col.1 = 1, col.2 = 2, col.score=3, echo=1, main="Connectivity plot", min.score=NA, n.top = NA, hex=FALSE) {
if (is.list(NET)) {
net.list <- NET;
} else {
print("Importing network from text file:");
net.list <- import.net(NET, Lowercase = Lowercase, col.1 = col.1, col.2 = col.2, col.score = col.score, min.score=min.score, n.top = n.top, echo = echo);
}
if (net.list$Ntotal < 300) {
print(paste("The network has ", net.list$Ntotal, " nodes which is too little to analyze...", sep=""));
stop(call. = T);
}
c0 <- unlist(sapply(net.list$links, length))
h1 <- hist(c0, breaks=seq(from=min(c0-1, na.rm=T), to=max(c0, na.rm=T)+1, by=1), plot=F)
Br = c(0, 2**seq(from=log2(min(c0, na.rm=T)), to=log2(max(c0, na.rm=T))+1, by=1));
t0 <- table(c0)
X = log2(as.numeric(names(t0)));
Y = log2(t0);
pick = which(!is.infinite(X) & !is.infinite(Y) & !is.na(X) & !is.na(Y));
X = X[pick];
Y = Y[pick];
if (hex) {
plot.new()
plot(hexbin(X, Y, xbins = 2*(log10(net.list$Ntotal) + 1) - 1), xlab="log2(edges per node)", ylab="log2(no. of nodes)", main=main)
} else {
plot(X, Y, log="", type="b", cex=0.5, pch=19, xaxt = "n", yaxt = "n", xlab="Edges per node", ylab="No. of nodes", main=main);
intX = seq(from=round(X[1]), to=round(X[length(X)]), by=1);
intY = 2**seq(from=round(Y[1]), to=round(Y[length(Y)]), by=-1);
axis(1, at=intX, labels=2**intX, col.axis="black", las=2);
axis(2, at=intY, labels=2**intY, col.axis="black", las=2);
}
Lm <- lm(Y ~ X, na.action=na.exclude);
abline(coef = coef(Lm), col="red", lty=2,untf=T);
legend("topright", legend=paste(c("Edge cut-off =", "Ntop =", "Nedges =", "Nnodes ="), c(min.score, n.top, net.list$Ntotal, length(net.list$links)), collapse="\n"), bty = "n")
}
