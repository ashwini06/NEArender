#' Higher order topology and correlation between node degrees 
#' 
#' Auxiliary function to estimate how non-random are associations of node degrees in network edges. 

#' @param NET Either a text network file or an R list imported with \code{import.net}
#' @param Lowercase Render gene/protein IDs lower-case. 
#' @param col.1 Number of the column in the input file that corresponds to node 1 (gene/protein ID).
#' @param col.2 Number of the column in the input file that corresponds to node 2 (gene/protein ID).
#' @param echo if execution progress should be reported. 
#' @param main The plot title. 
#' @seealso \code{\link{connectivity}}, \code{\link{benchmark}}
#' @examples
#' file <-  system.file("extdata", "Connectivity.FC1_full", package = "NEArender")
#' \donttest{
#'  topology2nd(NET=file)
#' }
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom MASS kde2d
#' @importFrom graphics axis
#' @export

topology2nd <- function (NET, Lowercase = 1, col.1 = 1, col.2 = 2, echo=1, main="Higher order topology") {
  
  Rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- Rf(32)
  if (is.list(NET)) {net.list <- NET} else {
    print("Importing network from text file:");
    net.list <- import.net(NET, Lowercase = Lowercase, col.1 = 1, col.2 = 2, echo = echo);
  } 
  c0 <- unlist(sapply(net.list$links, length));
  
  inew=1; 
  n2 <- matrix(NA, nrow=2 * net.list$Ntotal, ncol=2); 
  for (n1 in names(net.list$links)) {
    nnew = unlist(net.list$links[[n1]]); 
    n2[inew:(inew+length(nnew)-1),] <- cbind(n1, nnew); 
    inew=inew+length(nnew);
  }
  X = log10(c0[n2[,1]]+1);
  Y = log10(c0[n2[,2]]+1);
  k<-kde2d(X, Y, n=20);
  image(k, col=r, main=main, cex.main=1, xlab="First node degree", ylab="Second node degree", log="xy", , xaxt = "n", yaxt = "n"); 
  Shift = 0;
  intX = seq(from=round(min(X, na.rm=T)), to=round(max(X, na.rm=T)), by=1);
  intY = seq(from=round(min(Y, na.rm=T)), to=round(max(Y, na.rm=T)), by=1);
  axis(1, at=intX+Shift, labels=10**intX, col.axis="black", las=2);
  axis(2, at=intY+Shift, labels=10**intY, col.axis="black", las=2);
}