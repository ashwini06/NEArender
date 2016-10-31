#' Node degree distribution and network scale-freeness
#' 
#' Auxiliary function to estimate how much a certain network satisfy the  property of scale-freeness. 

#' @param NET Either a text network file or an R list imported with \code{import.net}
#' @param Lowercase Render gene/protein IDs lower-case. 
#' @param col.1 Number of the column in the input file that corresponds to node 1 (gene/protein ID).
#' @param col.2 Number of the column in the input file that corresponds to node 2 (gene/protein ID).
#' @param echo if execution progress should be reported. 
#' @param main The plot title. 
#' @seealso \code{\link{topology2nd}}, \code{\link{benchmark}}
#' @examples
#' data(net.kegg)
#' net <- net.kegg
#' connectivity(NET=net, main="Scale-free property of KEGG network")
#' @importFrom graphics axis
#' @export


connectivity <- function (NET, Lowercase = 1, col.1 = 1, col.2 = 2, echo=1, main="Connectivity plot") {

    if (class(NET)=="list") {net.list <- NET} else {
    print("Importing network from text file:");
    net.list <- import.net(NET, Lowercase = Lowercase, col.1 = col.1, col.2 = col.2, echo = echo);
  } 

  c0 <- unlist(sapply(net.list$links, length))
  # print(length(c0));
  h1 <- hist(c0, breaks=seq(from=min(c0-1, na.rm=T), to=max(c0, na.rm=T)+1, by=1), plot=F)
  
  Br = c(0, 2**seq(from=log2(min(c0, na.rm=T)), to=log2(max(c0, na.rm=T))+1, by=1));
  t0 <- table(cut(h1$counts, breaks=Br))
  X=log2(Br[2:length(Br)] + 0)
  Y=log2(t0 + 0)
  pick = which(!is.infinite(X) & !is.infinite(Y) & !is.na(X) & !is.na(Y));
  X = X[pick];
  Y = Y[pick];
  # print(t0[length(t0)]);
  plot(X, Y, log="", type="b", cex=0.5, pch=19, xaxt = "n", yaxt = "n", xlab="Edges per node", ylab="No. of nodes", main=main);
  intX = seq(from=round(X[1]), to=round(X[length(X)]), by=1);
  intY = 2**seq(from=round(Y[1]), to=round(Y[length(Y)]), by=-1);
  axis(1, at=intX, labels=2**intX, col.axis="black", las=2);
  axis(2, at=intY, labels=2**intY, col.axis="black", las=2);
  Lm <- lm(Y ~ X, na.action=na.exclude);
  abline(coef = coef(Lm), col="red", lty=2,untf=T);
}





