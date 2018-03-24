#' Benchmark networks using Network Enrichment Analysis (NEA)
#' 
#' Tests the ability of a given network to perform well in a network enrichment analysis. It executes a series of multiple individual tests: for each member gene of a pathway or another functional set calculates the network enrichment score against other members of the same gene set. This procedure gives true positive and false negative test results. In order to complement it with false positives and true negatives, the same is done for randomly picked genes (with matching node connectivity values) against the same functional sets. The two vectors allow plotting a ROC curve where at each sequential cutoff represents a ratio of true positive vs. false positive predictions. This approach (first presented in  \href{http://www.biomedcentral.com/1471-2105/15/308}{Merid et al. (2012)}) is an alternative to the trivial counting edges shared between different networks  and is superior to the latter because: 1) the analysis can be done without knowing the "true" reference network, 2) benchmarks can be context-dependent by using domain-specific test sets (e.g. cancer, diabetes etc.), 3) one can compare more than two networks at a time, and 4) given dense global networks and due to the use of multi-gene sets, presence or absence of particular links is unlikely to affect the overall result. 


#' @param NET A network to benchmark. See Details in \code{\link{nea.render}}. 
#' @param GS a test set, typically a set of pathways with known members. 
#' @param gs.gene.col number of the column containing GS genes (only needed if GS is submitted as a text file)
#' @param gs.group.col number of the column containing group IDs (only needed if GS is submitted as a text file)
#' @param net.gene1.col number of the column containing first nodes of each network edge (only needed if NET is submitted as a text file)
#' @param net.gene2.col number of the column containing second nodes of each network edge (only needed if NET is submitted as a text file)
#' @param mask when the test set contains various GSs, they can be used selectively by applying a mask. The mask follows the regular expression synthax, since \code{fixed=FALSE} in \code{\link{grep}}.
#' @param minN the minimal number of network edges that must connect a tested member with the GS genes for the test to be considered positive. (Default:0).
#' @param coff.z a parameter to \code{\link{roc}}.
#' @param coff.fdr to make significance levels comparable between different curves, the point where FDR=\code{coff.fdr} will be labeled with a circle (think of TP/FP ratio at this level).
#' @param echo if messages about execution progress should appear
#' @param Parallelize The number of CPU cores to be used for the step "Counting actual links" (while the other steps are sufficiently fast). The option is not supported in Windows.
#' @param graph Plot the ROC curve immediately. Alternatively, the returned list is plotted afterwards by \code{\link{roc}}. In the latter case, it could be a combined list of lists for multiple test sets and networks which are then plotted as separate curves (see Examples).
#' @param na.replace replace NA values. Default=0, i.e. do not replace.
#' 
#' @return An object, i.e. a list of three equal-length vectors from a \code{\link{prediction}} object of ROCR package (\code{prediction@@cutoffs}, \code{prediction@@fp}, \code{prediction@@tp}) and the point that matches \code{coff.fdr}. These are needed to plot a ROC curve for the given network and test set by using \code{\link{roc}}.
#' 
#' @details The function would either plot a ROC curve for the analyzed network, or return an object with the following slots from function prediction (package ROCR):\cr
#'  tp, vector of true positives;\cr
#'  fp, vector of false positives;\cr
#'  tn, vector of true negatives;\cr
#'  fn, vector of false negatives;\cr
#'  cutoffs, z-score cutoffs from \code{nea.render};\cr
#'  cross.z, a z-score value which corresponds to FDR=\code{coff.fdr} (will be denoted with a special marker at the curve);
#'
#' @seealso \code{\link{roc}}, \code{\link{nea.render}}
#' @references \url{http://www.biomedcentral.com/1471-2105/15/308}
#' @import parallel ROCR stats
#' @keywords ROC

#' @examples
#' data(can.sig.go);
#' fpath <- can.sig.go
#' gs.list <- import.gs(fpath, Lowercase = 1, col.gene = 2, col.set = 3);
#' data(net.kegg)
#' netpath <- net.kegg
#' net <- import.net(netpath)
#' \donttest{
#' b0 <- benchmark (NET = net,
#'  GS = gs.list, 
#'  echo=1, graph=TRUE, na.replace = 0, mask = ".", minN = 0,
#'  coff.z = 1.965, coff.fdr = 0.1, Parallelize=2);
#' }
#' \dontrun{
#' ## Benchmark a number of networks on GO terms and KEGG pathways separately, using masks:
#' b1 <- NULL;
#' for (mask in c("kegg_", "go_")) {
#' b1[[mask]] <- NULL;
#' for (file.net in c("netpath")) {
#' # a series of networks can be put here: c("netpath1", "netpath2", "netpath3")
#' net <- import.net(netpath, col.1 = 1, col.2 = 2, Lowercase = 1, echo = 1)
#' b1[[mask]][[file.net]] <- benchmark (NET = net, GS = gs.list, echo=1, 
#' graph=FALSE, na.replace = 0, mask = mask, minN = 0,  Parallelize=1);
#' }}
#' par(mfrow=c(2,1));
#' roc(b1[["kegg_"]], coff.z = 2.57,main="kegg_");
#' roc(b1[["go_"]], coff.z = 2.57,main="go_");
#' }
#' @export


benchmark <- function(NET, GS, gs.gene.col = 2, gs.group.col = 3, net.gene1.col = 1, net.gene2.col = 2, echo=1, graph=FALSE, 
na.replace = 0, mask = '.', minN = 0, coff.z = 1.965, coff.fdr = 0.1, 
Parallelize=1) {
  if (is.list(NET)) {net.list <- NET} else {net.list <- import.net(NET, Lowercase=1, col.1 = net.gene1.col, col.2 = net.gene2.col);}
  fgs.list <- as_genes_fgs(net.list);
  if (is.list(GS)) {
    ags.list <- GS;
  } else {
    ags.list <- import.gs (GS, Lowercase=1, col.gene = gs.gene.col, col.set = gs.group.col, gs.type = 'a');
  }
  n1  <- nea.render(AGS=ags.list, FGS=fgs.list, NET = net.list, Lowercase = 1, echo=echo, graph=FALSE, na.replace = na.replace, Parallelize=Parallelize);
  gc();
  l1 <- lapply(net.list$links, function (x) round(log2(length(x))));
  l0 <- unlist(l1);
  connectivity <- NULL;
  for (bin in unique(l0)) {
    All <- names(l0)[which(l0 == bin)]; 
    connectivity[[as.character(bin)]] <- All # [which(All %in% pw.members)]; 
  }
  zreal = znull =
    as.vector(rep(NA, times=length(which(unlist(ags.list) %in% rownames(n1$z)))))
  Agss = names(ags.list)[grep(mask, names(ags.list), fixed=FALSE)];
  i = 0; j = 0;
  print("Matching each positive test node with a negative one, regarding node connectivity: ");
  for (ags in Agss) {
    mem <- ags.list[[ags]][which(ags.list[[ags]] %in% rownames(n1$z))]
    
    if (length(mem) > 0) {
      vec <- n1$z[mem, ags];
      zreal[(j+1):(j+length(vec))] <- vec;
      j = j + length(vec);
      for (ge in mem) {
        i = i + 1;
        if (grepl("00$", as.character(i), fixed=F)) {print(i);}
        bin <- as.character(round(log2(length(net.list$links[[ge]]))))
        znull[i] <- n1$z[sample(connectivity[[bin]],1), ags]; 
      }}}
  n1$q[which(n1$z < 0)] = 1;
  cross.z = n1$z[which(abs(n1$q - coff.fdr) == min(abs(n1$q - coff.fdr), na.rm=TRUE))][1];
  p0 <- prediction(c(znull, zreal), c(rep(0, times=length(znull)), rep(1, times=length(zreal))));
  pred<-list(tp=p0@tp, fp=p0@fp, tn=p0@tn, fn=p0@fn, cutoffs=p0@cutoffs, cross.z=cross.z, nv=net.list$Ntotal, ne=length(net.list$links))
  if (graph) {
    roc(pred, coff.z = coff.z)#, coff.fdr = coff.fdr);
  }
  return(pred);
}