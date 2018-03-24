#' Network Enrichment Analysis (NEA)
#'
#' Given the altered gene sets (AGS), functional gene sets (FGS) and the network, calculates enrichment statistics based on the number of edges (network links) that connect individual genes in each AGS-FGS pair (edges confined within AGS only or FGS only are not taken into account). Returns relevant statistics in matrices of size \code{length(FGS)} x \code{length(AGS)} (see "Value"). Each of the first three parameters can be submitted as either a text file or as a list which has been preloaded with \code{\link{import.gs}} and \code{\link{import.net}}. The latter scenario could save much time in a batch mode(see Details).
#' @param AGS Either a text file or a list of members of each AGS (see Details). Group IDs should be found in \code{ags.group.col} and gene IDs should be found in \code{ags.gene.col}. As a special option, an AGS list can also be pre-created from raw data matrices of the format "genes X samples". Such matrices can contain information on gene copy number, gene/protein expression, gene methylation etc (\code{\link{samples2ags}}) or point mutation data (\code{\link{mutations2ags}}). These two functions return AGS lists of length=ncol(Nsamples).
#' @param FGS Either a text file or a list of members of each FGS (see Details). Group IDs should be found in \code{fgs.group.col} and gene IDs should be found in \code{fgs.gene.col}. As an  alternative, this function can use  single-node FGSs from the network nodes, which should be pre-created with \code{\link{as_genes_fgs}}. In this case each network node X becomes an FGS on its own. Respective names in the FGS list equal 'X' and the content becomes c('X'). The total length of the FGS list then equals the no. of nodes in the network. If fewer nodes are needed for such single-node analysis, then one can submit a text file in which \code{fgs.gene.col} and \code{fgs.group.col} contain identical gene IDs. This option uses \code{\link{as_genes_fgs}}.
#' @param NET The global network for the analysis (see Details) .
#' @param Lowercase render node and group IDs lower-case (default:1, i.e. 'yes').
#' @param ags.gene.col number of the column containing AGS genes (only needed if AGS is submitted as a text file).
#' @param ags.group.col number of the column containing group IDs (only needed if AGS is submitted as a text file).
#' @param fgs.gene.col number of the column containing FGS genes (only needed if FGS is submitted as a text file).
#' @param fgs.group.col number of the column containing group IDs (only needed if FGS is submitted as a text file).
#' @param net.gene1.col number of the column containing first nodes of each network edge (only needed if NET is submitted as a text file).
#' @param net.gene2.col number of the column containing second nodes of each network edge (only needed if NET is submitted as a text file).
#' @param echo if messages about execution progress should appear.
#' @param Parallelize The number of CPU cores to be used for the step "Counting actual links" (the other steps are sufficiently fast). The option is not supported on Windows.
#' @param graph Plot the heat map
#' @param na.replace replace NA values. Default=0, i.e. not to replace
#' @param digitalize If the node ID strings should be converted to internal integer ID, and then back to present the results (can speed up the computation). Since this procedure also takes some time, setting \emph{digitalize=TRUE} only makes sense for large computations, with a big network, many FGS and/or AGS.
#' @param members If matrices \code{members.fgs} and \code{members.fgs} should be created in addition (time- and memory-consuming). These matrices contain lists of genes that contributed to respective AGSxFGS enrichment statistic.


#' @return An object, i.e. a list of elements \code{n.actual}, \code{n.expected}, \code{chi}, \code{z}, \code{p}, \code{q}, each of which is a matrix of size \code{length(FGS)} x \code{length(AGS)}. The two former ones contain the number of network edges between any nodes of AGS and any nodes of FGS, respectively those observed in the actual network and expected by chance. \code{chi} are the original chi-squared network enrichment statistic values. \code{z} are respective z-scores which are normally distributed under null and are thus suitable as input to regression modelling and other parametric methods. \code{p} and \code{q} are p-values and respective FDR estimates from \code{p.adjust(p, method="BH")}.

#' @details both AGS and FGS can be either \cr
#' 1) a list preloaded from a text file using e.g. GS=import.gs("text_file.groups"); names of the list entries are gene set IDs and the entries contain gene/protein IDs that belong to the respective set or \cr
#' 2) name of the file "text_file.groups" to be read from the disk using \code{\link{import.gs}}. \cr
#' The TAB-delimited file "text_file.groups" should contain pre-compiled gene sets, so that gene set IDs will be found in column \code{ags.group.col} / \code{fgs.group.col} and gene IDs will be found in \code{ags.gene.col} / \code{fgs.gene.col}. Option (1) is much more efficient than (2) when \code{\link{nea.render}} has to be run multiple times. Similarly to AGS and FGS, NET could be submitted as either a list, pre-loaded with \code{\link{import.net}}, or a text TAB-delimited file where two columns \code{net.gene1.col} and \code{net.gene2.col} represent nodes of respective edges.

#' @seealso \code{\link{samples2ags}}
#' @references \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1534-y}
#'
#' @keywords AGS
#' @keywords FGS
#' @examples
#' ags.list <- samples2ags(fantom5.43samples, Ntop=20, method="topnorm")
#' data(can.sig.go)
#' fpath <- can.sig.go
#' fgs.list <- import.gs(fpath)
#' summary(fgs.list)
#' data(net.kegg)
#' netpath <- net.kegg
#' net <- import.net(netpath)
#' n1 <- nea.render(AGS=ags.list[1:10], FGS=fgs.list[1:10], NET=net, graph=FALSE)
#' hist(n1$chi, breaks=100)
#' hist(n1$z, breaks=100)
#' hist(n1$p, breaks=100)
#' hist(n1$q, breaks=100)
#' @export

nea.render <- function (AGS, FGS, NET,  Lowercase = 1, ags.gene.col = 2, ags.group.col = 3, fgs.gene.col = 2, fgs.group.col = 3, net.gene1.col = 1, net.gene2.col = 2, echo = 1, graph = FALSE, na.replace = 0, members=FALSE, digitalize = TRUE, Parallelize=1) {

if (!is.na(na.replace) & !is.numeric(na.replace)) {stop("Parameter 'na.replace' should contain a numeric value or NA...");}
if (echo>0) {print("Preparing input datasets:");}
if (is.list(NET)) {net.list <- NET} else {net.list <- import.net(NET, Lowercase=Lowercase, net.gene1.col, net.gene2.col)}
if (echo>0) {print(paste("Network: ",  length(net.list$links), " genes/proteins.", sep = ""));}
if (is.list(FGS)) {
fgs.list <- FGS;
} else {
if (FGS == "nw_genes") {
fgs.list <- as_genes_fgs(net.list);
if (echo>0) {print(paste("FGS: ",  length(unique(unlist(fgs.list))), " genes in as FGS groups.", sep = ""));}
} else {
fgs.list <- import.gs(FGS, Lowercase=Lowercase, fgs.gene.col, fgs.group.col, gs.type = 'f');
}}
if (echo>0) {
print(paste("FGS: ", length(unique(unlist(fgs.list))), " genes in ", length(fgs.list), " groups.", sep = ""));
}


if (is.list(AGS)) {
ags.list <- AGS;
} else {
ags.list <- import.gs(AGS, Lowercase=Lowercase, ags.gene.col, ags.group.col, gs.type = 'a');
}
if (echo>0) {
print(paste("AGS: ",  length(unique(unlist(ags.list))), " genes in ", length(ags.list), " groups...", sep = ""));
print("Calculating N links expected by chance...");
}

stats <- NULL;

if (members) {
Members.in <- function(i, gs1, gs2) {
return(sapply(gs1,
function (x) {
aa = unique(unlist(net.list$links[x]));
return(paste(gs2[[i]][which(gs2[[i]] %in% aa)], collapse=", "));
}));
}

if (echo>0) {print("Enumerating connected nodes:");} #CAN BE OPTIMIZED FURTHER FOR FGS="nw_genes"...
for (gg in c("ags", "fgs")) {
print(paste(toupper(gg), "...", sep=""));
if (gg == "fgs") {
list1 <- fgs.list; list2 <- ags.list;
} else {
list1 <- ags.list; list2 <- fgs.list;
}
if (Parallelize > 1)  {
print(system.time(l0 <- mclapply(names(list1), Members.in, list2, list1, mc.cores = Parallelize)));
} else {
print(system.time(l0 <- lapply(names(list1), Members.in, list2, list1)));
}
ele <- paste("members", gg, sep=".");
stats[[ele]] <- matrix(unlist(l0), nrow = length(list1), ncol = length(list2), dimnames = list(names(list1), names(list2)), byrow=T);
}
}

if (digitalize) {
print("Rendering integer IDs...");
print(system.time(mapped <- char2int.fast(net.list, ags.list, fgs.list, Parallelize=Parallelize)));
net.list$links <- mapped$net;
ags.list <- mapped$gs[["a"]];
fgs.list <- mapped$gs[["b"]];
}
net.list$cumulativeDegrees$AGS <- sapply(ags.list, FUN=function(x) length(unlist(net.list$links[unlist(x)])));
net.list$cumulativeDegrees$FGS <- sapply(fgs.list, FUN=function(x) length(unlist(net.list$links[unlist(x)])));
N.ags_fgs.expected <- outer(net.list$cumulativeDegrees$FGS, net.list$cumulativeDegrees$AGS , "*") / (2 * net.list$Ntotal);

Nin <- function(i, ags) {
return(lapply(ags, function (x) length(which(unlist(net.list$links[unlist(x)]) %in% unlist(unlist(fgs.list[[i]]))))));
}



if (echo>0) {print("Counting actual links...");}
if (Parallelize > 1)  {
print(system.time(l1 <- mclapply(names(fgs.list), Nin, ags.list, mc.cores = Parallelize)));
} else {
print(system.time(l1 <- lapply(names(fgs.list), Nin, ags.list)));
}

N.ags_fgs.actual <- matrix(unlist(l1), nrow = length(fgs.list), ncol = length(ags.list), dimnames = list(names(fgs.list), names(ags.list)), byrow=T);
N.ags_fgs.actual[which(N.ags_fgs.actual == 0)] = na.replace;
if (echo>0) {print("Calculating statistics...");}
chi = (
((N.ags_fgs.actual - N.ags_fgs.expected) ** 2) / N.ags_fgs.expected +
(((net.list$Ntotal - N.ags_fgs.actual) - (net.list$Ntotal - N.ags_fgs.expected)) ** 2) /
(net.list$Ntotal - N.ags_fgs.expected));
#############################################################
#############################################################
#############################################################
p.chi <- pchisq(chi, df=1, log.p = FALSE, lower.tail = FALSE);# - log(2)
Z <- qnorm(p.chi/2, lower.tail = FALSE)

Depleted <- which(sign(N.ags_fgs.actual  - N.ags_fgs.expected) < 0);
Z[Depleted] = -1 * abs(Z[Depleted]); #all depletion cases should produce negative Z values!
Z[which(is.nan(Z))] <- NA;
stats$cumulativeDegrees <- net.list$cumulativeDegrees;
stats$n.actual <- N.ags_fgs.actual;
stats$n.expected <- N.ags_fgs.expected;
stats$chi <- chi;
va = 'z'; stats[[va]] <- Z;
Min = min(stats[[va]][which(!is.infinite(stats[[va]]))], na.rm=T) - 1;
Max = max(stats[[va]][which(!is.infinite(stats[[va]]))], na.rm=T) + 1;
stats[[va]][which(is.infinite(stats[[va]]) & (stats[[va]] < 0))] = Min;
stats[[va]][which(is.infinite(stats[[va]]) & (stats[[va]] > 0))] = Max;
# stats$p <- 10 ^ p.chi;
stats$p <- p.chi;
stats$q <- matrix(
p.adjust(stats$p, method="BH"),
nrow = nrow(stats$p), ncol = ncol(stats$p),
dimnames = list(rownames(stats$p), colnames(stats$p)),
byrow=FALSE);
stats$q[which(stats$z < 0)] = 1;
if (graph) {set.heat(ags.list,fgs.list,stats$z)}
if (echo>0) {print("Done.");}
return(stats);
}
