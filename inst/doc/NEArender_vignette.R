## ----include=F, echo=F---------------------------------------------------
packages <- c("NEArender","data.table")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos = "http://cran.us.r-project.org")  
}

## ----include=T-----------------------------------------------------------
library(NEArender)
library(data.table)


## ------------------------------------------------------------------------
#input <- fread("http://research.scilifelab.se/andrej_alexeyenko/downloads/test_data/FANTOM5.43samples.txt",sep="\t",header=T,stringsAsFactors=FALSE,data.table = F)
##  Converting genenames as rownames
# rownames(input) <-input[,1]
# input <- as.matrix(input[,c(2:ncol(input))])

data("fantom5.43samples")
input <- fantom5.43samples
dim(input)
ags.list1 <- samples2ags(input, Ntop=20, method="topnorm")

## ------------------------------------------------------------------------
data("tcga.gbm",package="NEArender")
ags.list3 <- mutations2ags(tcga.gbm, col.mask="[-.]01$")

## ------------------------------------------------------------------------
ags.list2 <-import.gs(
"http://research.scilifelab.se/andrej_alexeyenko/downloads/test_data/cluster2_Downregulated_ags.txt", 
Lowercase = 1, col.gene = 1,col.set = 3, gs.type = 'ags')

## ------------------------------------------------------------------------
fgs.list <- import.gs(
"http://research.scilifelab.se/andrej_alexeyenko/downloads/test_data/can.sig.go.txt",
Lowercase = 1, col.gene = 2, col.set = 3, gs.type = 'fgs')

## ------------------------------------------------------------------------
data(can.sig.go)
fgs.list <- import.gs(can.sig.go)

## ------------------------------------------------------------------------
data(net.kegg)
net <- import.net(net.kegg)
print(paste(names(net$links)[10], net$links[[10]], sep=": "))

## ------------------------------------------------------------------------
net <- import.net("http://research.scilifelab.se/andrej_alexeyenko/downloads/test_data/net.kegg.txt")

## ------------------------------------------------------------------------
net.merged<-"http://research.scilifelab.se/andrej_alexeyenko/downloads/test_data/merged6_and_wir1_HC2"
net <- import.net(net.merged)

## ------------------------------------------------------------------------
data(net.kegg)
net <- import.net(net.kegg);
fgs.genes <- as_genes_fgs(net);
#save_gs_list(fgs.genes, File = "~/single_gene_ags.groups.tsv");

## ----include=TRUE--------------------------------------------------------
n1 <- nea.render(AGS=ags.list1, FGS=fgs.list, NET=net)

## ----fig.width=4, fig.height=4, fig.cap=c("n1$chi - chi-square estimate","n1$z- zscores", "NEA- pvalues","NEA-qvalues"),fig.show='asis'----
hist(n1$chi, breaks=100)    
hist(n1$z, breaks=100) 
hist(n1$p, breaks=100)  
hist(n1$q, breaks=100)  

## ----include=T, fig.cap=c("g1$estimate - an estimate of the odds ratio", "g1$n - number of ags genes that belongs to corresponding fgs","g1$p - the p-value of the fishers.exact test","g1$q - Adjusted p-values by \"BH-method\"")----
ags.list2 <- samples2ags(fantom5.43samples, Ntop=1000, method="topnorm")
g1 <- gsea.render(AGS=ags.list2, FGS=fgs.list, Lowercase = 1, 
ags.gene.col = 2, ags.group.col = 3, fgs.gene.col = 2, fgs.group.col = 3, 
echo=1, Ntotal = 20000, Parallelize=1)

hist(log(g1$estimate), breaks=100)  
hist(g1$n, breaks=100)  
hist(g1$p, breaks=100)  
hist(g1$q, breaks=100)  

## ----include=T, fig.cap=c("ROC curve evaluating KEGG network (net.kegg) for specific member term - \"kegg_04270_vascular_smooth_muscle_contraction\"")----
b0 <- benchmark (NET = net,
 GS = fgs.list[c("kegg_04270_vascular_smooth_muscle_contraction")], 
 echo=1, graph=TRUE, na.replace = 0, mask = ".", minN = 0,
 coff.z = 1.965, coff.fdr = 0.1, Parallelize=1);

## ----eval=FALSE----------------------------------------------------------
#  b1 <- NULL;
#  for (mask in c("kegg_", "go_")) {
#  b1[[mask]] <- NULL;
#  ref_list <- list(net.kegg=net.kegg,net.merged=net.merged)
#  for (file.net in c("net.kegg","net.merged")) {
#  # a series of networks can be put here: c("net.kegg1", "net.kegg2", "net.kegg3") in ref_list
#  net <- import.net(ref_list[[file.net]], col.1 = 1, col.2 = 2, Lowercase = 1, echo = 1)
#  b1[[mask]][[file.net]] <- benchmark (NET = net, GS = fgs.list, echo=1,
#  graph=FALSE, na.replace = 0, mask = mask, minN = 0,  Parallelize=1);
#  }}
#  

## ----eval=FALSE----------------------------------------------------------
#  roc(b1[["kegg_"]], coff.z = 2.57, coff.fdr = 0.01,main="kegg_");
#  roc(b1[["go_"]], coff.z = 2.57, coff.fdr = 0.01,main="go_");

## ----include=FALSE-------------------------------------------------------

library(RColorBrewer)
library(MASS)
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


## ----fig.cap=c("Node degree distribution of `net.merged`"),fig.show='asis'----
connectivity(NET="http://research.scilifelab.se/andrej_alexeyenko/downloads/test_data/merged6_and_wir1_HC2",
             Lowercase = 1, col.1 = 1, col.2 = 2, echo=1, main="Higher order topology")

## ----include=FALSE-------------------------------------------------------

topology2nd <- function (NET, Lowercase = 1, col.1 = 1, col.2 = 2, echo=1, main="Higher order topology") {
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
cr <- colorRampPalette( c('gray', "red"));

if (class(NET)=="list") {net.list <- NET} else {
  print("Importing network from text file:");
  net.list <- import.net(NET, Lowercase = Lowercase, col.1 = col.1, col.2 = col.2, echo = echo);
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
image(k, col=r, main=main, cex.main=1, xlab="First node degree", ylab="Second node degree", log="xy",xaxt = "n", yaxt = "n"); 
Shift = 0;
intX = seq(from=round(min(X, na.rm=T)), to=round(max(X, na.rm=T)), by=1);
intY = seq(from=round(min(Y, na.rm=T)), to=round(max(Y, na.rm=T)), by=1);
axis(1, at=intX+Shift, labels=10**intX, col.axis="black", las=2);
axis(2, at=intY+Shift, labels=10**intY, col.axis="black", las=2);
}



## ----include=T, fig.cap=c("Second order topology `net.merged`"), fig.show='asis'----
topology2nd(NET="http://research.scilifelab.se/andrej_alexeyenko/downloads/test_data/merged6_and_wir1_HC2", 
            Lowercase = 1, col.1 = 1, col.2 = 2, echo=1, main="Higher order topology")

