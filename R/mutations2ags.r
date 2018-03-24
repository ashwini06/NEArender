#' Create AGS from a mutation matrix
#'
#' Imports a TAB-delimited file with mutations.
#' This function creates a new list of AGSs from a table listing point (or otherwise qualitatively defined) mutations. Such a matrix M typically has size Ngenes x Nsamples, so that the current function returns a list of \code{length=ncol(M)}. For each of the Nsamples, AGSs are created as simple lists of all mutated genes G in a given sample S, i.e. any value X in the matrix M that satisfies condition \code{!is.na(X)} would be treated as a mutation. Eventual mutation types / categories are ignored. Wild type gene states in the original TAB-delimited file should be represented with NAs.

#' @param MUT Matrix of size \emph{Ngenes x Nsamples} (the both Ns are positive integers, depending on the screen scale).
#' @param col.mask To include only columns with IDs that contain the specified mask. This parameter is aware of regular expression syntax, i.e. uses \code{grep(..., fixed = FALSE)}.
#' @param namesFromColumn Number of the column (if any) that contains the gene/protein names. Note that it is only necessary if the latter are NOT the unique rownames of the matrix. This could be sometimes useful for processing redundant gene  profiles with one-to-many mapping etc. Otherwise (i.e. the default), rownames shall contain gene IDs.
#' @param permute If the list of AGSs should be created via random permutation of sample labels. This might be needed for testing the null hypothesis that mutated genes are randomly combined into individual genomes, while having the same frequency distribution as in the actual cohort. Since reproducing the original distribution of AGS sizes is a non-trivial set theoretical problem, the procedure is accompanied by plotting gene set sizes in actual vs. permuted AGS (the latter is usually smaller, which might be unavoidable without a sophisticated algortihm...).

#'
#' @examples
#' data("tcga.gbm",package="NEArender")
#' dim(tcga.gbm)
#' ags.list <- mutations2ags(tcga.gbm, col.mask="[-.]01$")
#' length(ags.list)
#' length(unique(unlist(ags.list)))
#' @export


mutations2ags  <- function(MUT, col.mask=NA, namesFromColumn=NA, permute=FALSE
# , Lowercase = 1
) {
if (is.null(MUT)) {stop("Not enough parameters...");}
mgs.list <- NULL;

if (is.na(namesFromColumn)) {m1 <- MUT;}
else {m1 <- MUT[,(namesFromColumn+1):ncol(MUT)];}
if (!is.na(col.mask)) {m1 <- m1[,colnames(m1)[grep(col.mask,colnames(m1))]];}
mgs.list <- apply(m1, 2, function (x) unique(tolower(names(x))[which(!is.na(x) )]));

if (permute) {mgs.list <- permute.gs(mgs.list);}

return(mgs.list);
}

permute.gs <- function (GS, Plot=FALSE) {
pmgs <- as.list(NULL); mmgs <- unlist(GS);
fmgs <- table(mmgs) / length(GS); # fmgs <- fmgs / sum(fmgs);
for (m in names(GS)) {
pmgs[[m]] <- sample(x = names(fmgs), size = length(GS[[m]]), replace = FALSE, prob = fmgs);
}
print("Gene set permutation done.");
if (Plot) {
plot(table(unlist(GS))[names(fmgs)], table(unlist(pmgs))[names(fmgs)], xlab="Original", ylab="Permuted", main="#Samples / gene");
abline(0,1,lty=2, col="grey");
plot(sapply(GS, length)[names(GS)], sapply(pmgs, length)[names(GS)], xlab="Original", ylab="Permuted", main="#Genes / sample")
abline(0,1,lty=2, col="grey");
}
return(pmgs);
}
