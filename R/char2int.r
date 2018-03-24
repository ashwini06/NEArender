#' char2int function
#' 
#' This function reads input list objects and maps unique numbers to each member in the network/ genesets. 
#' @seealso \code{\link{nea.render}}
#'
#' 
#' @param net.list  the global network object, pre-created using \code{\link{import.net}}
#' @param gs.list.1 AGS or FGS object that lists members of each individual AGS/FGS, pre-created using \code{\link{samples2ags}}
#' @param gs.list.2 FGS or AGS object that lists members of each individual FGS/AGS, pre-created using \code{\link{samples2ags}}
#' @param Parallelize The number of CPU cores to be used for the step "Counting actual network links" (the other steps are sufficiently fast). The option is not supported on Windows.
#' @keywords internal



char2int.fast  <- function (net.list, gs.list.1, gs.list.2 = NULL, Parallelize=1) {
all.names <- unique(c(names(net.list$links), unlist(net.list$links), unlist(gs.list.1)));
map.names <- 1:length(all.names);

Mapp <- function (x) {map.names[unlist(net.list$links[[x]])]}

names(map.names) <- all.names;
mapped <- NULL; mapped$net <- NULL; mapped$gs <- NULL;
# for (n in names(net.list$links)) {mapped$net[[map.names[n]]] <- as.list(map.names[net.list$links[[n]]]);}
if (Parallelize > 1)  {
mapped$net <- mclapply(names(net.list$links), Mapp, mc.cores = Parallelize); 
} else {
mapped$net <- lapply(net.list$links, function (x) {map.names[unlist(x)]});  
}
# names(mapped$net) <- map.names[names(net.list$links)]
# return(mapped);
for (i in c("a", "b")) {
if (i == "a") {gsl = gs.list.1;}
if (i == "b") { 
if (!is.null(gs.list.2)) {
gsl = gs.list.2;
} else {
gsl = NULL;
}}
if (!is.null(gsl)) {
for (n in names(gsl)) {
mapped$gs[[i]][[n]] <- as.list(map.names[gsl[[n]]]);
}
}
}
return(mapped);
}