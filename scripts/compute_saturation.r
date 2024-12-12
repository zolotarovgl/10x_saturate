options(warn=-1)
options(repos = c(CRAN = "https://cran.rstudio.com/"))
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel")
}
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
 
library(parallel)
library(argparse)
library(data.table)

.get_sat <- function(u, n = NULL) {
  # u - combination matrix 
  # n - number of rows to subsample the matrix to
  if (!is.null(n)) {
    i <- sample(1:nrow(u), n)
    #message('done sampling')
  } else {
    i <- 1:nrow(u)
    n <- nrow(u)
  }
  
  #message('computing with whole matrix')
  nuniq <- nrow(u[i, ][CB != '', .N, by = .(CB, UMI)])
  
  #message('done counting')
  ntot <- n
  return(data.frame(ntot = n, nuniq = nuniq, sat = round(1 - (nuniq / ntot),4)))
}

# Argument parsing
parser <- ArgumentParser(description = 'Compute saturation metrics')
parser$add_argument('matf', type = 'character', help = 'Path to the combination matrix file')
parser$add_argument('outfile', type = 'character', help = 'Path to the output file')
parser$add_argument('num_cells', type = 'integer', help = 'Number of cells')
parser$add_argument('--ncpu', type = 'integer', default = 10, help = 'Number of CPU cores to use')
parser$add_argument('--mapping_rate', type = 'numeric', help = 'Mapping rate')
parser$add_argument('--nstep', type = 'integer', default = 30, help = 'Number of plot steps')
args <- parser$parse_args()

# Main
ncpu <- args$ncpu

u <- fread(args$matf)
if(ncol(u) == 4){
	message('Found 4 columns: Assuming CB-tag-UMI-Gene structure')
	colnames(u) <- c('CB', 'tag', 'UMI', 'Gene')
}else if(ncol(u)==2){
	message('Found 2 columns: Assuming CB-UMI structure')
	colnames(u) <- c('CB', 'UMI')
}

n_input <- round(nrow(u) / args$mapping_rate)

ns <- seq(0.1, 1,length.out = as.integer(args$nstep)) * nrow(u)

o <- do.call(rbind,
  mclapply(ns, FUN = function(n) {
    .get_sat(u,n)
  }, mc.cores = ncpu)
)
o <- as.data.frame(o)
o$ntot <- as.integer(o$ntot)
o$nuniq <- as.integer(o$nuniq)
o$sat <- as.numeric(o$sat)
colnames(o)[colnames(o) == 'ntot'] <- 'nmap'
# Account for the mapping rates:
o$ninput <- round(o$nmap / args$mapping_rate)
o$mean_reads <- round(o$ninput / args$num_cells)
o <- o[, c('ninput', 'nmap', 'nuniq', 'sat', 'mean_reads')]
write.table(o, args$outfile, sep = '\t', quote = FALSE, row.names = FALSE)
message(args$outfile, ' has been written')
