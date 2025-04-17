library(minet)
library(Rgraphviz)

expr = t(read.delim('../fantom5/data/lognorm_dermal.fibroblast_iPSC.txt'))
date()
res <- minet(
    dataset = expr,
    disc = "equalwidth", nbins = round(sqrt(nrow(expr)),0),
    estimator = "mi.shrink",
    method = "mrnet"
)
date()

# Plots
tab <- validate(res, syn.net)
show.pr(tab)
g1 <- as(res, 'graphNEL')
plot(g1)