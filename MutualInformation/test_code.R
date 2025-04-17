library(minet)
library(Rgraphviz)
# 1. data discretization (DISC)
# 2. MIM computation (MIM)
# 3. network inference (NET)
# 4. normalization of the network (NORM)

# Example
data(syn.data)
expr = syn.data
res <- minet(
    dataset = expr,
    disc = "equalwidth", nbins = sqrt(nrow(expr)),
    estimator = "mi.shrink",
    method = "mrnet"
)

# The above is a compact way to execute the following sequence of runs:
# discretize with the binning algorithm (i.e. equal frequency or equal size interval) with number of bins 10
discdata <- discretize(expr, "equalwidth", 10)
# entropy estimator used for the computation of mutual information (empirical, Miller-Madow, shrink, Schurmannn-Grassberger)
mim <- build.mim(discdata, "mi.shrink")
# inference algorithm (ARACNE, CLR, MRNET etc)
net <- mrnet(mim)
# normalizes all the weights of the inferred adjancy matrix between 0 and 1.
res <- norm(net)

# Plots
tab <- validate(res, syn.net)
show.pr(tab)
g1 <- as(res, 'graphNEL')
plot(g1)
