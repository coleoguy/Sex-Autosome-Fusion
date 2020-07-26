# This analysis test for an excess of
# sex chromosome autosome fusions
# in the genus drosophila

library(ape)

# read in the tree and data
tree <- read.nexus("tree.nex")
dat <- read.csv("drosdata.csv", as.is=T, header = T)[,1:4]
colnames(dat) <- c("sp","hap","scs", "sim.state")
# trim the data and phylogeny to match
dat <- dat[dat$sp %in% tree$tip.label, ]
tree <- keep.tip(tree, dat$sp)

# lets reorder tip data to match tree order
new.dat <- dat
for(i in 1:nrow(dat)){
  hit <- which(dat$sp == tree$tip.label[i])
  new.dat[i, ] <- dat[hit, ]
}
dat <- new.dat
rm(new.dat)
plot(tree, cex=.3)
tiplabels(pch=16,cex=.4,col=c("red","blue")[as.factor(dat$scs)])

library(phytools)
# lets read in the transition matrix
# that describes our model of karyotype
# evolution
mat <- as.matrix(read.csv("transition.matrix.csv",header=F))
# not allowing for NeoXY->XY
# mat[mat==4]<-0
x <- matrix(0,  nrow=nrow(dat), 8)
rownames(x) <- dat$sp
for(i in 1:nrow(x)){
  x[i, dat$sim.state[i]] <- 1
}
colnames(x) <- 1:8
hists <- make.simmap(tree,
                     x=x,
                     model=mat,
                     pi=c(0,0,0,1,0,0,0,0),
                     nsim=1000)
counts <- describe.simmap(hists)$count
#colnames(counts)[c(8,15,22,36,43)]
colnames(counts)[c(9,17,25,33,41,49,57)]
#AAfusioncounts <- rowSums(counts[, c(8,15,22,36,43)])
AAfusioncounts <- rowSums(counts[, c(9,17,25,33,41,49,57)])
colnames(counts)[c(12, 20, 28)]
ASfusioncounts <- rowSums(counts[,c(12, 20, 28)])
obspropSA <- ASfusioncounts/AAfusioncounts

library(evobiR)
# diploid autosome count
# sim.state   #dipa count
#    1            4
#    2            6
#    3            8
#    4            10
#    5            4
#    6            6
#    7            7
#    8            10
expSA <- c()
for(i in 1:1000){
  times <- describe.simmap(hists[[i]])$times[2,]
  sa4 <- Pfsa(Da = 4, scs = "XY")
  sa6 <- Pfsa(Da = 6, scs = "XY")
  sa8 <- Pfsa(Da = 8, scs = "XY")
  sa10 <- Pfsa(Da = 10, scs = "XY")
  expSA[i] <- sum(sa4 * times[c(1, 5)],
                  sa6 * times[c(2, 6)],
                  sa8 * times[c(3, 7)],
                  sa10 * times[c(4, 8)])
}



plot(density(expSA, bw = .009),
     xlim = c(.1, .5), main = "",
     xlab = "Proportion sex-autosome fusion",
     cex.axis = .7, cex.lab = .7)
polygon(density(expSA, bw = .009),
        col = rgb(1, 0, 0, .3))
lines(density(obspropSA))
polygon(density(obspropSA),
        col = rgb(0, 0, 1, .3))
points(x = c(.1, .1),
       y = c(40, 35),
       pch = 15,
       col = c(rgb(1, 0, 0, .5),
               rgb(0, 0, 1, .5)))
text(x = c(.1, .1),
     y = c(40, 35),
     labels = c("Expected", "Inferred"),
     pos = 4, cex = .6)
