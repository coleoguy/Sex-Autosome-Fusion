######
# Pfsa <- function(Da, scs){
#   if(scs=="XO"){
#     Xs <- 1
#     Y <- 0
#     Ds <- Da + 1
#     Dd <- Da + 2
#   }
#   if(scs=="XY"){
#     Xs <- 1
#     Y <- 1
#     Ds <- Da + 2
#     Dd <- Da + 2
#   }
#   if(scs=="XYY"){
#     Xs <- 1
#     Y <- 2
#     Ds <- Da + 3
#     Dd <- Da + 2
#   }
#   if(scs=="XXY"){
#     Xs <- 2
#     Y <- 1
#     Ds <- Da + 3
#     Dd <- Da + 4
#   }
#   res <- 1 - ((Da*(Da-2)+2*Xs*(2*Xs-2))/(2*Dd*(Dd-2))) -
#              ((Da*(Da-2)+max(c(Xs,Y))*(max(c(Xs,Y))-1))/(2*Ds*(Ds-2)))
#   return(res)
# }
# SimpPfsa <- function(Da, scs){
#   if(scs=="XO"){
#     Xs <- 1
#     Y <- 0
#     Ds <- Da + 1
#     Dd <- Da + 2
#   }
#   if(scs=="XY"){
#     Xs <- 1
#     Y <- 1
#     Ds <- Da + 2
#     Dd <- Da + 2
#   }
#   if(scs=="XYY"){
#     Xs <- 1
#     Y <- 2
#     Ds <- Da + 3
#     Dd <- Da + 2
#   }
#   if(scs=="XXY"){
#     Xs <- 2
#     Y <- 1
#     Ds <- Da + 3
#     Dd <- Da + 4
#   }
#   res <- ((2*Da*Xs + 2*Da*Y) / (2*Ds^2 - 4*Ds)) +
#          ((4*Da*Xs) / (2*Dd^2 - 4*Dd))
#   return(res)
# }
#####
# function for P(SA)
Pfsa <- function(Da, scs, mud = 0.5){ 
  #splitting the scs string apart
  sexchroms <- toupper(strsplit(scs, split = "")[[1]])
  # determining sex chromosome system
  if("X" %in% sexchroms){
    # counting the number of sex chromosomes
    Xs <- sum(sexchroms == "X")
    Y <- sum(sexchroms == "Y")
    # calculating other parameters
    Ds <- Da + Xs + Y
    Dd <- Da + 2 * Xs
    # performing the calculation
    res <- 1 - mud * ((Da * (Da - 2) + 4 * Xs * (Xs - 1)) / (Dd * (Dd - 2))) - 
      (1 - mud) * (((Xs * (Xs - 1))/(Ds * (Da + Xs - 1))) + 
                     ((Y * (Y - 1))/(Ds * (Da + Y - 1))) + 
                     ((Da * (Da - 2))/(Ds * (Ds - 2))))
    # determining sex chromsome system
  }else if("Z" %in% sexchroms){
    # counting sex chromosome
    Zd <- sum(sexchroms == "Z")
    W <- sum(sexchroms == "W")
    # calculating other paramters
    Dd <- Da + Zd + W
    Ds <- Da + 2 * Zd
    mus <- 1-mud
    # performing the calculation
    res <- 1 - mus * ((Da * (Da - 2) + 4 * Zd * (Zd - 1)) / (Ds * (Ds - 2))) - 
      (1 - mus) * (((Zd * (Zd - 1))/(Dd * (Da + Zd - 1))) + 
                     ((W * (W - 1))/(Dd * (Da + W - 1))) + 
                     ((Da * (Da - 2))/(Dd * (Dd - 2))))
    # determining the scs
  }else if("U" %in% sexchroms){
    # counting sex chromosomes
    U <- sum(sexchroms == "U")
    V <- sum(sexchroms == "V")
    # ensuring that the UV probability will be accurate
    if(U != V)stop('UV scs results are only accurate if the number of U chromosomes 
                   equals the number of V')
    # all diploids are heterogametic so setting mud = 0 gets rid of that area of the equation
    if(mud != 0)stop('Because the diploid phase is always heterogametic, mud must be 
                     set to 0 for UV scs.')
    # calcuating other parameters
    Ds <- Dd <- U + V + Da
    # performing the calculation
    res <- 1 - mud * ((Da * (Da - 2) + 4 * U * (U - 1)) / (Dd * (Dd - 2))) - 
      (1 - mud) * (((U * (U - 1))/(Ds * (Da + U - 1))) + 
                     ((V * (V - 1))/(Ds * (Da + V - 1))) + 
                     ((Da * (Da - 2))/(Ds * (Ds - 2))))
    # if the scs us not some sort of XY ZW or UV
  }else stop('scs must contain an X, Z or U (not case sensitive).')
  
  return(res)
}







maxnum <- 60
XO <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XO")
XY <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XY")
XYY <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XYY")
XXY <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XXY")
X5Y5 <- Pfsa(Da=seq(from=2, to=maxnum, by=2), scs="XXXXXYYYYY")
rates <- c(XO,XY,XYY,XXY, X5Y5)
types <- rep(c("XO", "XY", "XYY","XXY","X5Y5"), each=length(XO))
autosomes <- rep(seq(from=2, to=maxnum, by=2), times=5)
res <- data.frame(rates, types, autosomes)
library(ggplot2)
ggplot(res, aes(y=rates, x=autosomes)) + 
  geom_point(aes(colour=types), stat="identity", position="identity", alpha=0.5, size=3) + 
  geom_line(aes(colour=types), stat="identity", position="identity", alpha=0.5) + 
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  scale_size(range=c(1, 3)) + 
  xlab("Diploid autosome count") + 
  ylab("Proportion of SA fusions")


# maxnum <- 60
# XO2 <- SimpPfsa(Da=seq(from=2, to=maxnum, by=2), scs="XO")
# XY2 <- SimpPfsa(Da=seq(from=2, to=maxnum, by=2), scs="XY")
# XYY2 <- SimpPfsa(Da=seq(from=2, to=maxnum, by=2), scs="XYY")
# XXY2 <- SimpPfsa(Da=seq(from=2, to=maxnum, by=2), scs="XXY")
# rates <- c(XO,XY,XYY,XXY)
# types <- rep(c("XO", "XY", "XYY","XXY"), each=length(XO))
# autosomes <- rep(seq(from=2, to=maxnum, by=2), times=4)
# res <- data.frame(rates, types, autosomes)
# 
# ggplot(res, aes(y=rates, x=autosomes)) + geom_point(aes(colour=types), stat="identity", position="identity", alpha=0.5, size=3) + geom_line(aes(colour=types), stat="identity", position="identity", alpha=0.5) + theme_bw() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 3)) + xlab("Diploid autosome count") + ylab("Proportion of fusions joining autosome and gonosome")
