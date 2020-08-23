# Nathan W. Anderson and Heath Blackmon
# 23 July 2020
# coleoguy@gmail.com

# Each of the three functions below contains three arguments
# Da: diploid count of autosomes
# scs: sex chromosome sytem, string of sex chromosomes
#      'XO', 'XXO', 'XY', 'XYY', 'ZO', 'ZW', 'ZZWW', 'UV' etc. etc
# mud: proportion of fusions that originate in females
#      should be 0.5 unless other data is avaialable
#      must be set to 0 for UV sex chromosome systems

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


# function for P(AA)
Pfaa <- function(Da, scs, mud = 0.5){
  # splitting the scs string
  sexchroms <- toupper(strsplit(scs, split = "")[[1]])
  # determing scs
  if("X" %in% sexchroms){
    # counting sex chromosomes
    Xs <- sum(sexchroms == "X")
    Y <- sum(sexchroms == "Y")
    # calculating other parameters
    Ds <- Da + Xs + Y
    Dd <- Da + 2 * Xs
    # performing the calculation
    res <- mud * ((Da * (Da - 2)) / (Dd * (Dd - 2))) + (1 - mud) * ((Da * (Da - 2))/(Ds * (Ds - 2)))
  }else if("Z" %in% sexchroms){
    # counting sex chromosmes
    Zd <- sum(sexchroms == "Z")
    W <- sum(sexchroms == "W")
    # calculating parameters
    Dd <- Da + Zd + W
    Ds <- Da + 2 * Zd
    mus <- 1-mud
    # performing the calculation
    res <- mus * ((Da * (Da - 2)) / (Ds * (Ds - 2))) + (1 - mus)((Da * (Da - 2))/(Dd * (Dd - 2)))
  }else if("U" %in% sexchroms){
    # counting sex chromosome
    U <- sum(sexchroms == "U")
    V <- sum(sexchroms == "V")
    # making sure result will be accurate
    if(U != V)stop('UV scs results are only accurate if the number of U chromosomes 
                   equals the number of V')
    # dropping the homogametic part of the equation
    if(mud != 0)stop('Because the diploid phase is always heterogametic, mud must be 
                     set to 0 for UV scs.')
    # calculating other parameters
    Ds <- Dd <- U + V + Da
    # performing the calculation
    res <- mud * ((Da * (Da - 2)) / (Dd * (Dd - 2))) + (1 - mud) * ((Da * (Da - 2))/(Ds * (Ds - 2)))
  }else stop('scs must contain an X, Z or U (not case sensitive).')
  return(res)
}

# function for P(SS)
Pfss <- function(Da, scs, mud = 0.5){
  # splitting the scs string
  sexchroms <- toupper(strsplit(scs, split = "")[[1]])
  # determining scs
  if("X" %in% sexchroms){
    # counting sex chromosomes
    Xs <- sum(sexchroms == "X")
    Y <- sum(sexchroms == "Y")
    # calculating other parameters
    Ds <- Da + Xs + Y
    Dd <- Da + 2 * Xs
    # performing calculation
    res <- mud * ((4 * Xs * (Xs - 1)) / (Dd * (Dd - 2))) + 
      (1 - mud) * (((Xs * (Xs - 1))/(Ds * (Da + Xs - 1))) + ((Y * (Y - 1))/(Ds * (Da + Y - 1))))
  }else if("Z" %in% sexchroms){
    # counting sex chromosomes
    Zd <- sum(sexchroms == "Z")
    W <- sum(sexchroms == "W")
    # calculating other parameters
    Dd <- Da + Zd + W
    Ds <- Da + 2 * Zd
    mus <- 1-mud
    # performing calculation
    res <- mus * ((4 * Zd * (Zd - 1)) / (Ds * (Ds - 2))) + 
      (1 - mus) * (((Zd * (Zd - 1))/(Dd * (Da + Zd - 1))) + ((W * (W - 1))/(Dd * (Da + W - 1))))
  }else if("U" %in% sexchroms){
    # counting sex chromosomes
    U <- sum(sexchroms == "U")
    V <- sum(sexchroms == "V")
    # ensuring anser will be accurate
    if(U != V)stop('UV scs results are only accurate if the number of U chromosomes 
                   equals the number of V')
    # getting rid of calculations for homogametic dipoloid
    if(mud != 0)stop('Because the diploid phase is always heterogametic, mud must be 
                     set to 0 for UV scs.')
    # calculating other parameters
    Ds <- Dd <- U + V + Da
    # perfomring calculation
    res <- mud * ((4 * V * (V - 1)) / (Dd * (Dd - 2))) + 
      (1 - mud) * (((V * (V - 1))/(Ds * (Da + V - 1))) + ((V * (V - 1))/(Ds * (Da + V - 1))))
  }else stop('scs must contain an X, Z or U (not case sensitive).')
  return(res)
}


# calculating pval from line 194
sum = 0
for(i in 8:10){
  for(j in 0:(10-i)){
    sum <- sum + (factorial(10)/(factorial(i)*factorial(j)*factorial(10-i-j))) *
      Pfsa(22, 'XXXY', .5)^i * Pfaa(22, 'XXXY', .5)^j * Pfss(22, 'XXXY', .5)^(10-i-j)
  }
}
p.val <- sum
p.val
