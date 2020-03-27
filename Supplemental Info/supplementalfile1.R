# function for P(SA)
Pfsa <- function(Da, scs, mud){
  if(scs=="XO" | scs=="XY" | scs=="XYY" | scs=="XXY"){
    if(scs=="XO"){
      Xs <- 1
      Y <- 0
      Ds <- Da + 1
      Dd <- Da + 2
    }
    if(scs=="XY"){
      Xs <- 1
      Y <- 1
      Ds <- Da + 2
      Dd <- Da + 2
    }
    if(scs=="XYY"){
      Xs <- 1
      Y <- 2
      Ds <- Da + 3
      Dd <- Da + 2
    }
    if(scs=="XXY"){
      Xs <- 2
      Y <- 1
      Ds <- Da + 3
      Dd <- Da + 4
    }
    res <- 1 - mud*((Da*(Da-2)+2*Xs*(2*Xs-2))/(Dd*(Dd-2))) - 
      (1 - mud)*((Da*(Da-2)+max(c(Xs,Y))*(max(c(Xs,Y))-1))/(Ds*(Ds-2)))
    return(res)
  }
  if(scs=="ZO" | scs=="ZW" | scs=="ZWW" | scs=="ZZW"){
    if(scs=="ZO"){
      Zd <- 1
      W <- 0
      Dd <- Da + 1
      Ds <- Da + 2
    }
    if(scs=="Zw"){
      Zd <- 1
      W <- 1
      Dd <- Da + 2
      Ds <- Da + 2
    }
    if(scs=="ZWW"){
      Zd <- 1
      W <- 2
      Dd <- Da + 3
      Ds <- Da + 2
    }
    if(scs=="ZZW"){
      Zd <- 2
      W <- 1
      Dd <- Da + 3
      Ds <- Da + 4
    }
    mus <- 1-mud
    res <- 1 - mus*((Da*(Da-2)+2*Zd*(2*Zd-2))/(2*Ds*(Ds-2))) - 
      (1 - mus)*((Da*(Da-2)+max(c(Zd,W))*(max(c(Zd,W))-1))/(2*Dd*(Dd-2)))
    return(res)
  }
}
# function for P(AA)
Pfaa <- function(Da, scs, mud){
  if(scs=="XO" | scs=="XY" | scs=="XYY" | scs=="XXY"){
    if(scs=="XO"){
      Xs <- 1
      Y <- 0
      Ds <- Da + 1
      Dd <- Da + 2
    }
    if(scs=="XY"){
      Xs <- 1
      Y <- 1
      Ds <- Da + 2
      Dd <- Da + 2
    }
    if(scs=="XYY"){
      Xs <- 1
      Y <- 2
      Ds <- Da + 3
      Dd <- Da + 2
    }
    if(scs=="XXY"){
      Xs <- 2
      Y <- 1
      Ds <- Da + 3
      Dd <- Da + 4
    }
    res <- mud * ((Da*(Da-2))/(Dd*(Dd-2))) + (1-mud)* ((Da*(Da-2))/(Ds*(Ds-2)))
    return(res)
  }
  if(scs=="ZO" | scs=="ZW" | scs=="ZWW" | scs=="ZZW"){
    if(scs=="ZO"){
      Zd <- 1
      W <- 0
      Dd <- Da + 1
      Ds <- Da + 2
    }
    if(scs=="Zw"){
      Zd <- 1
      W <- 1
      Dd <- Da + 2
      Ds <- Da + 2
    }
    if(scs=="ZWW"){
      Zd <- 1
      W <- 2
      Dd <- Da + 3
      Ds <- Da + 2
    }
    if(scs=="ZZW"){
      Zd <- 2
      W <- 1
      Dd <- Da + 3
      Ds <- Da + 4
    }
    mus <- 1-mud
    res <- mus * ((Da*(Da-2))/(Ds*(Ds-2))) + (1-mus)* ((Da*(Da-2))/(Dd*(Dd-2)))
    return(res)
  }
}
# function for P(SS)
Pfss <- function(Da, scs, mud){
  if(scs=="XO" | scs=="XY" | scs=="XYY" | scs=="XXY"){
    if(scs=="XO"){
      Xs <- 1
      Y <- 0
      Ds <- Da + 1
      Dd <- Da + 2
    }
    if(scs=="XY"){
      Xs <- 1
      Y <- 1
      Ds <- Da + 2
      Dd <- Da + 2
    }
    if(scs=="XYY"){
      Xs <- 1
      Y <- 2
      Ds <- Da + 3
      Dd <- Da + 2
    }
    if(scs=="XXY"){
      Xs <- 2
      Y <- 1
      Ds <- Da + 3
      Dd <- Da + 4
    }
    res <- mud * ((2*Xs*(2*Xs - 2))/(Dd*(Dd-2))) + (1-mud)*((max(Xs,Y)*(max(Xs,Y) - 1))/(Ds*(Ds-2)))
    return(res)
  }
  if(scs=="ZO" | scs=="ZW" | scs=="ZWW" | scs=="ZZW"){
    if(scs=="ZO"){
      Zd <- 1
      W <- 0
      Dd <- Da + 1
      Ds <- Da + 2
    }
    if(scs=="Zw"){
      Zd <- 1
      W <- 1
      Dd <- Da + 2
      Ds <- Da + 2
    }
    if(scs=="ZWW"){
      Zd <- 1
      W <- 2
      Dd <- Da + 3
      Ds <- Da + 2
    }
    if(scs=="ZZW"){
      Zd <- 2
      W <- 1
      Dd <- Da + 3
      Ds <- Da + 4
    }
    mus <- 1-mud
    res <- mus * ((2*Zd*(2*Zd - 2))/(Ds*(Ds-2))) + (1-mus)*((max(Zd,W)*(max(Zd,W) - 1))/(Dd*(Dd-2)))
    return(res)
  }
}
# calculating pval from line 194
sum = 0
for(i in 8:10){
  for(j in 0:(10-i)){
    sum <- sum + (factorial(10)/(factorial(i)*factorial(j)*factorial(10-i-j))) * 
      Pfsa(26, 'XXY', .5)^i * Pfaa(26, 'XXY', .5)^j * Pfss(26, 'XXY', .5)^(10-i-j)
  }
}
p.val <- sum
p.val
