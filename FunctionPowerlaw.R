# defining distribution parameters
param_rlk = numeric(3) # i parametri del modello power law

## defining function A(r) # è la funzione che nel power law definisce il cutoff 
# Nei data Tara oceans va da 0 a 10^-2 
A_func = function(param_rlk){
  out = 0
  for (i in 1 : length(vec)){
    if (vec[i] <= x_max){
      out = out + vec[i]
    }
  }
  out = -(out*param_rlk[1])
  return(out)
}

## defining function N(r,l,k)
N_func = function(param_rlk){
  out = 0
  for (x in 1 : x_max){
    out = out + exp(-param_rlk[1]*x)*exp(lgamma(x+param_rlk[2])-lgamma(x+param_rlk[3]+1))
  }
  return(out)
}

## defining function C(l)
C_func = function(param_rlk){
  out = 0
  for (i in 1 : length(vec)){
    if (vec[i] <= x_max){
      out = out + lgamma(vec[i]+param_rlk[2])
    }
  }
  return(out)
}   


## defining function D(k)
D_func = function(param_rlk){
  out = 0
  for (i in 1 : length(vec)){
    if (vec[i] <= x_max){
      out = out + lgamma(vec[i]+param_rlk[3]+1)
    }
  }
  out = -out
  return(out)
}   

## defining likelihood
L_func = function(param_rlk){
  out = 0
  out = A_func(param_rlk) + -n*log(N_func(param_rlk)) + C_func(param_rlk) + D_func(param_rlk)
  return(-out)
}