###################################################################
################ REDUCED QUANTUM GENETIC ALGORITHM ################
#################  by R. Gil-Merino, June, 2020	  #################
###################################################################
##  run code with: > source('RQGA.R')
##  prompt: even number of qubits (in laptops keep it <12)
##  uncomment lines 124 & 144 to send outputs to file 'RQGAout.txt'
###################################################################

####################### binary string TO decimal
#######################
bin2dec <- function(string_num) {
  return(strtoi(string_num,base=2))
}

####################### decimal TO binary vector
#######################
dec2bin <- function(x) {
  if (x == 1)
    1
  else if (x == 0)
    NULL
  else {
   mod <- x %% 2
   c(dec2bin((x-mod) %/% 2), mod)
  }
}

####################### decimal TO binary string
#######################
bin <- function(x) {
  result <- dec2bin(x)
  return(paste0(result, collapse = ""))
}

####################### decimal & length TO vector with 1
#######################
dec2vec <- function(dec,n) {
  vec <- matrix(0.0,2**n,1)
  vec[dec+1,1] <- 1;
  return(vec)
}

####################### quantum vector FROM string
#######################
psi <- function(string_num) {
  dec <- bin2dec(string_num)
  return (dec2vec(dec,nchar(string_num)))
}

####################### N-bit HADAMARD matrix
#######################
hadamard <- function(n) {
  r2 <- sqrt(2)
  H0 <- matrix(c(1/r2,1/r2,1/r2,-1/r2),2,2)
  if (n==1) {
    H <- H0
  }else{
    H = 1
    for (i in 1:n) {
      H <- kronecker(H,H0)
    }
  }
  return(H)
}

####################### FITNESS quantum gate
#######################
fitness_eval <- function(n) {
  fitness <- rep(0,2**n)
  for (i in 0:(2**n-1)) {
    ## string of length n with binary value of i
    binary_i <- paste(strrep('0',n-nchar(bin(i))), bin(i), sep='') 
    ## decimal value of x and y from string binary_i
    x <- strtoi(substr(binary_i,1,n/2),2)  
    y <- strtoi(substr(binary_i,(n/2+1),n),2)
    ## function to optimize
    fitness[i+1] <- x**2 + y**2
  }
  ## get the best qu-chromosome
  best_chrom <- 0
  fitness_max <- fitness[1]
  for (i in 1:2**n) {	
    if (fitness[i] >= fitness_max){
      fitness_max <- fitness[i]
      best_chrom <- i
    }
  }
  return(best_chrom)
}

####################### ORACLE: conditional inverter
#######################
oracle_matrix <- function(n) {
  matriz <- matrix(0,2**n,2**n)
  for (i in c(1:2**n)) {
    if ( i == fitness_eval(n) ) {
      oraclevalue <- 1
    }else{
      oraclevalue <- 0
    }
    matriz[i,i] <- (-1)**oraclevalue
  }
  return(matriz)
}

####################### GROVER's DIFUSSION
#######################
# inversion around average
avr_inversion <- function(n){
  matrizavr <- rep(2/2**n,2**n)
  matrizavr <- matrizavr - diag(2**n)
  return(matrizavr)
}

####################### GROVER's ITERATIONS
#######################
grover_iter_max <- function(n) {
  return( (pi/4)*sqrt(2**n) )
}

####################### REDUCED QUANTUM GENETIC ALGORITHM (rQGA)
#######################
rqga <- function(n) {
  ## uncomment next and last lines to send outputs to file: sink(...) --- sink()
  ##sink("RQGAout.txt") 
  string_num <- paste(replicate(n,'0'),collapse='')
  ## product psi from string with Hadamard
  psi_vec <- hadamard(n) %*% psi(string_num)
  cat("psi_vec inicial: ",'\n',psi_vec,'\n')
  ## inversion about the average
  diffusion <- avr_inversion(n)
  ## apply oracle and inversion about the average a number of Grover's iterations
  for (i in 1:trunc(grover_iter_max(n))) {
    oracle <- oracle_matrix(n)
    psi_vec <- diffusion %*% oracle %*% psi_vec
  }  
  ## output final psi, maximum amplitud and decimal values of x and y
  cat("final psi_vec amplitudes: ",'\n',t(psi_vec),'\n')
  cat("max. amplitude: ",max(psi_vec),'\n')
  base_state_max <- bin(which(psi_vec==max(psi_vec))-1)
  cat("base state with max. amplitude: ", base_state_max,'\n')
  cat("equivalent decimal values: ")
  cat("--> x=", strtoi(substr(base_state_max,1,nchar(base_state_max)/2),2))
  cat("  y= ", strtoi(substr(base_state_max,nchar(base_state_max)/2+1, nchar(base_state_max)),2), '\n')
  ##sink()
}

####################### MAIN PROGRAM
#######################
cat("****REDUCED QUANTUM GENETIC ALGORITHM****",'\n')
nqubits <- readline(prompt="Enter number of qubits (even): ")
rqga(as.integer(nqubits))








