###################################################################
################# BASIC QUANTUM GENETIC ALGORITHM #################
#################  by R. Gil-Merino, June, 2020	  #################
###################################################################
##  run code with: > source('BQGA.R')
##  parameters must be fixed before execution
##  outputs save in file 'BQGAevol.txt'
##  columns in file: 
##    generation, mean_fitness, best_fitness, best_chrom
###################################################################

## parameters:
##

popSize		 	<- 50	## population size 
quChromLength		<- 14	## quantum chromosome length
total_generations 	<- 300	## number of generations

quChromMutation		<- 0.2	## probability of chromosome mutation
qubitMutation		<- 0.02	## probability of qubit mutation

####################### INIT POPULATION
#######################
init_population <- function(){
  ## quantum chromosome for the whole population
  quChrom  <- array(0,c(popSize,quChromLength,2))
  ## Hadamard gate
  r2 <- as.double(sqrt(2))
  h <- matrix(c(1/r2,1/r2,1/r2,-1/r2),2,2)
  ## rotation gate
  theta <- 0
  rot <- matrix(0,2,2)
  ## qubit with basis states
  quBitZero <- c(1,0)
  quBitOne  <- c(0,1)
  alphaBeta <- matrix(0,1,2)
  ## Initial population array, quChrom 
  for (i in 1:popSize){
    for (j in 1:quChromLength) {
      theta <- runif(1,0,1) * 90
      theta <- theta * pi / 180
      rot[1,1] <- cos(theta); rot[1,2] <- -sin(theta)
      rot[2,1] <- sin(theta); rot[2,2] <-  cos(theta)
      alphaBeta[1] <- rot[1,1] * (h[1,1] * quBitZero[1]) + rot[1,2] * (h[1,2] * quBitZero[2])
      alphaBeta[2] <- rot[2,1] * (h[2,1] * quBitZero[1]) + rot[2,2] * (h[2,2] * quBitZero[2])
      ## alpha squared
      quChrom[i,j,1] <- round(2 * alphaBeta[1]**2,2)
      ## beta squared
      quChrom[i,j,2] <- round(2 * alphaBeta[2]**2,2)
    }
  }
  return(quChrom)
}

####################### MEASUREMENT: uses only alpha from quBit
#######################
measure <- function(alpha_prob,quChrom){
  chromosome <- array(0,c(popSize,quChromLength)) #classical chromosome!
  for (i in 1:popSize) {
    for (j in 1:quChromLength) {
      if (alpha_prob <= quChrom[i,j,1]) {
        chromosome[i,j] = 0
      }else{
        chromosome[i,j] = 1
      }
    }
  }
  return(chromosome)
}

####################### FITNESS EVALUATION
#######################
fitness_eval <- function(generation,chromosome) {
  fitness <- rep(0,popSize)
  ## loop to encode/decode problem
  for (i in 1:popSize) { 
    x<-0; y<-0
    for (j in 1:(quChromLength/2)) {
      ## binary to decimal
      x <- x + chromosome[i,j] * 2**(quChromLength/2-j)
      y <- y + chromosome[i,j+quChromLength/2] * 2**(quChromLength/2-j)
    }
    ## compute fitness value, function to optimize
    fitness[i] <- x**2 + y**2 
    cat("fitness ",i," = ",fitness[i],"\n")
  }
  ## best chromosome
  best_fitness = fitness[1]
  for (i in 1:popSize) {
    if (fitness[i] >= best_fitness) {
      best_fitness = fitness[i]
      best_chrom = i
    }
  }
  ## plot stats
  write(c(generation,mean(fitness),max(fitness),best_chrom), sep = ",", file="BQGAevol.dat", append=T)
  cat("Population Size = ",popSize,"\n")
  cat("mean fitness = ", round(mean(fitness)),"\n")
  cat("fitness sum = ",sum(fitness),"\n")
  cat("variance = ",round(var(fitness))," std.deviation = ",round(sd(fitness)),"\n")
  cat("best fitness = ",max(fitness),"\n")
  return(list(fitness,best_chrom))
}

####################### QUANTUM ROTATION GATE
#######################
rotation <- function(chromosome,quChrom,fitness,best_chrom){
  newquChrom <- array(0,c(popSize,quChromLength,2))
  rot <- matrix(0,2,2)
  for (i in 1:popSize) {
    for (j in 1:quChromLength) {
      if (fitness[i] < fitness[best_chrom]) { 
        if (chromosome[i,j] == 0 && chromosome[best_chrom,j] == 1) {
          ## define rotation angle
          ## delta_theta = 0.0785398163
          delta_theta = pi / 32
          rot[1,1] <- cos(delta_theta); rot[1,2] <- -sin(delta_theta)
          rot[2,1] <- sin(delta_theta); rot[2,2] <-  cos(delta_theta)
          newquChrom[i,j,1] <- (rot[1,1] * quChrom[i,j,1]) + (rot[1,2] * quChrom[i,j,2])
          newquChrom[i,j,2] <- (rot[2,1] * quChrom[i,j,1]) + (rot[2,2] * quChrom[i,j,2])
          quChrom[i,j,1]  <- round(newquChrom[i,j,1],2)
          quChrom[i,j,2]  <- round(1-newquChrom[i,j,1],2)
        }
        if (chromosome[i,j] == 1 && chromosome[best_chrom,j] == 0) {
          ## define rotation angle
          ## delta_theta = -0.0785398163
          delta_theta = -pi / 32
          rot[1,1] <- cos(delta_theta); rot[1,2] <- -sin(delta_theta)
          rot[2,1] <- sin(delta_theta); rot[2,2] <-  cos(delta_theta)
          newquChrom[i,j,1] <- (rot[1,1] * quChrom[i,j,1]) + (rot[1,2] * quChrom[i,j,2])
          newquChrom[i,j,2] <- (rot[2,1] * quChrom[i,j,1]) + (rot[2,2] * quChrom[i,j,2])
          quChrom[i,j,1]  <- round(newquChrom[i,j,1],2)
          quChrom[i,j,2]  <- round(1-newquChrom[i,j,1],2)
        }
      }
    }
  }
  return(quChrom)
}

####################### QUANTUM MUTATION GATE: X-PAULI
#######################
## pop_mutation_rate: whether a chromosome mutates
## mutation_rate: whether a bit mutates
mutation <- function(quChromMutation,qubitMutation,quChrom) {
  newquChrom <- array(0,c(popSize,quChromLength,2))
  for (i in 1:popSize) {
    if (runif(1,0,1) <= quChromMutation) {
      for (j in 1:quChromLength) {
        if (runif(1,0,1) <= qubitMutation) {
          newquChrom[i,j,1] <- quChrom[i,j,2]
          newquChrom[i,j,2] <- quChrom[i,j,1]
        }else{
          newquChrom[i,j,1] <- quChrom[i,j,1]
          newquChrom[i,j,2] <- quChrom[i,j,2]  
        }
      }
    }else{
      for (j in 1:quChromLength) {
        newquChrom[i,j,1] <- quChrom[i,j,1]
        newquChrom[i,j,2] <- quChrom[i,j,2]
      }
    }
  }
  for (i in 1:popSize) {
    for (j in 1:quChromLength) {
      quChrom[i,j,1] <- newquChrom[i,j,1]
      quChrom[i,j,2] <- newquChrom[i,j,2]  
    }
  }
 return(quChrom)
}

####################### BASIC QUANTUM GENETIC ALGORITHM (BQGA)
#######################
mainQGA <- function() {
  generation = 1
  cat("****GENERATION: ",generation,"\n")
  quChrom <- init_population()
  chromosome <- measure(0.5,quChrom)
  fitnessdata <- fitness_eval(generation,chromosome)
  fitness <- fitnessdata[[1]]; best_chrom <- fitnessdata[[2]]
  for (generation in 2:total_generations) {
    cat("****GENERATION: ",generation,"\n")
    quChrom <- rotation(chromosome,quChrom,fitness,best_chrom)
    ## mutation has chromosome mutation and bit mutation
    quChrom <- mutation(quChromMutation,qubitMutation,quChrom)
    ## measurement with prob. alpha = 0.5
    chromosome <- measure(0.5,quChrom)
    fitnessdata <- fitness_eval(generation,chromosome)
    fitness <- fitnessdata[[1]]; best_chrom <- fitnessdata[[2]]
    cat("best chromosome of generation ",generation,": ",best_chrom,"\n")
  }
}

####################### MAIN PROGRAM
#######################
cat("****QUANTUM GENETIC ALGORITHM****",'\n')
cat('population size: ', popSize,'\n')
cat('chromosome length (number of qubits): ',quChromLength,'\n')
cat('number of generations: ',total_generations,'\n')
cat('probability of chromosome mutation: ', chromMutation,'\n')
cat('probability of qubit mutation: ', qubitMutation,'\n')
readline(prompt="Press [enter] to continue")
cat(c('generation','mean_fitness','best_fitness','best_chrom\n'), sep=',', file="BQGAevol.dat", append=T)
mainQGA()

  
