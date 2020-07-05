###################################################################
################### CLASSICAL GENETIC ALGORITHM ###################
#################  by R. Gil-Merino, June, 2020	  #################
###################################################################
##  run code with: > source('SGA.R')
##  parameters must be fixed before execution
##  outputs save in file 'SGAevol.txt'
##  columns in file: 
##    generation, mean_fitness, best_fitness, best_chrom
###################################################################

## parameters:
##

popSize		 	<- 50	## population size 
chromLength		<- 14	## quantum chromosome length
total_generations 	<- 300	## number of generations

crossrate		<- 0.75 ## crossover probability
chromMutation		<- 0.2	## probability of chromosome mutation
bitMutation		<- 0.02	## probability of qubit mutation

###################### INIT POPULATION
######################
init_population <- function() {
  ## array with all individuals
  chromosome <- array(0,c(popSize,chromLength)) 
  for (i in 1:popSize) {
    for (j in 1:chromLength) {
      ## around half 0's and half 1's randomly
      if (runif(1,0,1) <= 0.5) {
        chromosome[i,j] <- 0
      }else{
        chromosome[i,j] <- 1
      }
    }
  }
  return(chromosome)
}

###################### FITNESS EVALUATION
######################
fitness_eval <- function(generation,chromosome) {
  fitness <- rep(0,popSize)
  ## loop to encode/decode problem
  for (i in 1:popSize) {  
    x<-0; y<-0
    for (j in 1:(chromLength/2)) {
      # binary to decimal values of x and y; uses two halves of chromosome
      x <- x + chromosome[i,j] * (2**(chromLength/2-j))
      y <- y + chromosome[i,j+chromLength/2] * (2**(chromLength/2-j))
    }
    ## compute fitness value, function to optimize
    fitness[i] <- x**2 + y**2
    cat("fitness ",i," = ",fitness[i],"\n")
  }
  ## best chromosome selection
  best_fitness = fitness[1]
  for (i in 1:popSize) {
    if (fitness[i] >= best_fitness) {
      best_fitness = fitness[i]
      best_chrom = i
    }
  }
  ## write stats in file SGAevol.txt
  write(c(generation,mean(fitness),max(fitness),best_chrom), sep = ",", file="SGAevol.dat", append=T)
  cat("Population Size = ",popSize,"\n")
  cat("mean fitness = ", mean(fitness),"\n")
  cat("fitness sum = ",sum(fitness),"\n")
  cat("variance = ",round(var(fitness))," std.deviation = ",round(sd(fitness)),"\n")
  cat("fitness max = ",max(fitness),"\n")
  return(list(fitness,best_chrom))
}

###################### TOURNAMENT SELECTION
######################
tournament <- function(fitness){
  u1 <- round(runif(1,1,popSize))		
  u2 <- round(runif(1,1,popSize))	
  if (fitness[u1] <= fitness[u2]){
    parent = u1
  }else{
    parent = u2
  }
  return(parent)
}

####################### CROSSOVER (1-point)
#######################
crossover <- function (crossrate,chromosome,fitness) {
  child1 <- array(0,c(popSize,chromLength))
  child2 <- array(0,c(popSize,chromLength))
  for (i in 1:popSize){
    papa <- tournament(fitness)
    mama <- tournament(fitness)
    if (runif(1,0,1) <= crossrate) {
      crosspoint <- round(runif(1,1,chromLength-1))
    }else{
      crosspoint <- 0
    }
    for (j in 1:chromLength) {
      if (j <= crosspoint) {
        child1[papa,j] <- chromosome[papa,j]
        child2[mama,j] <- chromosome[mama,j]
      }else{
        child1[papa,j] <- chromosome[mama,j]
        child2[mama,j] <- chromosome[papa,j]
      }
    }
    for (j in 1:chromLength) {
      chromosome[papa,j] <- child1[papa,j]
      chromosome[mama,j] <- child2[mama,j]
    }
  }
  return(chromosome)
}

####################### WHEEL SELECTION
#######################
wheelSelect <- function(chromosome,fitness){
  nchromosome <- array(0,c(popSize,chromLength))
  probability <- fitness / sum(fitness)
  cumulative  <- cumsum(probability)
  for (k in 1:popSize){
    random<-runif(1,0,1)
    for (i in 1:popSize){
      if (random < cumulative[i]) {
        parent <- i
        for (j in 1:chromLength) {
          nchromosome[k,j] <- chromosome[parent,j]
        }
        break
      }
    }
  }
  for (i in 1:popSize) {
    for (j in 1:chromLength) {
      chromosome[i,j] <- nchromosome[i,j]
    }
  }
  return(chromosome)
}

####################### MUTATION: chromosome and bit
#######################
mutation <- function (chromMutation,bitMutation,chromosome) {
  for (i in 1:popSize){
    if (runif(1,0,1) <= chromMutation) {
      for (j in 1:chromLength) {
        if (runif(1,0,1) <= bitMutation) {
          chromosome[i,j] <- 1 - chromosome[i,j]
        }
      }
    }
  }
  return(chromosome)
}

####################### SIMPLE GENETIC ALGORITHM (SGA)
#######################
mainSGA <- function() {
  generation = 1
  cat("****GENERATION: ",generation,"\n")
  chromosome <- init_population()
  fitnessdata <- fitness_eval(generation,chromosome)
  fitness <- fitnessdata[[1]]; best_chrom <- fitnessdata[[2]]
  for (generation in 2:total_generations) {
    cat("****GENERATION: ",generation,"\n")
    chromosome <- wheelSelect(chromosome,fitness)
    fitnessdata <- fitness_eval(generation,chromosome)
    fitness <- fitnessdata[[1]]; best_chrom <- fitnessdata[[2]]
    chromosome <- crossover(crossrate,chromosome,fitness)
    chromosome <- mutation(chromMutation, bitMutation,chromosome)
    cat("best chromosome of generation ",generation,": ",best_chrom,"\n")
  }
}

####################### MAIN PROGRAM
#######################
cat("****CLASSICAL GENETIC ALGORITHM****",'\n')
cat('population size: ', popSize,'\n')
cat('chromosome length: ',chromLength,'\n')
cat('number of generations: ',total_generations,'\n')
cat('probability of chromosome mutation: ', chromMutation,'\n')
cat('probability of qubit mutation: ', bitMutation,'\n')
readline(prompt="Press [enter] to continue")
cat(c('generation','mean_fitness','best_fitness','best_chrom\n'), sep=',', file="SGAevol.dat", append=T)
mainSGA()



