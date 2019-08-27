SimulatePassT <- function(inputMatrix, lrate, presRate, presTime,
                          runs = 100, n_output_units=round(ncol(inputMatrix)/2),
                          pulsesPerSecond = 1,
                          attentionFunction = rep("constant", length(presTime)),
                          interfUpdating = 0, attentionA = 3, attentionB = 10,
                          attentionC = .15, interfPres = 0,
                          random_order = T){
  # Simulates human frequency and time estimation with a competitive learning
  # network based on the PASS-family (Sedlmeier, 1999; 2002; Burkhardt, 2013;
  # Titz, 2014)
  #
  # Args:
  #   inputMatrix: matrix with input patterns that will be presented to the
  #                neural network, the row is one pattern
  #   lrate: determines how fast the neural network will learn, value
  #                   must be between 0-1, values too low and too high will
  #                   disable learning
  #   presRate: presentation rates for the patterns in the matrix
  #   presTime: presentation times for the patterns in the matrix
  #   runs: number of simulations to be run
  #   n_output_units: number of output units in the neural network
  #   pulsesPerSecond: determines how many time pulses one second has
  #   attentionFunction: determines how attention changes during presentation
  #                      of one stimulus
  #   interfUpdating: interference when updating winner-weights from normal
  #                   distribution with mean=0, sd=interfUpdating
  #   attentionA: determines attentionFunction
  #   attentionB: determines attentionFunction
  #   attentionC: determines attentionFunction
  #   interfPres: percent of stimuli that will be changed to a random pattern;
  #               the lower the attention the higher the probability that this
  #               will happen for a stimulus
  #
  # Returns:
  #   list with following values for the number of simulations specified above
  #   outputStrength: the sum of the activation strengths of the output units
  #                   for the input patterns
  #   sd: the standard deviation between output-units for every pattern
  #   weightMatrix: final weightMatrix
  #   presentationMatrix: input matrix

  outputStrength <- NULL
  outputVar <- NULL

  for (j in 1:runs) {
    n_input_units <- ncol(inputMatrix)
    weightMatrix <- InitializeWeightMatrix(n_input_units, n_output_units)
    pres <- MakePresMatrix(inputMatrix, lrate, presRate, presTime,
                           pulsesPerSecond=pulsesPerSecond,
                           attentionFunction=attentionFunction,
                           attentionA=attentionA,
                           attentionB=attentionB, attentionC=attentionC,
                           random_order = random_order)
    if (interfPres>0) {
      nInterfPres <- round(interfPres*nrow(pres$input))
      interfPresRows <- sample(1:nrow(pres$input), nInterfPres,
                               prob = pres$lrate)
      for (i in 1:nInterfPres) {
        interfPresCol <- sample(c(1, 0), n_input_units, replace = T)
        pres$input[interfPresRows[i]] <- interPresCol
      }
    }

    # present input i and update weights
    for (i in 1:nrow(pres$input)) {
      weightMatrix <- UpdateWinnerWeights(pres$input[i, ], weightMatrix,
                                          pres$lrate[i],
                                          interfUpdating)
    }

    # output for simulation j
    outputStrength <- rbind(outputStrength, OutputActivation(inputMatrix,
                                                             weightMatrix))
    outputVar <- rbind(outputVar, OutputVar(inputMatrix, weightMatrix))
  }
  finalList <- list("output"=outputStrength, "outputSd"=outputVar,
                    "weightMatrix"=weightMatrix,
                    "presentationMatrix"=pres)
  return(finalList)
}

MakePresMatrix <- function(inputMatrix, lrate, presRate, presTime,
                           pulsesPerSecond, attentionFunction,
                           attentionA, attentionB, attentionC,
                           random_order = T){
  # Makes a presentation matrix with learning weights for the function
  # SimulatePassT
  #
  # Args:
  #   inputMatrix: matrix with input patterns that will be presented to
  #                the neural network, the row is one pattern
  #   lrate: determines how fast the neural network will learn, value
  #                   must be between 0-1, values too low and too high will
  #                   disable learning
  #   presRate: presentation rates for the patterns in the matrix
  #   presTime: presentation times for the patterns in the matrix
  #   pulsesPerSecond: determines how many time pulses one second has
  #   attentionFunction: determines how attention changes during presentation
  #                      of one stimulus
  #   attentionA: determines attentionFunction
  #   attentionB: determines attentionFunction
  #   attentionC: determines attentionFunction
  #
  # Returns:
  #   list with two lists: the presentation matrix and the learning weights
  #   (both in one random order)

  attention <-  GetAttention(lrate, presTime,
                             pulsesPerSecond=pulsesPerSecond,
                             attentionFunction, attentionA=attentionA,
                             attentionB=attentionB, attentionC=attentionC)
  presList <- list()
  finalList <- list()
  finalResult <- list()
  # produce list (stimuli) of lists (presentation time) of smallest entity
  # (pulses*time)
  for (i in 1:nrow(inputMatrix)) {
    presTimePulse <- presTime[i]*pulsesPerSecond
    presList[[i]] <- list(cbind(matrix(rep(inputMatrix[i, ], presTimePulse),
                                       nrow=presTimePulse, byrow=T),
                                attention[[i]]))
    presList[[i]] <- lapply(presList[i], rep, presRate[i])
  }
  # unlist so that all presentations remain lists, apply unlist two times
  presList <- unlist(presList, recursive=F)
  presList <- unlist(presList, recursive=F)

  # randomize presentations
  if (random_order) {
    randomIndex <- sample(1:length(presList), length(presList))
  } else {
    randomIndex <- 1:length(presList)
  }
  presList <- presList[randomIndex]

  # dataFrame
  presDf <- do.call(rbind, presList)

  # produce result
  finalResult[["input"]] <- presDf[, 1:(ncol(presDf)-1)]
  finalResult[["lrate"]] <- presDf[, ncol(presDf)]
  return(finalResult)
}


InitializeWeightMatrix <- function(n_input_units, n_output_units){
  # Initializes a weight matrix for a competitive learning network algorithm
  #
  # Args:
  #   n_input_units: number of input units in the net
  #   n_output_units: number of output units in the net
  #
  # Returns:
  #   matrix with n_input_units rows and n_output_units columns, the sum of
  #   every column is 1

  nWeights <- n_output_units * n_input_units
  # initialize random weight matrix
  weightMatrix <- matrix(rnorm(nWeights, 0.5, 0.005), n_output_units, byrow = T)
  # weights for every output unit (row sum) must equal 1
  weightMatrix <- prop.table(weightMatrix, margin = 1)
  return(weightMatrix)
}

UpdateWinnerWeights <- function(input, weightMatrix, lrate,
                                interfUpdating = 0){
  # Updates weights for competitive learning network algorithm
  #
  # Args:
  #   input: input vector
  #   weightMatrix: weight matrix of network
  #   lrate: learning weight for algorithm
  #   interfUpdating: sd for a normal distribution with mean = 0, used for
  #   interference
  #
  # Returns:
  #   new weight matrix for step t + 1

  output <- weightMatrix %*% input
  winner <- which(output == max(output))
  if (length(winner) > 1) {
    print ("more than one winner, random selection")
    winner <- sample(winner, 1)
  }

  # update winner row
  # interfere with updating, interference from normal distribution with mean=0
  # and sd=interfUpdating
  # learning rule for competitive learning: Rummelhart & Zipser (1986)
  if (interfUpdating>0) {
    deltaW <- lrate*(input/sum(input)-weightMatrix[winner,]+
                                rnorm(length(input), 0, interfUpdating))
    weightMatrix[winner,] <- weightMatrix[winner,]+deltaW
    weightMatrix[winner,] <- prop.table(weightMatrix[winner,]+
                                         max(c(0, -min(weightMatrix[winner,]))))
  } else {
    deltaW <- lrate * (input / sum(input) - weightMatrix[winner, ])
    weightMatrix[winner, ] <- weightMatrix[winner, ] + deltaW
  }
  return(weightMatrix)
}

OutputVar <- function(inputs, weightMatrix){
    # calculates variance in output unit activation for one or several input
    # patterns
    #
    # Args:
    #   input: input vector or input matrix
    #   weightMatrix: weight matrix of network
    #
    # Returns:
    #   variance in output unit activation

  if (class(inputs)=="matrix") {
    return(apply(apply(inputs, 1, function(x) weightMatrix%*%x), 2, var))
  }
  if (is.vector(inputs)==T) {
    return(var(weightMatrix%*%inputs))
  }
}

OutputActivation <- function(inputs, weightMatrix){
  # calculates sum of activation in output units for one or several input
  # patterns
  #
  # Args:
  #   input: input vector or input matrix
  #   weightMatrix: weight matrix of network
  #
  # Returns:
  #   activation in output units

  if (class(inputs)=="matrix") {
    return(colSums(apply(inputs, 1, function(x) weightMatrix%*%x)))
  }
  if (is.vector(inputs)==T) {
    return(sum(weightMatrix%*%inputs))
  }
}


GetAttention <- function(lrate, presTime, pulsesPerSecond,
                         attention=rep("constant", length(presTime)),
                         attentionA=10, attentionB=2, attentionC=.15){
  # gives learning weight from an attention function
  #
  # Args:
  #
  #   attention: specifies how attention should evolve over time
  #   attentionA: modifies high attention function
  #   attentionB: modifies high attention function
  #   attentionC: percent value that attention sinks to after first time step
  #
  # Returns:
  #   attention for every time pulse

  duration <- presTime*pulsesPerSecond
  y <- list()
  x <- lapply(duration, function(x) seq(0, x-1))

  for (i in 1:length(presTime)) {
    switch(attention[i],
           high={
             y[[i]] <- c(lrate,
                         lrate*(exp(-attentionA*(x[[i]])+attentionB)/
                                           (1+exp(-attentionA*(x[[i]])+
                                                     attentionB)))[-1])
           },
           low={
             y[[i]] <- c(lrate*(exp(-1.3*(x[[i]]))))
           },
           max={
             y[[i]] <- rep(lrate, length(x[[i]]))
           },
           min={
             y[[i]] <- c(lrate, rep(0, length(x[[i]])-1))
           },
           constant={
             y[[i]] <- c(lrate, rep(lrate*attentionC,
                                             length(x[[i]])-1))
           })
  }
  return(y)
}