source("pass-t.R")
# very simple example with 9 orthogonal stimuli

# simulation parameter
inputMatrix <- diag(9) # create 9 orthogonal input patterns
presTime <- c(2, 4, 2, 8, 4, 2, 8, 4, 8)
presRate <- c(2, 2, 4, 2, 4, 8, 4, 8, 8)
learningWeight <- 0.1
runs <- 100
AttentionC <- .15 # learning weight falls to AttentionC*learningWeight

# start simulation
MySim <- SimulatePassT(inputMatrix = inputMatrix,
                       learningWeight = learningWeight,
                       presRate = presRate, presTime = presTime,
                       runs=runs, AttentionC=AttentionC, nOutputUnits = 9)

# extract overall activation strength for every input pattern
MySimStrength <- as.data.frame(MySim[[1]])

MySimStrengthSd <- apply(MySimStrength, 2, sd)
MySimStrengthMean <- apply(MySimStrength, 2, mean)

presTimeSum <- presTime*presRate

plot(presTimeSum, MySimStrengthMean, xlab="actual presentation time",
     ylab="Simulation: overall activation",
     main="Activation in simulation as a function of actual presentation time",
     pch=16, bty="l",)
text(presTimeSum, MySimStrengthMean, paste(presRate, "x", presTime, "s",
                                           sep=""), pos=4, cex=0.9)

plot(presRate, MySimStrengthMean, xlab="actual presentation frequency",
     ylab="Simulation: overall activation",
     main="Activation in simulation as a function of actual presentation
     frequency", pch=16, bty="l")
text(presRate, MySimStrengthMean, paste(presRate, "x", presTime, "s", sep=""),
     pos=4, cex=0.9)
