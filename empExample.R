source("pass-t.R")

# simulation parameters
presTime <- rep(c(2, 4, 2, 8, 4, 2, 8, 4, 8), 2)
presRate <-rep(c(2, 2, 4, 2, 4, 8, 4, 8, 8), 2)
learningWeight <- 0.01
pulsesPerSecond <- 1
attentionFreq <- c(rep("constant", length(presTime)))
attentionTime <- c(rep("constant", length(presTime)))
inputs <- diag(18)

# Simulation for words
attentionC <- .015

namesSimFreq <- SimulatePassT(inputMatrix = inputs,
                       learningWeight = learningWeight,
                       presRate = presRate, presTime = presTime,
                       runs=500, attentionFunction = attentionFreq,
                       attentionC=attentionC, nOutputUnits = 9)

namesSimFreqStrength <- as.data.frame(namesSimFreq$output)
namesSimFreqStrengthAgg <- aggregate(t(namesSimFreqStrength),
                                     by=list(presRate, presTime*presRate),
                                     mean)[,-c(1,2)]
namesSimFreqJofMean <- apply(namesSimFreqStrengthAgg, 1, mean)

# Simulation for pictures
attentionC <- .1

picsSimFreq <- SimulatePassT(inputMatrix = inputs,
                              learningWeight = learningWeight,
                              presRate = presRate, presTime = presTime,
                              runs=500, attentionFunction = attentionFreq,
                              attentionC=attentionC, nOutputUnits = 9)

picsSimFreqStrength <- as.data.frame(picsSimFreq$output)
picsSimFreqStrengthAgg <- aggregate(t(picsSimFreqStrength),
                                    by=list(presRate, presTime*presRate),
                                    mean)[,-c(1,2)]
picsSimFreqJofMean <- apply(picsSimFreqStrengthAgg, 1, mean)

simJof <- c(namesSimFreqJofMean, picsSimFreqJofMean)

# get empirical data from Isabells experiment 3
if(!exists("isabExp3Judgm.csv")) source("getDataFromIsab.R")
isabExp3 <- read.table("isabExp3Judgm.csv", sep=";", header=T)
# just take names and neutral pics
isabExp3 <- isabExp3[isabExp3$stimulusType %in% c("Namen", "neutrale Bilder"),]

cor(isabExp3$meanJof, simJof)
cor(isabExp3$meanJod, simJof)

# plot ####
mycolors.colorblindsafe <- c("#1b9e77", "#d95f02", "#7570b3")
mycolor <- c(rep(mycolors.colorblindsafe[1], 9), rep(mycolors.colorblindsafe[2], 9))
mypch <- c(rep(16, 9), rep(15, 9))
mycex <- 1.4

pdf("corrEmpSim.pdf", width=10, height=7)
  par(cex=mycex)
  plot(isabExp3$meanJof, simJof,
       xlim=c(min(isabExp3$meanJof), max(isabExp3$meanJof)+0.3),
       xlab="Empirie: Häufigkeitsschätzung",
       ylab="Simulation: Gesamtaktivierung",
       main="Vergleich Simulation und Empirie für Häufigkeitsschätzung",
       col=mycolor, pch=mypch, bty="l")
  # labels for data points
  text(isabExp3$meanJof, simJof, paste(presRate, "x", presTime, "s", sep=""),
       pos=4, cex=0.9)
  # description and correlation
  text(min(isabExp3$meanJof), max(simJof*0.98),
       paste("Zahlen an den Punkten: Häufigkeit x Dauer ",
             "(einzeln)\nKorrelation = ",
             round(cor.test(isabExp3$meanJof, simJof)$estimate,2), sep=""),
       adj=c(0.0, 0.7), offset=0)
  legend("bottomright",
         legend=c("niedrige Aufmerksamkeit", "hohe Aufmerksamkeit"),
         pch=unique(mypch), col=unique(mycolor), bty="n")
dev.off()