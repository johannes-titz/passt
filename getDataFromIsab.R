function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}

data <- read.table("isabExp3Compl.csv", sep=";", header=T)
empJofMean <- aggregate(data$avHF, by=list(data$Häuf, data$Zeit, data$Stimuli), mean)
empJodMean <- aggregate(data$avZeitG, by=list(data$Häuf, data$Zeit, data$Stimuli), mean)

isabExp3 <- RoundDf(cbind(empJofMean, empJodMean$x), 3)
names(isabExp3) <- c("frequency", "totalDuration", "stimulusType", "meanJof", "meanJod")
write.table(isabExp3, file="isabExp3Judgm.csv", col.names=T, row.names = F, sep=";")