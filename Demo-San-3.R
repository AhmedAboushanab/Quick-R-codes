# Demo-3 Binding Data for hRyIIA
# Odds ratio
# Same results with functions
data <- read.csv("Odds-Ratio-HfcyRIIA.csv")
Odds <- metabin(Ee, Ne, Ec, Nc, sm="OR", method = "I", data= data)
forest(Odds, comb.random = FALSE,hetstat = FALSE)

# Risk Difference
RiskD <- metabin(Ee, Ne, Ec, Nc, sm="RD", method = "I", data= data)
forest(RiskD, comb.random = FALSE,hetstat = FALSE)

# Risk Ratio
Risk <- metabin(Ee, Ne, Ec, Nc, sm="RR", method = "I", data= data)
forest(Risk, comb.random = FALSE,hetstat = FALSE)

# Arc Sinus Difference
ArcSD <- metabin(Ee, Ne, Ec, Nc, sm="ASD", method = "I", data= data)
forest(ArcSD, comb.random = FALSE,hetstat = FALSE)