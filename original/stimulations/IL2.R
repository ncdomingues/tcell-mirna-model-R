#IL2 stimulation (+anti CD28 CD3)


STATE0 <- c(1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

STATE0[23]=1 #Turn on BCM complex (BCL10 + CARD11 + MALT1)


nodes <- read.delim("nodes_values.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules <- read.delim("formulae.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)

model_new <- ginsimrun(STATE0, nodes, rules, 200, 200)
STATE0
model_new

model_new_155 <- model_new$evol[85,] #miR-155-5p
model_new_34c <- model_new$evol[84,] #miR-34c-5p
model_new_FOXO3 <- model_new$evol[46,] #FOXO3
model_new_TP53 <- model_new$evol[47,] #TP53
model_new_MDM2 <- model_new$evol[48,] #MDM2
model_new_GATA3 <- model_new$evol[72,] #GATA3
model_new_TBX21 <- model_new$evol[71,] #TBX21
model_new_IKB <- model_new$evol[25,] #NFKBIE
model_new_STAT3 <- model_new$evol[58,] #STAT3
model_new_RORC <- model_new$evol[73,] #RORC
model_new_FOS <- model_new$evol[17,]
model_new_JUN <- model_new$evol[16,] 
model_new_RAC1 <- model_new$evol[13,] 
model_new_FOXO1 <- model_new$evol[18,] 
model_new_SIRT1 <- model_new$evol[83,]
model_new_TGFBR <- model_new$evol[61,] 
model_new_FOXP3 <- model_new$evol[74,]
model_new_PRDM1_BLIMP1 <- model_new$evol[82,]
model_new_MYC <- model_new$evol[50,]
model_new_NFKB1 <- model_new$evol[24,] 
model_new_NFATC1 <- model_new$evol[22,]

#Other species

model_new_STAT1 <- model_new$evol[54,]
model_new_STAT4 <- model_new$evol[55,]
model_new_STAT5 <- model_new$evol[56,]
model_new_STAT6 <- model_new$evol[57,]
model_new_CDKN1A_p21 <- model_new$evol[12,]
model_new_BCL6 <- model_new$evol[75,]
model_new_AHR <- model_new$evol[78,]
model_new_FOXP1 <- model_new$evol[51,]
model_new_MAF <- model_new$evol[77,]
model_new_SPI1_PU1 <- model_new$evol[76,]
model_new_ZAP70 <- model_new$evol[7,]
model_new_AKT1 <- model_new$evol[75,]
model_new_ERK <- model_new$evol[28,]
model_new_JNK <- model_new$evol[27,]
model_new_CCNI <- model_new$evol[77,]



length(model_new_155)
behaviour_plot_IL2 <-plot(c(1:201), model_new_155, main="Species evolution with IL2", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), lwd=2, ylim=c(0,2), lwd=1)
lines(c(1:201), model_new_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), model_new_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)

legend(169.8,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend



lines(c(1:201), model_new_CDKN1A_p21, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), model_new_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), model_new_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), model_new_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)


legend(169.8,2,c(
    "CDKN1A_p21",
    "BCL6",
    "AHR",
    "FOXP1",
    "MAF",
    "SPI1_PU1",
    "STAT1",
    "STAT4",
    "STAT5",
    "STAT6",
    "ZAP70",
    "AKT1",
    "ERK",
    "JNK",
    "CCNI"
    
),lty=c(1,1),lwd=c(3,3),col = c(
    
    "blue",
    "cyan",
    "green",
    "yellow",
    "#234543",
    "#ab45ab",
    "black",
    "purple",
    "grey",
    "red",
    "#007788",
    "cyan",
    "green",
    "yellow",
    "#234543"
    
),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend

plot.new()

# - CYTOKYNES -

model_new_IL2 <- model_new$evol[34,] #miR-155-5p
model_new_IL6 <- model_new$evol[35,] #miR-34c-5p
model_new_IL12 <- model_new$evol[37,] #FOXO3
model_new_IL21 <- model_new$evol[38,] #TP53
model_new_IL17 <- model_new$evol[39,] #MDM2
model_new_IL23 <- model_new$evol[41,] #GATA3
model_new_IL4 <- model_new$evol[42,] #TBX21
model_new_IL10 <- model_new$evol[44,] #NFKBIE
model_new_IL2R <- model_new$evol[59,] #STAT3
model_new_IL6R <- model_new$evol[60,] #RORC
model_new_IL21R <- model_new$evol[63,]
model_new_IL23R <- model_new$evol[65,] 
model_new_IL4R <- model_new$evol[66,] 
model_new_IL10R <- model_new$evol[67,] 
model_new_IL22 <- model_new$evol[69,]
model_new_IL9 <- model_new$evol[70,] 
model_new_IL12R <- model_new$evol[79,]
model_new_IL1 <- model_new$evol[80,]
model_new_IL1R <- model_new$evol[81,]
model_new_TGFB <- model_new$evol[36,] 
model_new_TGFBR <- model_new$evol[61,]
model_new_TNF <- model_new$evol[40,] 
model_new_TNFR <- model_new$evol[64,]

cytokines_plot_IL2 <-plot(c(1:201), model_new_IL2, main="Cytokiner evolution with IL2", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), model_new_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2)
lines(c(1:201), model_new_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)

