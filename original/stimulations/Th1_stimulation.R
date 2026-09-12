#Th1 stimulation

STATE0 <- c(1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

STATE0[23]=1 #Turn on BCM complex (BCL10 + CARD11 + MALT1)
STATE0[37]=1 #Turn on Th1 stimulation with IL12 + APC + CD28 + CD3

Th1_stim <- ginsimrun(STATE0, nodes, rules, 200, 200)


Th1_stim_155 <- Th1_stim$evol[85,] #miR-155-5p
Th1_stim_34c <- Th1_stim$evol[84,] #miR-34c-5p
Th1_stim_FOXO3 <- Th1_stim$evol[46,] #FOXO3
Th1_stim_TP53 <- Th1_stim$evol[47,] #TP53
Th1_stim_MDM2 <- Th1_stim$evol[48,] #MDM2
Th1_stim_GATA3 <- Th1_stim$evol[72,] #GATA3
Th1_stim_TBX21 <- Th1_stim$evol[71,] #TBX21
Th1_stim_IKB <- Th1_stim$evol[25,] #NFKBIE
Th1_stim_STAT3 <- Th1_stim$evol[58,] #STAT3
Th1_stim_RORC <- Th1_stim$evol[73,] #RORC
Th1_stim_FOS <- Th1_stim$evol[17,]
Th1_stim_JUN <- Th1_stim$evol[16,] 
Th1_stim_RAC1 <- Th1_stim$evol[13,] 
Th1_stim_FOXO1 <- Th1_stim$evol[18,] 
Th1_stim_SIRT1 <- Th1_stim$evol[83,]
Th1_stim_TGFBR <- Th1_stim$evol[61,] 
Th1_stim_FOXP3 <- Th1_stim$evol[74,]
Th1_stim_PRDM1_BLIMP1 <- Th1_stim$evol[82,]
Th1_stim_MYC <- Th1_stim$evol[50,]
Th1_stim_NFKB1 <- Th1_stim$evol[24,] 
Th1_stim_NFATC1 <- Th1_stim$evol[22,]




length(Th1_stim_155)
behaviour_plot_IL2 <-plot(c(1:201), Th1_stim_155, main="Species evolution with IL2+IL4", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), Th1_stim_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)


legend(169.5,2,c("miR-155-5p","34c","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend


plot.new()


Th1_stim_IL2 <- Th1_stim$evol[34,] 
Th1_stim_IL6 <- Th1_stim$evol[35,] 
Th1_stim_IL12 <- Th1_stim$evol[37,] 
Th1_stim_IL21 <- Th1_stim$evol[38,]
Th1_stim_IL17 <- Th1_stim$evol[39,] 
Th1_stim_IL23 <- Th1_stim$evol[41,] 
Th1_stim_IL4 <- Th1_stim$evol[42,] 
Th1_stim_IL10 <- Th1_stim$evol[44,] 
Th1_stim_IL2R <- Th1_stim$evol[59,]
Th1_stim_IL6R <- Th1_stim$evol[60,]
Th1_stim_IL21R <- Th1_stim$evol[63,]
Th1_stim_IL23R <- Th1_stim$evol[65,] 
Th1_stim_IL4R <- Th1_stim$evol[66,] 
Th1_stim_IL10R <- Th1_stim$evol[67,] 
Th1_stim_IL22 <- Th1_stim$evol[69,]
Th1_stim_IL9 <- Th1_stim$evol[70,] 
Th1_stim_IL12R <- Th1_stim$evol[79,]
Th1_stim_IL1 <- Th1_stim$evol[80,]
Th1_stim_IL1R <- Th1_stim$evol[81,]
Th1_stim_TGFB <- Th1_stim$evol[36,] 
Th1_stim_TGFBR <- Th1_stim$evol[61,]
Th1_stim_TNF <- Th1_stim$evol[40,] 
Th1_stim_TNFR <- Th1_stim$evol[64,]

cytokines_plot_IL2 <-plot(c(1:201), Th1_stim_IL2,main="Cytokine and receptor evolution with IL2+IL4", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), Th1_stim_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), Th1_stim_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)

legend(169.8,2,c("IL2","IL6","IL12","IL21","IL17","IL23","IL4","IL10","IL2R","IL6R","IL21R","IL23R","IL4R","IL10R","IL22","IL9","IL12R","IL1","IL1R","TGFB","TGFBR","TNF","TNFR"
                 
),lty=c(1,1),lwd=c(3,3),col = c("#4d0000",
                                "#4d1300",
                                "#4d2600",
                                "#4d3900",
                                "#4d4d00",
                                "#394d00",
                                "#456789",
                                "black",
                                "purple",
                                "grey",
                                "red",
                                "#a2ebe8",
                                "pink",
                                "darkblue",
                                "#b5000d",
                                "pink",
                                "darkblue",
                                "orange",
                                "yellow",
                                "#aab500",
                                "#e7b9f8",
                                "#119955",
                                "black"
),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend

#Other species

Th1_stim_STAT1 <- Th1_stim$evol[54,]
Th1_stim_STAT4 <- Th1_stim$evol[55,]
Th1_stim_STAT5 <- Th1_stim$evol[56,]
Th1_stim_STAT6 <- Th1_stim$evol[57,]
Th1_stim_CDKN1A_p21 <- Th1_stim$evol[12,]
Th1_stim_BCL6 <- Th1_stim$evol[75,]
Th1_stim_AHR <- Th1_stim$evol[78,]
Th1_stim_FOXP1 <- Th1_stim$evol[51,]
Th1_stim_MAF <- Th1_stim$evol[77,]
Th1_stim_SPI1_PU1 <- Th1_stim$evol[76,]
Th1_stim_ZAP70 <- Th1_stim$evol[7,]
Th1_stim_AKT1 <- Th1_stim$evol[75,]
Th1_stim_ERK <- Th1_stim$evol[28,]
Th1_stim_JNK <- Th1_stim$evol[27,]
Th1_stim_CCNI <- Th1_stim$evol[77,]



behaviour_plot_Th1_others <-plot(c(1:201), Th1_stim_CDKN1A_p21, main="Other species evolution with IL2+IL4", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)

lines(c(1:201), Th1_stim_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), Th1_stim_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), Th1_stim_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th1_stim_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)

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
    
),cex=0.6)
