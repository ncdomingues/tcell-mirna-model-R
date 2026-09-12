#SIMULATIONS (GINsim + R)

getordervar=function(vars){
    #order variables according to string length
    #this avoids matching errors when using the gsub function  
    ncharvec=nchar(vars[,2])
    ordervar=order(ncharvec,decreasing=T)
    ordervar
}


ginsimnxt=function(state,vars,rules,ordervar){
    nvar=length(state)
    syncchange=rep(0,nvar)
    ruleset=rules[,4]
    #substitute var names in the rule strings by "0" or "1" according
    #to var values in the state vector
    for (i in 1:nvar){
        varvalue=state[ordervar[i]]
        varname=vars[ordervar[i],2]
        maxvalue=vars[ordervar[i],3]
        if (maxvalue==1){
            #if variable is boolean, it is enough to substitute var name
            #by its value in the state vector
            ruleset=gsub(varname,as.character(varvalue),ruleset)
        } else {
            #if variable has more levels, it is necessary to substitue 
            #first "varname:value" strings by "0" or "1" according to 
            #value in state vector
            aa=seq(maxvalue,1,-1)
            for (j in 1:(length(aa)-1)){
                tomatch=paste(varname,":",aa[j],sep="")
                if (varvalue>=aa[j]){
                    tosub="1"
                } else {
                    tosub="0"
                }
                ruleset=gsub(tomatch,tosub,ruleset)
            }
            #and then substitute varname by "0" or "1" according to 
            # value in state vector
            if (varvalue>=1){
                ruleset=gsub(varname,"1",ruleset)
            } else {
                ruleset=gsub(varname,"0",ruleset)
            }
        }
    }
    #solve the logic formulae in the rule strings
    solvedrules=vector()
    for (i in 1:length(ruleset)){
        solvedrules[i]=eval(parse(text=ruleset[i]))
    }
    #compute potential next state for each variable
    for (i in 1:nvar){
        #for each variable, select results of rules that target them
        TorF=solvedrules[rules[,1]==i]
        if (length(TorF)>0){
            #if there are rules that target this var, see target values
            rulesetlevel=rules[rules[,1]==i,3]
            if (sum(TorF)==0){
                #if there are no true rules next value is zero
                syncchange[i]=0
            } else {
                #if there are true rules, next value is the max target value
                # among true rules
                syncchange[i]=max(rulesetlevel[TorF])
            }
        } else {
            #if there are no rules that target this var, next state is 
            # the actual state (as for input nodes)
            syncchange[i]=state[i]
        }  
    }
    #find wich vars can change on the next state
    changed=which(state!=syncchange)
    if (length(changed)>0){
        #if there are vars to change
        if (length(changed)==1){
            # if there is only one var to change
            tochange=changed
        } else {
            # if there are more than one var to change, randomly pick one
            tochange=sample(changed,1)  
        }
        # change state of the chosen (or unique) variable
        if (syncchange[tochange]>state[tochange]){
            # if target value is greater than actual state, increase var
            #by 1
            state[tochange]=state[tochange]+1
        } else {
            #if target value is lower than actual state, decrease var 
            #by 1
            state[tochange]=state[tochange]-1
        }
    }
    #output is modified (or same) state vector
    state
}


ginsimrun=function(state, vars, rules,ntime,nrep){
    ordervar=getordervar(vars)
    nvar=length(state)
    statecum=matrix(0,nrow=nvar,ncol=ntime+1)
    steadymat=matrix(state,nrow=1,ncol=nvar)
    
    for (i in 1:nrep){
        statemat=matrix(0,nrow=nvar,ncol=ntime+1)
        statemat[,1]=state
        steady=0
        for (j in 1:ntime){
            if (steady==0){
                statemat[,j+1]=ginsimnxt(statemat[,j],vars,rules,ordervar)
                if (all(statemat[,j+1]==statemat[,j])){
                    steady=1
                    steadymat=rbind(steadymat,statemat[,j])
                }
            } else {
                statemat[,j+1]=statemat[,j]
            }
        }
        statecum=statecum+statemat
    }
    statecum=statecum/nrep
    list(evol=statecum, ss=steadymat)
    # output$evol is a matrix where each row is a variable and each
    #column is a time point. Values are averages of the nrep simulations
    #
    #output$ss is a matrix with found steady states - each row is a 
    #steady state and each column is a variable. The first row is the 
    #initial state of the simulation and not a steady state. There may 
    #be repeated steady states
}




nodes <- read.delim("nodes_values.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules <- read.delim("formulae.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)


STATE0 <- c(1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

STATE0[23]=1 #Turn on BCM complex (BCL10 + CARD11 + MALT1)
STATE0[37]=1 #Turn on Th1 stimulation with IL12
STATE0[42]=1 #Turn on Th2 stimulation with IL4

STATE0[41]=0 #Turn off IL23 stimulation
STATE0[80]=0 #Turn off IL1 stimulation

nodes <- read.delim("nodes_values.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules <- read.delim("formulae.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)

TFs_IL2_15_2_17 <- ginsimrun(STATE0, nodes, rules, 50, 200)



TFs_IL2_15_2_17_155 <- TFs_IL2_15_2_17$evol[85,] #miR-155-5p
TFs_IL2_15_2_17_34c <- TFs_IL2_15_2_17$evol[84,] #miR-34c-5p
TFs_IL2_15_2_17_FOXO3 <- TFs_IL2_15_2_17$evol[46,] #FOXO3
TFs_IL2_15_2_17_TP53 <- TFs_IL2_15_2_17$evol[47,] #TP53
TFs_IL2_15_2_17_MDM2 <- TFs_IL2_15_2_17$evol[48,] #MDM2
TFs_IL2_15_2_17_GATA3 <- TFs_IL2_15_2_17$evol[72,] #GATA3
TFs_IL2_15_2_17_TBX21 <- TFs_IL2_15_2_17$evol[71,] #TBX21
TFs_IL2_15_2_17_IKB <- TFs_IL2_15_2_17$evol[25,] #NFKBIE
TFs_IL2_15_2_17_STAT3 <- TFs_IL2_15_2_17$evol[58,] #STAT3
TFs_IL2_15_2_17_RORC <- TFs_IL2_15_2_17$evol[73,] #RORC
TFs_IL2_15_2_17_FOS <- TFs_IL2_15_2_17$evol[17,]
TFs_IL2_15_2_17_JUN <- TFs_IL2_15_2_17$evol[16,] 
TFs_IL2_15_2_17_RAC1 <- TFs_IL2_15_2_17$evol[13,] 
TFs_IL2_15_2_17_FOXO1 <- TFs_IL2_15_2_17$evol[18,] 
TFs_IL2_15_2_17_SIRT1 <- TFs_IL2_15_2_17$evol[83,]
TFs_IL2_15_2_17_TGFBR <- TFs_IL2_15_2_17$evol[61,] 
TFs_IL2_15_2_17_FOXP3 <- TFs_IL2_15_2_17$evol[74,]
TFs_IL2_15_2_17_PRDM1_BLIMP1 <- TFs_IL2_15_2_17$evol[82,]
TFs_IL2_15_2_17_MYC <- TFs_IL2_15_2_17$evol[50,]
TFs_IL2_15_2_17_NFKB1 <- TFs_IL2_15_2_17$evol[24,] 
TFs_IL2_15_2_17_NFATC1 <- TFs_IL2_15_2_17$evol[22,]

#Other species

TFs_IL2_15_2_17_STAT1 <- TFs_IL2_15_2_17$evol[54,]
TFs_IL2_15_2_17_STAT4 <- TFs_IL2_15_2_17$evol[55,]
TFs_IL2_15_2_17_STAT5 <- TFs_IL2_15_2_17$evol[56,]
TFs_IL2_15_2_17_STAT6 <- TFs_IL2_15_2_17$evol[57,]
TFs_IL2_15_2_17_CDKN1A_p21 <- TFs_IL2_15_2_17$evol[12,]
TFs_IL2_15_2_17_BCL6 <- TFs_IL2_15_2_17$evol[75,]
TFs_IL2_15_2_17_AHR <- TFs_IL2_15_2_17$evol[78,]
TFs_IL2_15_2_17_FOXP1 <- TFs_IL2_15_2_17$evol[51,]
TFs_IL2_15_2_17_MAF <- TFs_IL2_15_2_17$evol[77,]
TFs_IL2_15_2_17_SPI1_PU1 <- TFs_IL2_15_2_17$evol[76,]
TFs_IL2_15_2_17_ZAP70 <- TFs_IL2_15_2_17$evol[7,]
TFs_IL2_15_2_17_AKT1 <- TFs_IL2_15_2_17$evol[75,]
TFs_IL2_15_2_17_ERK <- TFs_IL2_15_2_17$evol[28,]
TFs_IL2_15_2_17_JNK <- TFs_IL2_15_2_17$evol[27,]
TFs_IL2_15_2_17_CCNI <- TFs_IL2_15_2_17$evol[77,]




plot_IL2_TF0 <-plot(c(1:201), TFs_IL2_15_2_17_155, main="Species evolution with IL2", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)

legend(169.8,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend



lines(c(1:201), TFs_IL2_15_2_17_CDKN1A_p21, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), TFs_IL2_15_2_17_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)


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

TFs_IL2_15_2_17_IL2 <- TFs_IL2_15_2_17$evol[34,] #miR-155-5p
TFs_IL2_15_2_17_IL6 <- TFs_IL2_15_2_17$evol[35,] #miR-34c-5p
TFs_IL2_15_2_17_IL12 <- TFs_IL2_15_2_17$evol[37,] #FOXO3
TFs_IL2_15_2_17_IL21 <- TFs_IL2_15_2_17$evol[38,] #TP53
TFs_IL2_15_2_17_IL17 <- TFs_IL2_15_2_17$evol[39,] #MDM2
TFs_IL2_15_2_17_IL23 <- TFs_IL2_15_2_17$evol[41,] #GATA3
TFs_IL2_15_2_17_IL4 <- TFs_IL2_15_2_17$evol[42,] #TBX21
TFs_IL2_15_2_17_IL10 <- TFs_IL2_15_2_17$evol[44,] #NFKBIE
TFs_IL2_15_2_17_IL2R <- TFs_IL2_15_2_17$evol[59,] #STAT3
TFs_IL2_15_2_17_IL6R <- TFs_IL2_15_2_17$evol[60,] #RORC
TFs_IL2_15_2_17_IL21R <- TFs_IL2_15_2_17$evol[63,]
TFs_IL2_15_2_17_IL23R <- TFs_IL2_15_2_17$evol[65,] 
TFs_IL2_15_2_17_IL4R <- TFs_IL2_15_2_17$evol[66,] 
TFs_IL2_15_2_17_IL10R <- TFs_IL2_15_2_17$evol[67,] 
TFs_IL2_15_2_17_IL22 <- TFs_IL2_15_2_17$evol[69,]
TFs_IL2_15_2_17_IL9 <- TFs_IL2_15_2_17$evol[70,] 
TFs_IL2_15_2_17_IL12R <- TFs_IL2_15_2_17$evol[79,]
TFs_IL2_15_2_17_IL1 <- TFs_IL2_15_2_17$evol[80,]
TFs_IL2_15_2_17_IL1R <- TFs_IL2_15_2_17$evol[81,]
TFs_IL2_15_2_17_TGFB <- TFs_IL2_15_2_17$evol[36,] 
TFs_IL2_15_2_17_TGFBR <- TFs_IL2_15_2_17$evol[61,]
TFs_IL2_15_2_17_TNF <- TFs_IL2_15_2_17$evol[40,] 
TFs_IL2_15_2_17_TNFR <- TFs_IL2_15_2_17$evol[64,]

cytokines_plot_IL2 <-plot(c(1:201), TFs_IL2_15_2_17_IL2, main="Cytokiner evolution with IL2", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), TFs_IL2_15_2_17_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)


#miR-34c TF combinations

STATE0 <- c(1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

STATE0[23]=1 #Turn on BCM complex (BCL10 + CARD11 + MALT1)
STATE0[37]=0 #Turn on/off Th1 stimulation with IL12
STATE0[42]=1 #Turn on Th2 timulation with IL4
STATE0[41]=0 #Turn off IL23 stimulation
STATE0[80]=0 #Turn off IL1 stimulation

nodes <- read.delim("nodes_values.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)


rules <- read.delim("formulae_16_2_17.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_MYC_TP53_FOXO3 <- read.delim("formulae_MYC_TP53_FOXO3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_GATA3_MYC_TP53_FOXO3 <- read.delim("formulae_GATA3_MYC_TP53_FOXO3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_GATA3_TP53_FOXO3 <- read.delim("formulae_GATA3_TP53_FOXO3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_GATA3_FOS_MYC <- read.delim("formulae_GATA3_FOS_MYC.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_TP53_FOXO3 <- read.delim("formulae_TP53_FOXO3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)


TF0 <- ginsimrun(STATE0, nodes, rules, 200, 200)

TF1 <- ginsimrun(STATE0, nodes, rules_MYC_TP53_FOXO3, 200, 200)

TF2 <- ginsimrun(STATE0, nodes, rules_GATA3_MYC_TP53_FOXO3, 200, 200)

TF3 <- ginsimrun(STATE0, nodes, rules_GATA3_TP53_FOXO3, 200, 200)

TF4 <- ginsimrun(STATE0, nodes, rules_GATA3_FOS_MYC, 200, 200)

TF5 <- ginsimrun(STATE0, nodes, rules_TP53_FOXO3, 200, 200)


#TF0 PLOT

TF0_155 <- TF0$evol[85,] #miR-155-5p
TF0_34c <- TF0$evol[84,] #miR-34c-5p
TF0_FOXO3 <- TF0$evol[46,] #FOXO3
TF0_TP53 <- TF0$evol[47,] #TP53
TF0_MDM2 <- TF0$evol[48,] #MDM2
TF0_GATA3 <- TF0$evol[72,] #GATA3
TF0_TBX21 <- TF0$evol[71,] #TBX21
TF0_IKB <- TF0$evol[25,] #NFKBIE
TF0_STAT3 <- TF0$evol[58,] #STAT3
TF0_RORC <- TF0$evol[73,] #RORC
TF0_FOS <- TF0$evol[17,]
TF0_JUN <- TF0$evol[16,] 
TF0_RAC1 <- TF0$evol[13,] 
TF0_FOXO1 <- TF0$evol[18,] 
TF0_SIRT1 <- TF0$evol[83,]
TF0_TGFBR <- TF0$evol[61,] 
TF0_FOXP3 <- TF0$evol[74,]
TF0_PRDM1_BLIMP1 <- TF0$evol[82,]
TF0_MYC <- TF0$evol[50,]
TF0_NFKB1 <- TF0$evol[24,] 
TF0_NFATC1 <- TF0$evol[22,]


TF0_plot_IL2 <-plot(c(1:201), TF0_155, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF0_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)

#Other species

TF0_STAT1 <- TF0$evol[54,]
TF0_STAT4 <- TF0$evol[55,]
TF0_STAT5 <- TF0$evol[56,]
TF0_STAT6 <- TF0$evol[57,]
TF0_CDKN1A_p21 <- TF0$evol[12,]
TF0_BCL6 <- TF0$evol[75,]
TF0_AHR <- TF0$evol[78,]
TF0_FOXP1 <- TF0$evol[51,]
TF0_MAF <- TF0$evol[77,]
TF0_SPI1_PU1 <- TF0$evol[76,]
TF0_ZAP70 <- TF0$evol[7,]
TF0_AKT1 <- TF0$evol[75,]
TF0_ERK <- TF0$evol[28,]
TF0_JNK <- TF0$evol[27,]
TF0_CCNI <- TF0$evol[77,]





legend(169.8,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend



plot_TF0_others <- plot(c(1:201), TF0_CDKN1A_p21, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF0_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), TF0_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)


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

TF0_IL2 <- TF0$evol[34,]
TF0_IL6 <- TF0$evol[35,]
TF0_IL12 <- TF0$evol[37,]
TF0_IL21 <- TF0$evol[38,] 
TF0_IL17 <- TF0$evol[39,] 
TF0_IL23 <- TF0$evol[41,] 
TF0_IL4 <- TF0$evol[42,] 
TF0_IL10 <- TF0$evol[44,] 
TF0_IL2R <- TF0$evol[59,] 
TF0_IL6R <- TF0$evol[60,] 
TF0_IL21R <- TF0$evol[63,]
TF0_IL23R <- TF0$evol[65,] 
TF0_IL4R <- TF0$evol[66,] 
TF0_IL10R <- TF0$evol[67,] 
TF0_IL22 <- TF0$evol[69,]
TF0_IL9 <- TF0$evol[70,] 
TF0_IL12R <- TF0$evol[79,]
TF0_IL1 <- TF0$evol[80,]
TF0_IL1R <- TF0$evol[81,]
TF0_TGFB <- TF0$evol[36,] 
TF0_TGFBR <- TF0$evol[61,]
TF0_TNF <- TF0$evol[40,] 
TF0_TNFR <- TF0$evol[64,]

plot_IL2_TF0_cytokynes <-plot(c(1:201), TF0_IL2, main="IL2 stimulation: cytokynes", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF0_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2)
lines(c(1:201), TF0_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), TF0_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)

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
),cex=0.6) 


#TF1 PLOT


TF1_155 <- TF1$evol[85,] #miR-155-5p
TF1_34c <- TF1$evol[84,] #miR-34c-5p
TF1_FOXO3 <- TF1$evol[46,] #FOXO3
TF1_TP53 <- TF1$evol[47,] #TP53
TF1_MDM2 <- TF1$evol[48,] #MDM2
TF1_GATA3 <- TF1$evol[72,] #GATA3
TF1_TBX21 <- TF1$evol[71,] #TBX21
TF1_IKB <- TF1$evol[25,] #NFKBIE
TF1_STAT3 <- TF1$evol[58,] #STAT3
TF1_RORC <- TF1$evol[73,] #RORC
TF1_FOS <- TF1$evol[17,]
TF1_JUN <- TF1$evol[16,] 
TF1_RAC1 <- TF1$evol[13,] 
TF1_FOXO1 <- TF1$evol[18,] 
TF1_SIRT1 <- TF1$evol[83,]
TF1_TGFBR <- TF1$evol[61,] 
TF1_FOXP3 <- TF1$evol[74,]
TF1_PRDM1_BLIMP1 <- TF1$evol[82,]
TF1_MYC <- TF1$evol[50,]
TF1_NFKB1 <- TF1$evol[24,] 
TF1_NFATC1 <- TF1$evol[22,]


TF1_plot_IL2 <-plot(c(1:201), TF1_155, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF1_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)

#Other species

TF1_STAT1 <- TF1$evol[54,]
TF1_STAT4 <- TF1$evol[55,]
TF1_STAT5 <- TF1$evol[56,]
TF1_STAT6 <- TF1$evol[57,]
TF1_CDKN1A_p21 <- TF1$evol[12,]
TF1_BCL6 <- TF1$evol[75,]
TF1_AHR <- TF1$evol[78,]
TF1_FOXP1 <- TF1$evol[51,]
TF1_MAF <- TF1$evol[77,]
TF1_SPI1_PU1 <- TF1$evol[76,]
TF1_ZAP70 <- TF1$evol[7,]
TF1_AKT1 <- TF1$evol[75,]
TF1_ERK <- TF1$evol[28,]
TF1_JNK <- TF1$evol[27,]
TF1_CCNI <- TF1$evol[77,]





legend(169.8,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend



plot_TF1_others <- plot(c(1:201), TF1_CDKN1A_p21, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF1_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), TF1_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)


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

TF1_IL2 <- TF1$evol[34,]
TF1_IL6 <- TF1$evol[35,]
TF1_IL12 <- TF1$evol[37,]
TF1_IL21 <- TF1$evol[38,] 
TF1_IL17 <- TF1$evol[39,] 
TF1_IL23 <- TF1$evol[41,] 
TF1_IL4 <- TF1$evol[42,] 
TF1_IL10 <- TF1$evol[44,] 
TF1_IL2R <- TF1$evol[59,] 
TF1_IL6R <- TF1$evol[60,] 
TF1_IL21R <- TF1$evol[63,]
TF1_IL23R <- TF1$evol[65,] 
TF1_IL4R <- TF1$evol[66,] 
TF1_IL10R <- TF1$evol[67,] 
TF1_IL22 <- TF1$evol[69,]
TF1_IL9 <- TF1$evol[70,] 
TF1_IL12R <- TF1$evol[79,]
TF1_IL1 <- TF1$evol[80,]
TF1_IL1R <- TF1$evol[81,]
TF1_TGFB <- TF1$evol[36,] 
TF1_TGFBR <- TF1$evol[61,]
TF1_TNF <- TF1$evol[40,] 
TF1_TNFR <- TF1$evol[64,]

plot_IL2_TF1_cytokynes <-plot(c(1:201), TF1_IL2, main="IL2 stimulation: cytokynes", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF1_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2)
lines(c(1:201), TF1_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), TF1_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)

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
),cex=0.6) 

#TF2 PLOT


TF2_155 <- TF2$evol[85,] #miR-155-5p
TF2_34c <- TF2$evol[84,] #miR-34c-5p
TF2_FOXO3 <- TF2$evol[46,] #FOXO3
TF2_TP53 <- TF2$evol[47,] #TP53
TF2_MDM2 <- TF2$evol[48,] #MDM2
TF2_GATA3 <- TF2$evol[72,] #GATA3
TF2_TBX21 <- TF2$evol[71,] #TBX21
TF2_IKB <- TF2$evol[25,] #NFKBIE
TF2_STAT3 <- TF2$evol[58,] #STAT3
TF2_RORC <- TF2$evol[73,] #RORC
TF2_FOS <- TF2$evol[17,]
TF2_JUN <- TF2$evol[16,] 
TF2_RAC1 <- TF2$evol[13,] 
TF2_FOXO1 <- TF2$evol[18,] 
TF2_SIRT1 <- TF2$evol[83,]
TF2_TGFBR <- TF2$evol[61,] 
TF2_FOXP3 <- TF2$evol[74,]
TF2_PRDM1_BLIMP1 <- TF2$evol[82,]
TF2_MYC <- TF2$evol[50,]
TF2_NFKB1 <- TF2$evol[24,] 
TF2_NFATC1 <- TF2$evol[22,]


TF2_plot_IL2 <-plot(c(1:201), TF2_155, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF2_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)

#Other species

TF2_STAT1 <- TF2$evol[54,]
TF2_STAT4 <- TF2$evol[55,]
TF2_STAT5 <- TF2$evol[56,]
TF2_STAT6 <- TF2$evol[57,]
TF2_CDKN1A_p21 <- TF2$evol[12,]
TF2_BCL6 <- TF2$evol[75,]
TF2_AHR <- TF2$evol[78,]
TF2_FOXP1 <- TF2$evol[51,]
TF2_MAF <- TF2$evol[77,]
TF2_SPI1_PU1 <- TF2$evol[76,]
TF2_ZAP70 <- TF2$evol[7,]
TF2_AKT1 <- TF2$evol[75,]
TF2_ERK <- TF2$evol[28,]
TF2_JNK <- TF2$evol[27,]
TF2_CCNI <- TF2$evol[77,]





legend(169.8,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend



plot_TF2_others <- plot(c(1:201), TF2_CDKN1A_p21, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF2_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), TF2_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)


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

TF2_IL2 <- TF2$evol[34,]
TF2_IL6 <- TF2$evol[35,]
TF2_IL12 <- TF2$evol[37,]
TF2_IL21 <- TF2$evol[38,] 
TF2_IL17 <- TF2$evol[39,] 
TF2_IL23 <- TF2$evol[41,] 
TF2_IL4 <- TF2$evol[42,] 
TF2_IL10 <- TF2$evol[44,] 
TF2_IL2R <- TF2$evol[59,] 
TF2_IL6R <- TF2$evol[60,] 
TF2_IL21R <- TF2$evol[63,]
TF2_IL23R <- TF2$evol[65,] 
TF2_IL4R <- TF2$evol[66,] 
TF2_IL10R <- TF2$evol[67,] 
TF2_IL22 <- TF2$evol[69,]
TF2_IL9 <- TF2$evol[70,] 
TF2_IL12R <- TF2$evol[79,]
TF2_IL1 <- TF2$evol[80,]
TF2_IL1R <- TF2$evol[81,]
TF2_TGFB <- TF2$evol[36,] 
TF2_TGFBR <- TF2$evol[61,]
TF2_TNF <- TF2$evol[40,] 
TF2_TNFR <- TF2$evol[64,]

plot_IL2_TF2_cytokynes <-plot(c(1:201), TF2_IL2, main="IL2 stimulation: cytokynes", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF2_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2)
lines(c(1:201), TF2_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), TF2_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)

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
),cex=0.6) 

#TF3


TF3_155 <- TF3$evol[85,] #miR-155-5p
TF3_34c <- TF3$evol[84,] #miR-34c-5p
TF3_FOXO3 <- TF3$evol[46,] #FOXO3
TF3_TP53 <- TF3$evol[47,] #TP53
TF3_MDM2 <- TF3$evol[48,] #MDM2
TF3_GATA3 <- TF3$evol[72,] #GATA3
TF3_TBX21 <- TF3$evol[71,] #TBX21
TF3_IKB <- TF3$evol[25,] #NFKBIE
TF3_STAT3 <- TF3$evol[58,] #STAT3
TF3_RORC <- TF3$evol[73,] #RORC
TF3_FOS <- TF3$evol[17,]
TF3_JUN <- TF3$evol[16,] 
TF3_RAC1 <- TF3$evol[13,] 
TF3_FOXO1 <- TF3$evol[18,] 
TF3_SIRT1 <- TF3$evol[83,]
TF3_TGFBR <- TF3$evol[61,] 
TF3_FOXP3 <- TF3$evol[74,]
TF3_PRDM1_BLIMP1 <- TF3$evol[82,]
TF3_MYC <- TF3$evol[50,]
TF3_NFKB1 <- TF3$evol[24,] 
TF3_NFATC1 <- TF3$evol[22,]


TF3_plot_IL2 <-plot(c(1:201), TF3_155, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF3_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)

#Other species

TF3_STAT1 <- TF3$evol[54,]
TF3_STAT4 <- TF3$evol[55,]
TF3_STAT5 <- TF3$evol[56,]
TF3_STAT6 <- TF3$evol[57,]
TF3_CDKN1A_p21 <- TF3$evol[12,]
TF3_BCL6 <- TF3$evol[75,]
TF3_AHR <- TF3$evol[78,]
TF3_FOXP1 <- TF3$evol[51,]
TF3_MAF <- TF3$evol[77,]
TF3_SPI1_PU1 <- TF3$evol[76,]
TF3_ZAP70 <- TF3$evol[7,]
TF3_AKT1 <- TF3$evol[75,]
TF3_ERK <- TF3$evol[28,]
TF3_JNK <- TF3$evol[27,]
TF3_CCNI <- TF3$evol[77,]





legend(169.8,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend



plot_TF3_others <- plot(c(1:201), TF3_CDKN1A_p21, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF3_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), TF3_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)


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

TF3_IL2 <- TF3$evol[34,]
TF3_IL6 <- TF3$evol[35,]
TF3_IL12 <- TF3$evol[37,]
TF3_IL21 <- TF3$evol[38,] 
TF3_IL17 <- TF3$evol[39,] 
TF3_IL23 <- TF3$evol[41,] 
TF3_IL4 <- TF3$evol[42,] 
TF3_IL10 <- TF3$evol[44,] 
TF3_IL2R <- TF3$evol[59,] 
TF3_IL6R <- TF3$evol[60,] 
TF3_IL21R <- TF3$evol[63,]
TF3_IL23R <- TF3$evol[65,] 
TF3_IL4R <- TF3$evol[66,] 
TF3_IL10R <- TF3$evol[67,] 
TF3_IL22 <- TF3$evol[69,]
TF3_IL9 <- TF3$evol[70,] 
TF3_IL12R <- TF3$evol[79,]
TF3_IL1 <- TF3$evol[80,]
TF3_IL1R <- TF3$evol[81,]
TF3_TGFB <- TF3$evol[36,] 
TF3_TGFBR <- TF3$evol[61,]
TF3_TNF <- TF3$evol[40,] 
TF3_TNFR <- TF3$evol[64,]

plot_IL2_TF3_cytokynes <-plot(c(1:201), TF3_IL2, main="IL2 stimulation: cytokynes", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF3_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2)
lines(c(1:201), TF3_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), TF3_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)

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
),cex=0.6) 

#TF4 PLOT


TF4_155 <- TF4$evol[85,] #miR-155-5p
TF4_34c <- TF4$evol[84,] #miR-34c-5p
TF4_FOXO3 <- TF4$evol[46,] #FOXO3
TF4_TP53 <- TF4$evol[47,] #TP53
TF4_MDM2 <- TF4$evol[48,] #MDM2
TF4_GATA3 <- TF4$evol[72,] #GATA3
TF4_TBX21 <- TF4$evol[71,] #TBX21
TF4_IKB <- TF4$evol[25,] #NFKBIE
TF4_STAT3 <- TF4$evol[58,] #STAT3
TF4_RORC <- TF4$evol[73,] #RORC
TF4_FOS <- TF4$evol[17,]
TF4_JUN <- TF4$evol[16,] 
TF4_RAC1 <- TF4$evol[13,] 
TF4_FOXO1 <- TF4$evol[18,] 
TF4_SIRT1 <- TF4$evol[83,]
TF4_TGFBR <- TF4$evol[61,] 
TF4_FOXP3 <- TF4$evol[74,]
TF4_PRDM1_BLIMP1 <- TF4$evol[82,]
TF4_MYC <- TF4$evol[50,]
TF4_NFKB1 <- TF4$evol[24,] 
TF4_NFATC1 <- TF4$evol[22,]


TF4_plot_IL2 <-plot(c(1:201), TF4_155, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF4_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)

#Other species

TF4_STAT1 <- TF4$evol[54,]
TF4_STAT4 <- TF4$evol[55,]
TF4_STAT5 <- TF4$evol[56,]
TF4_STAT6 <- TF4$evol[57,]
TF4_CDKN1A_p21 <- TF4$evol[12,]
TF4_BCL6 <- TF4$evol[75,]
TF4_AHR <- TF4$evol[78,]
TF4_FOXP1 <- TF4$evol[51,]
TF4_MAF <- TF4$evol[77,]
TF4_SPI1_PU1 <- TF4$evol[76,]
TF4_ZAP70 <- TF4$evol[7,]
TF4_AKT1 <- TF4$evol[75,]
TF4_ERK <- TF4$evol[28,]
TF4_JNK <- TF4$evol[27,]
TF4_CCNI <- TF4$evol[77,]





legend(169.8,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend



plot_TF4_others <- plot(c(1:201), TF4_CDKN1A_p21, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF4_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), TF4_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)


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

TF4_IL2 <- TF4$evol[34,]
TF4_IL6 <- TF4$evol[35,]
TF4_IL12 <- TF4$evol[37,]
TF4_IL21 <- TF4$evol[38,] 
TF4_IL17 <- TF4$evol[39,] 
TF4_IL23 <- TF4$evol[41,] 
TF4_IL4 <- TF4$evol[42,] 
TF4_IL10 <- TF4$evol[44,] 
TF4_IL2R <- TF4$evol[59,] 
TF4_IL6R <- TF4$evol[60,] 
TF4_IL21R <- TF4$evol[63,]
TF4_IL23R <- TF4$evol[65,] 
TF4_IL4R <- TF4$evol[66,] 
TF4_IL10R <- TF4$evol[67,] 
TF4_IL22 <- TF4$evol[69,]
TF4_IL9 <- TF4$evol[70,] 
TF4_IL12R <- TF4$evol[79,]
TF4_IL1 <- TF4$evol[80,]
TF4_IL1R <- TF4$evol[81,]
TF4_TGFB <- TF4$evol[36,] 
TF4_TGFBR <- TF4$evol[61,]
TF4_TNF <- TF4$evol[40,] 
TF4_TNFR <- TF4$evol[64,]

plot_IL2_TF4_cytokynes <-plot(c(1:201), TF4_IL2, main="IL2 stimulation: cytokynes", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF4_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2)
lines(c(1:201), TF4_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), TF4_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)

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
),cex=0.6) 

#TF5 PLOT


TF5_155 <- TF5$evol[85,] #miR-155-5p
TF5_34c <- TF5$evol[84,] #miR-34c-5p
TF5_FOXO3 <- TF5$evol[46,] #FOXO3
TF5_TP53 <- TF5$evol[47,] #TP53
TF5_MDM2 <- TF5$evol[48,] #MDM2
TF5_GATA3 <- TF5$evol[72,] #GATA3
TF5_TBX21 <- TF5$evol[71,] #TBX21
TF5_IKB <- TF5$evol[25,] #NFKBIE
TF5_STAT3 <- TF5$evol[58,] #STAT3
TF5_RORC <- TF5$evol[73,] #RORC
TF5_FOS <- TF5$evol[17,]
TF5_JUN <- TF5$evol[16,] 
TF5_RAC1 <- TF5$evol[13,] 
TF5_FOXO1 <- TF5$evol[18,] 
TF5_SIRT1 <- TF5$evol[83,]
TF5_TGFBR <- TF5$evol[61,] 
TF5_FOXP3 <- TF5$evol[74,]
TF5_PRDM1_BLIMP1 <- TF5$evol[82,]
TF5_MYC <- TF5$evol[50,]
TF5_NFKB1 <- TF5$evol[24,] 
TF5_NFATC1 <- TF5$evol[22,]


TF5_plot_IL2 <-plot(c(1:201), TF5_155, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF5_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)

#Other species

TF5_STAT1 <- TF5$evol[54,]
TF5_STAT4 <- TF5$evol[55,]
TF5_STAT5 <- TF5$evol[56,]
TF5_STAT6 <- TF5$evol[57,]
TF5_CDKN1A_p21 <- TF5$evol[12,]
TF5_BCL6 <- TF5$evol[75,]
TF5_AHR <- TF5$evol[78,]
TF5_FOXP1 <- TF5$evol[51,]
TF5_MAF <- TF5$evol[77,]
TF5_SPI1_PU1 <- TF5$evol[76,]
TF5_ZAP70 <- TF5$evol[7,]
TF5_AKT1 <- TF5$evol[75,]
TF5_ERK <- TF5$evol[28,]
TF5_JNK <- TF5$evol[27,]
TF5_CCNI <- TF5$evol[77,]





legend(169.8,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend



plot_TF5_others <- plot(c(1:201), TF5_CDKN1A_p21, main="IL2 stimulation", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF5_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), TF5_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)


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

TF5_IL2 <- TF5$evol[34,]
TF5_IL6 <- TF5$evol[35,]
TF5_IL12 <- TF5$evol[37,]
TF5_IL21 <- TF5$evol[38,] 
TF5_IL17 <- TF5$evol[39,] 
TF5_IL23 <- TF5$evol[41,] 
TF5_IL4 <- TF5$evol[42,] 
TF5_IL10 <- TF5$evol[44,] 
TF5_IL2R <- TF5$evol[59,] 
TF5_IL6R <- TF5$evol[60,] 
TF5_IL21R <- TF5$evol[63,]
TF5_IL23R <- TF5$evol[65,] 
TF5_IL4R <- TF5$evol[66,] 
TF5_IL10R <- TF5$evol[67,] 
TF5_IL22 <- TF5$evol[69,]
TF5_IL9 <- TF5$evol[70,] 
TF5_IL12R <- TF5$evol[79,]
TF5_IL1 <- TF5$evol[80,]
TF5_IL1R <- TF5$evol[81,]
TF5_TGFB <- TF5$evol[36,] 
TF5_TGFBR <- TF5$evol[61,]
TF5_TNF <- TF5$evol[40,] 
TF5_TNFR <- TF5$evol[64,]

plot_IL2_TF5_cytokynes <-plot(c(1:201), TF5_IL2, main="IL2 stimulation: cytokynes", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), TF5_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2)
lines(c(1:201), TF5_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), TF5_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)

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
),cex=0.6) 





