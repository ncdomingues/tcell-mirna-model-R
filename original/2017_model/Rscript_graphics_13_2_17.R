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





#IL2 stimulation (+anti CD28 CD3)


STATE0 <- c(1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

STATE0[23]=1 #Turn on BCM complex (BCL10 + CARD11 + MALT1)
STATE0[37]=1 #Turn on Th1 stimulation with IL12 + APC + CD28 + CD3
STATE0[42]=1 #Turn on Th2 timulation with IL4

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

legend(186,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend


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




length(model_new_155)
behaviour_plot_IL2 <-plot(c(1:201), model_new_155, main="Species evolution with IL2", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=1)
lines(c(1:201), model_new_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201),lwd=2)
lines(c(1:201), model_new_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201),lwd=2)


#Th2 stimulation


Th2_stim_155 <- Th2_stim$evol[85,] #miR-155-5p
Th2_stim_34c <- Th2_stim$evol[84,] #miR-34c-5p
Th2_stim_FOXO3 <- Th2_stim$evol[46,] #FOXO3
Th2_stim_TP53 <- Th2_stim$evol[47,] #TP53
Th2_stim_MDM2 <- Th2_stim$evol[48,] #MDM2
Th2_stim_GATA3 <- Th2_stim$evol[72,] #GATA3
Th2_stim_TBX21 <- Th2_stim$evol[71,] #TBX21
Th2_stim_IKB <- Th2_stim$evol[25,] #NFKBIE
Th2_stim_STAT3 <- Th2_stim$evol[58,] #STAT3
Th2_stim_RORC <- Th2_stim$evol[73,] #RORC
Th2_stim_FOS <- Th2_stim$evol[17,]
Th2_stim_JUN <- Th2_stim$evol[16,] 
Th2_stim_RAC1 <- Th2_stim$evol[13,] 
Th2_stim_FOXO1 <- Th2_stim$evol[18,] 
Th2_stim_SIRT1 <- Th2_stim$evol[83,]
Th2_stim_TGFBR <- Th2_stim$evol[61,] 
Th2_stim_FOXP3 <- Th2_stim$evol[74,]
Th2_stim_PRDM1_BLIMP1 <- Th2_stim$evol[82,]
Th2_stim_MYC <- Th2_stim$evol[50,]
Th2_stim_NFKB1 <- Th2_stim$evol[24,] 
Th2_stim_NFATC1 <- Th2_stim$evol[22,]




length(Th2_stim_155)
behaviour_plot_IL2 <-plot(c(1:201), Th2_stim_155, main="Species evolution with IL2+IL4", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), Th2_stim_34c, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IKB, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_FOXO3, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_TP53, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_TBX21, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_GATA3, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_STAT3, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_RORC, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_FOS, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_JUN, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_RAC1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_FOXO1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_SIRT1, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_FOXP3, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_PRDM1_BLIMP1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_MYC, xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_NFKB1, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_NFATC1, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)


legend(169.5,2,c("miR-155-5p","34c","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend


plot.new()


Th2_stim_IL2 <- Th2_stim$evol[34,] 
Th2_stim_IL6 <- Th2_stim$evol[35,] 
Th2_stim_IL12 <- Th2_stim$evol[37,] 
Th2_stim_IL21 <- Th2_stim$evol[38,]
Th2_stim_IL17 <- Th2_stim$evol[39,] 
Th2_stim_IL23 <- Th2_stim$evol[41,] 
Th2_stim_IL4 <- Th2_stim$evol[42,] 
Th2_stim_IL10 <- Th2_stim$evol[44,] 
Th2_stim_IL2R <- Th2_stim$evol[59,]
Th2_stim_IL6R <- Th2_stim$evol[60,]
Th2_stim_IL21R <- Th2_stim$evol[63,]
Th2_stim_IL23R <- Th2_stim$evol[65,] 
Th2_stim_IL4R <- Th2_stim$evol[66,] 
Th2_stim_IL10R <- Th2_stim$evol[67,] 
Th2_stim_IL22 <- Th2_stim$evol[69,]
Th2_stim_IL9 <- Th2_stim$evol[70,] 
Th2_stim_IL12R <- Th2_stim$evol[79,]
Th2_stim_IL1 <- Th2_stim$evol[80,]
Th2_stim_IL1R <- Th2_stim$evol[81,]
Th2_stim_TGFB <- Th2_stim$evol[36,] 
Th2_stim_TGFBR <- Th2_stim$evol[61,]
Th2_stim_TNF <- Th2_stim$evol[40,] 
Th2_stim_TNFR <- Th2_stim$evol[64,]

cytokines_plot_IL2 <-plot(c(1:201), Th2_stim_IL2,main="Cytokine evolution with IL2+IL4", xlab="Time", ylab="Relative quantity (a.u.)", col="#4d0000", xlim=c(0,201), ylim=c(0,2), lwd=2)
lines(c(1:201), Th2_stim_IL6, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d1300", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), Th2_stim_IL12, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d2600", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL21, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d3900", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL17, xlab="Time", ylab="Relative quantity (a.u.)", col="#4d4d00", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL23, xlab="Time", ylab="Relative quantity (a.u.)", col="#394d00", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL4, xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL10, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL2R, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL6R, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL21R, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL23R, xlab="Time", ylab="Relative quantity (a.u.)", col="#a2ebe8", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL4R, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL10R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL22, xlab="Time", ylab="Relative quantity (a.u.)", col="#b5000d", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL9, xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL12R, xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL1, xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_IL1R, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_TGFB, xlab="Time", ylab="Relative quantity (a.u.)", col="#aab500", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_TGFBR, xlab="Time", ylab="Relative quantity (a.u.)", col="#e7b9f8", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_TNF, xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_TNFR, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)

legend(169.8,2,c("IL2","IL6","IL12","IL21","IL17","IL23","IL4","IL10","IL2R","IL6R","IL21R","IL23R","IL4R","IL10R","IL22","IL9","IL12R","IL1","IL1R","TGFB","TGFBR","TNF","TNFR"
                 
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue","darkgreen","magenta","brown"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend


behaviour_plot_Th2_others <-plot(c(1:201), Th2_stim_CDKN1A_p21, main="Other species evolution with IL2+IL4", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,201), ylim=c(0,2), lwd=2)

lines(c(1:201), Th2_stim_BCL6, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_AHR, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_FOXP1, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_MAF, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_SPI1_PU1, xlab="Time", ylab="Relative quantity (a.u.)", col="#ab45ab", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_STAT1, xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_STAT4, xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_STAT5, xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_STAT6, xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,201), lwd=2)

lines(c(1:201), Th2_stim_ZAP70, xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,201), lwd=2, lwd=2)
lines(c(1:201), Th2_stim_AKT1, xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_ERK, xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_JNK, xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,201), lwd=2)
lines(c(1:201), Th2_stim_CCNI, xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,201), lwd=2)

legend(186,2,c("miR-155-5p","miR-34c-5p","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
