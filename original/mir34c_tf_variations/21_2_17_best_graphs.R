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

#-----------------

#miR-34c TF combinations

STATE1 <- c(1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

STATE0[23]=1 #Turn on BCM complex (BCL10 + CARD11 + MALT1)
STATE0[37]=1 #Turn on/off Th1 stimulation with IL12
STATE0[42]=1 #Turn on Th2 timulation with IL4
STATE0[41]=0 #Turn off IL23 stimulation
STATE0[80]=0 #Turn off IL1 stimulation

STATE1 <- c(1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

STATE1[23]=1 #Turn on BCM complex (BCL10 + CARD11 + MALT1)
STATE1[37]=0 #Turn on/off Th1 stimulation with IL12
STATE1[42]=1 #Turn on Th2 timulation with IL4
STATE1[41]=0 #Turn off IL23 stimulation
STATE1[80]=0 #Turn off IL1 stimulation

nodes <- read.delim("nodes_values.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)


rules <- read.delim("formulae_23_2_17.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_MYC_TP53_FOXO3 <- read.delim("formulae_MYC_TP53_FOXO3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_GATA3_MYC_TP53_FOXO3 <- read.delim("formulae_GATA3_MYC_TP53_FOXO3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_GATA3_TP53_FOXO3 <- read.delim("formulae_GATA3_TP53_FOXO3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_GATA3_FOS_MYC <- read.delim("formulae_GATA3_FOS_MYC.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_TP53_FOXO3 <- read.delim("formulae_TP53_FOXO3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_GATA3 <- read.delim("formulae_GATA3.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_MYC <- read.delim("formulae_MYC.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules_MYCandGATA3

IL2 <- ginsimrun(STATE1, nodes, rules, 400, 200)

IL2+IL4 <- ginsimrun(STATE0, nodes, rules, 200, 200)

TF2 <- ginsimrun(STATE0, nodes, rules_GATA3, 200, 200)

TF3 <- ginsimrun(STATE0, nodes, rules_GATA3_TP53_FOXO3, 200, 200)

TF4 <- ginsimrun(STATE0, nodes, rules_GATA3_FOS_MYC, 200, 200)

TF5 <- ginsimrun(STATE0, nodes, rules_TP53_FOXO3, 200, 200)

time <- c(1:201)
mir155 <- TF0$evol[85,] #miR-155-5p
mir34c <- TF0$evol[84,] #miR-34c-5p
FOXO3 <- TF0$evol[46,] #FOXO3
TP53 <- TF0$evol[47,] #TP53
MDM2 <- TF0$evol[48,] #MDM2
GATA3 <- TF0$evol[72,] #GATA3
TBX21 <- TF0$evol[71,] #TBX21
IKB <- TF0$evol[25,] #NFKBIE
RORC <- TF0$evol[73,] #RORC
FOS <- TF0$evol[17,]
JUN <- TF0$evol[16,] 
RAC1 <- TF0$evol[13,] 
FOXO1 <- TF0$evol[18,] 
SIRT1 <- TF0$evol[83,]
TGFBR <- TF0$evol[61,] 
FOXP3 <- TF0$evol[74,]
PRDM1_BLIMP1 <- TF0$evol[82,]
MYC <- TF0$evol[50,]
NFKB1 <- TF0$evol[24,] 
NFATC1 <- TF0$evol[22,]
STAT1 <- TF0$evol[54,]
STAT3 <- TF0$evol[58,]
STAT4 <- TF0$evol[55,]
STAT5 <- TF0$evol[56,]
STAT6 <- TF0$evol[57,]
CDKN1A_p21 <- TF0$evol[12,]
BCL6 <- TF0$evol[75,]
AHR <- TF0$evol[78,]
FOXP1 <- TF0$evol[51,]
MAF <- TF0$evol[77,]
SPI1_PU1 <- TF0$evol[76,]
ZAP70 <- TF0$evol[7,]
AKT1 <- TF0$evol[75,]
ERK <- TF0$evol[28,]
JNK <- TF0$evol[27,]
CCNI <- TF0$evol[77,]
IL2 <- TF0$evol[34,]
IL6 <- TF0$evol[35,]
IL12 <- TF0$evol[37,]
IL21 <- TF0$evol[38,] 
IL17 <- TF0$evol[39,] 
IL23 <- TF0$evol[41,] 
IL4 <- TF0$evol[42,] 
IL10 <- TF0$evol[44,] 
IL2R <- TF0$evol[59,] 
IL6R <- TF0$evol[60,] 
IL21R <- TF0$evol[63,]
IL23R <- TF0$evol[65,] 
IL4R <- TF0$evol[66,] 
IL10R <- TF0$evol[67,] 
IL22 <- TF0$evol[69,]
IL9 <- TF0$evol[70,] 
IL12R <- TF0$evol[79,]
IL1 <- TF0$evol[80,]
IL1R <- TF0$evol[81,]
TGFB <- TF0$evol[36,] 
TGFBR <- TF0$evol[61,]
TNF <- TF0$evol[40,] 
TNFR <- TF0$evol[64,]

#ggplot2
library(ggplot2)
best <- qplot(time, to_plot, data=, color=, shape=, size=, alpha=, geom=, method=, formula=, facets=, xlim=, ylim= xlab=, ylab=, main=, sub=)



#plot.ly

library(plotly)
packageVersion('plotly')

to_plot <- data.frame(time, mir155, mir34c, GATA3, MYC)

#The default order will be alphabetized unless specified as below:
data$month <- factor(data$month, levels = data[["month"]])

r <- plot_ly(to_plot, x = time, y = mir155, name = 'miR-155-5p', type = 'scatter', mode = 'lines',
             line = list(color = 'rgb(205, 12, 24)', width = 2)) %>%
    add_trace(y = mir34c, name = 'miR-34c-5p', line = list(color = '#ff0044', width = 2)) %>%
    add_trace(y = GATA3, name = 'GATA3', line = list(color = 'rgb(205, 12, 24)', width = 2, dash = 'dash')) %>%
    add_trace(y = MYC, name = 'MYC', line = list(color = 'rgb(22, 96, 167)', width = 2, dash = 'dash')) %>%
    add_trace(y = FOS, name = 'FOS', line = list(color = '#555555', width = 2, dash = 'dash')) %>%
    add_trace(y = NFATC1, name = 'NFATC1', line = list(color = 'rgb(77, 77, 66)', width = 2, dash = 'dash')) %>%
    
    add_trace(y = TP53, name = 'TP53', line = list(color = 'rgb(99, 99, 99)', width = 2, dash = 'dash')) %>%
    add_trace(y = FOXO3, name = 'FOXO3', line = list(color = 'rgb(11, 11, 11)', width = 2, dash = 'dash')) %>%
    layout(title = "IL2 stimulation",
           xaxis = list(title = "Time (a.u.)"),
           yaxis = list (title = "Evolution (a.u.)"))


s <- plot_ly(to_plot, x = time, y = STAT1, name = 'STAT1', type = 'scatter', mode = 'lines',
             line = list(color = 'rgb(205, 12, 24)', width = 2)) %>%
    add_trace(y = STAT3, name = 'STAT3', line = list(color = '#ff0044', width = 2)) %>%
    add_trace(y = STAT4, name = 'STAT4', line = list(color = 'rgb(205, 12, 24)', width = 2, dash = 'dash')) %>%
    add_trace(y = STAT5, name = 'STAT5', line = list(color = 'rgb(22, 96, 167)', width = 2, dash = 'dash')) %>%
    add_trace(y = STAT5, name = 'STAT6', line = list(color = 'rgb(77, 77, 66)', width = 2, dash = 'dash')) %>%
    
    add_trace(y = FOXP3, name = 'FOXP3', line = list(color = 'rgb(99, 99, 99)', width = 2, dash = 'dash')) %>%
    add_trace(y = BCL6, name = 'BCL6', line = list(color = 'rgb(11, 11, 11)', width = 2, dash = 'dash')) %>%
    layout(title = "IL2 stimulation",
           xaxis = list(title = "Time (a.u.)"),
           yaxis = list (title = "Evolution (a.u.)"))


cytokines <- plot_ly(to_plot, x = time, y = IL2, name = 'IL2', type = 'scatter', mode = 'lines',
             line = list(color = 'rgb(205, 12, 24)', width = 2)) %>%
    add_trace(y = IL4, name = 'IL4', line = list(color = '#ff0044', width = 2)) %>%
    add_trace(y = TGFB, name = 'TGFB', line = list(color = 'rgb(205, 12, 24)', width = 2, dash = 'dash')) %>%
    add_trace(y = IL21, name = 'IL21', line = list(color = 'rgb(22, 96, 167)', width = 2, dash = 'dash')) %>%
    add_trace(y = IL17, name = 'IL17', line = list(color = 'rgb(77, 77, 66)', width = 2, dash = 'dash')) %>%
    
    add_trace(y = IL6, name = 'IL6', line = list(color = 'rgb(99, 99, 99)', width = 2)) %>%
    add_trace(y = IL10, name = 'IL10', line = list(color = 'rgb(11, 11, 11)', width = 2, dash = 'dash')) %>%
    layout(title = "IL2 stimulation",
           xaxis = list(title = "Time (a.u.)"),
           yaxis = list (title = "Evolution (a.u.)"))
t



master_TFs <- plot_ly(to_plot, x = time, y = TBX21, name = 'TBX21', type = 'scatter', mode = 'lines',
             line = list(color = 'rgb(205, 12, 24)', width = 2)) %>%
    add_trace(y = GATA3, name = 'GATA3', line = list(color = '#ff0044', width = 2)) %>%
    add_trace(y = RORC, name = 'RORC', line = list(color = 'rgb(205, 12, 24)', width = 2, dash = 'dash')) %>%
    add_trace(y = FOXP3, name = 'FOXP3', line = list(color = 'rgb(22, 96, 167)', width = 2, dash = 'dash')) %>%
    add_trace(y = AHR, name = 'AHR', line = list(color = 'rgb(77, 77, 66)', width = 2, dash = 'dash')) %>%
    
    add_trace(y = SPI1_PU1, name = 'SPI1_PU1', line = list(color = 'rgb(99, 99, 99)', width = 2)) %>%
    add_trace(y = BCL6, name = 'BCL6', line = list(color = 'rgb(11, 11, 11)', width = 2, dash = 'dash')) %>%
    layout(title = "IL2 stimulation",
           xaxis = list(title = "Time (a.u.)"),
           yaxis = list (title = "Evolution (a.u.)"))
master_TFs
