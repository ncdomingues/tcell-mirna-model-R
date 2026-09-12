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


#______________________

STATE0 <- c(1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)

STATE0[23]=1 #Turn on BCM complex (BCL10 + CARD11 + MALT1)

STATE0[23]=1 #Turn on Th1 stimulation with IL12 + APC + CD28 + CD3
STATE0[42]=1 #Turn on Th2 timulation with IL4
nodes <- read.delim("nodes_values.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
rules <- read.delim("formulae.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)

Th2_stim <- ginsimrun(STATE0, nodes, rules, 100, 200)
STATE0
model_new
model_to_plot <- model_new$evol
model_to_plot
new_plot <- plot.data.frame(c(1:101), model_new$evol, main="Species evolution with IL2", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,101), ylim=c(0,2), lwd=1)
length(model_new$evol)
plot(model_new$evol, main="Species evolution with IL2", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,101), ylim=c(0,2), lwd=1)
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
behaviour_plot_IL2 <-plot(c(1:101), model_to_plot, main="Species evolution with IL2", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,101), ylim=c(0,2), lwd=1)
lines(c(1:101), model_new_34c, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="#007788", xlim=c(0,101))
lines(c(1:101), model_new_IKB, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,101))
lines(c(1:101), model_new_FOXO3, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,101))
lines(c(1:101), model_new_TP53, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,101))
lines(c(1:101), model_new_TBX21, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="#234543", xlim=c(0,101))
lines(c(1:101), model_new_GATA3, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="#456789", xlim=c(0,101))
lines(c(1:101), model_new_STAT3, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,101))
lines(c(1:101), model_new_RORC, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="purple", xlim=c(0,101))
points(c(1:101), model_new_FOS, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="grey", xlim=c(0,101))
lines(c(1:101), model_new_JUN, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="red", xlim=c(0,101))
lines(c(1:101), model_new_RAC1, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="orange", xlim=c(0,101))
lines(c(1:101), model_new_FOXO1, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,101))
lines(c(1:101), model_new_SIRT1, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,101))
lines(c(1:101), model_new_TGFBR, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="#006655", xlim=c(0,101))
lines(c(1:101), model_new_FOXP3, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="#119955", xlim=c(0,101))
lines(c(1:101), model_new_PRDM1_BLIMP1, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,101))
lines(c(1:101), model_new_MYC, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="#009955", xlim=c(0,101))
lines(c(1:101), model_new_NFKB1, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="pink", xlim=c(0,101))
lines(c(1:101), model_new_NFATC1, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="darkblue", xlim=c(0,101))

IL2
IL6
IL12
IL21
IL17
IL23
IL4
IL10
IL2R
IL6R
IL21R
IL23R
IL4R
IL10R
IL22
IL9
IL12R
IL1
IL1R
TGFB
TGFBR
TNF
TNFR

legend(86,2,c("miR-155-5p","34c","IKB","FOXO3","TP53","TBX21","GATA3","STAT3","RORC","FOS","JUN","RAC1","FOXO1","SIRT1","TGFBR","FOXP3","PRDM1_BLIMP","MYC","NFKB1","NFATC1"
),lty=c(1,1),lwd=c(3,3),col = c("blue", "#007788","cyan", "green", "yellow", "#234543", "#456789", "black", "purple", "grey", "red", "orange", "yellow", "#009955", "#006655", "#119955", "black", "#009955", "pink","darkblue"),cex=0.6) # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend


plot.new()
