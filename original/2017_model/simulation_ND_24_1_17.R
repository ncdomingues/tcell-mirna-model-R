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



#----------------------


#---------- SIMULATIONS WITH EXPANDED Naldi's MODEL ------------------------

Naldi_nodes <- read.delim("Naldi_2010_nodes_23_1_17.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
Naldi_rules <- read.delim("Naldi rules_no inputs_next.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)


Naldi_IL2 <- ginsimrun(state_IL2_Naldi, Naldi_nodes, Naldi_rules, 50, 1000)
Naldi_Th1 <- ginsimrun(state_Th1_Naldi, Naldi_nodes, Naldi_rules, 50, 100)
Naldi_Th2 <- ginsimrun(state_Th2_Naldi, Naldi_nodes, Naldi_rules, 50, 1000)
Naldi_Th17 <- ginsimrun(state_Th17_Naldi, Naldi_nodes, Naldi_rules, 50, 500)
Naldi_iTreg <- ginsimrun(state_iTreg_Naldi, Naldi_nodes, Naldi_rules, 50, 100)


# Naldi simulations

#IL2 stimulation
state_IL2_Naldi <- c(1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

length(state_IL2_Naldi)

#Th1 stimulation (IL2 + IL12 + APC)
state_Th1_Naldi <- c(1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Th2 stimulation (IL2 + IL4 + APC)
state_Th2_Naldi <- c(1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Th17 stimulation (IL2 + IL1 + IL6 + IL23 + TGFB + APC)
state_Th17_Naldi <- c(1,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#iTreg stimulation (IL2 + TGFB + APC)
state_iTreg_Naldi <- c(1,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#---------------

Naldi_Th17_155 <- Naldi_Th17$evol[67,] #miR-155-5p
Naldi_Th17_34c <- Naldi_Th17$evol[66,] #miR-34c-5p
Naldi_Th17_MDM2 <- Naldi_Th17$evol[80,] #MDM2
Naldi_Th17_TP53 <- Naldi_Th17$evol[71,] #TP53
Naldi_Th17_TBX21 <- Naldi_Th17$evol[51,] #TBX21
Naldi_Th17_STAT3 <- Naldi_Th17$evol[54,] #STAT3
Naldi_Th17_RORGT <- Naldi_Th17$evol[64,] #RORGT

behaviour_plot <- plot(c(1:51), Naldi_Th17_155, main="miR-155-5p ; miR-34c-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,50))
points(c(1:51), Naldi_Th17_34c, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,50))
points(c(1:51), Naldi_Th17_FOXO3, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,50))
points(c(1:51), Naldi_Th17_TP53, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="yellow", xlim=c(0,50))
points(c(1:51), Naldi_Th17_TBX21, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,50))
points(c(1:51), Naldi_Th17_GATA3, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,50))
points(c(1:51), Naldi_Th17_STAT3, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,50))
points(c(1:51), Naldi_Th17_RORGT, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="black", xlim=c(0,50))

#-------------------

Naldi_IL2_155 <- Naldi_IL2$evol[67,] #miR-155-5p
Naldi_IL2_34c <- Naldi_IL2$evol[66,] #miR-34c-5p
Naldi_IL2_FOXO3 <- Naldi_IL2$evol[70,] #FOXO3
Naldi_IL2_TP53 <- Naldi_IL2$evol[71,] #TP53
Naldi_IL2_MDM2 <- Naldi_IL2$evol[80,] #MDM2
Naldi_IL2_GATA3 <- Naldi_IL2$evol[50,] #GATA3
Naldi_IL2_TBX21 <- Naldi_IL2$evol[51,] #TBX21
Naldi_IL2_IKB <- Naldi_IL2$evol[63,] #TBX21


behaviour_plot_IL2 <- plot(c(1:51), Naldi_IL2_155, main="miR-155-5p ; miR-34c-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="blue", xlim=c(0,50))
points(c(1:51), Naldi_IL2_34c, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="green", xlim=c(0,50))
points(c(1:51), Naldi_IL2_IKB, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="cyan", xlim=c(0,50))

