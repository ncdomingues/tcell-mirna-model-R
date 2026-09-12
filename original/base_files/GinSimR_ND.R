
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

#---------- SIMULATIONS WITH EXPANDED Abou-Jaoudé's MODEL ------------------------


AJ_nodes <- read.delim("AJ_2015_nodes_19_1_17.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
AJ_rules <- read.delim("AJ rules_no inputs_19_1_17_use.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)

# AJ stimulations
#If necessary change length of vector and set values to 1 if not regulated


#IL2 stimulation

state_IL2_AJ <- c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

str(state_IL2_AJ)
length(state_IL2_AJ)
state_IL2_AJ[21]

#Th1 stimulation
state_Th1_AJ <- c(0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Th2 stimulation
state_Th2_AJ <- c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Th17 stimulation
state_Th17_AJ <- c(0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0)

#iTreg stimulation
state_iTreg_AJ <- c(0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0)


IL2 <- ginsimrun(state_IL2_AJ, AJ_nodes, AJ_rules, 100, 100)
Th1 <- ginsimrun(state_Th1_AJ, AJ_nodes, AJ_rules, 50, 100)
Th2 <- ginsimrun(state_Th2_AJ, AJ_nodes, AJ_rules, 50, 1000)
Th17 <- ginsimrun(state_Th17_AJ, AJ_nodes, AJ_rules, 1000, 10000)
iTreg <- ginsimrun(state_iTreg_AJ, AJ_nodes, AJ_rules, 50, 1000)



#---------- SIMULATIONS WITH EXPANDED Naldi's MODEL ------------------------

Naldi_nodes <- read.delim("Naldi_2010_nodes_19_1_17.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
Naldi_rules <- read.delim("Naldi rules_no inputs_19_1_17.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)


Naldi_IL2 <- ginsimrun(state_IL2_Naldi, Naldi_nodes, Naldi_rules, 200, 200)
Naldi_Th1 <- ginsimrun(state_Th1_Naldi, Naldi_nodes, Naldi_rules, 50, 1000)
Naldi_Th2 <- ginsimrun(state_Th2_Naldi, Naldi_nodes, Naldi_rules, 50, 1000)
Naldi_Th17 <- ginsimrun(state_Th17_Naldi, Naldi_nodes, Naldi_rules, 1000, 20000)
Naldi_iTreg <- ginsimrun(state_iTreg_Naldi, Naldi_nodes, Naldi_rules, 50, 1000)


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

#---------------------------------------------------------------------

#Data search

#Naldi

mir155_Naldi_IL2 <- Naldi_IL2$evol[67,]
mir155_Naldi_IL2_ss <- Naldi_IL2$ss[,67]


mir155_Naldi_Th1 <- Naldi_Th1$evol[67,]
mir155_Naldi_Th17 <- Naldi_Th17$evol[67,]
mir155_Naldi_Th17_ss <- Naldi_Th17$ss[,67]



mir155_AJ_IL2 <- IL2$evol[103,]
mir155_AJ_IL2_ss <- IL2$ss[,103]


mir155_AJ_Th1 <- Th1$evol[103,]
mir155_AJ_Th17 <- Th17$evol[103,]
mir155_AJ_Th17_ss <- Th17$ss[,103]



Naldi_Th17_b <- ginsimrun(state_Th17_Naldi, Naldi_nodes, Naldi_rules, 50, 1000)

Naldi_IL2_155 <- Naldi_IL2$evol[67,] #miR-155-5p
Naldi_IL2_34c <- Naldi_IL2$evol[66,] #miR-34c-5p
Naldi_IL2_FOXO3 <- Naldi_IL2$evol[70,] #FOXO3
Naldi_IL2_TP53 <- Naldi_IL2$evol[71,] #TP53
Naldi_IL2_MDM2 <- Naldi_IL2$evol[80,] #MDM2

Naldi_Th17_155 <- Naldi_Th17_b$evol[67,] #miR-155-5p
Naldi_Th17_34c <- Naldi_Th17_b$evol[66,] #miR-34c-5p
Naldi_Th17_MDM2 <- Naldi_Th17_b$evol[80,] #MDM2

Naldi_Th1_155 <- Naldi_Th1$evol[67,] #miR-155-5p
Naldi_Th1_34c <- Naldi_Th1$evol[66,] #miR-34c-5p

Naldi_Th2_155 <- Naldi_Th1$evol[67,] #miR-155-5p
Naldi_Th2_34c <- Naldi_Th1$evol[66,] #miR-34c-5p

Naldi_iTreg_155 <- Naldi_iTreg$evol[67,] #miR-155-5p
Naldi_iTreg_34c <- Naldi_iTreg$evol[66,] #miR-34c-5p

#PLOTS------------------------

time=c(1:51)
plot(Naldi_Th17_155, time, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="blue")
points(Naldi_Th17_MDM2, time, main="miR-155-5p", xlab="Time", ylab="Relative quantity (a.u.)", col="red")




