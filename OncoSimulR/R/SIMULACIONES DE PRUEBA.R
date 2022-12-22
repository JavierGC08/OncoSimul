
dfat <- data.frame(Genotype = c("WT", "B", "A", "B, A"),
                   Fitness = c("0*n_+1.5",
                               "1.2",
                               "1.2",
                               "1.0"))

afe2 <- allFitnessEffects(genotFitness = dfat,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

evalAllGenotypes(afe2,spPopSizes = c(900,33,33,33))


intervenciones <- list(
  list(ID="ANTIBIOTICO SOBRE WT",
       Trigger       = "TRUE",
       WhatHappens   = "n_ = (n_ -x*n_ -y*n_) -u*n_ + n_A*u*0.5  + n_B*u*0.5",
       Periodicity   = 1,
       Repetitions   = Inf),
 list(ID="ANTIBIOTICO SOBRE B",
      Trigger       = "TRUE",
      WhatHappens   = "n_B = (n_B - x*n_B-y*n_B*0.2) - u*n_B + n_*u*0.5  + n_A_B*u*0.5",
      Periodicity   = 1,
      Repetitions   = Inf),
 list(ID="ANTIBIOTICO SOBRE A",
      Trigger       = "TRUE",
      WhatHappens   = "n_A = (n_A -x*n_A*0.2-y*n_A) - u*n_A + n_*u*0.5  + n_A_B*u*0.5",
      Periodicity   = 1,
      Repetitions   = Inf),
 list(ID="ANTIBIOTICO SOBRE A,B",
      Trigger       = "TRUE",
      WhatHappens   = "n_A_B = (n_A_B - x*n_A_B*0.2-y*n_A_B*0.2) - u*n_A_B + n_B*u*0.5  + n_A_B*u*0.5",
      Periodicity   = 1,
      Repetitions   = Inf))



variables_modelo <- list(
  list(Name = "a",
       Value = 0.4
  ),
  list(Name = "b",
       Value       = 0.5
  ),
  list(Name = "x",
       Value = 0.4
  ),
  list(Name = "y",
       Value = 0.3),
  list(Name = "u",
       Value= 0.0001),
  list(Name = "user_var1",
       Value= 0.0001)
  )

rules <- list(
  list(ID = "rule_1",
       Condition = "(n_B + n_A_B) > n_",
       Action = "y = 0"
  ),list(ID = "rule_2",
         Condition = "(n_A + n_A_B) > n_",
         Action = "x = 0"
  ),list(ID = "rule_3",
         Condition = "(n_A + n_A_B) < n_",
         Action = "x = 0.4"
  ),list(ID = "rule_4",
         Condition = "(n_B + n_A_B) < n_",
         Action = "y = 0.3"
  )
)

rules <- createRules(rules, afe2)

v_Model <- createUserVars(variables_modelo)  
                        
inter2 <- createInterventions(intervenciones,afe2)

simu2 <- oncoSimulIndiv(afe2,
               initMutant = c("WT", "B", "A", "B, A"),
               initSize = c(900,300,300,300),
               finalTime = 50,
               interventions = inter2,
               mu=0.00000000001,
               userVars = v_Model,
               rules=rules,
               keepEvery = 1)
plot(simu2,show="genotypes")

pobs <- unlist(simu2$pops.by.time)[,2:5]
totalpob <- rowSums(unlist(simu2$pops.by.time)[,2:5])
freqs <- pobs/totalpob
time <- unlist(simu2$pops.by.time)[,1]
max <- totalpob/totalpob
plot(max~time,type="l", col="black",ylim=c(0,1.1),ylab="FREQUENCY")
lines(time,freqs[,1],type="l",col="#A6761D")
lines(time,freqs[,2],type="l",col="#666666")
lines(time,freqs[,3],type="l",col="#1B9E89")
lines(time,freqs[,4],type="l",col="red")

unlist(simu2$other)
a= unlist(simu2$other$userVarValues)
a
