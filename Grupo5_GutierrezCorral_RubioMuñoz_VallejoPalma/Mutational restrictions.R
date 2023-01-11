#Libraries

library("OncoSimulR")

#Aplication of the mutational restrictions to the fitness of each strain

dfat <- data.frame(Genotype = c("WT", "B", "A", "B, A"),
                   Fitness = c("1.1 - 0.01*n_/N + n_A*0.01*0.5/N + n_B*0.01*0.5/N",
                               "1.05 - 0.01*n_B/N + n_*0.01*0.5/N + n_A_B*0.01*0.5/N",
                               "1.05 - 0.01*n_A/N + n_*0.01*0.5/N + n_A_B*0.01*0.5/N",
                               "1.01 - 0.01*n_A_B/N + n_B*0.01*0.5/N + n_A*0.01*0.5/N"))

afe2 <- allFitnessEffects(genotFitness = dfat,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

#Variable creations. Treatment that will act as an antibiotic mixture

variables_modelo <- list(
  list(Name = "treatment_1",
       Value = 0
  )
)

v_Model <- createUserVars(variables_modelo) 

# List of rules that simulate adaptative therapy

rules <- list(
  list(ID = "rule_1",
       Condition = "n_ <= 500",
       Action = "treatment_1 = 0"
  ),list(ID = "rule_5",
         Condition = "n_ >= 1000",
         Action = "treatment_1 = 1"
  )
)

#Interventions simulating treatment use

intervenciones <- list(
  list(ID="ANTIBIOTICO 1 y 2 SOBRE WT",
       Trigger       = "treatment_1 == 1",
       WhatHappens   = "n_ = n_ - 0.3 * n_ - 0.3 * n_",
       Periodicity   = 0.07,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO 1 y 2 SOBRE A",
       Trigger       = "treatment_1 == 1",
       WhatHappens   = "n_A = n_A - 0.3* n_A * 0.2 - 0.3 * n_A",
       Periodicity   = 0.07,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO 1 y 2 SOBRE B",
       Trigger       = "treatment_1 == 1",
       WhatHappens   = "n_B = n_B - 0.3 * n_B - 0.3 * n_B * 0.2",
       Periodicity   = 0.07,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO 1 y 2 SOBRE A_B",
       Trigger       = "treatment_1 == 1",
       WhatHappens   = "n_A_B = n_A_B - 0.3 * n_A_B * 0.2 - 0.3 * n_A_B * 0.2",
       Periodicity   = 0.07,
       Repetitions   = Inf)
  
)
rules <- createRules(rules, afe2)

inter2 <- createInterventions(intervenciones, afe2)

#The initial mutational rate is a dummy mutational rate because we are
#implementing it in the fitness
evalAllGenotypes(afe2,spPopSizes = c(900, 33, 33, 33))
initSize = c(1000, 30, 30, 30)
simu2 <- oncoSimulIndiv(afe2,
                        initMutant = c("WT", "B", "A", "B, A"),
                        initSize = initSize,
                        finalTime = 1000,
                        mu=0.00000000001,
                        userVars = v_Model,
                        interventions = inter2,
                        rules=rules,
                        keepEvery = 1)
plot(simu2,show="genotypes")

#tests to check if the adaptative therapy is working as intended

pobs <- unlist(simu2$pops.by.time)[,2:5]
totalpob <- rowSums(unlist(simu2$pops.by.time)[,2:5])

#first test checking if any strain is growing out of control

stopifnot(totalpob<(length(initSize)*sum(initSize)))

#Second test checking if the original strain is still alive

stopifnot(pobs[, 1]>0)



