library(Sim.DiffProc)

for(demo in demo(package="Sim.DiffProc")$results[,"Item"]) 
  demo(demo, package="Sim.DiffProc", character.only=TRUE)
