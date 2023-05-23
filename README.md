# Transmission time

Code for: "Inferring viral transmission time from phylogenies for known transmission pairs"

* `trees/` : example tree file  
   `exampletree.ttn` from the coalescent simulator in [biophybreak](https://github.com/MolEvolEpid/biophybreak)
* `fit/` : fit the time-slice model  
  `python3 run_slice.py ../trees/` creates `trees/mk2/exampletree-mk2.csv`
* `acc/` : other tools
  - phylogical window  
   `Rscript run_window.R ../trees/` creates `trees/window.csv`
  - donor-recipient topology  
   `Rscript run_topology.R ../trees/` creates `trees/topology.csv`
  - node-based rules  
   `Rscript run_rules.R ../trees/` creates `trees/rules.csv`

## Legal

FCI Open Source Copyright Assertion: C20001

(C) 2021. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.  Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
