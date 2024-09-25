# libraries
import pandas as pd, numpy as np
import cobra
import os
from cobra import Model, Reaction, Metabolite

# Load the model - read_sbml, load_matlab_model, etc,.
model = cobra.io.load_matlab_model("/content/Actinomyces_odontolyticus_ATCC_17982.mat")

#diet
diet = pd.read_csv('/content/WesternDietAGORA2.txt', sep = '\t', header = None, names = ['exchange', 'lb'])

#check model details
print('Number of reactions:', len(model.reactions))
print('Number of metabolites:', len(model.metabolites))
print('Number of genes:', len(model.genes))
print('Objective function:', model.objective)
[(reaction.bounds, reaction) for reaction in model.reactions[:5]]

#objective function
for var in model.objective.variables:
    print(var, var.primal)
print(model.objective)

#find biomass rxn
for rxn in model.reactions:
    if "biomass" in rxn.id:
      print (rxn.id)
    #rxn.lower_bound = 0
    #rxn.upper_bound = 1000
print('Done!')

#model.use_diet
exch_rxns_dict = dict(zip(diet['exchange'], diet['lb']))
for rxn in model.reactions:
  if 'EX_' in rxn.id:
    if rxn.id in exch_rxns_dict.keys():
      #print('Reaction:', rxn.id)
      #print('Original lb:', rxn.lower_bound)
      #print('Diet lb:', exch_rxns_dict[rxn.id])
      rxn.lower_bound = exch_rxns_dict[rxn.id]
      #print('New lb:', rxn.lower_bound)
    else:
      rxn.lower_bound = 0
print('Done setting the diet!')

#change objective
model.objective = 'biomass525'
model.optimize()