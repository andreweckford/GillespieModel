# This example sets up a Gillespie simulation of a simple enzyme-substrate reaction.
# The reaction is: S+E <-> SE -> E+P
# where S is the substrate, E is the enzyme, SE is the substrate-enzyme intermediate, and P is the product
# The example produces two plots: one with a simulation up to time 400000, and a continuation up to time 800000 (arbitrary units)

import numpy as np
import matplotlib.pyplot as plt
from GillespieModel import *

##### parameters

# initial quantities for reactants
initial_substrate = 100
initial_enzyme = 10
initial_intermediate = 0
initial_product = 0

# rates
s_e_se_forward_rate = 0.1 # rate of the S+E->SE reaction
s_e_se_backward_rate = 0.1 # rate of the SE->S+E reaction
se_e_p_production_rate = 0.1 # rate of the SE->E+P reaction
s_creation_rate = 0 # rate of S creation, zero for now (doesn't happen) - set it to something other than zero to see how creation can affect the result
p_degradation_rate = 0 # rate of P degradation, zero for now (doesn't happen) - set it to something other than zero to see how degradation can affect the result

# other
vol = 1000 # volume of the container, used to calculate concentration

##### set up the simulator

# reactants - the string is the label
S = Reactant(initial_substrate,'Substrate')
E = Reactant(initial_enzyme,'Enzyme')
SE = Reactant(initial_intermediate,'Intermediate')
P = Reactant(initial_product,'Product')

reactants = [S, E, SE, P]

reactions = []

# forward reaction S+E->SE - parameters are explained in detail here
reactions.append(Reaction({
  #
  # Reactants involved in the reaction
  # Can be specified either using the reactant object itself, or by label
  # When specifying by label, you must provide the all_reactants parameter, below
  # In this case, the commented-out line should give the same behaviour, try it
  'reactants' : ['Substrate','Enzyme','Intermediate'],
  #'reactants' : [S,E,SE], 
  #
  # here, indicate the number of reactants involved in the reaction
  # this vector must be the same length as the 'reactants' vector, and the positions in the list correspond to those reactants
  # in this example, the reaction consumes one S and one E (reducing them by 1), and produces one SE (increasing it by 1)
  # Any number is possible here
  'stoichiometry matrix' : [-1,-1,1],
  #
  # rate constant - this is multiplied to any concentrations (see the sensitivity list) to obtain the reaction rate
  'rate' : s_e_se_forward_rate,
  #
  # volume - you need a volume to calculate concentration
  'volume' : vol
},all_reactants=reactants))

# For all the rest of the reactions, we specify the Reactant objects rather than labels
# In these cases, you don't need to provide the all_reactants parameter to Reaction
# backward reaction SE->S+E
reactions.append(Reaction({
  'reactants' : [S,E,SE],
  'stoichiometry matrix' : [1,1,-1],
  'rate' : s_e_se_backward_rate,
  'volume' : vol
}))

# production reaction SE->E+P
reactions.append(Reaction({
  'reactants' : [SE,E,P],
  'stoichiometry matrix' : [-1,1,1],
  'rate' : se_e_p_production_rate,
  'volume' : vol
}))

# substrate production
reactions.append(Reaction({
  'reactants' : [S],
  'stoichiometry matrix' : [1],
  'rate' : s_creation_rate,
  'volume' : vol
}))

# product degradation
reactions.append(Reaction({
  'reactants' : [P],
  'stoichiometry matrix' : [-1],
  'rate' : p_degradation_rate,
  'volume' : vol
}))


g = GillespieModel({'reactants':reactants,'reactions':reactions})

g.simulateToTime(4e6)

g.plot()
plt.show()

# continue the simulation to time 8e6 and plot, 
# preserving the state of the simulator from the first simulation
g.simulateToTime(8e6)

# if you want to START OVER and simulate to 8e6, pass reset=True as a parameter to the function
# as in the line below, commented out. If you uncomment the line below to try this,
# You can leave the g.SimulateToTime(8e6) line above, that simulation is overwritten with the version below
#g.simulateToTime(8e6,reset=True)

g.plot()
plt.show()
