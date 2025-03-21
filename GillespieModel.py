import numpy as np
import matplotlib.pyplot as plt

class Reactant:
  
  def __init__(self,n,label=''):
    self.n = n
    self.initial_n = n
    self.label = label

  def reset(self):
    self.n = self.initial_n

  def set_n(self,n):
    self.n = n
    
  def delta_n(self,d):
    self.n += d
    
class Reaction:
  
  def __init__(self,params,all_reactants=None):
    # assemble the reactant list - can either be Reactant objects or labels
    # if the label can't be found, a ValueError is raised
    self.reactants = []
    for i in range(len(params['reactants'])):
      if (isinstance(params['reactants'][i],Reactant)):
        self.reactants.append(params['reactants'][i])
      else:
        reactant_not_found = True # checks if we didn't find the reactant
        if all_reactants is None:
          raise ValueError('When specifying reactants by label, all_reactants must be given')
        for r in all_reactants:
          if r.label == params['reactants'][i]:
            reactant_not_found = False
            self.reactants.append(r)
        if reactant_not_found:
          raise ValueError('Reactant not found')

    self.stoichiometry_matrix = params['stoichiometry matrix'] # vector, same length as reactants
    self.rate = params['rate']
    self.volume = params['volume']
    # if the sensitivity list isn't provided, create it from the stoichiometry matrix
    if 'sensitivity list' in params:
      self.sensitivity_list = params['sensitivity list']
    else:
      self.sensitivity_list = []
      for i in range(len(params['stoichiometry matrix'])):
        if params['stoichiometry matrix'][i] < 0:
          self.sensitivity_list += abs(params['stoichiometry matrix'][i]) * [i]
    
  def get_rate(self):
    r = self.rate
    for i in range(len(self.sensitivity_list)):
      reactant_index = self.sensitivity_list[i]
      r *= self.reactants[reactant_index].n / self.volume
    return r

  def set_rate_parameter(self,rate):
    self.rate = rate
  
  def do_reaction(self):
    for i in range(len(self.stoichiometry_matrix)):
      self.reactants[i].delta_n(self.stoichiometry_matrix[i])
    
class GillespieModel:
  
  def __init__(self,params):
    self.reactants = params['reactants']
    self.reactions = params['reactions']
    self.t = 0
    self.history = []
    self.add_to_history()

  def reset(self):
    for r in self.reactants:
      r.reset()
    self.t = 0
    self.history = []
    self.add_to_history()
    
  def add_to_history(self):
    self.history.append(self.get_state_and_time())
    
  def get_state_and_time(self):
    return {'s':np.array([r.n for r in self.reactants]),'t':self.t}
  
  def iterate(self):
    rate_v = [r.get_rate() for r in self.reactions]
    r0 = np.sum(rate_v)
    if r0 > 0:
      # perform the Gillespie simulation
      (v1,v2) = np.random.rand(2)
      tau = 1/r0 * np.log(1/v1) # next reaction occurs after a wait of tau
      react = np.sum(np.cumsum(rate_v)/r0 <= v2) # identity of the reaction (starting at 0)
      self.reactions[react].do_reaction()
      self.t += tau      
      self.add_to_history()
      return True
    return False
  
  def simulate(self,max_iter,reset=False):
    if reset:
      self.reset()
    i = 0
    while self.iterate() and i < max_iter:
      i += 1
    return i

  # Simulate up to the given time
  def simulateToTime(self,end_time,reset=False):
    if reset:
      self.reset()
    at_least_one_iteration = False
    keep_going = True
    while (self.t <= end_time) and keep_going:
      at_least_one_iteration = True
      last_state = self.get_state_and_time()['s']
      keep_going = self.iterate()

    # if at least one iteration happened, rewind
    # this depends on the only modifiable state being recorded in last_state,
    # and also the reactants being in the same order as in last_state (should be fine in the list comprehension)
    if at_least_one_iteration:
      for i in range(len(last_state)):
        self.reactants[i].n = last_state[i]
      self.history = self.history[:-1]
      
    self.t = end_time # no matter what, we end at the end_time
    self.add_to_history() # make sure the history ends at the correct time - the state won't change but the time will
    return self.t

  # returns all event times as a numpy array
  def getTimeVector(self):
    return np.array([h['t'] for h in self.history])
  
  # returns the state at each event as a numpy array
  # if you don't specify an index, returns everything; 
  # if you specify an index, returns the state for only that reactant
  def getStateVector(self,index=None):
    if index is None:
      allRows = []
      for i in range(len(self.history[0]['s'])):
        allRows.append(np.array([h['s'][i] for h in self.history]))
      return np.vstack(allRows)
    return np.array([h['s'][index] for h in self.history])   
  
  # Plots all (or specified by index) reactants
  # You can add custom colors if you want
  # Legend comes from the labels associated with each reactant
  # todo: specify reactants by label
  def plot(self,reactants=None,colors=None,legend=True):
    t = self.getTimeVector()
    m = self.getStateVector()
    if reactants is None:
      reactants = range(np.shape(m)[0])
    for i in reactants:
      if colors is None:
        plt.plot(t,m[i,:],label=self.reactants[i].label)
      else:
        plt.plot(t,m[i,:],color=colors[i],label=self.reactants[i].label)
    if legend:
      plt.legend()
    ax = plt.gca()
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of molecules')

