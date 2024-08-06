
# how to write the propensity matrix
#propensities = [
#    lambda mrna, protein, k1 : mrna+protein*k1,
#    lambda mrna, protein, k2 : protein + ...,
#    ...
#] 

import inspect
import numpy as np
import matplotlib.pyplot as plt
import math


class CRN():

    def __init__(self, stoichiometric_matrix, species_names, parameters_names, reaction_names, propensities):
        """
        Create a CRN object

        :param stoichiometric_matrix: numpy array
        :param species_names: list of strings
        :param parameters_names: list of strings
        :param reaction_names: list of strings
        :param propensities:
        """
        self.stoichiometric_matrix = np.array(stoichiometric_matrix) # numpy object
        
        self.species_names = species_names # list of strings
        self.reaction_names = reaction_names # list of strings
        # save the ordering for further use
        self.species_ordering = { self.species_names[i]:i for i in range(len(species_names)) }
        self.reaction_ordering = { self.reaction_names[i]: i for i in range(len(reaction_names))}

        # TODO update s.m. so that self.species_names labels the rows
        self.parameters_names = parameters_names # list of strings
        self.propensities = propensities # list of functions
        self.propensities_dependencies = [ inspect.getfullargspec(p).args for p in self.propensities ] # [[strings]]

        # define state and parameters 
        self.state = np.zeros( [len(self.species_names)] ) # the state is a numpy array
        self.parameters = {} # the set of parameters is just a dictionary


    # incomplete stub for calling the propensity functions on the right parametrs
    # self.propensities[i](*[  self.propensities_dependencies ])

    def _eval(self):
        """
        Evaluate propensity functions on the current state
        """
        out = []
        for prop in self.propensities:
            args = inspect.getfullargspec(prop).args
            input_args = {} # input args of this specific propensity
            for a in args:
                try:
                    input_args[a] = self.parameters[a]
                except KeyError:
                    try:
                        input_args[a] = self.state[self.species_ordering[a]]
                    except KeyError:
                        raise Exception(f'argument {a} is nor a parameter or a species')
            out.append(prop(**input_args))
        return np.array(out)

    def set_state(self, state):
        """Set the state to a specific value

        Args:
            state (Dictionary): an (exhaustive!) dictionary setting the values for each element of the state 
        """
        for a in self.species_names:
            self.state[self.species_ordering[a]] = state[a]

    def set_parameters(self, parameters):
        """Set the value of the parameters

        Args:
            parameters (Dictionary): an (exhaustive!) dictionary setting the values for each parameter
        """
        self.parameters = {}
        for p in self.parameters_names:
            self.parameters[p] = parameters[p]

    def update(self, reaction_idx):
        """Update the state accordingly to a given reaction

        Args:
            reaction_idx (int): the index of the firing reaction
        """
        self.state += self.stoichiometric_matrix[:, reaction_idx]

    def get_number_of_reactions(self):
        return len(self.reaction_names)

    def get_number_of_species(self):
        return len(self.species_names)

    def get_number_of_parameters(self):
        return len(self.parameters_names)

    def SSA(self, state, parameters, T0, Tf):
        """

        :param state: a dictionary
        :param parameters: a disctionary
        :param T0: the inital time
        :param Tf: final time
        :return: time, state
        """
        # Initialize time and state variables
        t = T0
        self.set_parameters(parameters)
        self.set_state(state)

        # Store initial state in the output
        time_output = [t]
        state_output = [self.state]

        while t < Tf:
            # Calculate propensities
            a = self._eval()

            # Generate two random numbers
            r1 = np.random.rand()
            r2 = np.random.rand()

            # Calculate the time until the next reaction occurs
            alpha = np.sum(a)
            tau = (1/alpha) * np.log(1/r1)

            # Choose the next reaction
            mu = np.argmax(np.cumsum(a) >= r2*alpha)

            # Update the time and state
            t = t + tau
            if t > Tf:
                break
            self.state = self.state + self.stoichiometric_matrix[:,mu]

            # Store the new state in the output
            time_output.append(t)
            state_output.append(self.state)

        # the last update
        time_output.append(Tf)
        state_output.append(self.state)

        return time_output, state_output

    # extract the dynamics of particular species
    def extract_trajectory(self, time_list, state_list, species_name_list):
        # initialization
        species_ordering = {species: count for species, count in zip(species_name_list, range(len(species_name_list)))}
        species_index = [self.species_ordering[species] for species in species_name_list]
        states = np.array(state_list)[:, species_index]
        time_new = [time_list[0]]
        state_new = [states[0,:]]

        for i in range(len(time_list)-2):
            dX = states[i+1,:] - states[i,:]
            if np.any(dX != 0):
                time_new.append(time_list[i+1])
                state_new.append(states[i+1,:])

        # record the final distribution
        time_new.append(time_list[-1])
        state_new.append( [states[-1, :]] )

        return time_new, state_new, species_ordering


    # plot the result
    def plot_trajectories(self, time_list, state_output):
        state_output = np.array(state_output)
        fig, axs = plt.subplots(math.ceil(self.get_number_of_species()/2), 2, figsize=(5,2))
        fig.subplots_adjust(wspace=0.5, hspace=1)

        for i, ax in enumerate(axs.flatten()):
            if i >= self.get_number_of_species():
                break
            ax.step(time_list, state_output[:, i], where='post')
            ax.set_xlabel('Time')
            ax.set_ylabel('Copy number of ' + self.species_names[i])

if __name__ == '__main__':
    
    # central dogma

    propensities = [
        lambda kr : kr,
        lambda kp, mrna : kp*mrna,
        lambda gr, mrna : gr*mrna,
        lambda gp, protein : gp*protein  
    ]

    stoichiometric_matrix = np.array([
        [1, 0, -1, 0],
        [0, 1, 0, -1]
    ])

    rn = CRN(
        stoichiometric_matrix=stoichiometric_matrix,
        species_names = ['mrna', 'protein'],
        parameters_names = ['kr', 'kp', 'gr', 'gp'],
        reaction_names= ['mRNA birth', 'Protein birth', 'mRNA degradation', 'Protein degradation'],
        propensities = propensities )

    rn.set_state({'mrna' : 10, 'protein' : 3})
    rn.set_parameters({'kr': 1, 'gr': 2, 'kp' : 1, 'gp' : 1})

    print(rn.state)
    print(rn._eval())
    rn.update(3)
    print(rn.state)
    print(rn._eval())
    rn.update(2)
    print(rn.state)
    print(rn._eval())

    print('')

    # test SSA
    time_list, state_list = rn.SSA({'mrna' : 10, 'protein' : 3}, {'kr': 1, 'gr': 2, 'kp' : 1, 'gp' : 1}, 0, 1)
    print(time_list)
    print(state_list)

    rn.plot_trajectories(time_list, state_list)
