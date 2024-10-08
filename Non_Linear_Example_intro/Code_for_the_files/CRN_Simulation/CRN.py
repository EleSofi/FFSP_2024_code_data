
# how to write the propensity matrix
#propensities = [
#    lambda mrna, protein, k1 : mrna+protein*k1,
#    lambda mrna, protein, k2 : protein + ...,
#    ...
#] 

import inspect

from scipy import sparse

from .Data import *
from .DistributionOfSystems import DistributionOfSystems
from .MarginalDistribution import MarginalDistribution
from .MatrixExponentialKrylov import MatrixExponentialKrylov


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
        #self.derivatives = derivatives
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
            if alpha == 0:
                break # no more reactions can occur
            else:
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
            state_output.append(self.state.copy())

        # the last update
        time_output.append(Tf)
        state_output.append(self.state.copy())

        return time_output.copy(), state_output.copy()

    def SSA_rich_data_structure(self, state, parameters, T0, Tf):  # TODO this version of the SSA algorithm was changed to create the Trajectory object, check its details as well
        """
        This procedure return a Trajectory object that tracks
        time, state, firing reactions and propensities of a SSA simulation

        Consider that the interval stored contains the initial time and the final time
        (with no associated firing reactions/propensities)

        :param state: a dictionary
        :param parameters: a disctionary
        :param T0: the inital time
        :param Tf: final time
        :return: Trajectory
        """
        # Initialize time and state variables
        t = T0
        self.set_parameters(parameters)
        self.set_state(state)

        # Store initial state in the output

        trj = Trajectory()
        trj.add_time(t)
        trj.add_state(self.state.copy())
        trj.add_propensity(self._eval())
        trj.add_firing_reaction(-1)

        while t < Tf:
            # Calculate propensities
            a = self._eval()
            # save the indexes of firing reactions

            # Generate two random numbers
            r1 = np.random.rand()
            r2 = np.random.rand()

            # Calculate the time until the next reaction occurs
            alpha = np.sum(a)
            tau = (1 / alpha) * np.log(1 / r1)

            # Choose the next reaction
            mu = np.argmax(np.cumsum(a) >= r2 * alpha)

            # Update the time and state
            t = t + tau
            if t > Tf:
                break

            trj.add_firing_reaction(mu)  # save the index of the firing reaction
            trj.add_propensity(a)  # save the propensities of all reactions
            self.state = self.state + self.stoichiometric_matrix[:, mu]

            # Store the new state in the output
            trj.add_time(t)
            trj.add_state(self.state.copy())

        # the last update
        trj.add_firing_reaction(mu)
        trj.add_propensity(a)
        trj.add_time(Tf)
        trj.add_state(self.state.copy())
        return trj

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


    # marginal time average distribution for a given species
    def marginal_time_average_distribution(self, time_list, state_list, species_name):
        if species_name not in self.species_names:
            raise Exception(f'Species {species_name} is not in the model')


        species_index = self.species_ordering[species_name]
        state_list = np.array(state_list)
        max_state = int(np.max(state_list[:,species_index]))
        marginal_time_average = np.zeros(max_state+1)
        for i in range(len(time_list)-1):
            marginal_time_average[int(state_list[i,species_index])] += time_list[i+1] - time_list[i]

        return marginal_time_average / (time_list[-1] - time_list[0])

    ###############################################
    # FSP algorithm
    ###############################################

    def FSP(self, Initial_marginal_distributions, range_of_spcies, parameters, T0, Tf, normalization=False):
        """

        :param Initial_distributions: A list of marginal distributions for each species (assume that the species are independent
        :param range_of_spcies: a df with the range of each species
        :param parameters: a dictionary of parameters
        :param T0: initial time
        :param Tf: final time
        :return:
        """

        # Prepare the initial distribution
        species_ordering, states = self.FSP_prepare_the_state_space(range_of_spcies)
        distribution = DistributionOfSystems(states, species_ordering)
        initial_joint_distribution = self.generate_joint_distribution_from_marginal_distributions(Initial_marginal_distributions, distribution)
        initial_joint_distribution = initial_joint_distribution.reshape(-1,1)
        distribution.extend_distributions([0],[initial_joint_distribution])

        # Construct the A matrix
        A = self.constract_A_matrix_for_CME(distribution, parameters)

        # Solve the CME
        time_list, distribution_list = \
            MatrixExponentialKrylov.exp_AT_x(A, T0, Tf, initial_joint_distribution)

        # normalize the distribution
        if normalization == True:
            distribution_list = [propability / np.sum(propability) for propability in distribution_list]

        # save the result in the distribution
        distribution.extend_distributions(time_list, distribution_list)

        return distribution

    def FSP_prepare_the_state_space(self, range_of_species):
        """

        :param range_of_species: a df with the range of each species
        :return:
        """
        coords_species = [
            np.arange(range_of_species.loc[species, 'min'], range_of_species.loc[species, 'max'] + 1) \
            for species in self.species_ordering.keys()]

        size_of_state_space = np.prod([len(coords) for coords in coords_species])
        if size_of_state_space > 1e6:
            raise Exception(f'The size of state space is {size_of_state_space}, which is too large for FSP algorithm')

        # construct the state
        meshes = np.meshgrid(*coords_species, indexing='ij')
        matrix = np.stack(meshes, axis=-1)
        states = matrix.reshape(-1, matrix.shape[-1])  # each row is a state of hidden species and species

        return self.species_ordering, states

    # generate uniform marginal distributions for every species and parameters
    def generate_uniform_marginal_distributions_via_speceis_range(self, range_of_species):
        Marginal_distributions = {}
        # transverse all species
        for species in self.species_names:
            states = list(range(range_of_species.loc[species, 'min'], range_of_species.loc[species, 'max'] + 1))
            uniform_distribution = np.ones(len(states)) / len(states)
            marginal_uniform_distribution = MarginalDistribution(species, states, uniform_distribution)
            Marginal_distributions.update({species: marginal_uniform_distribution})
        return Marginal_distributions

    def generate_joint_distribution_from_marginal_distributions(self, Marginal_distributions, distribution):
        """

        :param Marginal_distributions: a dictionary of marginal distributions
        :param distribution: an object of DistributionOfSystems
        :return:
        """
        joint_distribution = np.ones(distribution.states.shape[0])
        for count, state in enumerate(distribution.states):
            for species in self.species_names:
                marginal_distribution = Marginal_distributions[species]
                state_index =   marginal_distribution.states.index(state[self.species_ordering[species]])
                joint_distribution[count] *= marginal_distribution.distribution[state_index]

        return joint_distribution

    def constract_A_matrix_for_CME(self, distribution, parameters):
        """

        :param distribution: an object of DistributionOfSystems
        :param parameters:
        :return:
        """
        # return row_index, column_index, reaction_list, state_column, sign
        column = []
        row = []
        value = []
        self.set_parameters(parameters)

        for state_index, state in enumerate(distribution.states): # transverse all states where the probability outflows
            state_dic = {species: state[distribution.species_ordering[species]] for species in distribution.species_ordering.keys()}
            self.set_state(state_dic)
            propensities = self._eval() # a list of propensities

            # probability outflows
            column.append(state_index)
            row.append(state_index)
            value.append(-np.sum(propensities))

            # probability inflows to other states
            for reaction_index, propensity in enumerate(propensities):
                state_after_reaction = state + self.stoichiometric_matrix[:, reaction_index]
                state_after_reaction_index = np.where(np.all(distribution.states == state_after_reaction, axis=1))[0]
                if len(state_after_reaction_index) > 0:
                    column.append(state_index)
                    row.append(state_after_reaction_index[0])
                    value.append(propensity)

        number_of_states = distribution.states.shape[0]
        A = sparse.coo_matrix((value, (row, column)), shape=(number_of_states, number_of_states)).tocsr()

        return A







    # plot the result
    def plot_trajectories(self, time_list, state_output):
        state_output = np.array(state_output)
        rows = math.ceil(self.get_number_of_species()/2)
        columns = 2
        fig, axs = plt.subplots(rows, columns, figsize=(columns*3,rows*3))
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
