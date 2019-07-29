# Introduction
This text explains how functions of each script work, what parameters require and why, and what outputs they produce.
If you're interested in just one function or script, use Table of Content to directly find it: it follows the same order of the scripts.

# Table of contents (scripts and their functions)
* [Main.py](#main)
* [Main_analysis.py](#main_analysis)
* [Evolutions.py](#evolutions)
  * [Network generation functions](#network_generation_dar-and-network_generation_tgrg)
  * [Degree functions](#degree-functions)
  * [Communicability and rankings](#communicability-and-rankings)
* [Propagation SI.py](#propagation_si)
  * [Easing simulation functions](#easing-simulation)
  * [Propagation](#propagation)
  * [Epidemic scores functions](#epidemic-scores)
* [Propagation LTM](#propagation_ltm)
* [Saves.py](#saves)
  * Networks save/load
  * Epidemics save/load

# Main
        From this script, user can generate temporal networks and perform epidemic upon them.
        System parameters must be set in "inputs.ini" file.

Here's a list of the changeable parameters you can set in "inputs.ini" file:
* **N, T**: nodes and duration of temporal network;
* **isDAR**, **isDIRECTED**: tempnet evolution's law and symmetry;
* **net_REAL**,**K**: number of tempnet realizations and epidemic simulations; generating several tempnet is up to you, but iterating propagation for nodes is "compulsory", since scores MUST be averages; so net_REAL can be 1, K must be > 1 (some assertions make sure of this);
* **net_name**: identificative name you give to network (make sure to variate: networks with same N,T,type,iteration and identification may be overwritten;
* **P**,**alpha** and **xi** are defined, for DAR networks, in explanation, as well as, **phi0,phi1,epsilon** for TGRG;
* **beta**: infection rate, that defines epidemic virulence.

To create a new network or epidemic, just add a new section in the input file: this script performs a simulation for each section. Similarly, remove a section if it's not more of your interest, otherwise it will be always used for simulations.

The script produces, and saves, a **temporal_network** and a **label** structure, whose syntax is: label\[index_case]\[iteration]\[time_step]\[node]; so, it is a list of N dictionaries (one for each node-set as index case), containing a list of the results K iterations (one for each same epidemic simulation), which in turn are dictionary of dictionaries, describing evolution of nodes states (as you will see in *propagation function*);

# Main analysis
         From this script, user can perform structural and epidemiological measures upon networks.
         System parameters must be set in "inputs.ini" file

After loading both network and epidemic, centrality and virulence measures are performed, and these data strucrures are used to store informations:
* **spec_radius** and **Q** are the inverse of maximum spectral radius of all adjacencies and Communicability matrix;
* **B/R Centrality** and **virulence** scores are computed in namesake lists.

Script ends with printing of structural and virulence results.

# Evolutions
    Functions in this script work in Pyhon3, may require numpy (v1.16) and allow to:
    1) generate and update a temporal network, according to the dar(p) or trgr laws of evolution.
    2) perform structural analysis, such as degree evolution and centrality measures (BC/RC, AD, BD)

### network_generation_dar and network_generation_tgrg
These functions generate DAR(P)s and TGRGs, according to definion you can find in [Explanation](/explanation.md). Each update is a stochastic process, so it's difficult to get twice the same results, and you should call these functions multiple times to infer some statistically relevant properties of the system.

Let's have a look at what parameters these function require:

    network_generation_dar(alpha,xi,P=1, T=100, directed = False)

    network_generation_tgrg(phi0,phi1,sigma,T=100, directed = False)

Compulsory parameters are those defining networks' nature (alpha, xi; phi0,phi1,sigma), while one can change simuation' time steps (default is 100), whether to generate a directed or undirected graph (default is undirected), and, only for DAR(P), the number P of past step to consider (default is 1).
The output is np.array with shape T\*N\*N, i.e. an ordered sequence of T square(N) matrices, each of which is the adjacency at that time step.

Inputs undergo following assertions (mostly imported from Test_suite script), so be sure to respect them: 
* matrices/vectors provided must be np.arrays with square or vectorial shape and proper length;
* alpha and xi values must be probabilities, i.e., all their values should be included within 0 and 1; 
* if all these are verified, functions extract the information about the number of nodes (contained in matrices/vector "length") and verify it is < 1000 (it's improbable to work with more than 1000 nodes, so  if it's the case, it may be due to human error); 
* the provided T and P (for DAR(P)) must be natural numbers, and also T must be < 1000;
* P needs to be < T "for definition"

Strucure of both codes follow these ideas:
* to easily switch for undirected to directed case, they firstly generate an upper triangoular matrix, and then, for the lower one, they just copy the result or keep on extracting values, accordingly;
* they both exploit "np.random.choice", that performs an extractions between 2 numbers according to a given probability for each one; probability conservation (i.e. sum = 1) is guaranteed by initial assertions (for DAR) or by mathematical properties of distribuion (for TGRG).

This is what the function actually do:
* In DAR, firstly, the initial adjacency is extracted (with 50% probability of having or not a link for each initial couple) and the whole T\*N\*N tensor is generated, copying in the the initial state. So, everything is set to continue evolution. Both initial state and evolutions are handled with two nested comprehensions and for cicles (that distinguish triangles by setting i> or <j), the inner one containing "choice". For the initial states, choice is between 0 and 1 with probability 0.5 for each, while for following moments, they are took from the previous P steps with probability alpha, or extracted by a tossing coin (whose "inner" probability is xi), just as DAR(P) says.

      #EVOLUTION
      #Initial adiacency, as tossing simple coin (0.5); if undirected, it's made symmetric
      initial_state = [[np.random.choice((0,1),p=(0.5,1-0.5)) if j>i else 0 for j in range(N)] for i in range(N)] #upper triangle first (j>i)
      if directed == False: #setting lower triangle as directed value
          initial_state = [[initial_state[j][i] if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]
      else:
          initial_state = [[np.random.choice((0,1),p=(0.5,1-0.5)) if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]

      #Temporal network initilialization
      temporal_network = np.zeros((T,N,N))
      temporal_network[0] = initial_state

      #Network evolution
      for t in range(1,T):
          #UPPER TRIANGLE (j>i):
          temporal_network[t] = [[np.random.choice((temporal_network[t-P,i,j],np.random.choice((1,0), p=(xi[i,j],1-xi[i,j]))),
                          p=(alpha[i,j],1-alpha[i,j])) if j>i else 0 for j in range(N)] for i in range(N)] #DAR(P) generation
          #LOWER TRIANGLE (j<i)
          if directed == False:    
              temporal_network[t] = [[temporal_network[t,j,i] if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
              #copy upper triangle, if undirected
          else:
              temporal_network[t] = [[np.random.choice((temporal_network[t-P,i,j],np.random.choice((1,0), p=(xi[i,j],1-xi[i,j]))),
                          p=(alpha[i,j],1-alpha[i,j])) if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
              #follow the same rule as upper, if directed
      return temporal_network

          Examples
          --------
        >>> network_generation_dar(0.6*np.ones((3,3)),0.5*np.ones((3,3)),P=1, T=3, directed = False)
        array([[[0., 1., 0.],
        [1., 0., 1.],
        [0., 1., 0.]],
        [[0., 1., 0.],
        [1., 0., 1.],
        [0., 1., 0.]],
        [[0., 0., 0.],
        [0., 0., 0.],
        [0., 0., 0.]]])

* In TGRG, firstly, all fitnesses evolutions are extracted (to see how this works, check "expalanation.md"), for each node, after initializing them to be equal to the respective "phi0". Then, T\*N\*N tensor is generated and, making use of nested for-cicles, and distinguishing triangles by setting i> or <j, it is computed the probability of link for each time step and couple of nodes, and it is put in "choice".

        #FITNESSES EVOLUTION
      theta = np.zeros((N,T)) #definition
      theta[:,0] = phi0
      for t in range(1,T):
          theta[:,t] = phi0 + phi1*theta[:,t-1] + np.random.normal(0,sigma,size=N) #evolution for each node and each time

      #ADIACENCIES COMPUTATION
      temporal_network = np.zeros((T,N,N)) #tensor definition
      for t in range(T):
          for i in range(N):
              for j in range(N):
                  #TGRG FOR UPPER TRIANGLE
                  temporal_network[t] = [[np.random.choice((1,0), p=(np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])),
                                  1-np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])))) if j>i else 0 for j in range(N)] for i in range(N)]
                  #LOWER TRIANGLE (j<i)
                  if directed == False:    
                      temporal_network[t] = [[temporal_network[t,j,i] if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
                      #copy upper triangle, if undirected
                  else:
                      temporal_network[t] = [[np.random.choice((1,0), p=(np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])),
                                      1-np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])))) if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
      return temporal_network, theta
      
        Examples
        --------
        >>> network_generation_tgrg(0.6*np.ones(3),0.4*np.ones(3),0.1*np.ones(3),T=3, directed = False)
        ( array([[[0., 1., 0.],
         [1., 0., 1.],
         [0., 1., 0.]],
        [[0., 1., 1.],
         [1., 0., 1.],
         [1., 1., 0.]],
        [[0., 1., 1.],
         [1., 0., 1.],
         [1., 1., 0.]]]), 
        array([[0.6       , 0.99664444, 1.09950828],
        [0.6       , 0.9285223 , 1.00334695],
        [0.6       , 0.73700156, 0.85323848]]) )
                                      
DAR gives back the network only, TGRG also the matrix of fitnesses.

### degree functions
These functions can compute degree for a single node or for a whole temporal network, and will be very useful in further analysis. 
They just require the adjacency/cies, from which extracting informations, and the specification of computing out- or in-going degree. Both functions will be useful in further development of this project

Before operating both functions check that adjacencies ar NN or NNT np.array, and have null diagonal; if node index is required, it is checked whether it exists in input network.

* **degree_node** requires 1 adjacency and returns out or in-going degree of the provided node, just by summing (respectively) column of line of node in provided network, and returns it as a natural number;

      def degree_node(network,node,out = True):
           #FUNCTION
           if out:
              return sum(network[node])
               else:
                   return sum(network[:,node])
       
         Examples
         --------
        >>> degree_node(np.array([[0,0,1],[1,0,1],[1,0,0]]),0)
        1
        
        >>> degree_node(np.array([[0,0,1],[1,0,1],[1,0,0]]),1)
        2

* **degree_mean** computes mean (out or in)-degree of all nodes in a temporal network of any duration (so, even for T = 1). Output is always a list, with T entries. To use the same code in case of both 2- or 3-dimensional arrays, function first checks dimensions, and, if bidimensional, copies adjacency in a 1NN array.

      def degree_mean(tempnet, out = True):
        assert isinstance(tempnet,np.ndarray)

        N = tempnet.shape[1] #this is true anyway
        #Infer temporal duration and perform assertions
        if len(tempnet.shape) == 3:
            T = tempnet.shape[0]
            [Test_suite.assert_square(tempnet[t]) for t in range(T)] #check square for each step
            [Test_suite.assert_nulldiagonal(tempnet[t]) for t in range(T)] #check null diagonal for each step
            network = tempnet #same new for 3- or 2-dimensional arrays
        elif len(tempnet.shape) == 2:
            T = 1
            Test_suite.assert_square(tempnet) #check square for each step
            Test_suite.assert_nulldiagonal(tempnet) #check null diagonal for each step
            network = np.zeros((T,N,N))
            network[0] = tempnet #same new for 3- or 2-dimensional arrays
        else:
            raise AssertionError ("You must provide an adjacency or a sequence of adjacencies")

        #FUNCTION
        d_seq = []
        for t in range(T):
            degrees = []
            for node in range(N):
                if out:
                    degrees.append(sum(network[t][node]))
                else:
                    degrees.append(sum(network[t][:,node]))
            d_seq.append(np.average(degrees))
        return d_seq #output is always a list
      
      Examples
      --------
        >>> degree_mean(np.array([[0,0,1],[1,0,1],[1,0,0]]))
        [1.3333333333333333]
        
        >>> degree_mean(np.array([[[0,0,1],[1,0,1],[1,0,0]],
        [[0,1,1],[1,0,1],[1,0,0]]]))
        [1.3333333333333333, 1.6666666666666667]


### communicability and rankings
These functions compute Communicability matrix for a provided temporal network, perform Broadcast and Receive Centrality and rank nodes accordingly.
The only assertions they perform are about shape (sequence of squadre matrices) and diagonal nullity.

* **communicability** extracts the number of nodes and evolution duration by the shape of the temporal network, then computes for each adjacency the maximum spectral radius (each of which makes use of np.linalg.eigvals) and selects the maximum from these. A quarter of the highest one is the factor applied to each adjacency in Communicability updating (as in Chen paper), but in further developments it could be different, or it could be given the possibility to choose among different ones.
Before operating, input is checked to be a temporal network, with square adjacencies and null diagonal, and each adjacency is checked not to be a 0 matrix: since we're working with inverse matrices, this operation must be well defined.
Communicability is initialized as a normalized (np.linalg.norm) identity matrix, and multiplied, at each time step, by the normalized inverse (np.linalg.inv) of the difference between identity (np.identity) and the the weighted adjacency; normalizations take place step by step, multiplications are performed by np.matmul.

    def communicability(temporal): 
         #FUNCTION
      T = temporal.shape[0]
      N = temporal.shape[1]
      #Find max spectral radius:
      spec = []
      for t in range(T):
          spec.append(np.real(max(np.linalg.eigvals(temporal[t])))) #find eigenval with max real part for each adiacency
      rec_maxradius = 1/max(spec) #reciprocal of the maximum eigenvalue
      #Communicability builing:
      Q = np.identity(N)/np.linalg.norm(np.identity(N)) #initialization (and normalization)
      for t in range(T):
          inv = np.linalg.inv(np.identity(N)-0.25*rec_maxradius*temporal[t]) #new factor for that time step
          Q = np.matmul(Q,inv)/np.linalg.norm(np.matmul(Q,inv)) #Q updating and normalizing
      return(rec_maxradius,Q)
    
     Examples
     --------
        >>> communicability(np.array([[[0,0,1],[1,0,1],[1,0,0]],
        [[0,1,1],[1,0,1],[1,0,0]]]))
        ( 0.6180339887498951, 
        array([[0.54377529, 0.0840179 , 0.17483766],
        [0.20185156, 0.52294101, 0.20185156],
        [0.16411784, 0.0253576 , 0.53305547]]))

* **broadcast_ranking/receive_ranking** require such a communicability matrix, and returns a list with nodes BCs/RCs by summing over lines/columns; then it generates a list, lines_sum, containing nodes indeces ordered by their scores, the firsts being those with highest centralities, using np.argsort (that return indeces but from the worst to the best) and np.flip (that invertes lists); note that "lines_sum" doesn't follow such a order;

      def broadcast_ranking(Q):
        #FUNCTION    
        lines_sum = np.sum(Q, axis = 1) #Broadcast -> sum over lines:
        rank = np.flip(np.argsort(lines_sum)) #argsort -> increasing score; flip -> decreasing
        return(lines_sum,rank)
      
       Examples
       --------
        >>> broadcast_ranking(np.array([[0.54377529, 0.0840179 , 0.17483766],
        [0.20185156, 0.52294101, 0.20185156],
        [0.16411784, 0.0253576 , 0.53305547]]))
        (array([0.80263085, 0.92664413, 0.72253091]), 
        array([1, 0, 2], 
        dtype=int64))
      

# Propagation_SI
    Functions in this script work in Pyhon3, may require numpy (v1.16) and allow to:
    1) spread an epidemic, in SI mode, over a temporal network
    2) measure virulence of each node

### easing simulation
These functions allow [propagation](#propagation) to perform well and in an efficient way, solving simple tasks as analyzing adjacencies and node epidemic states.
* **neighbourhood**: 

       def neighbourhood(adjacency,node):
       "Extracts the subset of nodes that reach the given one at a certain time. So, it doesn't care about network's directness."
            #FUNCTION
           neigh = {i for i in range(len(adjacency)) if adjacency[i,node]==1}
           return neigh
       
       Examples
       -------
        >>> neighbourhood(np.array([[0,1],[0,0]]), 1)
        {0}
        
        >>> neighbourhood(np.array([[0,1],[0,0]]), 0)
        set()

* **onlyzeros**: extracts the subset of susceptible nodes at that time, from a general set, list, etc, by checking the provided states dict; set is create through a comprehension;

      def onlyzeros(nodes_set,states_dict):
          selected = {node for node in nodes_set if states_dict[node]==0} #0 means susceptible
          return selected
       
       Examples
       -------
        >>> onlyzeros([0,1,2], {0:1,1:0,2:0})
        {1,2}

* **contact_lasting**: computes the number of time steps an S and a I node were connected, starting from a provided time. This is accomplished by checking backwards, time by time, the existence of the link and the state of the I node, increasing the value of a counter variable as long these conditions are satisfied; so, the loop interrupts if it has gone so back that the link wasn't  yet present, or the I-node wasn't yet infected. Clearly, the higher the output is, the higher the probability of having an infection (the value is checked in the already mentioned probabilities dictionary).

      def contact_lasting(adjacencies,state,t,infected_node,susceptible_node):
          #FUNCTION
         counter = 0
         for instant in range(t+1):
             if (tempnet[t-instant,infected_node,susceptible_node] == 1 and states_sequence[t-instant][infected_node] != states_sequence[t-instant][susceptible_node]):
                 counter +=1
             else:
                 break
          return counter
         
         Examples
         -------
        >>> contact_lasting(np.ones((2,3,3)) - np.identity(3),
        {0:{0:1,1:0,2:0}, 1:{0:1,1:0,2:0}}, 1, 0, 1)
        1
        
        >>> contact_lasting(np.ones((2,3,3)) - np.identity(3),
        {0:{0:1,1:0,2:0}, 1:{0:1,1:0,2:0}}, 1, 0, 1)
        2

* **Poisson_Probability** reproduces the Poisson PDF (Probability Density Function), whose integral, defined over a temporal interval, is the probability of having an infection, for a SI couple, within that time; as stated in [Explanation](/explanation.md), since several integrals are often re-computed, they should be computed once for all, by *quad* function, and results are stored in **probabilities** dictionary; moreover, assertions ensure the input beta to be a probability
       
       def poisson_probability(t,beta): 
         lamda = -np.log(1-beta) # Chen, Benzi use 60
         return(lamda*np.exp(-lamda*t))

### propagation
This function requires a temporal network, an index case and the probabilities dictionary, and performs an iteration of a whole epidemic propagation in SI mode, from the first to last snapshot of the tempnet; the result is a dictionary of dictionaries, describing nodes' states for each time step.

First operations are, as usual, definitions: T and N are extracted from the tempnet, and then is defined an internal function, *set_infected*, that sets - from a given time to the end - its input node states equals to 1 (once infected, always infected) and will be used on every new infected; then, output is initialized by setting every node susceptible at every time, with the expection of the index case who is always infected (set_infected since t=0). Finally, infected and susceptibles nodes sets are initialized: they will be updated every time a new infection occurs. 
Following Chen approach, node states are updated time by time by processing each susceptible node (at previous time step), finding its infective neighbourhood at the previous time step (i.e. infected nodes whose links "point" to the susceptible one), and performing the *infection extraction* for each of them: if infection occurs, infective state is set "1" at present time step, S/I sets are updated and the program jumps directly to next temporal step (this is a conscious choice; future developments may adopt different solutions), by for...else loop (check [here](http://book.pythontips.com/en/latest/for_-_else.html) if you want to learn more). Infection extraction compares the integral of Poisson distribution with a random uniform (from 0 to 1) extraction, performed by *np.random.uniform* function. This is iterated until the end of network evolution.

    def propagation(tempnet,index_case,probabilities):
        #FUNCTION
        def set_infected(node,t): #once one is infected,it stays infeced
            for instant in range(t,T):
                states_sequence[instant][node] = 1
            return

        #Output initialization:
        states_sequence = dict()
        for t in range(T):
            states_sequence[t] = dict.fromkeys(range(N),0)
        set_infected(index_case,0)

        #Sets initialization
        susceptibles = onlyzeros(range(N),states_sequence[0]) #"targets"
        infecteds = {index_case} #they will be decisive to change target state

        for t in range(1,T):
            for s in susceptibles.copy(): #copy avoids rising an error when the iteration set changes
                infectneighbourhood = neighbourhood(tempnet[t-1],s).intersection(infecteds)
                for i in infectneighbourhood.copy(): 
                    if probabilities[contact_lasting(tempnet,states_sequence,t-1,i,s)]>np.random.uniform(0,1): #rand extraction
                            set_infected(s,t) #if successful, change the state of the node, at next t
                            susceptibles.remove(s)
                            infecteds.add(s)
                            break
                else:
                    continue # only executed if the inner loop did NOT break
                break  # only executed if the inner loop DID break
        return(states_sequence)
        
        Examples
        -------
        >>> propagation(np.ones((2,3,3)) - np.identity(3), 0, {0:0,1:0.999})
        {0: {0: 1, 1: 0, 2: 0}, 1: {0: 1, 1: 1, 2: 0}}


### epidemic scores
They are computed making use of two functions:
* **infected_counter**: takes the states dictionary at a certain time, and counts of many nodes are infected;

      def infected_counter(set_of_nodes):
          counter = 0
          for i in range(len(set_of_nodes)):
              if set_of_nodes[i]==1:
                  counter+=1
          return counter

* **time_score**: uses infected_counter for each time step of the propagation, and saves the first step when the fraction of infected nodes was higher than a given value. If this never happens, it returns the total duration of the propagation.

      def time_score(scores,fraction):
           #FUNCTION
       time_spent = T-1 #initialized as the final temporal step
       for t in range(T):
           if infected_counter(scores_evolution[t])>=fraction*N:
               time_spent = t
               break
       return time_spent
       
       Examples
       -------
    
        >>> time_score({0:{0:0,1:1,2:1,3:0,4:0,5:0},1:{0:1,1:1,2:1,3:0,4:0,5:0},2:{0:1,1:1,2:1,3:1,4:1,5:1}},0.5)
        1
        
        >>> time_score({0:{0:0,1:1,2:1,3:0,4:0,5:0},1:{0:1,1:1,2:1,3:0,4:0,5:0},2:{0:1,1:1,2:1,3:1,4:1,5:0}},0.8)
        2
        
        >>> time_score({0:{0:0,1:1,2:1,3:0,4:0,5:0},1:{0:1,1:1,2:1,3:0,4:0,5:0},2:{0:1,1:1,2:1,3:1,4:1,5:0}},0.9)
        3

# Propagation_LTM
This section will be deepened in further developments.

# Saves
     Functions in this script work in Pyhon3, require os, pickle, and allow to save results.
Results are saved, or load (you may need to load a previously generated network) by these functions:
 * **network_save**: serializes data structures in binary protocol, using **pickle** (if you don't know what does it mean, check [here](https://docs.python.org/3/library/pickle.html)), following the foldering/naming rule expressed in Readme (so, distinguishing networks according to their parameters and realizations), allowing user to add a particoular identification name to the file, and using **os** libary to generate folders. Os and pickle belong to Python standard library.
 * **network_load**: loads pickle files, checking them by parameters and identification name.
 
 * **infection_save/load** work in a very similar way; the first just doesn't make preliminary check, since at that point of the script everything should have worked as expected. They require to underline parameter "beta", which is used in output file name.
