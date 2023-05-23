from math import log, exp
import sys
import numpy as np

from TreeStruct import Nstates

def NegLogL(var_rates, root, root_prior, rate_arrange, \
             fixed_rate, which_fixed, time_slice):
    '''
    Runs the machinery for calculating the Mk likelihood.
    Returns the negative log-likelihood of the tree and states.
    (If negative parameter values are passed, a very large positive value is 
      returned.  This prevents negative parameter values from being considered
      in the minimization.)
    '''

    # Can't have negative rate values.
    for r in var_rates:
        if r < 0:
            return np.inf

    # If some rates are being fixed, arrange them
    if rate_arrange == "fix":
        rates = list(var_rates)
        rates.insert(which_fixed, fixed_rate)
        rates = np.array(rates)
    elif rate_arrange == "equal":
        rates = list(var_rates) * Nstates
        rates = np.array(rates)
    else:
        rates = var_rates

    # Assign conditional likelihoods to each node in the tree
    if time_slice == None:
        _GetTreeCLs(root, rates)
    else:
        _GetTreeCLsSlice(root, rates, time_slice)

    # Compute the likelihood by combining the cl's at the root.
    loglike = _CombineAtRoot(rates, root, root_prior)

    return -loglike

#--------------------------------------------------
# Functions for the guts of the likelihood calculation
#-------------------------------------------------- 

def _GetTreeCLs(node, rates):
    '''
    Get the conditional likelihoods (likelihood of being in each state)
      for each node in the tree.
    '''

    ### get somewhere in the tree
    if node.daughters != None:
        for d in node.daughters:
            _GetTreeCLs(d, rates)

    ### compute the cl for this node

    # if it's a tip, the cl's are determined by the observed state
    if node.daughters == None:
        try:
            for i in range(Nstates):
                node.cl[i] = 0
            node.cl[node.state] = 1
        except TypeError:
            print("ERROR: Tip state not specified.  Aborting in Mk2Like...")
            sys.exit()

    # otherwise, the cl's are computed from the cl's of the daughters
    else:
        _ComputeCL(node, rates)

def _GetTreeCLsSlice(node, rates, time_slice):
    '''
    Get the conditional likelihoods (likelihood of being in each state)
      for each node in the tree, for the special case of transitions allowed
      only within time_slice.
    '''

    ### get somewhere in the tree
    if node.daughters != None:
        for d in node.daughters:
            _GetTreeCLsSlice(d, rates, time_slice)

    ### if it's a tip, the cl's are the observed state
    if node.daughters == None:
        try:
            for i in range(Nstates):
                node.cl[i] = 0
            node.cl[node.state] = 1
        except TypeError:
            print("ERROR: Tip state not specified.  Aborting in Mk2Like...")
            sys.exit()

    ### if the node is outside the time slice, don't allow transitions
    elif node.time < time_slice[0] or node.time >= time_slice[1]:
        _ComputeCL(node, (0, 0))

    ### otherwise, allow transitions normally
    else:
        _ComputeCL(node, rates)

def _ComputeCL(node, rates):
    '''
    Get the conditional likelihood for a node, computed from the cl's of 
      its daughters.
    '''

    # For each daughter, compute the transition probabilities
    #  from daughter state to parent state.
    # Multiply the daughter probabilities together, applying log-compensation.
    # Note: Because of log-comp, the cl value at a node is not its full CL.
    #       log(node's full CL) = log(node.cl) + node.lq
    #       But cl is still on the regular, not-log scale.

    CLresults = [None] * len(node.daughters)
    # one element for each daughter
    # each element will have length Nstates, for the possible parent states

    for i_d, d in enumerate(node.daughters):
        # compute the value for the parent state (s) considering all possible
        #   states of the daughter (s_d)
        cl_temp = [None] * Nstates
        for s in range(Nstates):
            d_sum = 0
            for s_d in range(Nstates):
                # transition prob from daughter state (s_d) to parent state (s)
                #   (d.cl already indicates if daughter state is fixed)
                d_sum += _TransitionProb(s, s_d, d.length, rates) * d.cl[s_d]
            cl_temp[s] = d_sum
        CLresults[i_d] = cl_temp

    # If it's a true node (not a dummy), do log-compensation
    if len(node.daughters) > 1:

        # Log-compensation (to avoid underflow in cl)
        CLresults = np.array(CLresults)
        # divide the cl's for each daughter by the q = sum(cl_i) for that daughter
        q = np.sum(CLresults, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            # (if some q = 0, nan will be returned, so okay to hide the warning)
            # (don't use /= because might be converting int to float)
            CLresults = CLresults / q[:,None]
            # lq for this node = sum of lq for each daughter now (branch bases), and
            #   of each daughter node's lq (branch inits)
            lq = np.sum(np.log(q))
        for d in node.daughters:
            lq += d.lq
        node.lq = lq

        # cl[i] for the node contains the product of the daughter cl[i]'s
        node.cl = np.prod(CLresults, axis=0).tolist()

    else:
        node.cl = CLresults[0]
        node.lq = node.daughters[0].lq

def _TransitionProb(to_state, from_state, time, rates):
    '''
    From transition rates and branch length, calculate transition probability.
    Currently, only valid for Nstates = 2.
    '''

    if all([r == 0 for r in rates]) or time == 0:
        if to_state == from_state:
            ans = 1
        else:
            ans = 0
    else:
        ans = ( rates[(from_state+1)%2] + pow(-1, to_state-from_state) * \
                (rates[to_state] * exp(-(rates[0]+rates[1])*time)) \
               ) / (rates[0]+rates[1]);

    return ans

def _CombineAtRoot(rates, root, root_prior, separate=False):
    '''
    Invoke the assumption about the root prior here.
    Log-compensation is restored.
    Get the total likelihood (default) or a list of likelihoods for each 
       state (separate=True)
    '''

    like = root.cl.copy()

    ### calculate the appropriate weights for each root state

    root_p = [None] * Nstates

    # stationary distribution
    if root_prior == "stationary":
        assert Nstates == 2
        root_p[0] = rates[1] / sum(rates)
        root_p[1] = 1 - root_p[0]

    # equal weights for each state
    elif root_prior == "uniform":
        root_p =  [1./Nstates] * Nstates

    # weight by the data itself
    elif root_prior == "condlike":
        for i in range(Nstates):
            try:
                root_p[i] = root.cl[i] / sum(root.cl)
            except ZeroDivisionError:
                root_p[i] = 0
                # then later get like = 0 and return -inf

    # arbitrary root prior
    elif type(root_prior) == list and len(root_prior) == Nstates:
        root_p = [ float(p) for p in root_prior ]

    else:
        print("ERROR: invalid root_prior specified.")
        sys.exit()        # more graceful exit?

    # apply the root state weightings
    for i in range(Nstates):
        like[i] *= root_p[i]

    ### return the log-likelihood

    if separate:
        ans = [-np.inf] * Nstates
        for i, x in enumerate(like):
            try:
                ans[i] = log(x) + root.lq
            except ValueError:
                pass
    else:
        try:
            ans = log(sum(like)) + root.lq
        except ValueError:
            ans = -np.inf

    return ans
