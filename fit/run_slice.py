### Compare Mk time-slice models via direct integration to get marginal likelihood.
### Slices have fixed width (1 month), as many as needed for the tree.
### Slice 1 starts at the latest-sampled time.

# Input: directory of treefiles
# Output: files of marginal likelihood values
#         one file per tree, with time-homogeneous + all time slices
#
# No parallelization.  (Can run as an array job on a cluster.)

import Newick, TreeExtra, Mk2Like
import sys, os, glob
import numpy as np
import scipy.integrate as integrate

#--------------------------------------------------
# Posterior = Likelihood * Prior
#--------------------------------------------------

def logprior(theta, prior_rate):
    lp = 0
    for param_rate in theta:
        # exponential prior (prior_mean = 1 / prior_rate)
        lp += np.log(prior_rate) - param_rate * prior_rate
    return lp

def loglike(theta, root, which_fixed, time_slice):
    # which_fixed = 0 for R -> D
    # which_fixed = 1 for D -> R
    return -Mk2Like.NegLogL(theta, root=root, root_prior="condlike", \
            rate_arrange="fix", fixed_rate=0, \
            which_fixed=which_fixed, time_slice=time_slice)

def post(q, root, which_fixed, time_slice, prior_rate):
    theta = [q]
    return np.exp(logprior(theta, prior_rate) + loglike(theta, root, which_fixed, time_slice))

def runme(treefile):

    ### Get the tree ###

    tree = Newick.ReadFromFileTTN(treefile)
    TreeExtra.AssignNodeTimes(tree)

    prior_rate = 1 # reconsider if not one-month slices on few-year-old tree

    outfile = os.path.join(os.path.dirname(treefile), "mk2",
                           os.path.basename(treefile).rsplit(".", 1)[0] + "-mk2.csv")
    with open(outfile, "w") as ofp:
        ofp.write("slice,direction,marglik\n")

    ### Time-homogeneous ###

    ans01 = integrate.quad(post, 0, np.inf, args=(tree, 1, None, prior_rate))
    ans10 = integrate.quad(post, 0, np.inf, args=(tree, 0, None, prior_rate))

    with open(outfile, "a") as ofp:
        ofp.write("ns,01," + str(np.log(ans01[0])) + "\n")
        ofp.write("ns,10," + str(np.log(ans10[0])) + "\n")

    print("done with time-homogeneous for", treefile)

    ### Time slices ###

    # time increases from root to tips
    # root time might not be zero

    # first slice (slice 1) ends at the latest-sampled tip
    # slices work back from then in 1 month increments
    tip_time = tree.time + tree.Age() # root time + max root-to-tip distance
    all_slice_times = np.arange(tip_time, tree.time, -1/12) # fixed slice size
    all_slice_times = np.append(all_slice_times, tree.time) # put in root time
    num_slices = len(all_slice_times) - 1

    for n in range(num_slices):

        t_slice = [all_slice_times[n+1], all_slice_times[n]] # note the flip
        tree = Newick.ReadFromFileTTN(treefile) # fresh copy, so slice can be inserted
        TreeExtra.AssignNodeTimes(tree)
        for t in t_slice:
            TreeExtra.InsertNodesSlice(tree, t)

        ans01 = integrate.quad(post, 0, np.inf, args=(tree, 1, t_slice, prior_rate))
        ans10 = integrate.quad(post, 0, np.inf, args=(tree, 0, t_slice, prior_rate))

        s = "s" + str(n+1).zfill(2)
        with open(outfile, "a") as ofp:
            ofp.write(s + ",01," + str(np.log(ans01[0])) + "\n")
            ofp.write(s + ",10," + str(np.log(ans10[0])) + "\n")

        print("done with slice", n+1, "of", num_slices, "for", treefile)

if __name__ == '__main__':

    wd = sys.argv[1]
    os.makedirs(os.path.join(wd, "mk2"), exist_ok=True)

    ttnfiles = glob.glob(wd + "*.ttn")
    ttnfiles.sort()

    for ttn in ttnfiles:
        runme(ttn)
