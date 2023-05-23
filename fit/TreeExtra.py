import TreeStruct

def AssignNodeTimes(node, root_time=0):
    '''
    Use given branch lengths to assign node times.
    '''

    if node.parent == None:
        node.time = root_time
    else:
        node.time = node.parent.time + node.length

    if node.daughters != None:
        for d in node.daughters:
            AssignNodeTimes(d, root_time)

def AssignBranchLengths(node, backwards=False):
    '''
    Use given node times to assign branch lengths. 
    Normally, tip times > root times.  Set backwards=True if otherwise.
    '''

    if node.parent == None:
        node.length = 0
    else:
        if backwards:
            node.length = node.parent.time - node.time

        else:
            node.length = node.time - node.parent.time

    if node.daughters != None:
        for d in node.daughters:
            AssignBranchLengths(d, backwards)

def GetTipDistances(node, distances):
    '''
    Return a list of distances from tip to root.
    Probably only useful for non-ultrametric tree.
    Probably want to call after defining distances = [].
    '''

    if node.daughters == None:

        x = node.length
        p = node.parent

        while p != None:
            x += p.length
            p = p.parent

        distances += [x]

    else:
        for d in node.daughters:
            GetTipDistances(d, distances)

######################################################

def InsertNodesSlice(root, time):
    '''
    Along each branch that spans "time", insert a node at "time".
    Return the list of new nodes.
    '''

    if root.time == None:
        print("ERROR: node times are needed by InsertNodesSlice()")
        return None

    # first, create a list of relevant branches

    branch_list = []
    _SliceList(root, time, branch_list)

    # then, insert the new nodes
    # (note: inserting nodes while traversing the tree causes trouble)

    new_nodes = []

    for node in branch_list:

        TreeStruct._nodenum += 1
        new = TreeStruct.TreeNode(label = "n%d" % TreeStruct._nodenum, time=time, parent=node.parent)
        node.parent = new

        new.daughters = [node]
        new.parent.daughters.remove(node)
        new.parent.daughters.append(new)

        new_nodes.append(new)

    AssignBranchLengths(root)

    return new_nodes

def _SliceList(node, time, branch_list):
    '''
    Create a list of branches that span "time".
    '''

    if node.daughters != None:
        for d in node.daughters:
            _SliceList(d, time, branch_list)

    if node.time > time and node.parent.time < time:
        branch_list.append(node)
