from TreeStruct import TreeNode
import sys
import re

class NewickError(Exception):
    def __init__(self, value):
        self.value = "ERROR: " + value
    def __str__(self):
        return repr(self.value)

#--------------------------------------------------
# for reading just Newick strings
#-------------------------------------------------- 

def Read(tree_string):
    '''
    Take a Newick string and return the root of the tree built from it.
    (note on tip/node labels: they are honored, and even spaces are okay)
    '''

    ### check if the Newick string seems well-formed
    try:
        Check(tree_string)
    except (NewickError, error):
        print(error.value)
        print("Bad Newick string.  Aborting in Newick.Read()...")
        sys.exit(1)

    # the first live char must be a (, so skip to there
    index = tree_string.find('(')
    root = TreeNode(parent=None)
    current_node = root

    ### now start building the tree

    index += 1
    next_char = tree_string[index]
    while next_char != ';':

        if next_char == ' ':
            pass

        # add an internal node, or a tip and its label
        elif next_char == '(':
            new_node = TreeNode(parent=current_node)
            if current_node.daughters == None:
                current_node.daughters = [new_node]
            else:
                current_node.daughters.append(new_node)
            current_node = new_node

        # back up to parent so a new sister can be added
        elif next_char == ',':
            current_node = current_node.parent

        # back up to parent, and possibly get its label
        elif next_char == ')':
            current_node = current_node.parent
            if tree_string[index+1] not in [',',')','(',';',':']:
                gotcha = _ToNextSymbol(tree_string, index+1)
                index = gotcha[1]
                current_node.label = gotcha[0]

        # get branch length
        elif next_char == ':':
            gotcha = _ToNextSymbol(tree_string, index+1)
            index = gotcha[1]
            current_node.length = float(gotcha[0])

        # at the start of a label, so add that node
        else:
            gotcha = _ToNextSymbol(tree_string, index)
            index = gotcha[1]
            new_node = TreeNode(parent=current_node, label=gotcha[0])
            if current_node.daughters == None:
                current_node.daughters = [new_node]
            else:
                current_node.daughters.append(new_node)
            current_node = new_node

        index += 1
        next_char = tree_string[index]

    # when the tree_string has surrounding (...), the current root has the 
    #   real root as its only daughter
    if len(root.daughters) == 1:
        root = root.daughters[0]
        root.parent = None
    if root.length == None:
        root.length = 0.0
    return root

def _ToNextSymbol(tree_string, start_index):

    stop_index = start_index+1
    next_char = tree_string[stop_index]
    while next_char not in [',',')','(',';',':']:
        stop_index += 1
        next_char = tree_string[stop_index]
    str = tree_string[start_index:stop_index]
    return [str, stop_index-1]

def Check(nstr):
    '''
    Check for various malformations of a Newick string.
    Raise a NewickError if there is a problem.
    '''
    nstr = nstr.strip()
    if nstr[-1] != ";":
        raise NewickError("Newick strings must end with a ;")
    if nstr.count(";") != 1:
        raise NewickError("too many ; in Newick string")
    if nstr.count("(") != nstr.count(")"):
        raise NewickError("mismatched ( ) in Newick string")
    if nstr.find(")(") != -1:
        raise NewickError("shouldn't have )( in Newick string")
    if nstr.find("()") != -1:
        raise NewickError("shouldn't have () in Newick string")

def ReadFromFile(filename):
    '''
    From the specified file, try to read in the first line as a Newick string.
    Blank lines and ones beginning with # or [ are skipped.
    Returns the root of the tree, or None if no tree was formed.
    '''

    try:
        infile = open(filename, "r")
    except IOError:
        print("ERROR: problem opening the file %s" % (filename))
        return None

    # get the first non-blank, non-comment line in the file
    while 1:
        try:
            line = infile.readline().strip()
            if line and line[0]!= "#" and line[0]!="[":
                break
                # found a useful line, so stop looking
        except StopIteration:
            print("ERROR: end of file reached before finding a possible tree description.")
            line = None
            break

    infile.close()

    # if a good line was found, try to read it in as a tree, and return the root
    if line != None:
        try:
            root = Read(line)
            return root
        except NewickError:
            pass
    else:
        return None

#--------------------------------------------------
# for reading my .ttn files
#-------------------------------------------------- 

def ReadFromFileTTN(filename, floatstates=False):
    '''
    Read in one of my .ttn files, with relaxed assumptions:
        The tree string is on a single line.
        The tip and node states come afterwards, in any order.
        The states can have integer or floating (floatstates=True) values,
           or can be NA or None.
    The comment character # is respected, and blank lines are skipped.
    '''

    # First, form a tree from the Newick string.
    root = ReadFromFile(filename)

    if root == None:
        nstates = 0

    # Then, deal with the state information if a tree was made.
    else:
        state_dict = _MakeStateDict(filename)
        nstates = len(set(state_dict.values()))
        _PutStates(root, state_dict, nstates, floatstates)

    return root

def _MakeStateDict(filename):
    '''
    Convert the states in a .ttn file into a dictionary.
    '''

    # already did error-checking in ReadFromFile()
    infile = open(filename, "r")
    while 1:
        line = infile.readline().strip()
        if line and line[0]!= "#" and line[0]!="[":
            break
            # now at the tree string

    # put the state values in a dictionary
    state_dict = {}
    while 1:
        try:
            line = next(infile).strip()
            line = line.partition("#")[0].strip()
            if line and line[0]!= "#" and line[0]!="[":
                (name, state) = line.split(None, 1)
                if name in state_dict:
                    print("WARNING: label %s is used more than once" % (name))
                state_dict[name] = state

        except ValueError:
            infile.close()
            raise NewickError("Problem reading character states.  " + \
                    "Something is wrong with:\n   " + line + \
                    "\nProper format is:\n   label1   state1")

        except StopIteration:
            break

    infile.close()
    return state_dict

def _PutStates(node, state_dict, nstates, floatstates):
    '''
    Give tips their states, from a dictionary.
    '''

    try:
        state = state_dict[node.label]

        if len(state) != 1:
            raise NewickError('too many states for tip "' \
                                + str(node.label) + '"')

        if floatstates:
            state = float(state[0])
        else:
            state = int(state[0])
            if state not in range(nstates):
                raise NewickError('invalid state "' + str(state) + '" for tip "' \
                                    + str(node.label) + '"')

        node.state = state

    except KeyError:
        # raise NewickError("can't find a state for " + node.label)
        # or just print a message about number of states/nodes found
        pass

    if node.daughters != None:
        for d in node.daughters:
            _PutStates(d, state_dict, nstates, floatstates)
