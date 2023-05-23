from io import StringIO
import TreeExtra

Nstates = 2
_nodenum = 1

class TreeNode:
    '''
        TreeNode contains the properties of a single node (or tip) in a 
        phylogenetic tree
           label: a name or number for the node
           time: the time (on the tree) of the node
           state: the state of a particular character
           parent: the node that is this node's ancestor 
                      (None for the root)
           daughters: a list of this node's immediate descendants
                         (None for a tip)
           length: the time from this node to its ancestor
                      (computed automatically if times are specified)
    '''
    def __init__(self, label=None, time=None, length=None, state=None, \
                 parent=None, daughters=None):
        self.label = label
        self.time = time
        self.state = state
        self.parent = parent
        self.daughters = daughters
        # use node times if given; otherwise, use branch lengths if given
        if self.parent!=None and self.time!=None and self.parent.time!=None:
            self.length = self.time - self.parent.time
        else:
            self.length = length
        self.cl = [0]*Nstates    # will hold conditional likelihood values
        self.lq = 0              # will hold log-compensation for all branches above this node
        self.fixed = False       # will note if the node is fixed

    def PrintNode(self):
        ''' prints information about this node '''
        print(self.label, end=" ")
        print(":", end=" ")
        if self.time != None:
            print("t = %2.4f," % (self.time), end=" ")
        if self.length != None:
            print("l = %2.4f," % (self.length), end=" ")
        if self.state != None:
            print("s = %s," % (str(self.state)), end=" ")
        if self.cl != None:
            print("cl = [%1.4f,%1.4f]," % (self.cl[0], self.cl[1]), end=" ")
        if self.lq != None:
            print("lq = %2.4f," % (self.lq), end=" ")
        print("p =", end=" "),
        if self.parent != None:
            print(self.parent.label, end=" ")
        else:
            print("--", end=" ")
        print(", d =", end=" ")
        if self.daughters != None:
            for d in self.daughters:
                print(d.label, end=" ")
        else:
            print("--", end=" ")
        print("")

    def PrintTree(self, indent=2):
        ''' prints a list of all descendants from this node '''
        print(" "*indent, end=" ")
        self.PrintNode()
        if self.daughters != None:
            for d in self.daughters:
                d.PrintTree(indent+2)

    def NewickString(self, outparen=True):
        ''' returns all descendants from this node as a Newick string '''
        nstr = StringIO()   # allows for more efficient string concatenation
        if outparen:
            nstr.write("(")
        self._NewickStringGuts(nstr)
        if outparen:
            nstr.write(")")
        nstr.write(";")
        returnme = nstr.getvalue()
        nstr.close()
        return returnme

    def _NewickStringGuts(self, nstr):
        if self.daughters == None:
            nstr.write( str(self.label) )
            if self.length != None:
                nstr.write( ":%f" % (self.length) )
        else:
            nstr.write("(")
            for d in self.daughters:
                d._NewickStringGuts(nstr)
                nstr.write(",")
            nstr.seek(nstr.tell()-1, 0)    # remove that last comma
            nstr.write(")")
            nstr.write( str(self.label) )
            if self.length != None:
                nstr.write( ":%f" % (self.length) )


    def TipStates(self):
        ''' returns a list of tip labels and trait values in left-to-right order '''
        tiplist = []
        self._TipStatesGuts(tiplist)
        return tiplist

    def _TipStatesGuts(self, tiplist):
        if self.daughters == None:
            tiplist.append([self.label, self.state])
        else:
            for d in self.daughters:
                d._TipStatesGuts(tiplist)

    def Age(self, ultrametric=False):
        '''
        returns the greatest distance between the root and any tip
        (take a shortcut if all tips are known to be the same age)
        '''

        if self.time == None:
            print("ERROR: node times are needed by Age()")
            return None

        if ultrametric:
            node = self
            while node.daughters != None:
                node = node.daughters[0]
            max_time = node.time
        else:
            distances = []
            TreeExtra.GetTipDistances(self, distances)
            max_time = max(distances)

        return abs(max_time - self.time)
