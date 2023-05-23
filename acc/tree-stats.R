library(ape)

### CODE FROM ETHAN 2019-09-23 ###

# The attached function takes a tree (in ape format--you need to have the ape library loaded) with tip labels given by the args lab1 and lab2 and returns a list with the root state, the number of monophyletic clades in lab1 and lab2, the topological class, and some other stuff. All the args after lab2 can be ignored.

#-----------------------------------------------
# Calculates the 
# 1) # of mutations in the tree
# 2) root label 
# 3) topological class
# 4) and number of monophyletic clades of class A and B

num.intro <- function(tr, lab1, lab2, mu=0.80e-2/365, sig=0.14e-2/365, seq.len=346, gen.len=2){
  if(!is.rooted(tr)){print("Not Rooted YO"); return(NULL)}
  if(!all(tr$tip.label%in%c(lab1, lab2))){print("Tip labels not covered"); return(NULL)}
  root = Ntip(tr)+1
  ind1 = which(tr$tip.label==lab1)
  ind2 = which(tr$tip.label==lab2)
  #if(is.monophyletic(tr, ind1) & is.monophyletic(tr, ind2)) return(list(anc=NULL, nintro=0))
  lab.vect = c(tr$tip.label, rep(NA, Nnode(tr)))
  anc.vect = c(tr$tip.label, rep(NA, Nnode(tr)))
  names(anc.vect) <- c(tr$tip.label, tr$node.label)
  for(x in (Ntip(tr)+Nnode(tr)):root){
    kids.id = tr$edge[tr$edge[,1]==x,2]
    kids.lab = lab.vect[kids.id]    
    kids.anc = anc.vect[kids.id]
    #--------------------Propagation for number monophyletic clades
    # lables are the same, propagate up
    if(kids.lab[1]==kids.lab[2]){
      lab.vect[x] = kids.lab[1]
      lab.vect[kids.id] = '?'
    }else{
      # labels are A? in either order 
      if(setequal(kids.lab, c(lab1, "?"))){lab.vect[x] = '?'}
      # labels are B? in either order 
      if(setequal(kids.lab, c(lab2, "?"))){lab.vect[x] = '?'}
      # labels are AB in either order 
      if(setequal(kids.lab, c(lab1, lab2))){lab.vect[x] = "?"}
    }
    #--------------------Propagation for root label
    # lables are the same, propagate up
    if(kids.anc[1]==kids.anc[2]){
      anc.vect[x] = kids.anc[1]
    }else{
      # labels are A? in either order 
      if(setequal(kids.anc, c(lab1, "?"))){anc.vect[x] = lab1}
      # labels are B? in either order 
      if(setequal(kids.anc, c(lab2, "?"))){anc.vect[x] = lab2}
      # labels are AB in either order 
      if(setequal(kids.anc, c(lab1, lab2))){anc.vect[x] = "?"}
    }
  }
  # All trees have one or more AB coal, if more than one AB coal->PP topo, if root AB-> MM topo else PM topo 
  if(sum(anc.vect=='?')>1){
    topo='PP'
    nd = node.depth.edgelength(tr)
    nd = nd/max(nd)
    tIA = nd[which(lab.vect==lab1)]
    tIB = nd[which(lab.vect==lab2)]
    tIA.mean = mean(tIA)
    tIB.mean = mean(tIB)
    tIA.sd = sd(tIA)
    tIB.sd = sd(tIB)
  }else{
    nd = node.depth.edgelength(tr)
    nd = nd/max(nd)
    tIA = nd[which(lab.vect==lab1)]
    tIB = nd[which(lab.vect==lab2)]
    tIA.mean = mean(tIA)
    tIB.mean = mean(tIB)
    tIA.sd = 0
    tIB.sd = 0
    if(anc.vect[root]!='?')topo='PM' else topo = "MM" 
  }
  theta = sig^2/mu
  k = mu^2/sig^2 
  evo.rate = rgamma(1, k, scale=theta)
  t.dist=sum(tr$edge.length*gen.len*evo.rate*seq.len)
  
  d = cophenetic(tr)
  ind.A = which(rownames(d)=="MP2")
  ind.B = which(rownames(d)=="MP3")
  tmp.AB = d[ind.A,ind.B]
  tmp.A = d[ind.A,ind.A]
  tmp.B = d[ind.B,ind.B]
  pwd.AB = mean(tmp.AB[upper.tri(tmp.AB)]*gen.len*seq.len*evo.rate)
  pwd.A = mean(tmp.A[upper.tri(tmp.A)]*gen.len*seq.len*evo.rate)
  pwd.B = mean(tmp.B[upper.tri(tmp.B)]*gen.len*seq.len*evo.rate)
  
  
  if(any(tr$edge.length<0))print("Negative branch lengths in sims")
  
  return(list(anc=anc.vect[root], IA=unname(table(lab.vect)[lab1]), IB=unname(table(lab.vect)[lab2]), topo=topo, 
              t.dist=t.dist, pwd.AB=pwd.AB, pwd.A=pwd.A, pwd.B=pwd.B, tIA.mean=tIA.mean, tIA.sd=tIA.sd, tIB.mean=tIB.mean, tIB.sd=tIB.sd, anc.vect=anc.vect))
}
