#!/usr/bin/python
import sys,os
from optparse import OptionParser
import copy
'''
Thomas Harrison
Motif Finder
'''
l=6
d=1

class lmerNode:
    '''Graph node Represented as an adjacency list'''
    def __init__(self,_lmer,_readnum):
        self.lmer=_lmer
        self.readnum=_readnum
        self.edges=[]
        self.color="w"
        self.l=[]

    def __repr__(self):
        return(self.lmer+" "+str(self.readnum) )

                   

def getReads(inp):
    '''Reads in a fasta file and stores reads in a list'''
    fh = open(inp,"r")
    reads = []
    for line in fh:
       
        if ">" in line:
             reads.append("")
        else:
             reads[len(reads)-1]+= (line.rstrip()) 
    return reads         
    
    
    
def lmerfy(reads,l):
    '''Converts reads into lmers'''
    lmers = []
    i=0
    for read in reads:
        beg = 0
        end=l
        while not(end==len(read)+1):
            lmers.append(lmerNode(read[beg:end],i))
            beg+=1
            end+=1
        i+=1    
    return lmers       
    
def constructGraph(lmers,d):
    '''Construct lmer graph O(n^2)'''
    
    for lmer_lhs in lmers:
         
        for lmer_rhs in lmers:
            # condition lmer is hamming distance 2d away
            if not (lmer_rhs.readnum==lmer_lhs.readnum):
                diffs = 0
                for ch1, ch2 in zip(lmer_lhs.lmer, lmer_rhs.lmer):
                    if ch1 != ch2:
                        diffs += 1
                if diffs <= 2*d:        
                    lmer_rhs.edges.append(lmer_lhs)

        

def findCliques(lmers,r,d,l):
    '''Get starting vertices and feed to DP algorithm'''
    i=0
    for item in lmers:
        if not item.readnum == 0:
            break
        newAlgo(item,r,d,l)
        i+=1
        
        print "Explored "+str(i)+" Nodes"
        #print motifList



def newAlgo(lmer,r,d,l):
    '''DP algorithm'''
    ###construct###
    s=[]
    for i in range(0,r):
        s.append([])
    s[0].append(lmer)
    for item in lmer.edges:
        s[item.readnum].append(item)
    
    
    ###################traverse############################
    s[0][0].l=[]
    for item in s[1]:
        item.l =[s[0][0]]
    for i in range(2,r):
        
        for j in range(0,len(s[i])):
            s[i][j].l=[]
            for k in range(0,len(s[i-1])-1):
            
                if findHammingPlain(s[i][j].lmer,s[i-1][k].lmer)<=2*d:
                     #union of s[i][j] and s[i-1][k]
                     temp= list(set(s[i][j].l) | set(s[i-1][k].l))
                     s[i][j].l=temp
                     s[i][j].l.append(s[i-1][k])
                     if len(s[i][j].l)>=r-1 and i == r-1:
                         
                         temp=copy.copy(s[i][j].l)
                         temp.append(copy.copy(s[i][j]))
                         # case multiple instances of the same level in l[i][j]
                         bins=nFilter(temp,r)       
                         if not bins == None:
                             widdleBins(bins,d)    
                             searchMotif(bins,0,[None]*len(bins),d,l)
                             break
                             
                             





def widdleBins(bins,d):
    '''Given the l[-1] element find the right vertices'''
    temp=copy.copy(bins)
    i=0
    j=0
    for i in range(len(bins)):
        for j in range(len(bins[j])):
            
            if findHammingPlain(temp[i][j].lmer,bins[-1][0].lmer) > 2*d:
                bins[i].remove(bins[i][j])










def findHammingPlain(s1,s2):
    diffs = 0
    for ch1, ch2 in zip(s1, s2):
        if ch1 != ch2:
            diffs += 1
    return diffs 


                
    
    
        
                
                
      
        

    
def nFilter(greys,n):
    '''Creates bins of certain depths'''
    binVector=[False]*n
    for item in greys:
        if binVector[item.readnum]== False:
            binVector[item.readnum]=[]
        binVector[item.readnum].append(item)

    if not False in binVector:
        return binVector
    return None        





def searchMotif(bins,depth,current,d,l):
   '''Recursively goes through bins combination'''
   for item in bins[depth]:
       current[depth]=item
       if depth == len(bins)-1:
           findD(current,d,l)
       else:
           searchMotif(bins,depth+1,current,d,l)    

def findD(current,d,l):
    '''Finds consensus'''
    consensus=""
    i=0
    for i in range(l):
        A=0
        G=0
        C=0
        T=0
        N=0
        for item in current:
             
             if item.lmer[i]=='A':
                 A+=1
             if item.lmer[i]=='G':
                 G+=1
             if item.lmer[i]=='C':
                 C+=1
             if item.lmer[i]=='T':
                 T+=1 
             if item.lmer[i]=='N':
                 N+=1
        #print "calMax"         
        maxCount=max([A,G,C,T,N])
        if (maxCount==A):
            consensus+="A"
            continue 
        if (maxCount==G):
            consensus+="G"
            continue 
        if (maxCount==C):
            consensus+="C"
            continue 
        if (maxCount==T):
            consensus+="T"
            continue
        if (maxCount==N):
            consensus+="N"
        i=0    
    findHamming(consensus,current,d)                          
                                    
def findHamming(consensus,current,d):
    '''get hamming distance of the consensus with all others'''
    for item in current:
        diffs = 0
        for ch1, ch2 in zip(item.lmer, consensus):
            if ch1 != ch2:
                diffs += 1
            if diffs > d:
               
                return False
                break
    #report
    f=open("output_"+str(l)+"_d_"+str(d)+".txt","a")
   
    f.write(consensus+"\t"+str(current)+"\n")
    f.close
    return True

      
            
                
         
        
  



if __name__ == "__main__":
   

    
    inp = None
    for i in range(0,len(sys.argv)):
        if sys.argv[i]=="-input":
            inp= sys.argv[i+1]
        if sys.argv[i]=="-l":
            l= int(sys.argv[i+1])
        if sys.argv[i]=="-d":
            d= int(sys.argv[i+1])
    
   
   
    
  
    ################## build graph
    
    reads = getReads(inp)
    # l-merfy
    lmers=lmerfy(reads,l)
    print ("Connecting "+str(len(lmers))+" Nodes")
    constructGraph(lmers,d)
    
    print "Looking for Cliques"
    findCliques(lmers,len(reads),d,l)

   
  

   
    
        
    


        
        










    
    
    
    
    

