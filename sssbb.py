import numpy as np
import random
from pylab import pcolor, show, colorbar, xticks, yticks
import matplotlib as ml
import matplotlib.pyplot as plt
import time

# SSSBB Project : Samiran Rohit Jyotsana Priyanka


# States

# 1 : Receptor Present only
# 2 : Ligand Present only
# 3 : Ligand and receptor are both present but unbound
# 4 : Complex Formation
# 5 : Dimer formation





OldM=np.empty((100, 100))
OldM.fill(0)



Contacts=0
ComplexFormations=0
ComplexDissiciations=0
DimerFormations=0




def Populate():

    M.fill(0)


    for i in R:
        Rx=i[0]
        Ry=i[1]
        M[Rx,Ry]=1



    for j in L:
        Lx=j[0]
        Ly=j[1]
        M[Lx,Ly]=2


    for i in R:
        for j in L:
            if i==j:
               Rx=i[0]
               Ry=i[1]
               M[Rx,Ry]=3
               
    for i in C:
        Cx=i[0]
        Cy=i[1]
        M[Cx,Cy]=4

    for i in D:
        Dx=i[0]
        Dy=i[1]
        M[Dx,Dy]=5



















# Function to return neighbors of a point

def GetNeighbors(a,b):
    Neighbors=[[a+1,b],[a,b+1],[a+1,b+1],[a-1,b],[a-1,b-1],[a,b-1],[a+1,b-1],[a-1,b+1]]
    ToRemove=[]
    for Ns in Neighbors:
        if (Ns[0]<0 or Ns[0]>99 or Ns[1]<0 or Ns[1]>99):
            ToRemove.append(Ns)
    return [x for x in Neighbors if x not in ToRemove]



# Store Random Numbers


R=[]
for i in range(1000):
    R.append(random.random())






# Initialize Master Array


M = np.empty((100, 100))
M=M.astype(int)
M.fill(0)





# Generate Receptors and populate the master matrix

n=0
R=[]
while (n<200):
    x=random.randint(0,99)
    y=random.randint(0,99)

    if [x,y] not in R:
        R.append([x,y])
        n=n+1



np.bincount(M.flatten())


# Generate Ligands and

n=0
L=[]

while (n<100):
    x=random.randint(0,99)
    y=random.randint(0,99)
    if [x,y] not in L:
        L.append([x,y])
        n=n+1




#Store Complex

C=[]





#Store Dimers

D=[]

#fig = plt.figure(figsize=(7, 7))

#plt.imshow(M)
#plt.ion()
#plt.show(block=False)

# Populate the master matrix





AssociationCounter=[]
DissociationCounter=[]
ContactCounter=[]
BoundAntigenCounter=[]
DimerCounter=[]

Populate()


pcolor(M)

show()


for i in range(10000):
    AssociationCounter.append(ComplexFormations)
    DissociationCounter.append(ComplexDissiciations)
    ContactCounter.append(Contacts)
    BoundAntigenCounter.append(len(C))
    DimerCounter.append(len(D))




   



    # Diffusion Step Pdiff=1

    for receptor in R:
       Neighbors=GetNeighbors(receptor[0],receptor[1])
       Slot=random.choice(Neighbors)
       if (M[Slot[0]][Slot[1]]==0 or M[Slot[0]][Slot[1]]==2):
           R.remove(receptor)
           R.append(Slot)
           
    for ligand in L:
       Neighbors=GetNeighbors(ligand[0],ligand[1])
       Slot=random.choice(Neighbors)
       if (M[Slot[0]][Slot[1]]==0 or M[Slot[0]][Slot[1]]==1):
           L.remove(ligand)
           L.append(Slot)

    


    # Calculate Contacts
    for k in R:
        for p in L:
            if k==p:
               Rx=p[0]
               Ry=p[1]
               if not (M[Rx,Ry]==3):
                   Contacts=Contacts+1




    Populate()

    temp=[]


    #Dissociation Step
    for c1 in C:
        if (random.random()<0.001): #Poff
            Cx=c1[0]
            Cy=c1[1]
            M[Cx,Cy]=3
            C.remove(c1)
            R.append(c1)
            L.append(c1)
            temp.append(c1)
            ComplexDissiciations=ComplexDissiciations+1












    #Complex forming step


    
    


    for i1 in R:
        for j1 in L:
            if i1==j1:
                if i1 not in temp: # Important step: the Complexes which have just dissociated, must not form complexes again in the same step
                    if (random.random()<0.1): #Pon=0.1
                        Rx=i1[0]
                        Ry=i1[1]
                        M[Rx,Ry]=4
                        L.remove(j1)
                        R.remove(i1)
                        C.append(j1)
                        ComplexFormations=ComplexFormations+1





    #Dimer forming step
    for c2 in C:
        Neighbors=GetNeighbors(c2[0],c2[1])
        dimer=False
        for N in Neighbors:
            if (M[N[0],N[1]]==4):
                dimer=True
                
                M[N[0],N[1]]=5
                DimerFormations=DimerFormations+1
                D.append(N)
                C.remove(N)
        if (dimer==True):
            M[c2[0],c2[1]]=5
            C.remove(c2)
            D.append(c2)
            DimerFormations=DimerFormations+1

    
                
    Populate()        
   
print "Steps : 10000"
print "Pon: 0.1"
print "Poff: 1"
print "Pdiff: 1"
print "BCR: ",len(R)
print "Antigens: ",len(L)

print "Contacts: ", Contacts
print "Complexes formed: ",ComplexFormations
print "Complex Dissiciations: ",ComplexDissiciations
print "Number of Bound Antigen: ",len(C)
print "Dimer Formations: ",len(D)
print 

print



from pylab import *
plot(AssociationCounter,'-b',label='Bindings')
plot(DissociationCounter,'-r',label='Unbindings')

legend(loc='upper right')
xlabel("Steps")

title("BCR-Antigen Bindings")
show()



plot(BoundAntigenCounter,'-b',label='Bound Antigens')


legend(loc='upper right')
xlabel("steps")

title("Number of Bound Antigen")
show()



plot(DimerCounter,'-b',label='Dimer Formations')


legend(loc='upper right')
xlabel("steps")

title("Number of Dimers")
show()










pcolor(M)

show()








