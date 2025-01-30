import matplotlib.pyplot as plt
import numpy as np
import random
import math

##################################################################################### Fonctions #################################################################################

def distance_euclidienne(x1,y1,x2,y2):
  return np.sqrt(((x2-x1)**2)+((y2-y1)**2))

def projete(Cylindre, i, j, k):
    vectDirect = np.array([Cylindre[0][j] - Cylindre[0][i], Cylindre[1][j] - Cylindre[1][i]])
    
    # Check if vectDirect is a zero vector (i.e., points i and j are the same)
    if np.all(vectDirect == 0):
        return Cylindre[0][i], Cylindre[1][i]  # Return point i as the projection (since i and j are the same)
    
    vectProj = np.array([Cylindre[0][k] - Cylindre[0][i], Cylindre[1][k] - Cylindre[1][i]])
    dot_product = np.dot(vectProj, vectDirect)
    norm_square = np.dot(vectDirect, vectDirect)
    
    # Avoid division by zero
    if norm_square == 0:
        return Cylindre[0][i], Cylindre[1][i]  # Return point i as the projection (same as above)
    
    proj = (dot_product / norm_square) * vectDirect
    point_proj = (Cylindre[0][i], Cylindre[1][i]) + proj
    return point_proj

def entre(xa,ya,xc,yc,i,j,k):
  a=distance_euclidienne(xa,ya,projete(Cylindre,i,j,k)[0],projete(Cylindre,i,j,k)[1])
  b=distance_euclidienne(projete(Cylindre,i,j,k)[0],projete(Cylindre,i,j,k)[1],xc,yc)
  c=distance_euclidienne(xa,ya,xc,yc)
  return a+b-c<0.001

def calcul_angle(x,y,x_obj,y_obj):
  global orientation_actuelle
  if y_obj-y>0 and x_obj-x>0:
    theta= np.arctan((y_obj-y)/(x_obj-x))*180/np.pi
  if y_obj-y<0 and x_obj-x>0:
    theta= -np.arctan(-(y_obj-y)/(x_obj-x))*180/np.pi
  if y_obj-y>0 and x_obj-x<0:
    theta= 180-np.arctan(-(y_obj-y)/(x_obj-x))*180/np.pi
  if y_obj-y<0 and x_obj-x<0:
    theta= 180+np.arctan((y_obj-y)/(x_obj-x))*180/np.pi

  orientation=theta-orientation_actuelle
  if orientation>180:
    orientation=orientation-360
  if orientation<-180:
    orientation= orientation+360
  orientation_actuelle=orientation_actuelle+orientation
  return str('TURN ')+str(orientation)

def somme_chemin(liste, M):
    m = 0
    for i in range(len(liste) - 1):
      m += M[liste[i]][liste[i + 1]]
    m += M[liste[-1]][liste[0]]
    return m

def Temperature(t):
    return 10/np.log(t)

def Pck_cest_notre_Pb(X,M,N):
    Xmin=X
    Valeurmin=somme_chemin(X, M)
    Valeuractuelle=somme_chemin(X, M)
    for i in range(N):
        Xp=X.copy()
        random.shuffle(Xp)
        XP=X.copy()
        XP[Xp[0]],XP[Xp[1]]=XP[Xp[1]],XP[Xp[0]]
        Vtransition=somme_chemin(XP, M)
        rot=np.exp(1/Temperature(i+2)*(Valeuractuelle-Vtransition))
        U=random.uniform(0,1)
        if(rot>=1 or U<rot):
            X=XP.copy()
            Valeuractuelle=Vtransition
            if(Valeuractuelle<Valeurmin):
                Xmin=X
                Valeurmin=Valeuractuelle
    return [X,Valeurmin,Xmin]


##################################################################################### Création cylindres#########################################################################

donnees=np.loadtxt("C:\CHALLENGE\donnees-map.txt")

Cylindre=[[],[],[]]

for i in range(len(donnees)):
  Cylindre[0].append(donnees[i][0])
  Cylindre[1].append(donnees[i][1])
  Cylindre[2].append(donnees[i][2])


valeur1=[[],[]]
valeur2=[[],[]]
valeur3=[[],[]]
X1=[]
X2=[]

for i in range(len(Cylindre[0])):
  if Cylindre[2][i]==1.0:
    valeur1[0].append(Cylindre[0][i])
    valeur1[1].append(Cylindre[1][i])
    X1.append(i)
  if Cylindre[2][i]==2.0:
    valeur2[0].append(Cylindre[0][i])
    valeur2[1].append(Cylindre[1][i])
    X2.append(i)
  if Cylindre[2][i]==3.0:
    valeur3[0].append(Cylindre[0][i])
    valeur3[1].append(Cylindre[1][i])
    X2.append(i)

x=Cylindre[0]
y=Cylindre[1]
valeur=Cylindre[2]


#################################################################################################################################################################################



M1=np.array([np.zeros(len(valeur1[0])) for i in range(len(valeur1[0]))])
for i in range(len(valeur1[0])):
  for j in range(len(valeur1[0])):
    M1[i][j]=distance_euclidienne(valeur1[0][i],valeur1[1][i],valeur1[0][j],valeur1[1][j])


M=np.array([np.zeros(20) for i in range(20)])
for i in range(len(x)):
  for j in range(len(y)):
    M[i][j]=distance_euclidienne(Cylindre[0][i],Cylindre[1][i],Cylindre[0][j],Cylindre[1][j])

M2=np.array([np.zeros(len(X2)) for i in range(len(X2))])
for i in range(len(X2)):
  for j in range(len(X2)):
    M2[i][j]=distance_euclidienne(Cylindre[0][X2[i]],Cylindre[1][X2[i]],Cylindre[0][X2[j]],Cylindre[1][X2[j]])


P=np.array([np.zeros(20) for i in range(20)])

c=0
for i in range(len(Cylindre[0])):
  for j in range(len(Cylindre[0])):
    P[i][j]=1
    if i==j:
      P[i][j]=0
    for k in range(len(Cylindre[0])):
      if k!=i and k!=j:
        if entre(Cylindre[0][i],Cylindre[1][i],Cylindre[0][j],Cylindre[1][j],i,j,k):
          c=c+1
          if distance_euclidienne(Cylindre[0][k],Cylindre[1][k],projete(Cylindre,i,j,k)[0],projete(Cylindre,i,j,k)[1])<0.5:
            P[i][j]=0


orientation_actuelle=0

L=[0, 1, 2, 3, 7, 6, 10, 9, 13, 14, 11, 15, 19, 18, 17, 16, 12, 8, 4, 5]

Turn=[calcul_angle(0,0,Cylindre[0][0],Cylindre[1][0])]

distance=[str('GO ')+str(distance_euclidienne(0,0,Cylindre[0][0],Cylindre[1][0]))]

for i in range(19):
  Turn.append(calcul_angle(Cylindre[0][L[i]],Cylindre[1][L[i]],Cylindre[0][L[i+1]],Cylindre[1][L[i+1]]))
  if i==18:
    distance.append(str('GO ')+str(distance_euclidienne(Cylindre[0][L[i]],Cylindre[1][L[i]],Cylindre[0][L[i+1]],Cylindre[1][L[i+1]])-0.49))
  else:
    distance.append(str('GO ')+str(distance_euclidienne(Cylindre[0][L[i]],Cylindre[1][L[i]],Cylindre[0][L[i+1]],Cylindre[1][L[i+1]])))
    
  

X=[i for i in range(20)]
X=Pck_cest_notre_Pb(X,M,1000000)[2]
X1_1=[i for i in range(len(X1))]
X2_1=[i for i in range(len(X2))]
X1_1=Pck_cest_notre_Pb(X1_1,M1,1000000)[2]
X2_1=Pck_cest_notre_Pb(X2_1,M2,1000000)[2]

X1_2=[i for i in range(len(X1))]
c=0
for i in X1_1:
  X1_2[c]=X1[i]
  c=c+1

X2_2=[i for i in range(len(X2))]
c=0
for i in X2_1:
  X2_2[c]=X2[i]
  c=c+1

print(Turn)
print(distance)

with open ("C:\CHALLENGE\script.txt", "w") as script:
    for i in range (len(Turn)):
        script.write(Turn[i])
        script.write("\n")
        script.write(distance[i])
        script.write("\n")
    script.write(str('FINISH'))