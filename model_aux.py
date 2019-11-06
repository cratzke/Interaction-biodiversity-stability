from numpy import *
from scipy import *
from numpy.random import random_sample
from scipy.integrate import odeint


def gaus(p,po,sig):
     return exp(-(p-po)**2/2./sig**2)
    
def Ma(P,t,po,pc,cp,K):
    kgrowth=10
    kdeath=10
    return list(r_[P[:-1]*(1-P[:-1])*(heaviside(-abs(P[-1]-po)+pc,0.5)*(kgrowth+kdeath)-kdeath),   
                    dot(K*cp,P[:-1])])

def Ma_gauss(P,t,po,pc,cp,K):
    kgrowth=10
    kdeath=10
    return list(r_[P[:-1]*(1-P[:-1])*(gaus(P[-1],po,pc)*(kgrowth+kdeath)-kdeath),   
                    dot(K*cp,P[:-1])])


def simul_gauss(NrSp,interaction_strength):
    #parameters
    K=1e10

    #start values
    P0=[10**-2]*NrSp
    P0.append(7)
    P0=array(P0).astype('float')

    #interaction and reaction parameters
    po=random_sample(size=NrSp)*5.+4.5
    pc=(array([2.5]*NrSp)).astype('float')
    cp=2*(random_sample(size=NrSp)-0.5)*interaction_strength 


    t=linspace(0,1,10000)

    MA_all=[]
    start=copy(P0)
    MA_final=[start]


    for i in range(80):
        MA = odeint(Ma_gauss, P0, t, args=(po,pc,cp,K))
        MA_all.extend(MA)
        MA_final.append(MA[-1])

        P0[:-1]=MA[-1,:-1]/10.     # dilute 
        P0[P0<10**(-9)]=0
        P0[-1]=0.9*7+0.1*MA[-1,-1]
    return array(MA_all)


def simul(NrSp,interaction_strength):
    #parameters
    K=1e10

    #start values
    P0=[10**-2]*NrSp
    P0.append(7)
    P0=array(P0).astype('float')

    #interaction and reaction parameters
    po=random_sample(size=NrSp)*5.+4.5
    pc=(array([2.5]*NrSp)).astype('float')
    cp=2*(random_sample(size=NrSp)-0.5)*interaction_strength 


    t=linspace(0,1,10000)

    MA_all=[]
    start=copy(P0)
    MA_final=[start]


    for i in range(80):
        MA = odeint(Ma, P0, t, args=(po,pc,cp,K))
        MA_all.extend(MA)
        MA_final.append(MA[-1])

        P0[:-1]=MA[-1,:-1]/10.     # dilute 
        P0[P0<10**(-9)]=0
        P0[-1]=0.9*7+0.1*MA[-1,-1]
    return array(MA_all)

def simulbis(P00,po,pc,cp):
    #parameters
    K=1e10
    NrSp=len(P00)
    P0=hstack([P00,7])
    
    t=linspace(0,1,100000)

    MA_all=[]
    start=copy(P0)
    MA_final=[start]


    for i in range(80):
        MA = odeint(Ma, P0, t, args=(po,pc,cp,K))
        #MA_all.extend(MA)
        #MA_final.append(MA[-1])
        P0[:-1]=MA[-1,:-1]/10.     # dilute 
        #P0[P0<0.1]=0
        P0[-1]=0.9*7+0.1*MA[-1,-1]
    return array(MA[-1,:-1])

def simulbis2(P00,po,pc,cp):
    #parameters
    K=1e10
    NrSp=len(P00)
    P0=hstack([P00,7])
    
    t=linspace(0,1,100000)

    MA_all=[]
    start=copy(P0)
    MA_final=[start]
    daily_fin=[]

    for i in range(80):
        MA = odeint(Ma, P0, t, args=(po,pc,cp,K))
        #MA_all.extend(MA)
        #MA_final.append(MA[-1])
        P0[:-1]=MA[-1,:-1]/10.     # dilute 
        #P0[P0<0.1]=0
        P0[-1]=0.9*7+0.1*MA[-1,-1]
        daily_fin.append(MA[-1,:-1])
    return array(daily_fin)



def simulbisc(interaction_strength,P00,po,pc,cp):
    #parameters
    K=1e10

    NrSp=len(P00)
    P0=hstack([P00,7])
    
  


    t=linspace(0,1,10000)

    MA_all=[]
    start=copy(P0)
    MA_final=[start]

    daily_fin=[]
    for i in range(80):
        MA = odeint(Ma, P0, t, args=(po,pc,cp,K))

        P0[:-1]=MA[-1,:-1]/5.     # dilute 
        #P0[P0<0.1]=0
        P0[-1]=0.8*7+0.2*MA[-1,-1]
        daily_fin.append(MA[-1,:-1])
    return array(daily_fin)

def diversity(x):
    #x[0] and x[1] are fractions of species 1 and 2
    # if both species are present
    if all(x>0):
        return exp(-sum(x*log(x)))
    # if one species dominates
    if  any(x!=0):
        return 1 #exp(0)
    if all(x==0):
        return 0
    
def diversity_com(x):
 
    # if at least two species are present
    if sum(x!=0)>1:
        ## take care here to not take log of 0 !!
        
        return exp(-sum(x[x>0]*log(x[x>0])))
    # if one species dominates
    if  sum(x!=0)==1:
        return 1 
    #no species survives
    if sum(x!=0)==0:
        return 0
