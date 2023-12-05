#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 22:08:23 2023

@author: isaias
"""

ReverseFault
i=58

ReverseFault["Base"]
ReverseFault["EjeH1"]
ReverseFault["EjeH2"]
CanalH1=np.int(np.float(RecuperaOrientacion(ReverseFault["Base"][i])[8]))
CanalH2=np.int(np.float(RecuperaOrientacion(ReverseFault["Base"][i])[9]))
    
    
(X1,NS1,f1,Y1,VelMuest1)=fourier(ReverseFault["Base"][i], CanalH1) #Usar registrocrudo cuando no quiero normalizado y centrado
(X2,NS2,f2,Y2,VelMuest2)=fourier(ReverseFault["Base"][i], CanalH2) #Usar registrocrudo cuando no quiero normalizado y centrado




plt.plot(f1[(f1<15)*(f1>0.1)],Y1[(f1<15)*(f1>0.1)], label=r"$A_{h_1}$",alpha=0.5)
plt.plot(f2[(f1<15)*(f1>0.1)],Y2[(f1<15)*(f1>0.1)], label=r"$A_{h_2}$", alpha=0.5)
plt.legend( prop={'size': 15})
#plt.title("Registro: "+ReverseFault["Base"][i]+" R="+str(round(ReverseFault["R"][i],2))+r" $M_w$="+str(ReverseFault["GMW"][i]))
plt.xlabel("$\mathfrak{f} \quad $(Hz)")
plt.ylabel(r"$A(\mathfrak{f}) \quad (\frac{cm}{s})$")




plt.plot(np.log10(f1[(f1<15)*(f1>0.1)]),np.log10(Y1[(f1<15)*(f1>0.1)]), label=r"$A_{h_1}$",alpha=0.5)
plt.plot(np.log10(f2[(f1<15)*(f1>0.1)]),np.log10(Y2[(f1<15)*(f1>0.1)]), label=r"$A_{h_2}$", alpha=0.5)
#plt.legend( prop={'size': 15})
#plt.title("Registro: "+ReverseFault["Base"][i]+" R="+str(round(ReverseFault["R"][i],2))+r" $M_w$="+str(ReverseFault["GMW"][i]))
plt.xlabel(r"$\log_{10}(\mathfrak{f})$")
plt.ylabel(r"$\log_{10}(A(\mathfrak{f}))$")




ReverseFault["GMW"][i]



plt.plot(np.log10(f1[(f1<15)*(f1>0.1)]),np.log10(Y1[(f1<15)*(f1>0.1)]), label=r"$A_{h_1}$",alpha=0.5)
plt.plot(np.log10(f2[(f1<15)*(f1>0.1)]),np.log10(Y2[(f1<15)*(f1>0.1)]), label=r"$A_{h_2}$", alpha=0.5)
plt.plot(np.log10(f1[(f1<15)*(f1>0.1)]),np.log10(np.sqrt((Y1[(f1<15)*(f1>0.1)])**2/2+(Y2[(f1<15)*(f1>0.1)])**2/2)), label=r"$A_{h_1}$")


from csaps import csaps

soporte=np.log10(f1[(f1<15)*(f1>0.1)])
respuesta=np.log10(np.sqrt((Y1[(f1<15)*(f1>0.1)])**2/2+(Y2[(f1<15)*(f1>0.1)])**2/2))


plt.plot(np.log10(f1[(f1<15)*(f1>0.1)]),np.log10(np.sqrt((Y1[(f1<15)*(f1>0.1)])**2/2+(Y2[(f1<15)*(f1>0.1)])**2/2)), label=r"$A_{h}$")
plt.plot(soporte, csaps(soporte,respuesta ,soporte,  smooth=0.8),label=r"$A^s_{h}$")
plt.xlabel(r"$\log_{10}(\mathfrak{f})$")
plt.ylabel(r"$\log_{10}(A(\mathfrak{f}))$")
plt.legend( prop={'size': 15})


######################################



fmax=15
Muestras=20
smooth=0.8
####################################


ReverseFault["Base"]
ReverseFault["EjeH1"]
ReverseFault["EjeH2"]


flag=0
for i,j,m in zip(ReverseFault["Base"],ReverseFault["EjeH1"],ReverseFault["EjeH2"]):
    if flag==0:
        (X,NS,f,Y, VelMuest)=fourier(i,int(j))
        (X,NS,f,Y2, VelMuest)=fourier(i,int(m))
        MCMCBase=[np.sqrt(Y*Y2),f]
        flag=1
    else :
        (X,NS,f,Y, VelMuest)=fourier(i,int(j))
        (X,NS,f,Y2, VelMuest)=fourier(i,int(m))
        MCMCBase.extend([np.sqrt(Y*Y2),f])
        #MCMCBase.extend([Y,f])

###########################################################

# plt.plot(MCMCBase[1][:],MCMCBase[0][:])


frec=10**(np.linspace(np.log10(0.1),np.log10(fmax),Muestras))



AShBase=np.ones(0)
for i in range(len(ReverseFault)):
    # f=interp1d(MCMCBase[2*i+1],MCMCBase[2*i] )
    # MCMCBase[i]=f(frec)
    soporte=MCMCBase[2*i+1]
    respuesta=MCMCBase[2*i]
    respuesta=respuesta[(soporte<15)*(soporte>=0.09)]
    soporte=soporte[(soporte<15)*(soporte>=0.09)]
    ASh=csaps(np.log10(soporte),np.log10(respuesta) ,np.log10(frec),  smooth=0.8)

    if len(AShBase)==0:
        AShBase=ASh
    else :
        AShBase=np.vstack((AShBase,ASh))

        
# for i in range(len(ReverseFault)):
#     plt.plot(np.log10(frec),AShBase[i])        
###########################################################
        




def Q(f,ca,cb):
    return ca*f**cb

# def G(R,Rx,alphacuad):
#     return (1/Rx)*(R/Rx)**alphacuad

def G(R,alphacuad):
    Rx=1
    return (1/Rx)*(R/Rx)**alphacuad


def S(f,M,beta,deltasigma):
#    M0=np.exp(3/2*(M+10.7))
    M0=10**(3/2*(M+10.7))
    fc=4.91*10**6*beta*(deltasigma/M0)**(1/3)
#    fc=4.91*10**6*(deltasigma/M0)**(1/3)
    return M0*f**2/(1+(f/fc)**2)

def exponencial(f,R,beta,ca,cb):
    return np.exp(-np.pi*f*R/(beta*Q(f,ca,cb)))

# def A(f,M,R,Rho,beta,P,Rtf,Fs,Rx,deltasigma,alphacuad):
#     Const=((2*np.pi)**2*Rtf*Fs*P)/(4*np.pi*Rho*beta**3)
#     return Const*S(f,M,beta,deltasigma)*exponencial(f,R,beta)*G(R,Rx,alphacuad)/(10**20)

def A(f,M,R,beta,deltasigma,alphacuad,Const,ca,cb):
    return Const/beta**3*S(f,M,beta,deltasigma)*exponencial(f,R,beta,ca,cb)*G(R,alphacuad)/(10**20)





def EnSupp(Theta):
    Verificacion=all(0 < Theta) and Theta[0]<100 #Const
    Verificacion=Verificacion and Theta[1]<4 and Theta[1]>2.8 #beta 
    Verificacion=Verificacion and Theta[2]<1.1 and Theta[2]>0.4 #Alphacuad
    Verificacion=Verificacion and (Theta[3]>50 and Theta[3]<200)  #deltasigma
    Verificacion=Verificacion and Theta[4]<10  #Var
    Verificacion=Verificacion and (Theta[5]>100 and Theta[5]<500)  #ca
    Verificacion=Verificacion and Theta[6]<2  #cb

    return Verificacion



def Energia(Theta):
    Const,beta,alphacuad,deltasigma, Var,ca,cb=Theta
    Error=0    
    Sumlogdet=0
    #Mataux=np.copy(Mat)

    
    
    for k in range(len(ReverseFault)):
        
        
        idx=ReverseFault.index[k]
        l=ReverseFault.loc[idx]["GMW"]
        n=ReverseFault.loc[idx]["R"]
        
        muf=np.log10(np.array(A(frec,l,n,beta,deltasigma,-alphacuad,Const,ca,cb),dtype=float))
        Obs=AShBase[k]
        Diferencia=np.array(Obs-muf)
        Diferencia.shape=(len(Obs),1)


        # Error=Error+np.transpose(Diferencia)@np.linalg.inv(Mat*abs(Amp)*Var)@Diferencia/(2)   
        # Error=Error+np.transpose(Diferencia)@np.diag(1/(abs(Amp)*Var))@Diferencia/(2)   
        Error=Error+np.transpose(Diferencia)@np.diag(1/(np.ones(len(muf))*Var))@Diferencia/(2)   
        # Sumlogdet=Sumlogdet+np.linalg.det(Mat*abs(Amp)*Var)
        # Sumlogdet=Sumlogdet+np.prod(abs(Amp)*Var)
        Sumlogdet=Sumlogdet+np.prod(np.ones(len(muf))*Var)

                
    #plt.plot(Diferencia,alpha=0.1)
    Energia=-(
    -Error#[0,0]
    -len(Obs)/2*Sumlogdet
     +(200-1)*np.log(Const)-4*Const
     +(320-1)*np.log(beta)-100*beta
     +(100-1)*np.log(deltasigma)-1*deltasigma
     +(8-1)*np.log(alphacuad)-10*(alphacuad)
     +(10-1)*np.log(Var)-10*Var     
     +(270-1)*np.log(ca)-1*ca
     +(1-1)*cb-1*cb
     # +(6-1)*cb-10*cb

     #np.sum(np.abs(b))/.1
    )
    #print(Error, Energia,Var)  
    return Energia



Const=45
beta=3.2
deltasigma=100
alphacuad=1/2
Var=0.01
ca=273
cb=0.66

Theta=np.array([Const,beta,alphacuad,deltasigma,Var,ca,cb])


np.random.seed(1)
TwEnergia = pytwalk.pytwalk( n=len(Theta), U=Energia, Supp=EnSupp)
TwEnergia.Run( T=1000000, x0=Theta+0.05, xp0=Theta)
TwEnergia.Ana()
#1000000



# np.savetxt(r'/Users/isaias/Desktop/Corrida2.csv',TwEnergia.Output,delimiter=',', fmt=('%s'))
Corrida=np.genfromtxt(r'/Users/isaias/Desktop/Corrida2.csv', delimiter=',')
Corrida=TwEnergia.Output





# plt.plot(-Corrida[2000:,-1][0::200 ] )

Efectiva=Corrida[0:,][20000::1000 ] 
plt.plot(-Efectiva[:,-1] )
len(Efectiva)

#########
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

# sm.graphics.tsa.plot_acf(Corrida[2000:,5], lags=50)
sm.graphics.tsa.plot_acf(Efectiva[:,5], lags=10)
###############



# plt.plot(-Corrida[:,-1])

# Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]


# Const,beta,alphacuad,deltasigma, Var,ca,cb=Theta

    
# k=0
# for k in range(len(ReverseFault)):
#     plt.plot()
#     idx=ReverseFault.index[k]
#     l=ReverseFault.loc[idx]["GMW"]
#     n=ReverseFault.loc[idx]["R"]
    
#     muf=np.log10(np.array(A(frec,l,n,beta,deltasigma,-alphacuad,Const,ca,cb),dtype=float))
    
    
#     plt.plot(log10(frec),muf)
#     plt.plot(log10(frec),AShBase[k])
#     plt.show()
    
def Prediccion(k, Theta):
    Const,beta,alphacuad,deltasigma, Var,ca,cb=Theta
    idx=ReverseFault.index[k]
    l=ReverseFault.loc[idx]["GMW"]
    n=ReverseFault.loc[idx]["R"]
    muf=np.log10(np.array(A(frec,l,n,beta,deltasigma,-alphacuad,Const,ca,cb),dtype=float))
    return muf    




###############################################


def Ajustes(k):

    Predicciones=np.ones(len(frec))
    for i in range(len(Efectiva)):
        muf=Prediccion(k, Efectiva[i,:-1])
        Predicciones=np.vstack((Predicciones,muf))
    
    return Predicciones[1:,]

# Predicho=Ajustes(0)

def Intervalos(Predicho,a):
    U=np.apply_along_axis(lambda x: np.quantile(x, q=1-a/2), 0, Predicho)
    L=np.apply_along_axis(lambda x: np.quantile(x, q=a/2), 0, Predicho)
    return(U,L)



def Ruido(Predicho,Efectiva,a):
    PredichoN=np.copy(Predicho)
    for i in range(len(Predicho)):
        PredichoN[i,:]=Predicho[i,:]+np.random.normal(size=len(frec), loc=0, scale=np.sqrt(Efectiva[i,-2]))
    U=np.apply_along_axis(lambda x: np.quantile(x, q=1-a/2), 0, PredichoN)
    L=np.apply_along_axis(lambda x: np.quantile(x, q=a/2), 0, PredichoN)
    return(U,L,PredichoN)





graf=0
a=0.2
Predicho=Ajustes(k)
L,U=Intervalos(Predicho, a)
np.random.seed(1)
LR,UR,PredichoR=Ruido(Predicho, Efectiva, a)
plt.plot()
# plt.plot(log10(frec),Prediccion(graf, Efectiva[np.argsort(Efectiva[:,-1])[0],:-1]), label="MAP")
#plt.plot(log10(frec),np.apply_along_axis(np.median, 0, PredichoR), label="Median",color="red")
plt.plot(log10(frec),np.apply_along_axis(np.median, 0, Predicho), label="Median",color="red",linewidth=0.4)
plt.fill_between(log10(frec), L, U, color='orange', alpha=.35)
plt.fill_between(log10(frec), LR, UR, color='b', alpha=.15)
plt.plot(log10(frec),AShBase[graf],label="Observado")
plt.legend()
plt.xlabel(r"$\mathfrak{f}\quad (Hz)$")
plt.ylabel(r"$A(\mathfrak{f})\quad \frac{cm}{s}$")
plt.title(ReverseFault["Base"][ReverseFault.index[graf]]+ " R="+str(np.round(ReverseFault["R"][ReverseFault.index[graf]],2))+" Mw="+str(np.round(ReverseFault["GMW"][ReverseFault.index[graf]],2)))
plt.show()        




################################################



Efectiva.shape

from scipy.stats import gamma


def GraficaPriorPost(indice,SopL,SopU,ParA,ParB):
    plt.hist(Efectiva[:,indice],density=True)
    L=np.min(Efectiva[:,indice])
    U=np.max(Efectiva[:,indice])
    SopGraf = np.linspace(L,U, 100)
    Densidad=gamma(a=ParA,scale=1/ParB)
    plt.plot(SopGraf, Densidad.pdf(SopGraf)/(Densidad.cdf(SopU)-Densidad.cdf(SopL)))

GraficaPriorPost(0,0,100,200,4)   #const
GraficaPriorPost(1,2.8,4,320,100) #beta
GraficaPriorPost(2,0.4,1.1,8,10)  #alphacuad
GraficaPriorPost(3,50,200,100,1)  #deltasigma
GraficaPriorPost(4,0,10,10,10)    #Var
GraficaPriorPost(5,100,500,270,1) #ca
GraficaPriorPost(6,0,2,1,1)       #cb



#Const,beta,alphacuad,deltasigma, Var,ca,cb=Theta




##############################################################################



Densidad=gamma(a=320,scale=1/100)
SopGraf = np.linspace(2.8,4, 100)
plt.plot(SopGraf, Densidad.pdf(SopGraf))



##############################################################################
##############################################################################
##############################################################################
##############################################################################




















##############################################################################
##############################################################################
##############################################################################
##############################################################################



BasesFourier=10
fmin=0.1





def EnSupp(Theta):
    Verificacion=all(0 < Theta[:7]) and Theta[0]<100 #Const
    Verificacion=Verificacion and Theta[1]<4 and Theta[1]>2.8 #beta 
    Verificacion=Verificacion and Theta[2]<1.1 and Theta[2]>0.4 #Alphacuad
    Verificacion=Verificacion and (Theta[3]>50 and Theta[3]<200)  #deltasigma
    Verificacion=Verificacion and Theta[4]<10  #Var
    Verificacion=Verificacion and (Theta[5]>100 and Theta[5]<500)  #ca
    Verificacion=Verificacion and Theta[6]<2  #cb
    Verificacion=Verificacion and all(np.abs(Theta[7:])<200 ) #cb
    return Verificacion



def Energia(Theta):
    Const,beta,alphacuad,deltasigma, Var,ca,cb=Theta[:7]
    ParamBases=Theta[7:]
    
    Error=0    
    Sumlogdet=0
    #Mataux=np.copy(Mat)
    Expansion=0

    for i in range(BasesFourier):
        Expansion+=ParamBases[i]*(2/(fmax-fmin))**(1/2)*np.cos(2*np.pi*(i+1)*(frec-fmax)/(fmax-fmin))

#    plt.plot(frec,Expansion)
    
    
    for k in range(len(ReverseFault)):
        
        
        idx=ReverseFault.index[k]
        l=ReverseFault.loc[idx]["GMW"]
        n=ReverseFault.loc[idx]["R"]
        
        muf=np.log10(np.array(A(frec,l,n,beta,deltasigma,-alphacuad,Const,ca,cb),dtype=float))
        muf+=Expansion
        Obs=AShBase[k]
        Diferencia=np.array(Obs-muf)
        Diferencia.shape=(len(Obs),1)


        # Error=Error+np.transpose(Diferencia)@np.linalg.inv(Mat*abs(Amp)*Var)@Diferencia/(2)   
        # Error=Error+np.transpose(Diferencia)@np.diag(1/(abs(Amp)*Var))@Diferencia/(2)   
        Error=Error+np.transpose(Diferencia)@np.diag(1/(np.ones(len(muf))*Var))@Diferencia/(2)   
        # Sumlogdet=Sumlogdet+np.linalg.det(Mat*abs(Amp)*Var)
        # Sumlogdet=Sumlogdet+np.prod(abs(Amp)*Var)
        Sumlogdet=Sumlogdet+np.prod(np.ones(len(muf))*Var)

                
    #plt.plot(Diferencia,alpha=0.1)
    Energia=-(
    -Error#[0,0]
    -len(Obs)/2*Sumlogdet
     +(200-1)*np.log(Const)-4*Const
     +(320-1)*np.log(beta)-100*beta
     +(100-1)*np.log(deltasigma)-1*deltasigma
     +(8-1)*np.log(alphacuad)-10*(alphacuad)
     +(10-1)*np.log(Var)-10*Var     
     +(270-1)*np.log(ca)-1*ca
     +(1-1)*cb-1*cb
     -np.sum(np.abs(ParamBases)*(2)**(1+np.arange(BasesFourier)) )
     # +(6-1)*cb-10*cb
     #np.sum(np.abs(b))/.1
    )
    #print(Error, Energia,Var)  
    return Energia





Const=45
beta=3.2
deltasigma=100
alphacuad=1/2
Var=0.01
ca=273
cb=0.66

Theta=np.array([Const,beta,alphacuad,deltasigma,Var,ca,cb])

Theta=np.hstack((Theta,np.ones(BasesFourier)))





np.random.seed(1)
TwEnergia = pytwalk.pytwalk( n=len(Theta), U=Energia, Supp=EnSupp)
TwEnergia.Run( T=1000000, x0=Theta+0.05, xp0=Theta)

# np.savetxt(r'/Users/isaias/Desktop/CorridaNoparam.csv',TwEnergia.Output,delimiter=',', fmt=('%s'))
Corrida=np.genfromtxt(r'/Users/isaias/Desktop/CorridaNoparam.csv', delimiter=',')
Efectiva=Corrida[2000:,][0::2000 ] 
plt.plot(-Efectiva[:,-1] )


# sm.graphics.tsa.plot_acf(Corrida[2000:,5], lags=50)
sm.graphics.tsa.plot_acf(Efectiva[:,5], lags=10)
###############

Corrida=TwEnergia.Output

TwEnergia.Ana()


Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]


Const,beta,alphacuad,deltasigma, Var,ca,cb=Theta[:7]
ParamBases=Theta[7:]

    
# k=0
# for k in range(len(ReverseFault)):
#     plt.plot()
#     idx=ReverseFault.index[k]
#     l=ReverseFault.loc[idx]["GMW"]
#     n=ReverseFault.loc[idx]["R"]
    
#     muf=np.log10(np.array(A(frec,l,n,beta,deltasigma,-alphacuad,Const,ca,cb),dtype=float))
        
#     Expansion=0
    
#     for i in range(BasesFourier):
#         Expansion+=ParamBases[i]*(2/(fmax-fmin))**(1/2)*np.cos(2*np.pi*(i+1)*(frec-fmax)/(fmax-fmin))

#     muf+=Expansion    
#     plt.plot(log10(frec),muf)
#     plt.plot(log10(frec),AShBase[k])
#     plt.show()


########################################
def Prediccion(k, Theta):
    Const,beta,alphacuad,deltasigma, Var,ca,cb=Theta[:7]
    ParamBases=Theta[7:]
    #Mataux=np.copy(Mat)
    Expansion=0
    for i in range(BasesFourier):
        Expansion+=ParamBases[i]*(2/(fmax-fmin))**(1/2)*np.cos(2*np.pi*(i+1)*(frec-fmax)/(fmax-fmin))

    idx=ReverseFault.index[k]
    l=ReverseFault.loc[idx]["GMW"]
    n=ReverseFault.loc[idx]["R"]
    muf=np.log10(np.array(A(frec,l,n,beta,deltasigma,-alphacuad,Const,ca,cb),dtype=float))+Expansion
    return muf    




###############################################


def Ajustes(k):

    Predicciones=np.ones(len(frec))
    for i in range(len(Efectiva)):
        muf=Prediccion(k, Efectiva[i,:-1])
        Predicciones=np.vstack((Predicciones,muf))
    
    return Predicciones[1:,]

# Predicho=Ajustes(0)

def Intervalos(Predicho,a):
    U=np.apply_along_axis(lambda x: np.quantile(x, q=1-a/2), 0, Predicho)
    L=np.apply_along_axis(lambda x: np.quantile(x, q=a/2), 0, Predicho)
    return(U,L)



def Ruido(Predicho,Efectiva,a):
    PredichoN=np.copy(Predicho)
    for i in range(len(Predicho)):
        PredichoN[i,:]=Predicho[i,:]+np.random.normal(size=len(frec), loc=0, scale=np.sqrt(Efectiva[i,4]))
    U=np.apply_along_axis(lambda x: np.quantile(x, q=1-a/2), 0, PredichoN)
    L=np.apply_along_axis(lambda x: np.quantile(x, q=a/2), 0, PredichoN)
    return(U,L,PredichoN)





graf=0

for graf in range(len(ReverseFault)):
    a=0.2
    Predicho=Ajustes(graf)
    L,U=Intervalos(Predicho, a)
    np.random.seed(1)
    LR,UR,PredichoR=Ruido(Predicho, Efectiva, a)
    plt.plot()
    # plt.plot(log10(frec),Prediccion(graf, Efectiva[np.argsort(Efectiva[:,-1])[0],:-1]), label="MAP")
    #plt.plot(log10(frec),np.apply_along_axis(np.median, 0, PredichoR), label="Median",color="red")
    plt.plot(log10(frec),np.apply_along_axis(np.median, 0, Predicho), label="Median",color="red",linewidth=0.4)
    plt.fill_between(log10(frec), L, U, color='orange', alpha=.35)
    plt.fill_between(log10(frec), LR, UR, color='b', alpha=.15)
    plt.plot(log10(frec),AShBase[graf],label="Observado")
    plt.title(ReverseFault["Base"][ReverseFault.index[graf]]+ " R="+str(np.round(ReverseFault["R"][ReverseFault.index[graf]],2))+" Mw="+str(np.round(ReverseFault["GMW"][ReverseFault.index[graf]],2)))
    plt.legend()
    plt.show()        



########################################


# plt.plot(-Corrida[2000:,-1][0::200 ] )



from scipy.stats import gamma


def GraficaPriorPost(indice,SopL,SopU,ParA,ParB):
    plt.hist(Efectiva[:,indice],density=True)
    L=np.min(Efectiva[:,indice])
    U=np.max(Efectiva[:,indice])
    SopGraf = np.linspace(L,U, 100)
    Densidad=gamma(a=ParA,scale=1/ParB)
    plt.plot(SopGraf, Densidad.pdf(SopGraf)/(Densidad.cdf(SopU)-Densidad.cdf(SopL)))

GraficaPriorPost(0,0,100,200,4)   #const
GraficaPriorPost(1,2.8,4,320,100) #beta
GraficaPriorPost(2,0.4,1.1,8,10)  #alphacuad
GraficaPriorPost(3,50,200,100,1)  #deltasigma
GraficaPriorPost(4,0,10,10,10)    #Var
GraficaPriorPost(5,100,500,270,1) #ca
GraficaPriorPost(6,0,2,1,1)       #cb





###################################

def Expans(k, Theta):
    # Const,beta,alphacuad,deltasigma, Var,ca,cb=Theta[:7]
    ParamBases=Theta[7:]
    #Mataux=np.copy(Mat)
    Expansion=0
    for i in range(BasesFourier):
        Expansion+=ParamBases[i]*(2/(fmax-fmin))**(1/2)*np.cos(2*np.pi*(i+1)*(frec-fmax)/(fmax-fmin))

    # idx=ReverseFault.index[k]
    # l=ReverseFault.loc[idx]["GMW"]
    # n=ReverseFault.loc[idx]["R"]
    # muf=np.log10(np.array(A(frec,l,n,beta,deltasigma,-alphacuad,Const,ca,cb),dtype=float))+Expansion
    return Expansion    




Expansion=np.ones(len(frec))
for i in range(len(Efectiva)):
    muf=Expans(k, Efectiva[i,:-1])
    Expansion=np.vstack((Expansion,muf))    


U=np.apply_along_axis(lambda x: np.quantile(x, q=1-a/2), 0, Expansion)
L=np.apply_along_axis(lambda x: np.quantile(x, q=a/2), 0, Expansion)

plt.plot(log10(frec),np.apply_along_axis(np.median, 0, Expansion))
plt.fill_between(log10(frec), L, U, color='orange', alpha=.35)

for graf in range(len(ReverseFault)):
    a=0.2
    Predicho=Ajustes(graf)
    plt.plot(log10(frec),np.apply_along_axis(np.median, 0, Predicho), label="Median",color="red",linewidth=0.4)
    # plt.fill_between(log10(frec), L, U, color='orange', alpha=.35)
    plt.legend()
    plt.show()        




##################################

































