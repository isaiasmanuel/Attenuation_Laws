#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 19:53:54 2021

@author: isaias
"""


from scipy.interpolate import interp1d

from scipy.stats import gamma
import plotly.graph_objects as go
from numpy import genfromtxt

from csaps import csaps

######## Correr antes Funciones y SeleccionFinal


frec=np.log10(10**(np.linspace(np.log10(0.1),np.log10(10),15)))


# hist=np.ones(len(Muestra2))
# for i in range(len(Muestra2)):
#     hist[i]=np.max(MCMCBase[2*i+1])
               
# plt.hist(hist)    
# plt.ylabel("Número de registros")    
# plt.xlabel("Frecuencia máxima")    

# Muestra["Base"].loc[Muestra["Base"]=="INMD8509.191"]
# Muestra.loc[143]
# np.sum(Muestra.index<143)
i=57


fasex=np.arange(-np.pi,np.pi,0.1)
for i in ReverseFault.index:
    plt.figure()
    Base=ReverseFault["Base"][i]
    #Base="INMD8509.191"
    R=ReverseFault["R"][ReverseFault["Base"]==Base].values[0]    
    Mag=ReverseFault["GMW"][ReverseFault["Base"]==Base].values[0]    
    j=ReverseFault["EjeH1"][ReverseFault["Base"]==Base].values[0]
    #Muestra["EjeH1"][Muestra["Base"]=="INMD8509.191"]    
    (X,NS,f,Amp, VelMuest)=fourier(Base,int(j))


    len(X)
    
    np.random.seed(2)    

    # NS=np.hstack((0,np.cumsum(np.random.normal(size=len(X)))))
    # NS=NS-np.arange(len(NS))/(len(NS))*NS 
    # NS=np.sin(np.arange(len(X))/len(X)*np.pi)+np.random.normal(size=len(X))/100
    # plt.plot(np.arange(len(NS))/(len(NS)),NS)    
    # # plt.plot(np.arange(len(X))/len(X)*np.pi,NS)
    # plt.title("PuenteBrowniano")

    Y=sp.fftpack.fft(NS)
    len(Y)
    N=Y.size
    Y=Y[:N//2]
    Y=Y[1:]
    len(Y)
    
    #plt.plot(np.arctan(np.real(Y)[:100]/np.imag(Y)[:100]))
    
    Fases=np.angle(Y)
    # #plt.plot(Fases[:100],"-o")
    # plt.plot(Fases,"o",markersize=0.8)
    # #plt.title(Base+" Distancia "+ str(np.round(R,2))+ " Magnitud " +str(np.round(Mag,2)))    
    # plt.xlabel("$n$")
    # plt.ylabel("$\phi_{n}$")

    # plt.hist(Fases)
    # plt.xlabel("$\phi_n$")
    # plt.ylabel("Frec")
    # plt.title(Base+" Distancia "+ str(np.round(R,2))+ " Magnitud " +str(np.round(Mag,2)))    
    
    #plt.plot(Fases[:-1],Fases[1:])
    plt.scatter(Fases[:-1],Fases[1:],s=5)
    plt.plot(fasex[fasex<0],(fasex+np.pi)[fasex<0],color="orange")
    plt.plot(fasex[fasex>0],(fasex-np.pi)[fasex>0],color="orange")
    plt.xlabel("$\phi_n$")
    plt.ylabel("$\phi_{n+1}$")
    plt.title(Base+", R="+ str(np.round(R,2))+ ", $M_w$=" +str(np.round(Mag,2)))    
    plt.show()
    #plt.plot(Fases[:-1],Fases[1:])
    
    
plt.plot(f[1:][f[1:]<15],np.angle(Y)[f[1:]<15], linewidth=0.5)
plt.xlabel("$\mathfrak{f}\quad (Hz)$")
plt.ylabel("$\phi_{\mathfrak{f}}$")
plt.title(Base+", R="+ str(np.round(R,2))+ ", $M_w$=" +str(np.round(Mag,2)))    

############################    
plt.hist(np.angle(Y), linewidth=0.5)
# plt.xlabel("$\mathfrak{f}\quad (Hz)$")
plt.xlabel("$\phi_{\mathfrak{f}}$")
plt.title(Base+", R="+ str(np.round(R,2))+ ", $M_w$=" +str(np.round(Mag,2)))    
############################
Fases=np.angle(Y)

np.random.seed(1)
Fases2=np.random.uniform(low=-np.pi, high=np.pi,size=len(Y))

plt.plot(X[1::2][1:],np.real(ifft(np.abs(Y)*np.exp(Fases*1j))))
plt.plot(X,NS)

plt.plot(X[1::2][1:],np.real(ifft(np.abs(Y)*np.exp(Fases2*1j))))
plt.xlabel(r"$s$")
plt.ylabel(r"$\frac{cm}{s^2}$")


################
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "browser"
df = px.data.iris()
fig = px.line_3d( z=np.arange(len(Fases)-1)[:100], y=Fases[:-1][:100], x=Fases[1:][:100])
fig.show()

################ Simulacion de seno con ruido, puente browniano y excursion

# np.random.seed(1)    
# NS=np.sin(np.arange(len(X))/len(X)*np.pi)+np.random.normal(size=len(X))/100
# NS[0]=0
# NS[-1]=0
# plt.plot(np.arange(len(NS))/(len(NS)),NS)    


# Soporte=np.arange(len(NS))/(len(NS))  
# # np.random.seed(1)    

# NS=np.hstack((0,np.cumsum(np.random.normal(size=len(X)))))
# NS=NS-np.arange(len(NS))/(len(NS))*NS 
# # plt.plot(np.arange(len(NS))/(len(NS)),NS)    
# ###############Excursion

# minimo= np.min(NS)

# LocMin=np.where(NS==minimo)[0][0]

# NT=((Soporte[LocMin]+Soporte)%1 * len(NS))//1
# NT=NT.astype(int)
# NS=NS[NT]-minimo

# # plt.plot(np.arange(len(NT))/len(NT), NS)


##############


Y=sp.fftpack.fft(NS)
len(Y)
N=Y.size
Y=Y[:N//2]
Y=Y[1:]
len(Y)

#plt.plot(np.arctan(np.real(Y)[:100]/np.imag(Y)[:100]))

Fases=np.angle(Y)
# #plt.plot(Fases[:100],"-o")
# plt.plot(Fases,"o",markersize=0.8)
# #plt.title(Base+" Distancia "+ str(np.round(R,2))+ " Magnitud " +str(np.round(Mag,2)))    
# plt.xlabel("$n$")
# plt.ylabel("$\phi_{n}$")

# plt.hist(Fases)
# plt.xlabel("$\phi_n$")
# plt.ylabel("Frec")
# plt.title(Base+" Distancia "+ str(np.round(R,2))+ " Magnitud " +str(np.round(Mag,2)))    
plt.scatter(Fases[:-1][:1000],Fases[1:][:1000])
plt.xlabel("$\phi_n$")
plt.ylabel("$\phi_{n+1}$")
plt.title(Base+" Distancia "+ str(np.round(R,2))+ " Magnitud " +str(np.round(Mag,2)))    
# plt.savefig(Base[:4]+".png",dpi=300)


# plt.hist(Fases)

















