#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:32:29 2023

@author: isaias
"""
#######################Anadir seleccion Final

###############################################
import os
from scipy.fft import irfft,rfft,fft, ifft,fftfreq
import seaborn as sns
# Get the current working directory
os.chdir('/Users/isaias/Desktop/Atenuacion')

####Restroingi a 18.94,13.06, -105.1,-93
from geopy import distance
import os
import pytwalk

from io import StringIO 
from numpy import loadtxt 
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.fftpack
import re
import pandas as pd
from scipy import signal
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
from plotly.graph_objs import *
import plotly.express as px
from datetime import date
from scipy import integrate

from pandas import read_csv

from scipy.stats import gamma, lognorm,laplace_asymmetric
# a=1;b=2;s2=1

import multiprocessing 

from multiprocessing import Pool
import multiprocessing as mp


from scipy.stats import gamma, lognorm,laplace_asymmetric

    
from scipy.ndimage import gaussian_filter

from multiprocessing import Pool
from multiprocessing import get_context
import multiprocessing 

############################

def RespuestaSimulada(PSDF):
    Angle=np.random.uniform(0,2*np.pi,size=PSDF.shape[0])
    Proceso=np.ones(PSDF.shape[1])
    
    PSDF.shape
    t2.shape
    f2.shape
    
    
    # for i in range(PSDF.shape[1]):
    #     Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(2*PSDF.clip(0)[:,i]*np.pi*f2[1])*np.cos(2*np.pi*np.arange(len(f2))*f2[1]*t2[i]+Angle))
    
    # for i in range(PSDF.shape[1]):
    #     #Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(2*PSDF.clip(0)[:,i]*f2[1])*np.cos(np.arange(len(f2))*f[1]*t2[i]+Angle))
    #     Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(PSDF.clip(0)[:,i]*f2[1])*np.cos(2*np.pi*np.arange(len(f2))*f2[1]*t2[i]+Angle))
    # from datetime import datetime
    # now = datetime.now()   
    # for i in range(PSDF.shape[1]):
    #     #Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(2*PSDF.clip(0)[:,i]*f2[1])*np.cos(np.arange(len(f2))*f[1]*t2[i]+Angle))
    #     Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(PSDF.clip(0)[:,i]*f2[1])*np.cos(2*np.pi*np.arange(len(f2))*f2[1]*t2[i]+Angle))
    # print(datetime.now()-now)
    
    # now = datetime.now()    
    R=(np.reshape(2*np.pi*f2[1]*np.arange(len(f2)),(len(f2),1))@np.reshape(t2,(1,len(t2))))
    R+=np.reshape(Angle,(Angle.shape[0],1))
    R=np.cos(R)
    R=np.einsum('ij,ji->i', np.sqrt(2)*(np.sqrt(PSDF.clip(0)[:,:]*f2[1])).T,R)
    # plt.plot(Proceso-R)
    # plt.plot(R)    
    # print(datetime.now()-now)

    Proceso=np.copy(R)
    
    # plt.plot(X,NS,alpha=0.6)
    # plt.plot(t2,Proceso,alpha=0.6)
    
    # print(integrate.cumtrapz(NS**2,X , initial=0)[-1],integrate.cumtrapz(Proceso**2,X , initial=0)[-1])
    
    # (fsim,tsim,Zxxsim)=sp.signal.stft(Proceso,fs=nperseg,nperseg=batch,window=Window,noverlap=over)#los indices estan volteados respecto al articulo
    
    
    
    # trace1 = go.Surface(x=t,y=f,z=Z, colorscale='Viridis')
    # trace1 = go.Surface(x=tsim,y=fsim,z=np.abs(Zxxsim)**2, colorscale='Viridis')
        
    # iplot([trace1])
    
    
    
    
    # Sintetico=Proceso
    # Registro=np.vstack((t2,Sintetico))
    # Registro=np.vstack((Registro,Sintetico))
    # Registro=np.vstack((Registro,Sintetico))
    # Registro.shape
    # ############### Sismos registrados lejanos
    
    # # Muestra.sort_values(by=["R"])[-10:]
    # # STFT(59,10,90) #Canditato
    # # STFT(52,10,90) #Canditato
    # # STFT(116,10,90) #Canditato
    # # STFT(168,10,90) #Canditato
    # # Muestra.sort_values(by=["R"])[-30:-20]
    
    
    # RegistroChristen=Registro.T
    # np.savetxt(r'./RegistroChristen2'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv',RegistroChristen,delimiter=',', fmt=('%s'))
    # # np.savetxt(r'./RegistroChristen2.csv',RegistroChristen,delimiter=',', fmt=('%s'))
    # tst = Acelerogram( r'./RegistroChristen2'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv', Dt=-1, skiprows=1)
    
    return Proceso













def SePSDF(i,Canal,b):
    Muestra=ReverseFault
    idsismo=np.copy(i)
    if Canal=="V":
        CanalV=np.int(np.float(RecuperaOrientacion(Muestra["Base"][i])[7]))
        Canal=CanalV
    if Canal=="H1":
        CanalH1=np.int(np.float(RecuperaOrientacion(Muestra["Base"][i])[8]))
        Canal=CanalH1
    if Canal=="H2":
        CanalH2=np.int(np.float(RecuperaOrientacion(Muestra["Base"][i])[9]))
        Canal=CanalH2
    (X,NS,f,Y,VelMuest)=fourier(Muestra["Base"][i], Canal) #Usar registrocrudo cuando no quiero normalizado y centrado
    
    
    LenMCMC=10000
    Hilos=46 ###Elegir 12 si  es local y 48 en el servidor
    # fmax=15
    Trayectorias=20
    tmax=max(X)
    tmin=00
    nu=30        #Grados de libertad
    cuantil=1 #La muestra se corta en este cuantil de la marginal en tiempo
    Reduccion=int(1*VelMuest )#Numero de muestras que se juntan, no hay interseccion
    
    # Dominio=X>Muestra["Inicio"][idsismo]
    # X=X[Dominio]
    # NS=NS[Dominio]
    
    # plt.plot(X,NS,linewidth=.5)
    # print(VelMuest)
    
    ###############Creamos registro
    
    
    
    # np.savetxt(r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv',RegistroChristen,delimiter=',', fmt=('%s'))
    
    
    ####### Calculamos STFT
    NS=signal.detrend(NS[X>ReverseFault["Inicio"][i]])
    NS=signal.tukey(len(NS),alpha=0.05)*NS
    Y=sp.fftpack.fft(NS)
    N=Y.size
    Y=np.abs(Y[:N//2])
    #X=np.linspace(0,Duracion, int((Duracion-1)*VelMuest))
    #A veces marca error
    X=np.arange(len(NS))/VelMuest
    T=X[1]-X[0]
    
    f = np.linspace(0, 1 / T, N)    
    f = f[:N//2]    
    
    nperseg=muestseg=VelMuest
    batch=np.sum(X<b) #50 es el bueno
    ####### Calculamos STFT
    # plt.plot()
    # plt.plot(X,NS)
    # plt.show()
    # NS=NS[X>tmin]
    # X=X[X>tmin]
    
    # NS=NS[X<tmax]
    # X=X[X>tmax]
    
    
    
    #Window=np.sqrt(scipy.signal.windows.gaussian(M=batch,std=0.25*batch))
    Window=np.sqrt(scipy.signal.windows.gaussian(M=batch,std=1*batch))
    Const=integrate.cumtrapz(Window**2,X[:batch])[-1]
    Window=Window/np.sqrt(Const)
    #Window=Window/10000
    
    over=batch*.9
    # plt.subplots()
    (f,t,Zxx)=sp.signal.stft(NS,fs=nperseg,nperseg=batch,window=Window,noverlap=over)#los indices estan volteados respecto al articulo
    
    Z=np.abs(Zxx)**2
    return (f,t,Z)
    


def Graficas(idsismo,b,log):
    
    
    
    f,t,Z=SePSDF(idsismo,"H1",b)                
    
    Marginalf=np.ones(len(f))
    
    for i in np.arange(len(f)):
        Marginalf[i]=integrate.cumtrapz(Z[i,:],t , initial=0)[-1]
        
    Const=integrate.cumtrapz(Marginalf,f,initial=0)[-1]    
        
    Marginalf=Marginalf/integrate.cumtrapz(Marginalf,f,initial=0)[-1]
    
    Marginalt=np.ones(len(t))
    
    for i in np.arange(len(t)):
        Marginalt[i]=integrate.cumtrapz(Z[:,i],f , initial=0)[-1]
        
    
    Marginalt=Marginalt/integrate.cumtrapz(Marginalt,t,initial=0)[-1]
    
    
    Znorm=Z/Const
    
    CopulaDensidad=np.zeros(Znorm.shape)
    
    for i in range(Znorm.shape[0]):
        for j in range(Znorm.shape[1]):
            CopulaDensidad[i,j]=Znorm[i,j]/(Marginalt[j]*Marginalf[i])
    
    
    DistMargf=integrate.cumtrapz(Marginalf,f,initial=0)
    DistMargt=integrate.cumtrapz(Marginalt,t,initial=0)
    
    
    # plt.pcolormesh(t[t<tmax], f[f<fmax], Z[f<fmax,:][:,t<tmax], shading='gouraud', rasterized=True)
    
    plt.plot()
    plt.pcolormesh(DistMargt[t<tmax], DistMargf[f<fmax], (CopulaDensidad[f<fmax,:][:,t<tmax]), shading='gouraud', rasterized=True)
    plt.colorbar()
    plt.show()


######################################

    
      
    ConjuntaProducto=np.zeros(Znorm.shape)
    
    for i in range(Znorm.shape[0]):
        for j in range(Znorm.shape[1]):
            ConjuntaProducto[i,j]=(Marginalt[j]*Marginalf[i])
    
    
    if log=="log":
        plt.plot()
        plt.title("Independiente")
        plt.pcolormesh(t[t<tmax], f[f<fmax], np.log10(ConjuntaProducto[f<fmax,:][:,t<tmax]), shading='gouraud', rasterized=True)
        plt.colorbar()
        plt.show()
        
        
        plt.plot()
        plt.title("ePSDF")
        plt.pcolormesh(t[t<tmax], f[f<fmax], np.log10(Znorm[f<fmax,:][:,t<tmax]), shading='gouraud', rasterized=True)
        plt.colorbar()
        plt.show()
    else:
        plt.plot()
        plt.title("Independiente")
        plt.pcolormesh(t[t<tmax], f[f<fmax], ConjuntaProducto[f<fmax,:][:,t<tmax], shading='gouraud', rasterized=True)
        plt.colorbar()
        plt.show()
        
        
        plt.plot()
        plt.title("ePSDF")
        plt.pcolormesh(t[t<tmax], f[f<fmax], Znorm[f<fmax,:][:,t<tmax], shading='gouraud', rasterized=True)
        plt.colorbar()
        plt.show()

        
    Const=integrate.cumtrapz(NS[X>ReverseFault["Inicio"][idsismo]]**2,X[X>ReverseFault["Inicio"][idsismo]])[-1]
    return ConjuntaProducto,Const,f,t


#############################


def Buscar(Busqueda,Base):
        P=open(Extension+"/"+Base, errors="ignore")
        Datos=P.read()
        
        # if Datos.find("\x1a")!=-1:
        #     Datos=Datos.replace("\x1a","")

        CE=Datos.find(Busqueda)
        if CE!=-1:
            Inicio=Datos[CE:].find(":")            
            if Busqueda.find("COORDENADAS")!=0:
                FINAL=Datos[CE:].find("\n")
                #print(Datos[(CE+Inicio+2):(CE+FINAL)])
                Cadena=Datos[(CE+Inicio+2):(CE+FINAL)]
            else:
                PREFINAL=Datos[CE+Inicio:].find("\n")
                FINAL=Datos[(CE+Inicio+PREFINAL+1):].find("\n")
                #    CoordenadasEstacion=np.hstack((Sismo2,Sismo[i]+" "+Datos[(CE+Inicio+2):(CE+FINAL)]))                
                Precadena=Datos[(CE+Inicio+2):(CE+Inicio+PREFINAL+FINAL+1)]                                
                #Cadena=Precadena[:Precadena.find("\n")]+" "+Precadena[Precadena.find(": ")+2:]
                Cadena=Precadena[:Precadena.find("\n")]+" "+Precadena[Precadena.find(":")+1:]
        else:
            Cadena="-"
        return Cadena
    
def fourier(Base,canal):
    #'ZIG39609.221' ZIG39606.051 No se puede arreglar
    Encontrado=Buscar("VEL. DE MUESTREO", Base)    
    Encontrado=Encontrado.replace(" ","")

    if len(Encontrado)>2 and Encontrado.find("--")==-1:
        Pos=Posiciones(Encontrado, "/")
        if canal<2:
            VelMuest=float(Encontrado[(Pos[canal]+1):Pos[canal+1]])
        else:
            VelMuest=float(Encontrado[Pos[canal]+1:])
    if len(Encontrado)<2  or Encontrado.find("--")!=-1:
        Encontrado=Buscar("INTERVALO DE MUESTREO", Base) 
        Pos=Posiciones(Encontrado, "/")        

        if canal<2:
            VelMuest=1/float(Encontrado[(Pos[canal]+1):Pos[canal+1]])
        else:
            VelMuest=1/float(Encontrado[Pos[canal]+1:])


    ###############
    
    #ZIG39609.221 ZIG39606.051 ZIG39606.032 ZIG39610.033
    #ZIG39408.281 ZIG39409.221 ZIG69603.091 ZIG69606.031
    #ZIG39409.031 ZIG69608.081
    
    #IN188510.29A Deja de medir

    Datos=rengloninicio(Base)
    Datos=Datos.replace("0-", "0 -").replace("1-", "1 -").replace("2-", "2 -")
    Datos=Datos.replace("3-", "3 -").replace("4-", "4 -").replace("5-", "5 -").replace("6-", "6 -")
    Datos=Datos.replace("7-", "7 -").replace("8-", "8 -").replace("9-", "9 -").replace("..", ".")
    


    if Datos.find("\x1a")!=-1:
        c = StringIO(Datos.replace("\x1a",""))
        Mediciones=np.loadtxt(c)
    
        # b=np.genfromtxt(c,delimiter='\n',dtype=str)
        # b[3501]

    else:
        c = StringIO(Datos)
        Mediciones=np.loadtxt(c)
        #Mediciones=loadtxt(Extension+"/"+Base, skiprows=rengloninicio(Base))

    NS=Mediciones[:,canal]
    #Correcion de media y tendencia
    NS=signal.detrend(NS)
    NS=signal.tukey(len(NS),alpha=0.05)*NS
    Y=sp.fftpack.fft(NS)
    N=Y.size
    Y=np.abs(Y[:N//2])
    #X=np.linspace(0,Duracion, int((Duracion-1)*VelMuest))
    # A veces marca error
    #X=np.linspace(0,len(NS)/VelMuest, len(NS))    
    #X=np.linspace(0,len(NS)/VelMuest, len(NS)+1)    
    X=np.arange(len(NS))/VelMuest
    T=X[1]-X[0]

    f = np.linspace(0, 1 / T, N)    
    f = f[:N//2]
    # plt.subplot(2, 1, 1)
    # plt.plot(X,NS)
    # plt.subplot(2, 1, 2)
    # plt.plot(f[:N//2],Y[:N//2])
    # plt.xlim(0.1,10)
    return (X,NS,f,Y,VelMuest)



def RecuperaOrientacion(Base):
    Orient=Buscar("ORIENTACION",Base)
    Orient=Orient.replace(" ","")
    #print(Orient)
    aux=np.array(["H1","H2",0.,0.,0.,"90","00",0.,0.,0.])
    eje=0
    horizontal=0    
    #print(np.sort(Separacadenas(Tamanio[i][0], Posiciones(Tamanio[i][0],"/"))))
    for j in Separacadenas(Orient, Posiciones(Orient,"/")):        
        if j.find("V")!=-1:
            if j.find("-")!=-1:
                print("HAY UN V NEGATIVO")
                aux[2]=-1
            else:
                aux[2]=1
            aux[7]=eje
    #    if horizontal<2:       
        
        if j.find("V")==-1 and horizontal<2:
            if j.find("-")!=-1:
                #print("HAY UN SIGNO NEGATIVO")
                if j.find("--")!=-1:
                    aux[3+horizontal]=1
                else :
                    aux[3+horizontal]=-1
                                

            else:
                aux[3+horizontal]=1
            
            aux[8+horizontal]=eje
            
            if j.find(";")==-1:
                aux[0+horizontal]=j[0]+j[-1]
                aux[5+horizontal]=j[1]+j[-2]
            else:
                aux[0+horizontal]=j[0]+j[:j.find(";")][-1]
                aux[5+horizontal]=j[1]+j[:j.find(";")][-2]
            horizontal+=1
        eje+=1
    if aux[7]=="0.0":     
        aux=np.array(["Error",Bases[i],0.0,0.,0.,"90","00",0.,0.,0.])
    return aux
    



def rengloninicio(Base):
    P=open(Extension+"/"+Base, errors="ignore")
    Datos=P.read()
    
    Datos
    
    Intento=Datos.find("---------+")
    
    Intento2=Datos[Intento:].find("\n")
    
    
    Intento3=Datos[Intento+Intento2:].find("---------+")
    
    Intento4=Datos[Intento+Intento2+Intento3:].find("\n")
    return Datos[Intento+Intento2+Intento3+Intento4+1:]
    


def Separacadenas(Cadena, lugares):

    return (Cadena[(lugares[0]+1):lugares[1]],Cadena[(lugares[1]+1):lugares[2]],Cadena[(lugares[2]+1):])

def Posiciones(Cadena,Caracter):
    Posicion=0
    Conteo=0
    for i in Cadena: 
        if i==Caracter:
            Posicion=np.hstack((Posicion,Conteo))
        #print(i,Conteo)
        Conteo+=1
    Posicion=Posicion[1:]
    return Posicion





###############################################
import os
from scipy.fft import irfft,rfft,fft, ifft,fftfreq
import seaborn as sns
# Get the current working directory
os.chdir('/Users/isaias/Desktop/Atenuacion')

####Restroingi a 18.94,13.06, -105.1,-93
from geopy import distance
import os
from io import StringIO 
from numpy import loadtxt 
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.fftpack
import re
import pandas as pd
from scipy import signal
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
from plotly.graph_objs import *
import plotly.express as px
from datetime import date
from scipy import integrate

import plotly.graph_objects as go
import plotly.express as px
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot

from pandas import read_csv

######
Extension='./RegistrosCrudos/'
Bases=os.listdir(Extension)
try :
    Bases.remove(".DS_Store")
    print("Se completo")
except:
    print("Se completo")
#Base=Bases[0]
###############################################



###############################################


GlobalCMTBase=pd.read_csv('/Users/isaias/Desktop/Atenuacion/BaseGCMT.csv')
GlobalCMTBase=GlobalCMTBase.rename(columns={"NOMBRE": "IDSismo"}, errors="raise")

FinalGCMT=pd.read_csv('/Users/isaias/Desktop/Atenuacion/FinalGCMT.csv')



FinalGCMT[FinalGCMT["GMW"]>0]
Muestra=pd.merge(FinalGCMT, GlobalCMTBase, on=['IDSismo'], how="inner", indicator=True)
Muestra=Muestra.drop(columns=Muestra.columns[-1])
Muestra["LonEstacion"]=-Muestra["LonEstacion"]
Muestra.columns



from geopy import distance

    
R=np.ones(len(Muestra))
k=0
for i in Muestra.index:
    R[k]=distance.distance((Muestra["LatEstacion"][i],Muestra["LonEstacion"][i]), (Muestra["LatSismo"][i],Muestra["LonSismo"][i])).km
    k+=1

Muestra["R"]=R

##############################################



indices=np.zeros(0)
for i in Muestra.index:
    if Muestra["Angle1"][i]<0:
        ang1=Muestra["Angle1"][i]+360
    else :
        ang1=Muestra["Angle1"][i]

    if Muestra["Angle2"][i]<0:
        ang2=Muestra["Angle2"][i]+360
    else :
        ang2=Muestra["Angle2"][i]
    if ang1<135 and ang1>45 and ang2<135 and ang2>45 :
        indices=np.append(indices,i)



ReverseFault=Muestra.loc[indices]
#############################

from obspy.taup import TauPyModel
from datetime import datetime

Sismos= pd.read_csv('Elegidos.csv', delimiter=',')
# ReverseFault

ReverseFault=pd.merge(Sismos, ReverseFault, on=['Base'], how="inner", indicator=True)

###########################


ReverseFault=ReverseFault[~ReverseFault["IDEstacion"].isin(["SCRU","MIHL","NILT","VHSA","CHPA","COMA","COLL","MANZ","TOTO","PZPU","SXPU","PBPZ","PHPU","LMPP","SODO","OZST","THEZ","OXJM","HUAM","XALA","PBP2"])]
ReverseFault=ReverseFault[ReverseFault["Profundidad_y"]>=20]
ReverseFault=ReverseFault[ReverseFault["LONGITUD"]<-95] ##################### Hay que quitar placa del caribe y placa rivera, ReverseFault[ReverseFault["LONGITUD"]<-95]
ReverseFault=ReverseFault[ReverseFault["LatEstacion"]<19.331]

ReverseFault=ReverseFault[ReverseFault["R"]>100] ##################### Hay que quitar placa del caribe y placa rivera, ReverseFault[ReverseFault["LONGITUD"]<-95]



ReverseFault=ReverseFault.groupby("IDSismo").filter(lambda x: len(x) >= 3)
# ReverseFault=ReverseFault.groupby("IDEstacion").filter(lambda x: len(x) > 1)
######################################


ReverseFault["ET"]=0

for idsismo in ReverseFault.index:
    CanalH1=np.int(np.float(RecuperaOrientacion(ReverseFault["Base"][idsismo])[8]))
    (X,NS,f,Y,VelMuest)=fourier(ReverseFault["Base"][idsismo], CanalH1)    
    Accel=NS[X>ReverseFault["Inicio"][idsismo]]
    Sup=X[X>ReverseFault["Inicio"][idsismo]]
    ReverseFault["ET"].loc[idsismo]=integrate.cumtrapz(Accel**2,Sup)[-1]

ReverseFault.columns


ReverseFault=ReverseFault[ReverseFault["IDEstacion"]!="SJLL"]
#######################################



def AjustaMarginalFrec(i):
    theta=np.array((2,0.01,0.01))
    print(i)
    np.random.seed(1)
    # theta=np.array((2,2,2,3))
        
    def ExpSupp(theta):
        Valor= all(theta>0)
        Valor= Valor and (theta[0]>1)
        Valor= Valor and (theta[1]<100)
        Valor= Valor and (theta[2]<100)
        return Valor
    
    
    def Veros(theta):
        a,b,s2=theta
        # a,b,c,s2=theta
        y=gamma.pdf(Directorio[i][0][Directorio[i][0]<fmax], a,scale=1/b)
        const=gamma.cdf(fmax,a,scale=1/b)
        # y=lognorm.pdf(Directorio[i][0][Directorio[i][0]<fmax], s=a,loc=b)
        # y=laplace_asymmetric.pdf(Directorio[i][0][Directorio[i][0]<fmax],loc=a,scale=b,kappa=c)
        # return np.sum((Directorio[i][1][Directorio[i][0]<fmax]-y)**2)/(2*s2)+np.sum([Directorio[i][0]<fmax])/2*np.log(s2)
        return np.sum((Directorio[i][1][Directorio[i][0]<fmax]-y/const)**2)/(2*s2)+np.sum([Directorio[i][0]<fmax])/2*np.log(s2)
    

    MCMC= pytwalk.pytwalk( n=len(theta), U=Veros, Supp=ExpSupp)
    MCMC.Run( T=100000, x0=theta, xp0=theta+0.1)
    
    Corrida=MCMC.Output
    
    # plt.plot(-Corrida[:,-1])
    # plt.plot(-Corrida[10000:,-1])
    Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]
    a,b,s2=Theta
    # a,b,c,s2=Theta

    # y=laplace_asymmetric.pdf(Directorio[i][0][Directorio[i][0]<fmax],loc=a,scale=b,kappa=c)
    # # plt.plot(Directorio[i][0][Directorio[i][0]<fmax],y)
    
    # plt.figure()
    # print(a,b,c)
    # Efectiva=(Corrida[10000:,:3])[0::500]
    # for Efe in range(len(Efectiva)):
    #     a,b,c=Efectiva[Efe,:]
    #     y=laplace_asymmetric.pdf(Directorio[Efe][0][Directorio[Efe][0]<fmax],loc=a,scale=b,kappa=c)
    #     plt.plot(Directorio[Efe][0][Directorio[Efe][0]<fmax],y,alpha=0.1,color="red")
    # plt.plot(Directorio[i][0][Directorio[i][0]<fmax],Directorio[i][1][Directorio[i][0]<fmax],alpha=0.5)
    # plt.title((ReverseFault["Base"][i])+" R="+str(np.round(ReverseFault["R"][i],2)))
    # plt.show()
    np.savetxt(r'./Corrida/CorridaFrecuencia'+ (ReverseFault["Base"][ReverseFault.index[i]]).replace(".","")  +".csv",Corrida,delimiter=',', fmt=('%s'))    
    
        





def AjustaMarginalTiempo(i):
    theta=np.array((2,2,3))
    print(i)
    np.random.seed(1)
    # theta=np.array((2,2,2,3))
        
    def ExpSupp(theta):
        Momentos=gamma.stats(a=theta[0],scale=1/theta[1],moments='mv')
        Valor= all(theta>0)
        Valor= Valor and (theta[0]<50)
        Valor= Valor and (theta[1]<50)
        Valor= Valor and (theta[0]>1)
        Valor= Valor and (Momentos[1]<150) and (Momentos[1]>2)
        # Valor= Valor and (theta[1]<100)
        # Valor= Valor and (theta[2]<10)
        return Valor
    
    
    def Veros(theta):
        a,b,s2=theta
        # a,b,c,s2=theta
        y=gamma.pdf(Directorio2[i][0][Directorio2[i][0]<tmax], a,scale=1/b)
        # y=lognorm.pdf(Directorio[i][0][Directorio[i][0]<tmax], s=a,loc=b)
        # y=laplace_asymmetric.pdf(Directorio[i][0][Directorio[i][0]<fmax],loc=a,scale=b,kappa=c)
        return np.sum((Directorio2[i][1][Directorio2[i][0]<tmax]-y)**2)/(2*s2)+np.sum([Directorio2[i][0]<tmax])/2*np.log(s2)
    
    theta=np.array((2,.2,3))
    MCMC= pytwalk.pytwalk( n=len(theta), U=Veros, Supp=ExpSupp)    
    MCMC.Run( T=100000, x0=theta, xp0=theta+0.1)
    
    Corrida=MCMC.Output
    
    # plt.figure()
    # Efectiva=(Corrida[20000:,:3])[0::500]
    # plt.plot(Directorio2[i][0][Directorio2[i][0]<tmax],Directorio2[i][1][Directorio2[i][0]<tmax],alpha=0.5)
    # for Efe in range(len(Efectiva)):
    #     a,b,c=Efectiva[Efe,:]
    #     y=gamma.pdf(Directorio2[i][0][Directorio2[i][0]<tmax], a,scale=1/b)
    #     plt.plot(Directorio2[i][0][Directorio2[i][0]<tmax],y,alpha=0.1,color="red")
    # plt.plot(Directorio2[i][0][Directorio2[i][0]<tmax],Directorio2[i][1][Directorio2[i][0]<tmax],alpha=0.5)
    # plt.title((ReverseFault["Base"][i])+" R="+str(np.round(ReverseFault["R"][i],2)))
    # plt.show()
    
    np.savetxt(r'./Corrida/CorridaTiempo'+ (ReverseFault["Base"][ReverseFault.index[i]]).replace(".","")  +".csv",Corrida,delimiter=',', fmt=('%s'))    
    









#######################################



import os

from numpy import loadtxt, floor, zeros, pi, exp, log, sin, cos, arange, ones, mean
from numpy import array, argmax, angle, linspace, log10, append, where, ndarray, cumsum
from numpy import abs as npabs
from numpy import round as npround
from scipy.integrate import solve_ivp
from scipy.signal import detrend, tukey
from scipy.fftpack import fft, ifft, irfft
from scipy.stats import uniform, norm
from scipy.stats import t as t_student
from scipy.interpolate import InterpolatedUnivariateSpline
# from matplotlib.pylab import plot, hist, subplots, close, ion, ioff
from pandas import read_csv

#from FitK import FitK

def Mod(t):
    return t-2*pi*floor(t/(2*pi))

def Q( f, c1, c2):
    return c1*f**c2

def Af( f, R, M0, Const=10, c1=273, c2=0.66, beta=3.3, fc=1, alpha=-1, Rx=1.0):
    a = 1/(beta**3)
    b = M0*f**2 / (1+(f/fc)**2)
    c = exp(-pi * f * R / (beta*Q(f,c1,c2)))
    d = (R/Rx)**alpha
    return Const*a*b*c*d

loge_10 = 1/log(10)
def log10A( f, R, M0, Const=10, c1=273, c2=0.66, beta=3.3, fc=1, alpha=-1, Rx=1.0):
    a = -3*log10(beta)
    b = log10(M0) + 2*log10(f) - log10(1+(f/fc)**2)
    c = (-pi * f * R / (beta*Q(f,c1,c2)))/loge_10
    d = alpha*log10(R) - log10(Rx)
    return log10(Const) + a + b + c + d




class Acelerogram:
    
    def __init__( self, fnam, Dt, th=0.0, dir=[1,1,1], skiprows=1, delimiter=','):
        """Reads an acelerogram from from file fnam, in 4 columns,
           col 0 index and three channels.
           Col 1 is V, Col 2 H1 at th degrees from N and col 3, H2 perpendicular to H1.
           dir: Optional multipliers for the directions in each channel.
           If th=0 and dir=[1,1,1], V is up, H1 is N and H2 is E.
           The H and V chanels are automatically corrected to the above default.

           Dt: The time step, for the three channels.
           
           The channels are also detrend and applied a windoing (tukey) filter. 
        """
        self.data = loadtxt( fnam, skiprows=skiprows, delimiter=',')
        if self.data.shape[0] % 2 == 1:
            self.data = self.data[:-1,:] #Force to have even number of elements, remove last line if odd
        # zeros_index = where(npabs(self.data[:,2]) > 0.02)[0][1]
        zeros_index = where(npabs(self.data[:,2]) >= 0)[0][1]        
        self.data = self.data[zeros_index:,:] #Remove initial zeros
        ### Apply rotation matrix on -th
        tmp = self.data[:,2].copy()
        self.data[:,2] = dir[2]*self.data[:,3]*sin(-th) + dir[1]*self.data[:,2]*cos(-th)
        self.data[:,3] = dir[2]*self.data[:,3]*cos(-th) - dir[1]*tmp*sin(-th)
        self.N = self.data.shape[0]
        if Dt == -1: #The column 0 is time, otherwise is an index
            self.t = self.data[:,0]-self.data[0,0] #start at zero
            self.Dt = self.t[1]-self.t[0]
        else:
            self.Dt = Dt # Delta t
            self.t = arange(self.N)*self.Dt
        ### Detrend and basic filtering
        for i in [1,2,3]:
            self.data[:,i] = detrend(self.data[:,i])
            self.data[:,i] = tukey( self.N, alpha=0.05)*self.data[:,i]
            ### To see that the consecutive angles get distorted
            #self.data[:,i] += 20 * self.t
            ### sum(tst.data[:,2]) is big now

        if 1/self.Dt > 110:
            Dt = 1/100 ##Max 1/Dt Hz
            N = int(self.N*(self.Dt/Dt))
            t = linspace( 0, self.Dt*self.N, num=N)
            data = zeros((N,4))
            data[:,0] = t
            for i in [1,2,3]:
                spl = InterpolatedUnivariateSpline( self.t, self.data[:,i], k=3)
                data[:,i] = spl(t)
            del self.data
            self.data = data
            self.N = N
            self.Dt = Dt
            self.t = t

    
    def PlotChannels( self, ax=None):

        if ax == None:
            fig, ax = subplots(nrows=3, ncols=1, sharex=True)

        ax[0].plot( self.t, self.data[:,1], 'b.')
        ax[0].set_ylabel("V")
        ax[1].plot( self.t, self.data[:,2], 'b.')
        ax[1].set_ylabel("H1")
        ax[2].plot( self.t, self.data[:,3], 'b.')
        ax[2].set_ylabel("H2")
        
        return ax
    
    def __transform_hoz__( self, th):
        """ Make an horizontal component at th degrees,
            **assuming** a pi/2 (right angle) in the two horizontal components.
            (Uses a rotation matrix)
            
            cos(th)  -sin(th)   |  h2
                                |
            sin(th)   cos(th)   |  h1
        """
        ### "y", N component
        self.hoz =  self.data[:,3]*sin(th) + self.data[:,2]*cos(th)
    
    def PlotHoz( self, th=0.0, ax=None):
        if ax == None:
            fig, ax = subplots()
        ### Make self.hoz
        self.__transform_hoz__(th)
        ax.plot( self.t, self.hoz)
        ax.set_ylabel("H, %f degrees, gal" % (th))
        ax.grid(color='gray', linestyle='--', linewidth=0.5)
        
        return ax

    def Acc( self, t):
        """Return the horizontal aceleration at time t,
           with linear interpolation.
        """
        n = int(floor(t/self.SA_Dt))
        if (n >= self.hoz.size-1) or (n < 0):
            #print("Warning, acceleration asked beyong limits")
            return 0.0 #acelerograms start anf finish at 0.
        return self.hoz[n] + (t-self.SA_t[n])*((self.hoz[n+1]-self.hoz[n])/self.SA_Dt)

    def Acc_v( self, vt):
        """Vectorized version of Acc."""
        return array([self.Acc(t) for t in vt])
        

    def SA_init( self, T, th=0.0, hoz=None, Dt=None, xi=0.05):
        """Initialize the Spectral Response analysis.
           T: Natural period of vibration, in s
           th: angle to transform the horizontal components, rad (=0.0)
               if None, used the provided hoz and Dt as data
           xi: Damping factor (=0.05)
           
           If th == None, then use hoz as time series with time interval Dt.
        """
        ### Make self.hoz
        if th == None:
            self.hoz = hoz
            self.SA_Dt = Dt
            self.SA_t = arange(self.hoz.size)*self.SA_Dt
        else:
            self.__transform_hoz__(th)
            self.SA_Dt = self.Dt
            self.SA_t = self.t
        self.T = T
        self.w0 = 2*pi/T #natural frequancy, in rad/s
        self.xi = xi
        self.par1 = -2*xi*self.w0
        self.par2 = -1*self.w0**2
        self.sa_jac = zeros((2,2))
        self.sa_jac[0,0] = self.par1
        self.sa_jac[0,1] = self.par2
        self.sa_jac[1,0] = 1.0

    def SA_rhs( self, t, X):
        return self.par1*X[0] + self.par2*X[1] + self.Acc(t), X[0]

    def SA_jac( self, t, X):
        return self.sa_jac

    def SA_solve(self, f_factor=1.5, method='LSODA'):
        self.rt = solve_ivp( fun=self.SA_rhs, t_span=(self.t[0],self.t[-1]*f_factor), y0=[0,0],\
          method=method, vectorized=False, jac=self.SA_jac)#, t_eval=tst.t)
        self.SA = self.par1*self.rt['y'][0,:] + self.par2*self.rt['y'][1,:] + self.Acc_v(self.rt['t'])
        self.SA_max_indx = argmax(abs(self.SA))
        self.SA_max = abs(self.SA[self.SA_max_indx])
        self.SD_max_indx = argmax(abs(self.rt['y'][1,:]))
        self.SD_max = max(abs(self.rt['y'][1,:]))

    def SA_plot(self, ax=None, color='blue'):        
        if not(isinstance( ax, ndarray)):
            fig, ax = subplots(nrows=4, ncols=1, sharex=True)

        hoz = self.Acc_v(self.rt['t']) 
        ax[0].plot( self.rt['t'], hoz, '-', color=color, alpha=0.5)
        ax[0].grid(color='gray', linestyle='--', linewidth=0.5)
        ax[0].set_ylabel(r"$a_b$, Gal")
        ax[0].set_title("T=%4.2f: max $a$=%6.1f Gal, PSV (SA): %6.1f Gal" % (self.T,self.SA_max,self.SD_max * self.w0**2))
        
        ax[1].plot( self.rt['t'], self.SA, '-', color=color, alpha=0.5)
        ax[1].set_ylabel(r"$a$, Gal")
        ax[1].grid(color='gray', linestyle='--', linewidth=0.5)
        ax[1].plot( self.rt['t'][self.SA_max_indx], self.SA[self.SA_max_indx], 'k.')
        ax[1].plot( self.rt['t'][self.SD_max_indx], self.SD_max * self.w0**2, '.', color=color)
        ax[1].text( self.rt['t'][self.SD_max_indx], self.SD_max * self.w0**2,\
              "%6.1f" % (self.SD_max * self.w0**2), fontsize=10, color=color, ha='right')
        ax[2].plot( self.rt['t'], self.rt['y'][0,:], '-', color=color, alpha=0.5)
        ax[2].set_ylabel(r"$v, cm/s$")
        ax[2].grid(color='gray', linestyle='--', linewidth=0.5)
        ax[3].plot( self.rt['t'], self.rt['y'][1,:], '-', color=color, alpha=0.5)
        ax[3].set_ylabel(r"$u, cm$")
        ax[3].grid(color='gray', linestyle='--', linewidth=0.5)
        ax[3].set_xlabel("$t, s$")
        
        for i in range(4):
            yl = max([-ax[i].get_ylim()[0], ax[i].get_ylim()[1]])
            ax[i].set_ylim((-yl,yl))
        return ax

    def SA_T_plot( self, T, th=0.0, ax=None):
        self.SA_init(T, th=0.0)
        self.SA_solve()
        ax = self.SA_plot(ax=ax)
        return ax

    """ Old:
    def Fourier( self):
        ###Apply fft transform on H components.
        tmp = fft(self.data[:,2])
        #The input is real, so the result is symetric, we keep only the first part
        #tmp.size was forced to be even
        self.A_size = tmp.size//2
        self.A2 = npabs(tmp[1:self.A_size])
        self.A2_th = angle(tmp[1:self.A_size])
        tmp = fft(self.data[:,3])
        self.A3 = npabs(tmp[1:self.A_size])
        self.A3_th = angle(tmp[1:self.A_size])
        ### Frequencies
        self.freq = linspace( 0, 1 / self.Dt, self.A_size-1)
        #self.freq = fftfreq( 2*self.A_size, d=self.Dt)
    """

    def Fourier(self, th=0.2):
        """Apply fft transform on H components, with angle th."""
        ### Creates the hoz variable
        self.__transform_hoz__(th)
        self.A = fft(self.hoz)
        #The input is no longer real, the result is not symetric,
        #but we neverthelss keep only the first part in moduls and angle
        #self.size was forced to be even
        self.A_size = self.A.size//2
        self.A_r = npabs(self.A[1:self.A_size])
        self.A_th = angle(self.A[1:self.A_size])
        ### Frequencies
        self.freq = linspace( 0, 1 / self.Dt, self.A_size-1)

    def A_plot(self, xlim=(-1,1.5), ax=None):        
        """Plot A(f)."""
        if ax==None:
            fig, ax = subplots()

        self.Fourier()
        ax.plot( log10(tst.freq[1:]), 0.5*log10(tst.A2[1:])+0.5*log10(tst.A3[1:]), '-')
        ax.grid(color='gray', linestyle='--', linewidth=0.5)
        ax.set_ylabel(r"$\log ~ A(f)$")
        ax.set_xlabel(r"log Hz")
        ax.set_xlim(xlim)






from numpy.random import permutation 
from scipy.ndimage import uniform_filter1d
from geopy.distance import distance

IDSismos = read_csv('./FinalGCMT.csv')
IDSismos = IDSismos.loc[IDSismos["GMW"]!=0] #Quitamos los sismos que no tenemos magnitud de momentos
IDSismos = IDSismos.loc[~IDSismos["Base"].isin(["PNIG9701.211","TNLP9509.142","MEZC9509.141","COPL9509.142"])] #Quitamos sismos con problemas de formato
#l_sismo = arange(242)[permutation(242)[:50]] # random selection of 50 registers
#l_sismo.sort()
#l_sismo = arange(242)
l_sismo = [238] #241
    
sismo_i=l_sismo[0]

### Sismogramas con problemas
#l_sismo = [13,191,81,103,29,66,45,200,190,66]

T_list = [ 0.5, 1, 1.5, 2, 3, 4]
lam = 9
tM=12





f_cut = 50
N = 50
### Select N frequncies to do the model fit, +- uniform in log10 and orig sacale
#points_f = append(10**linspace( log10(0.05), log10(5), num=20, endpoint=False),\
#                  10**linspace( log10(3), log10(f_cut), num=15))
points_f = linspace( 0.05, f_cut, num=N)

points_f.sort()

id_sismo = "Sismo"+IDSismos.loc[l_sismo]['IDSismo']
fnam_csv = IDSismos.loc[l_sismo]['Base'].replace(".","")+".csv"
fnam = "Sismos/" + id_sismo + "/" + fnam_csv

# tst = Acelerogram( "/Users/isaias/Desktop/SMR29603131.csv", Dt=-1, skiprows=0)
# tst = Acelerogram( "/Users/isaias/Desktop/RegistroChristen2.csv", Dt=-1, skiprows=1)
    #h_angle = 0.0 # Angle for the horizontal component, 0 E, pi/2 N

def RS(tst):
    for h_angle,h_name in [[pi/2,"NN"]]: #[0.0,"EE"],[pi/2,"NN"] [[0.0,"EE"],[pi/4,"NE"],[pi/2,"NN"],[3*pi/4,"NW"]]: #,[pi,"WW"],[5*pi/4,"SW"],[3*pi/2,"SS"],[7*pi/4,"SE"]]
        tst.Fourier(th=h_angle)
        K=argmax(tst.freq[tst.freq < f_cut]) #tst.A_size-1
        ### Simulate the synthetic phases
        t_l = t_student(df=1.5)
        th = zeros(K)
        for i in range(1,K):
            th[i] = th[i-1] - 2*pi*(tM/tst.t[-1]) + (pi/lam)*t_l.rvs()
        #U = uniform(loc=-1,scale=2) # U(-1,1)
        #th = th + (pi/2)/lam*U.rvs(size=K)
        #th[1:K] = 2*pi*(tst.freq[1:K]*tM - floor(tst.freq[1:K]*tM)) #- (pi/2)/lam*t_l.rvs(size=K-1)
        th = Mod(th)-pi
        ### Filter the norms of Fourire transfm, ie denosie
        Afil = uniform_filter1d( tst.A_r, size=N)
        spl = InterpolatedUnivariateSpline( tst.freq, Afil, k=3)
        
        points = spl(points_f)
        ### Fit the gamma type model to the 50 points
        P = N*(points_f[1]-points_f[0])
        #fk = FitK( points_f, points, P=N*(points_f[1]-points_f[0]), p_a=3, p_b=3*20**2)
        #fk.RunTwalk(T=100000, burnin=2000)
        #E, k, la, f0, p = fk.map
        ### Evaluate the fit in all frequencies
        #A_r = fk.Fit(tst.freq) #*(2*npround(uniform.rvs(size=tst.freq.size))-1)
    
        ### Create synthetic Fourier, with smoothed norms and
        ### synthetic phases
        #A_synth = A_r[:K]*(cos(th) + sin(th)*(1.0j))
        #A_synth = append( [0.0], A_synth) #Add the first elemnt, K+1 length
        ### Obtain the sismogram
        #res2 = irfft(A_synth)
        ### Force it to start and end in zero ???
        #res2 = tukey( 2*K, alpha=0.05)*res2
        ### And its Dt
        Dt_synth = tst.N*tst.Dt/(2*K)
        t = arange(2*K)*(Dt_synth)
    
        """
        for T in [1,2,4]:
            tst.SA_init( T=T )
            tst.SA_solve()
            ax1 = tst.SA_plot(color='red')
        
            tst.SA_init( T=T, th=None, hoz=res2, Dt=Dt_u)
            tst.SA_solve()
            tst.SA_plot(ax=ax1, color='blue')
        """
        
        fig, ax = subplots( nrows=2, ncols=1, figsize=(3*3.23,5.51))
    
        # ax[0,0].plot( log10(tst.freq[1:]), log10(tst.A_r[1:]), 'k-')
        # ax[0,0].plot( log10(tst.freq[1:K]), log10(tst.A_r[1:K]), 'b-')
        # ax[0,0].plot( log10(tst.freq[1:K]), log10(Afil[1:K]), '-', color='lightblue')
        # ax[0,0].plot( log10(points_f), log10(points), '.', color='lightblue')
        # # ax[0,0].plot( log10(tst.freq[1:K]), log10(A_r[1:K]), '-', color='darksalmon')
        # ax[0,0].set_ylim((0.0,4.5))
        # ax[0,0].set_ylabel(r"$\log ~ A(f)$")
        # ax[0,0].set_xlabel(r"Hz")
        # ax[0,0].set_xticks([-1,0,1,log10(1/tst.Dt)])
        # ax[0,0].set_xticklabels([r"$10^{-1}$",r"$10^{0}$",r"$10^{1}$",r"$\Delta_t^{-1} = %2.0f$" % (1/tst.Dt)])
        # ax[0,0].grid(color='gray', linestyle='--', linewidth=0.5)
         
        # ax[1,0].plot( tst.A_th[:(K-1)], tst.A_th[1:K], 'b.')
        # ax[1,0].set_xticks([-pi,-pi/2,0,pi/2,pi])
        # ax[1,0].set_xticklabels([r"$-\pi$",r"$-\pi/2$",r"0",r"$\pi/2$",r"$\pi$"])
        # ax[1,0].set_xlabel(r"$\theta_n$")
        # ax[1,0].set_yticks([-pi,-pi/2,0,pi/2,pi])
        # ax[1,0].set_yticklabels([r"$-\pi$",r"$-\pi/2$",r"0",r"$\pi/2$",r"$\pi$"])
        # ax[1,0].set_ylabel(r"$\theta_{n+1}$")
        # ax[1,0].grid(color='gray', linestyle='--', linewidth=0.5)
    
        ax[0].plot( tst.t, tst.hoz, 'b-')
        ylim = ax[0].get_ylim()
        #ax[1,1].plot( t, res2, 'r-', alpha=0.5)
        yl = max([-ylim[0], ylim[1]])
        ax[0].set_ylim((-yl,yl))
        ax[0].set_ylabel(r"$\frac{cm}{s^2}$, gal")
        ax[0].set_xlabel(r"s")
        ax[0].grid(color='gray', linestyle='--', linewidth=0.5)
        
    
        SR = []
        SR_synth = []
        for T in T_list:
            tst.SA_init( T=T, th=h_angle )
            tst.SA_solve()
            ax[1].plot( T, tst.SD_max * tst.w0**2, 'bo')
            SR += [tst.SD_max * tst.w0**2]
        ax[1].set_xlim((0, 1.05*T_list[-1]))
        ax[1].set_xticks(T_list)
        ax[1].set_ylim((0,max(0.1,1.05*max(SR+SR_synth))))
        ax[1].set_xlabel(r"$T$")
        ax[1].set_ylabel(r"$\frac{cm}{s^2}$, gal")
        ax[1].grid(color='gray', linestyle='--', linewidth=0.5)
    
    
        R = distance((IDSismos.loc[sismo_i]["LatEstacion"],-IDSismos.loc[sismo_i]["LonEstacion"]), (IDSismos.loc[sismo_i]["LatSismo"],IDSismos.loc[sismo_i]["LonSismo"])).km
        title_w = "%d, %s, %s, Mw: %3.1f, R= %6.1f km, %s" %\
            ( sismo_i, id_sismo + "/" + fnam_csv,\
             IDSismos.loc[sismo_i]['Fecha'], IDSismos.loc[sismo_i]['GMW'],\
             R, h_name)
        fig.suptitle( title_w, fontsize=14)
        fig.tight_layout()
        #close(fig)




def Respuesta(tst):
    for h_angle,h_name in [[pi/2,"NN"]]: #[0.0,"EE"],[pi/2,"NN"] [[0.0,"EE"],[pi/4,"NE"],[pi/2,"NN"],[3*pi/4,"NW"]]: #,[pi,"WW"],[5*pi/4,"SW"],[3*pi/2,"SS"],[7*pi/4,"SE"]]
        tst.Fourier(th=h_angle)
        K=argmax(tst.freq[tst.freq < f_cut]) #tst.A_size-1

        SR = []
        SR_synth = []
        for T in T_list:
            tst.SA_init( T=T, th=h_angle )
            tst.SA_solve()
            #plt.plot( T, tst.SD_max * tst.w0**2, 'bo')
            SR += [tst.SD_max * tst.w0**2]
    return T_list,SR
        #close(fig)








########################




def RespuestaSimuladaNueva(PSDF,t2,f2,idsismo):
    Angle=np.random.uniform(0,2*np.pi,size=PSDF.shape[0])
    Proceso=np.ones(PSDF.shape[1])
        
    
    # for i in range(PSDF.shape[1]):
    #     Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(2*PSDF.clip(0)[:,i]*np.pi*f2[1])*np.cos(2*np.pi*np.arange(len(f2))*f2[1]*t2[i]+Angle))
    
    # for i in range(PSDF.shape[1]):
    #     #Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(2*PSDF.clip(0)[:,i]*f2[1])*np.cos(np.arange(len(f2))*f[1]*t2[i]+Angle))
    #     Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(PSDF.clip(0)[:,i]*f2[1])*np.cos(2*np.pi*np.arange(len(f2))*f2[1]*t2[i]+Angle))
    # from datetime import datetime
    # now = datetime.now()   
    # for i in range(PSDF.shape[1]):
    #     #Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(2*PSDF.clip(0)[:,i]*f2[1])*np.cos(np.arange(len(f2))*f[1]*t2[i]+Angle))
    #     Proceso[i]=np.sqrt(2)*np.sum(np.sqrt(PSDF.clip(0)[:,i]*f2[1])*np.cos(2*np.pi*np.arange(len(f2))*f2[1]*t2[i]+Angle))
    # print(datetime.now()-now)
    
    # now = datetime.now()    
    R=(np.reshape(2*np.pi*f2[1]*np.arange(len(f2)),(len(f2),1))@np.reshape(t2,(1,len(t2))))
    R+=np.reshape(Angle,(Angle.shape[0],1))
    R=np.cos(R)
    R=np.einsum('ij,ji->i', np.sqrt(2)*(np.sqrt(PSDF.clip(0)[:,:]*f2[1])).T,R)
    # plt.plot(Proceso-R)
    # plt.plot(R)    
    # print(datetime.now()-now)

    Proceso=np.copy(R)
    
    # plt.plot(X,NS,alpha=0.6)
    # plt.plot(t2,Proceso,alpha=0.6)
    
    # print(integrate.cumtrapz(NS**2,X , initial=0)[-1],integrate.cumtrapz(Proceso**2,X , initial=0)[-1])
    
    # (fsim,tsim,Zxxsim)=sp.signal.stft(Proceso,fs=nperseg,nperseg=batch,window=Window,noverlap=over)#los indices estan volteados respecto al articulo
    
    
    
    # trace1 = go.Surface(x=t,y=f,z=Z, colorscale='Viridis')
    # trace1 = go.Surface(x=tsim,y=fsim,z=np.abs(Zxxsim)**2, colorscale='Viridis')
        
    # iplot([trace1])
    
    
    
    
    Sintetico=Proceso
    Registro=np.vstack((t2,Sintetico))
    Registro=np.vstack((Registro,Sintetico))
    Registro=np.vstack((Registro,Sintetico))
    Registro.shape
    # ############### Sismos registrados lejanos
    
    # # Muestra.sort_values(by=["R"])[-10:]
    # # STFT(59,10,90) #Canditato
    # # STFT(52,10,90) #Canditato
    # # STFT(116,10,90) #Canditato
    # # STFT(168,10,90) #Canditato
    # # Muestra.sort_values(by=["R"])[-30:-20]
    
    
    RegistroChristen=Registro.T
    np.savetxt(r'./RegistroChristen2'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv',RegistroChristen,delimiter=',', fmt=('%s'))
    # # np.savetxt(r'./RegistroChristen2.csv',RegistroChristen,delimiter=',', fmt=('%s'))
    tst = Acelerogram( r'./RegistroChristen2'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv', Dt=-1, skiprows=1)
    os.remove(r'./RegistroChristen2'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv')

    return Proceso, np.array(Respuesta(tst))[1]




##########################################################








def MapaCompleto(ReverseFault):
    contador=0
    
    for sismo in np.unique(ReverseFault["IDSismo"]):
        
        CoordEpic=np.array(ReverseFault[["LatSismo","LonSismo"]][ReverseFault["IDSismo"]==sismo].iloc[0])
        CoordEstacion=(ReverseFault[["LatEstacion","LonEstacion"]][ReverseFault["IDSismo"]==sismo])
        CoordEstacion=CoordEstacion.to_numpy()
        
        longs=np.ones(0)
        lats=np.ones(0)
        for i in range(len(CoordEstacion)):
            longs=np.hstack((longs,np.hstack((CoordEpic[1],CoordEstacion[i][1]))))
            lats=np.hstack((lats,np.hstack((CoordEpic[0],CoordEstacion[i][0]))))
        
        if contador==0:
            fig3 = go.Figure(go.Scattermapbox(
                mode = "markers+lines",
                lon = longs,
                lat = lats,
                marker = {'size': 10}))
            fig3.update_layout(mapbox_style="open-street-map")
            fig3.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
            contador+=1
        
        else :
        
            fig3.add_trace(go.Scattermapbox(
                mode = "markers+lines",
                lon = longs,
                lat = lats,
                marker = {'size': 10}))
    
    
    plot(fig3)
    
    fig3.write_html("Prueba.html")

#####################################################





def GaussianaAjustada(idsismo,b):
    f,t,Z=SePSDF(idsismo,"H1",4)#, Muestra=ReverseFault)                
    
    
    ##############Esto trunca las frecuencias
    Z=Z[f<fmax,:]
    f=f[f<fmax]
    
    Marginalf=np.ones(len(f))
    
    for i in np.arange(len(f)):
        Marginalf[i]=integrate.cumtrapz(Z[i,:],t , initial=0)[-1]
        
    Const=integrate.cumtrapz(Marginalf,f,initial=0)[-1]    
        
    Marginalf=Marginalf/integrate.cumtrapz(Marginalf,f,initial=0)[-1]
    
    # plt.plot(f,Marginalf)
    
    Marginalt=np.ones(len(t))
    
    for i in np.arange(len(t)):
        Marginalt[i]=integrate.cumtrapz(Z[:,i],f , initial=0)[-1]
        
    
    Marginalt=Marginalt/integrate.cumtrapz(Marginalt,t,initial=0)[-1]
    
    
    Znorm=Z/Const
        

    DistMarginalt=integrate.cumtrapz(Marginalt,t)
    DistMarginalf=integrate.cumtrapz(Marginalf,f)
    
    
    # b=-0.9
    
    Sigma=np.array([[1,b],[b,1]])
    # Sigma=np.array([[a,b],[b,c]])
    SigmaInv=np.linalg.inv(Sigma)
    SigmaDet=np.linalg.det(Sigma)
    
    
    
    InvNormalf=sp.stats.norm.ppf(DistMarginalf) 
    InvNormalt=sp.stats.norm.ppf(DistMarginalt) 
    
    global CopulaDensidad
    
    def CopulaDensidad(X):
        i=X[0]
        j=X[1]
        Vec=np.array([InvNormalf[i],InvNormalt[j]])
        # Vec=np.array([DistMarginalf[i],DistMarginalt[j]])
        Vec.shape=(2,1)
        Valor=0
        
        if  all(~np.isinf(Vec)) :
            Valor= np.exp(-Vec.T@(SigmaInv-np.identity(2))@Vec/2)*SigmaDet**(-1/2)* (Marginalt[j] * Marginalf[i]) 
        if np.isinf(Valor) or np.isnan(Valor) :
            Valor=0            
        return Valor         
    
    from scipy.ndimage import gaussian_filter
    
    from multiprocessing import Pool
    from multiprocessing import get_context
    import multiprocessing 
    
    x=np.arange(len(f)-1)
    y=np.arange(len(t)-1)    
    Malla=pd.DataFrame({'x':np.repeat(x,y.shape[0]),
                  'y':np.tile(y,x.shape[0])})
    Malla=Malla.to_numpy()
    
    
    Hilos=multiprocessing.cpu_count()
    p = get_context("fork").Pool(Hilos)
    Valores=p.map(CopulaDensidad,Malla)  
    p.close()    
    
    
    Conjunta=np.zeros((len(f)-1,len(t)-1))
    
    for i in range(Malla.shape[0]):    
        loc=Malla[i,:]
        Conjunta[loc[0],loc[1]]=Valores[i]
    
    return t[:-1],f[:-1][f[:-1]<fmax],Conjunta[f[:-1]<fmax,:]




#######################################




################################################## Recupera momentos



########################################Frecuencia

from scipy.stats import gamma
from scipy.optimize import minimize


def abFrecPredicho(Dis):
    muF,s2F=PrediceFrecMedia(Dis),PrediceFrecVar(Dis)
    
    def Objetivo(X):
        alpha,beta=X
        mu=(alpha/beta)*gamma.cdf(fmax, a=alpha+1,scale=1/beta)/gamma.cdf(fmax, a=alpha,scale=1/beta)
        mu2=((alpha+1)*alpha/beta**2)*gamma.cdf(fmax, a=alpha+2,scale=1/beta)/gamma.cdf(fmax, a=alpha,scale=1/beta)
        sigma2=mu2-mu**2
        return (muF-mu)**2+(s2F-sigma2)**2
    
    
    Minimiza=minimize(Objetivo, [2,2], method='Nelder-Mead', tol=1e-6, bounds=((1, None), (0, None)))
    
    alpha,beta=Minimiza["x"]
    return alpha,beta
    
    # MuestraGamma=gamma.rvs(size=10000, a=alpha,scale=1/beta)
    # (np.mean(MuestraGamma[MuestraGamma<fmax]),np.var(MuestraGamma[MuestraGamma<fmax]))
    # (muF,s2F)

######################################## Tiempo



# def abTimePredicho(Dis):
#     mu=PrediceAmpMedia(Dis)
#     s=PrediceAmpSD(Dis)**2

#     alpha=mu**2/s**2
#     beta=mu/s**2    
#     return alpha,beta



def abTimePredicho(Dis):
    mu=PrediceAmpMedia(Dis)
    s2=PrediceAmpSD(Dis)**2

    def Objetivo(X):
        alpha,beta=X
        muE=alpha/beta
        s2E=alpha/beta**2
        return (muE-mu)**2+(s2E-s2)**2
    Minimiza=minimize(Objetivo, [2,2], method='Nelder-Mead', tol=1e-6, bounds=((1, None), (0, None)))
    alpha,beta=Minimiza["x"]
    return alpha,beta





def Predicciones(r,afp,bfp,atp,btp,PredET,PredRho):

    global CopulaDensidad
    idsismo=R2.index[r]
    CanalH1=np.int(np.float(RecuperaOrientacion(Muestra["Base"][idsismo])[8]))

    (X,NS,f,Y,VelMuest)=fourier(ReverseFault["Base"][idsismo], CanalH1)   
    NS=NS[X>ReverseFault["Inicio"][idsismo]]
    X=X[X>ReverseFault["Inicio"][idsismo]]
    
    Dominiof=np.arange(0,fmax,VelocidadMuestreoFrec)
    Dominio=np.arange(np.min(Directorio2[r][0,:]),np.max(Directorio2[r][0,:]),VelocidadMuestreo)
    
    #######################################################
    
    
    DensMarginalf=gamma.pdf(Dominiof, afp,scale=1/bfp)/gamma.cdf(fmax, afp,scale=1/bfp)###############################Marginalf
    DensMarginalt=gamma.pdf(Dominio, atp,scale=1/btp)###############################Marginalt
    
    
    
    DistMarginalt=integrate.cumtrapz(DensMarginalt,Dominio)
    DistMarginalf=integrate.cumtrapz(DensMarginalf,Dominiof)
    
    
    ##########################3
    
    b=PredRho ###############################RhoCopula
    
    Sigma=np.array([[1,b],[b,1]])
    # Sigma=np.array([[a,b],[b,c]])
    SigmaInv=np.linalg.inv(Sigma)
    SigmaDet=np.linalg.det(Sigma)
    
    
    
    InvNormalf=sp.stats.norm.ppf(DistMarginalf) 
    InvNormalt=sp.stats.norm.ppf(DistMarginalt) 
    
    
    
    def CopulaDensidad(X):
        i=X[0]
        j=X[1]
        Vec=np.array([InvNormalf[i],InvNormalt[j]])
        # Vec=np.array([DistMarginalf[i],DistMarginalt[j]])
        Vec.shape=(2,1)
        Valor=0
        
        if  all(~np.isinf(Vec)) :
            Valor= np.exp(-Vec.T@(SigmaInv-np.identity(2))@Vec/2)*SigmaDet**(-1/2)* (DensMarginalt[j] * DensMarginalf[i]) 
        if np.isinf(Valor) or np.isnan(Valor) :
            Valor=0            
        return Valor         
    
    from scipy.ndimage import gaussian_filter
    
    from multiprocessing import Pool
    from multiprocessing import get_context
    import multiprocessing 
    
    x=np.arange(len(Dominiof)-1)
    y=np.arange(len(Dominio)-1)    
    Malla=pd.DataFrame({'x':np.repeat(x,y.shape[0]),
                  'y':np.tile(y,x.shape[0])})
    Malla=Malla.to_numpy()
    
    
    Hilos=multiprocessing.cpu_count()
    p = get_context("fork").Pool(Hilos)
    Valores=p.map(CopulaDensidad,Malla)  
    p.close()    
    
    
    Conjunta=np.zeros((len(Dominiof)-1,len(Dominio)-1))
    
    for i in range(Malla.shape[0]):    
        loc=Malla[i,:]
        Conjunta[loc[0],loc[1]]=Valores[i]
    
    
    Sintetico=NS
    Registro=np.vstack((X,NS))
    Registro=np.vstack((Registro,Sintetico))
    Registro=np.vstack((Registro,Sintetico))
    Registro.shape
    
    RegistroChristen=Registro.T
    np.savetxt(r'./RegistroChristen'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv',RegistroChristen,delimiter=',', fmt=('%s'))
    tst = Acelerogram( r'./RegistroChristen'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv', Dt=-1, skiprows=0)    
    Real=np.array(Respuesta(tst))
    os.remove(r'./RegistroChristen'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv')
    
    
    
    
    
    
    
    Const=PredET###############################Et
    # Const=EnTot
    PSDFEst=Conjunta*Const
    
    
    
    
    Proceso,Simulado=RespuestaSimuladaNueva(PSDFEst,Dominio[:-1],Dominiof[:-1],idsismo)
    for i in range(Trayectorias):
        print(i+1)
        np.random.seed(i+1)
        Proceso,Nuevo=RespuestaSimuladaNueva(PSDFEst,Dominio[:-1],Dominiof[:-1],idsismo)
        Simulado=np.vstack((Simulado,Nuevo))    
    
    Simulado=np.array(Simulado)
    
    #######################################3
    
    
    f,t,Z=SePSDF(idsismo,"H1",Ventana)       

    return [Proceso,Real,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]

    
    
    
    
    




















































