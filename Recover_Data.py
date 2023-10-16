#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 15:53:43 2023

@author: isaias
"""

############### Hablar de la base de datos, red acelerografica, Kramer, Espacial, Multisitio-Multiestacion

###############################################
import os
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

############################################


### Lee las fechas
# fechas = read_csv("fechas.csv", sep='\n').loc[:].values
# fechas = [f[0] for f in fechas]
mes_texto = ['ENE', 'FEB', 'MAR', 'ABR', 'MAY', 'JUN', 'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DIC']

def ConvierteFecha(f):
    ### La pasamos a mayúsculas y le quitamos espacios al principio y al final
    ### Sustituimos las diagonales por espacios
    f = f.upper().strip().replace(  '/', ' ')

    ### Si los cuatro primeros son dígitos, fecha al reves eg. 1999/09/12
    if f[0:4].isdigit():
        y = int(f[0:4])
        m = int(f[5:7])
        d = int(f[8:])
        return date( y, m, d)
    
    ### los dos primeros son dígitos, con el dia, fecha convencional
    if f[0:2].strip().isdigit():
        ### Dia
        d = int(f[0:2].strip())
        
        ### Año, los dos último dígitos son la fecha
        if f[-2:].isdigit():
            y = int(f[-2:])
            if y > 50: #Cambiar 70 por 50
                y += 1900
            else:
                y += 2000
        else:
            print("Año no reconocido.")
            return None
        
        ### Mes
        m = 0
        for i,mes in enumerate(mes_texto):
            if f.find(mes) > 0:
                m += i+1
                break
        if m == 0: #Mes no está en texto, solo queda que sea el numero de enmedio
            if f[3:5].isdigit():
                m = int(f[3:5])
            else:
                print("Mes no reconocido:")
                return None
        return date( y, m, d)
    else:
        print("Formato no reconocido:")
        return None


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

#Recupera las posiciones de un caracter
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

#Recibe una base y devuelve el dominio la aceleracion las frecuencias y la FFT
# Bases[1]
# Base="PNIG9701.211"
# canal=1
# Base=Bases[1]






def Recuperacoordenada(Base,Tipo=1):    
    #El 1 recupera Epicentro
    aux=np.array(["N","W"])
    if Tipo==1:
        Coordenada=Buscar("COORDENADAS DEL EPICENTRO", Base)
    else:
        Coordenada=Buscar("COORDENADAS DE LA ESTACION", Base)
    Coordenada=Coordenada.upper()    
    Coordenada=Coordenada.replace("LANT","LAT").replace("LONG","LON").replace("LOG","LON").replace(":","") .replace("LST","LAT")
    if Coordenada.find("LAT")!=-1:
            Coordenada=Coordenada.replace(" ","")    
            Lat=Coordenada[:Coordenada.find("LAT")]
            Lon=Coordenada[Coordenada.find("LAT")+5:Coordenada.find("LON")]
            Lon = re.sub(r"[a-z]", "", Lon, flags=re.I).replace(" ","").replace("'","")

    
    else:

        if len(Coordenada.replace(" ",""))<4:
            Lat=""
            Lon=""
    
        elif Coordenada.find(" ")!=-1:
            Lat=Coordenada[:Coordenada.find(" ")]
            Lon=Coordenada[Coordenada.find(" ")+1:]


    if Lat!="" and len(Posiciones(Lat,"."))>1:
        Lat=Lon[:(Posiciones(Lat,".")[0]+1)]+Lon[(Posiciones(Lat,".")[0]+1):].replace(".","") 

    if Lon!="" and len(Posiciones(Lon,"."))>1:
        Lon=Lon[:(Posiciones(Lon,".")[0]+1)]+Lon[(Posiciones(Lon,".")[0]+1):].replace(".","") 


    
    if Coordenada.find("S")!=-1:
        aux[0]="S"
    
    if Coordenada.find("E")!=-1:
        aux[1]="E"
    aux=np.hstack((aux,(Lat,Lon)))
    return aux    
    



def Separacadenas(Cadena, lugares):

    return (Cadena[(lugares[0]+1):lugares[1]],Cadena[(lugares[1]+1):lugares[2]],Cadena[(lugares[2]+1):])

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
    
############################################


########################################### Borrador base

# Coordenadaestacion=np.ones(4)


BaseAuxiliar=np.ones(25)


for i in range(len(Bases)):         
    Base=Bases[i]
    Profundidad=Buscar("PROFUNDIDAD FOCAL", Base).replace(" ","").replace("-","")
    censurado=np.sum(Profundidad.find("<")!=-1)
    Profundidad=Profundidad.replace("<","")
    X=np.append(Base,np.concatenate((RecuperaOrientacion(Base),Recuperacoordenada(Base,2),Recuperacoordenada(Base))))
    X=np.append(X,ConvierteFecha(Buscar("FECHA DEL SISMO", Base)))
    X=np.append(X,Buscar("INSTITUCION RESPONSABLE",Base).replace(",", " "))
    X=np.append(X,Buscar("CLAVE DE LA ESTACION",Base).replace(",", " "))    
    X=np.append(X,Buscar("HORA EPICENTRO",Base).replace(",","."))
    X=np.append(X,Buscar("HORA DE LA PRIMERA MUESTRA",Base).replace(",","."))
    X=np.append(X,Profundidad)
    BaseAuxiliar=np.vstack((BaseAuxiliar,X))
    #print(Buscar("COORDENADAS DE LA ESTACION", Bases[i]), Recuperacoordenada(Bases[i],2))    
     
BaseAuxiliar[0,:]=np.array(["Base","OrientacionH1","OrientacionH2","SignoV","SignoH1","SignoH2",
                            "AnguloH1","AnguloH2","EjeV","EjeH1","EjeH2","Norte","Oeste","LatEstacion",
                            "LonEstacion","Norte","Oeste","LatSismo","LonSismo","Fecha","Institucion",
                            "IDEstacion","HoraEpicentro","HoraPrimerMuestra","Profundidad"])

np.savetxt(r'/Users/isaias/Desktop/Atenuacion/BaseAuxiliar.csv',BaseAuxiliar,delimiter=',', fmt=('%s'))



BaseAuxiliar=pd.read_csv('/Users/isaias/Desktop/Atenuacion/BaseAuxiliar.csv')

BaseAuxiliar=BaseAuxiliar.loc[BaseAuxiliar['OrientacionH1'] != 'Error']
BaseAuxiliar=BaseAuxiliar.loc[BaseAuxiliar['Fecha'] != 'None']
BaseAuxiliar=BaseAuxiliar.loc[~np.isnan(BaseAuxiliar['LatSismo'])] 
# BaseAuxiliar["ConteosEstacion"]=BaseAuxiliar.groupby(["LatEstacion","LonEstacion"])["LatEstacion"].transform("count")
# BaseAuxiliar["ConteosSismo"]=BaseAuxiliar.groupby(["LatSismo","LonSismo"])["LatEstacion"].transform("count")
BaseAuxiliar["LonEstacion"]=-BaseAuxiliar["LonEstacion"]
BaseAuxiliar["LonSismo"]=-BaseAuxiliar["LonSismo"]
BaseAuxiliar=BaseAuxiliar.loc[BaseAuxiliar['Fecha'] != 'None']
BaseAuxiliar=BaseAuxiliar.loc[BaseAuxiliar['Profundidad'] != 'None']
#BaseAuxiliar=BaseAuxiliar.loc[~np.isnan(BaseAuxiliar['Profundidad'])] 



###############

VecMagn=np.array(["LatSismo","LonSismo","Mw","Mb","Ms","Mc","Ma","Me","Md","Ml","M"])
BaseMagnitudes=np.array(["LatSismo","LonSismo","Mw","Mb","Ms","Mc","Ma","Me","Md","Ml","M"])
MagnIndiv=np.zeros(len(VecMagn))

SismosCord=BaseAuxiliar.pivot_table(index=['LatSismo','LonSismo'], aggfunc='size')
#SismosCord[50]
#2049
for i in range(len(SismosCord)):
    Sis=BaseAuxiliar.loc[BaseAuxiliar["LatSismo"]==SismosCord.index[i][0]]
    ObservacionesSismo=Sis.loc[Sis["LonSismo"]==SismosCord.index[i][1]]["Base"]
    MagnIndiv=np.zeros(len(VecMagn))
    MagnIndiv[[0,1]]=SismosCord.index[i][:]
    #ObservacionesSismo.index[0]
    
    for j in ObservacionesSismo.index:
        Magnitudes=Buscar("MAGNITU",Bases[j])
        Magnitudes=Magnitudes.replace(" ", "").replace(",","/")
        if Magnitudes.find("/")==-1:
            Magnitudes="/"+Magnitudes
        if len(Magnitudes)<2:
            print(" ")

        else: 
            Pos=Posiciones(Magnitudes,"/")

            if len(Pos)==1:
                print(Magnitudes[Pos[0]+1:])
                Partir=Magnitudes[Pos[0]+1:]
                if Partir.find("=")!=-1:
                    Igual=Partir.find("=")
                    MagnIndiv[VecMagn==Partir[:Igual]]=Partir[Igual+1:]
                else: 
                    Igual=2
                    MagnIndiv[VecMagn==Partir[:Igual]]=Partir[Igual:]                    

                    

            else:
                for l in range(len(Pos)-1):
                    print(Magnitudes[Pos[l]+1:Pos[l+1]])
                    Partir=Magnitudes[Pos[l]+1:Pos[l+1]]
                    if Partir.find("=")!=-1:
                        Igual=Partir.find("=")
                        MagnIndiv[VecMagn==Partir[:Igual]]=Partir[Igual+1:]
                    else: 
                        Igual=2
                        MagnIndiv[VecMagn==Partir[:Igual]]=Partir[Igual:]                    

                
                print(Magnitudes[Pos[l+1]+1:])
                Partir=Magnitudes[Pos[l+1]+1:]
                if Partir.find("=")!=-1:
                    Igual=Partir.find("=")
                    MagnIndiv[VecMagn==Partir[:Igual]]=Partir[Igual+1:]
                else: 
                    Igual=2
                    MagnIndiv[VecMagn==Partir[:Igual]]=Partir[Igual:]                    
        
        
        print(i,j,Bases[j],Magnitudes)
    
    print(MagnIndiv)
    BaseMagnitudes=np.vstack((BaseMagnitudes,MagnIndiv))


np.savetxt(r'/Users/isaias/Desktop/Atenuacion/BaseMagnitudes.csv',BaseMagnitudes,delimiter=',', fmt=('%s'))



BaseMag=pd.read_csv('/Users/isaias/Desktop/Atenuacion/BaseMagnitudes.csv')


#Redondeo por si acaso
BaseMag["LatSismo"]=np.round(BaseMag["LatSismo"],3)
BaseMag["LonSismo"]=np.round(BaseMag["LonSismo"],3)


BaseAuxiliar["LatSismo"]=np.round(BaseAuxiliar["LatSismo"],3)
BaseAuxiliar["LonSismo"]=np.round(BaseAuxiliar["LonSismo"],3)



Final=pd.merge(BaseAuxiliar,BaseMag,how="left",left_on = ['LatSismo', "LonSismo"], right_on = ['LatSismo', "LonSismo"])

Final.to_csv(r'/Users/isaias/Desktop/Atenuacion/BaseFinal.csv',index=False)



BaseFinal=pd.read_csv('/Users/isaias/Desktop/Atenuacion/BaseFinal.csv')



BaseMag2=BaseMag.loc[BaseMag["Ms"]!=0]
len(BaseMag2)
BaseMag2=BaseMag2.loc[BaseMag["Mb"]!=0]


############################################


###########################################Aqui empieza GCMT
######Encuentra registro
def GCMTRegistro(NombreSismo):
    P=open('/Users/isaias/Desktop/Atenuacion/Catalogo2.txt', errors="ignore")
    Datos=P.read()
    
    LocIdSismo=Datos.find(NombreSismo)
    
    InicioRegistro=Posiciones(Datos[:LocIdSismo],"\n")[-2]
    Cadena=Datos[(InicioRegistro+1):]
    Caracter="\n"
    Posicion=0
    Conteo=0
    for i in Cadena: 
        if i==Caracter:
            Posicion=np.hstack((Posicion,Conteo))
        #print(i,Conteo)
        Conteo+=1
        
        if(np.size(Posicion)>6):
            break
    
    Posicion=Posicion[1:]
    
    Registro=Datos[(InicioRegistro+1):(InicioRegistro+Posicion[-2]+1)]
    
    return Registro


def GCMTDatos(NombreSismo):
    Registro=GCMTRegistro(NombreSismo).replace("/", "-")
    
    while Registro!=Registro.replace("  ", " "):
        Registro=Registro.replace("  ", " ")
        
    
    Espacios=Posiciones(Registro," ")

    Fecha=Registro[Espacios[0]+1:Espacios[1]]
    
    Hora=Registro[Espacios[1]+1:Espacios[2]]
    
    Latitud=Registro[Espacios[2]+1:Espacios[3]]

    Longitud=Registro[Espacios[3]+1:Espacios[4]]
    
    CambioLinea=Posiciones(Registro,"\n")
    
    Registro[CambioLinea[3]:]
    
    Magnitud=float(Registro[CambioLinea[-2]+1:Espacios[Espacios>CambioLinea[-2]][0]])
    
    
    M0=float(Registro[Espacios[-7]+1:Espacios[-6]])
    
    
    Mw=np.round((2/3)*(np.log10(M0*(10**Magnitud)) - 16.1),1)

    Vector=Registro[(CambioLinea[2]+1):CambioLinea[3]]
    Vector=np.fromstring(Vector, dtype=float, sep=" ")[2*np.arange(6)+1]
    v1,v2,v3,v4,v5,v6=Vector[0],Vector[1],Vector[2],Vector[3],Vector[4],Vector[5]
    angle1=Registro[Espacios[-4]:Espacios[-3]]
    angle2=Registro[Espacios[-1]:]
    Profundidad=Registro[Espacios[Espacios>CambioLinea[1]][6]+1:Espacios[Espacios>CambioLinea[1]][7]]
    print(NombreSismo)
    return np.array([NombreSismo,Fecha,Hora,Latitud, Longitud, Mw,v1,v2,v3,v4,v5,v6,angle1,angle2, Profundidad])

NombreSismo="B082077A"
print(GCMTRegistro("201710211659A"))   #### GCMT presenta diferentes valores en hora si se consulta el ASCII a si se consulta en la web, ver registro B082077A
GCMTDatos("201710211659A")


ListaGCMTSismos=pd.read_csv('/Users/isaias/Desktop/Atenuacion/NombresBaseNueva.csv')
# GlobalCMTInformacion=(ListaGCMTSismos["Nombre"]).apply(GCMTDatos)

# import swifter
# GlobalCMTInformacion =(ListaGCMTSismos["Nombre"]).swifter.allow_dask_on_strings(enable=True).apply(GCMTDatos)


from pandarallel import pandarallel 
import multiprocessing

pandarallel.initialize(nb_workers=multiprocessing.cpu_count())

GlobalCMTInformacion=(ListaGCMTSismos["Nombre"]).parallel_apply(GCMTDatos)


len(GlobalCMTInformacion)
for i in range(len(GlobalCMTInformacion)):
    if i==0:
        GlobalCMTBase=pd.DataFrame(GlobalCMTInformacion[0]).transpose()
    if i!=0:
        GlobalCMTBase=pd.concat([GlobalCMTBase,pd.DataFrame(GlobalCMTInformacion[i]).transpose()])

GlobalCMTBase.columns=["NOMBRE","FECHA","HORA","LATITUD","LONGITUD","Mw","V1","V2","V3","V4","V5","V6","Angle1","Angle2","Profundidad"]

GlobalCMTBase=GlobalCMTBase.set_index(np.linspace(0,len(GlobalCMTInformacion)-1,len(GlobalCMTInformacion)).astype(int))

GlobalCMTBase.to_csv(r'/Users/isaias/Desktop/Atenuacion/BaseGCMT.csv',index=False)
###################################

GlobalCMTBase=pd.read_csv('/Users/isaias/Desktop/Atenuacion/BaseGCMT.csv')




BaseAuxiliar=pd.read_csv('/Users/isaias/Desktop/Atenuacion/BaseFinal.csv')
BaseAuxiliar=BaseAuxiliar.loc[BaseAuxiliar['OrientacionH1'] != 'Error']
BaseAuxiliar=BaseAuxiliar.loc[BaseAuxiliar['Fecha'] != 'None']
BaseAuxiliar=BaseAuxiliar.loc[BaseAuxiliar['HoraEpicentro'] != 'None']
BaseAuxiliar=BaseAuxiliar.loc[pd.notnull((BaseAuxiliar["HoraEpicentro"]))]
BaseAuxiliar=BaseAuxiliar.loc[~np.isnan(BaseAuxiliar['LatSismo'])] 
# BaseAuxiliar["ConteosEstacion"]=BaseAuxiliar.groupby(["LatEstacion","LonEstacion"])["LatEstacion"].transform("count")
# BaseAuxiliar["ConteosSismo"]=BaseAuxiliar.groupby(["LatSismo","LonSismo"])["LatEstacion"].transform("count")
BaseAuxiliar["LonEstacion"]=-BaseAuxiliar["LonEstacion"]
#BaseAuxiliar["LonSismo"]=-BaseAuxiliar["LonSismo"]


##############################################


import datetime

#for i in BaseAuxiliar.index[:]:
indice=0

#BaseAuxiliar.loc[BaseAuxiliar["Base"]==indice]["HoraEpicentro"].values[0]
def FechaSismo(indice):
    hora=BaseAuxiliar["HoraEpicentro"][indice].replace(" ","")
    hora=hora.replace(";",":")
    hora=hora.replace(".",":").replace("O","0")
        
    if len(Posiciones(hora,":"))==2:
        if hora.find(".")==-1:
           hora=hora+".0"
    elif len(Posiciones(hora,":"))==3:
        if(len(hora)>Posiciones(hora,":")[-1]+1):
            hora=hora[:Posiciones(hora,":")[-1]]+"."+hora[Posiciones(hora,":")[-1]+1]
        else:
            hora=hora[:Posiciones(hora,":")[-1]]+".0"
        
    elif len(Posiciones(hora,":"))==1:
         hora=hora+":00.0"
    hora=hora+"0"
   
    date_time_str=BaseAuxiliar["Fecha"][indice]+" "+hora
    return (datetime.datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S.%f'))




for i in range(len(GlobalCMTBase)):
    FechaGCMT=GlobalCMTBase["FECHA"][i]+" "+GlobalCMTBase["HORA"][i]
    if i==0:
        BaseFechas=np.array([GlobalCMTBase["NOMBRE"][i],datetime.datetime.strptime(FechaGCMT, '%Y-%m-%d %H:%M:%S.%f')])
    else :
        BaseFechas=np.vstack([BaseFechas,np.array([GlobalCMTBase["NOMBRE"][i],datetime.datetime.strptime(FechaGCMT, '%Y-%m-%d %H:%M:%S.%f')])])

Registro=15
TiempUmbral=20
def GCMTMagnitud(Registro,TiempUmbral,DistUmbral):
    if np.min(np.abs(BaseFechas[:,1]-FechaSismo(Registro)))<datetime.timedelta(seconds=TiempUmbral):
    
        Candidato=np.argmin(np.abs(BaseFechas[:,1]-FechaSismo(Registro)))
        
        coord1=(BaseAuxiliar["LatSismo"][Registro],BaseAuxiliar["LonSismo"][Registro])
        coord2=(GlobalCMTBase.loc[Candidato]["LATITUD"],GlobalCMTBase.loc[Candidato]["LONGITUD"])
        
        if distance.distance(coord1, coord2).km<DistUmbral:
            return GlobalCMTBase.loc[Candidato][["NOMBRE","Mw"]]
        else:
            return [0,0]
    else:
        return [0,0]

#######################################



GCMTMagnitud(0,300,20) #Tiempo es en segundos #Distancia en km


for i in BaseAuxiliar.index:
    if i==BaseAuxiliar.index[0]:
        Consulta=GCMTMagnitud(i,200,200)
        GMW=Consulta[1]
        IDSismo=Consulta[0]

    else:
        GMW=np.hstack((GMW,GCMTMagnitud(i,200,200)[1]))
        IDSismo=np.hstack((IDSismo,GCMTMagnitud(i,200,200)[0]))
        
        
BaseAuxiliar["GMW"]=GMW
BaseAuxiliar["IDSismo"]=IDSismo

# BaseAuxiliar.toPandas()
# BaseAuxiliar2=pd.DataFrame(BaseAuxiliar)
# GlobalCMTBase=GlobalCMTBase.rename(columns={"IDSismo":"NOMBRE"}  )




BaseAuxiliar.to_csv(r'/Users/isaias/Desktop/Atenuacion/FinalGCMT.csv',index=False)





###########################################################################


###########################################################################


########################################################################### Analisis





GlobalCMTBase=pd.read_csv('/Users/isaias/Desktop/Atenuacion/BaseGCMT.csv')
GlobalCMTBase=GlobalCMTBase.rename(columns={"NOMBRE": "IDSismo"}, errors="raise")

FinalGCMT=pd.read_csv('/Users/isaias/Desktop/Atenuacion/FinalGCMT.csv')



FinalGCMT[FinalGCMT["GMW"]>0]
Muestra=pd.merge(FinalGCMT, GlobalCMTBase, on=['IDSismo'], how="inner", indicator=True)
Muestra=Muestra.drop(columns=Muestra.columns[-1])
Muestra.columns

plt.scatter(Muestra["Profundidad_x"],Muestra["Profundidad_y"])
plt.xlabel("Profundidad RAII")
plt.ylabel("Profundidad GCMT")


plt.scatter(Muestra["LatSismo"],Muestra["LATITUD"])
plt.scatter(Muestra["LonSismo"],Muestra["LONGITUD"])
plt.xlabel("Profundidad RAII")
plt.ylabel("Profundidad GCMT")


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
ReverseFault=ReverseFault[ReverseFault["Profundidad_y"]>=20]
ReverseFault=ReverseFault[ReverseFault["LONGITUD"]<-95] ##################### Hay que quitar placa del caribe y placa rivera, ReverseFault[ReverseFault["LONGITUD"]<-95]
######################################

from obspy.imaging.beachball import beachball
for i in ReverseFault.index:
    beachball(ReverseFault.loc[i][-9:-3],outfile='/Users/isaias/Desktop/Atenuacion/Playeras/'+ReverseFault["IDSismo"][i]+'.png')

###################################



import plotly.graph_objects as go
import plotly.express as px
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot


plt.hist(ReverseFault["Mw_y"],density=True)
plt.hist(GlobalCMTBase["Mw"],density=True)



fig2 = px.scatter_mapbox(ReverseFault, 
    lat="LATITUD", lon="LONGITUD", hover_name="IDSismo", hover_data=["FECHA","HORA","Mw_y","Profundidad_y"],
    opacity=1, zoom=3, height=900)

fig2.update_layout(mapbox_style="open-street-map")
fig2.update_layout(margin={"r":0,"t":0,"l":0,"b":0})

plot(fig2)

len(np.unique(Muestra["IDSismo"]))



# ReverseFault.to_csv(r'/Users/isaias/Desktop/Atenuacion/SeleccionSismos.csv',index=False)
#####################################

























