#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 17:29:47 2023

@author: isaias
"""

CU=[19.4,-99.184008]
CHIL=[17.502333,-99.509542]
CHIL=[17.0,-99.509542]
# CHIL=[15.22333,-99.509542]
OMET=[16.688475,-98.379861]
PINO=[16.338578,-98.031422]
TE85=[18.0,-104.0]
FM1=[17.4678,-102.8287]
FM2=[15.5169,-97.7521]


m1=(CU[1]-OMET[1])/(CU[0]-OMET[0])
m1=(CU[1]-PINO[1])/(CU[0]-PINO[0])

m2=(CU[1]-TE85[1])/(CU[0]-TE85[0])
m3=(FM1[1]-FM2[1])/(FM1[0]-FM2[0])

###No hay diferencia si me tomo a OMET o a PINO
def AristaDerecha(x):
    return m1*(x-CU[0])+CU[1]

def AristaIzquierda(x):
    return m2*(x-CU[0])+CU[1]

def AristaInferior(x):
    return m3*(x-CHIL[0])+CHIL[1]


# FinalGCMT["LonEstacion"]=-FinalGCMT["LonEstacion"]


# ReverseFault

# ReverseFault2=ReverseFault[ReverseFault["IDSismo"]=="091985B"]
# ReverseFault2=FinalGCMT[FinalGCMT["IDSismo"]=="091985B"]
# ReverseFault2=FinalGCMT[FinalGCMT["LonEstacion"]<AristaDerecha(FinalGCMT["LatEstacion"])]

ReverseFault=ReverseFault[ReverseFault["LonEstacion"]<AristaDerecha(ReverseFault["LatEstacion"])]
ReverseFault=ReverseFault[ReverseFault["LonEstacion"]>AristaInferior(ReverseFault["LatEstacion"])]
ReverseFault=ReverseFault[ReverseFault["LonEstacion"]>AristaIzquierda(ReverseFault["LatEstacion"])]
ReverseFault=ReverseFault[ReverseFault["LonSismo"]<AristaDerecha(ReverseFault["LatSismo"])]
ReverseFault=ReverseFault[ReverseFault["LonSismo"]>AristaIzquierda(ReverseFault["LatSismo"])]

print(len(ReverseFault),len(np.unique(ReverseFault["IDSismo"])),len(np.unique(ReverseFault["IDEstacion"])))

# MapaCompleto(ReverseFault)
MapaCompleto(ReverseFault)
# ReverseFault2.columns

# (ReverseFault[["LATITUD","LONGITUD"]][ReverseFault["IDSismo"]==sismo].iloc[0])
# (ReverseFault[["LatSismo","LonSismo"]][ReverseFault["IDSismo"]==sismo].iloc[0])

ReverseFault.columns
np.min(ReverseFault["Profundidad_x"])
np.min(ReverseFault["Profundidad_y"])

# ReverseFault2=FinalGCMT[FinalGCMT["IDSismo"]=="091985B"]
# fourier(FinalGCMT["Base"][i], 0) 

# ReverseFault2=ReverseFault2[ReverseFault2["LonEstacion"]<AristaDerecha(ReverseFault2["LatEstacion"])]
# ReverseFault2=ReverseFault2[ReverseFault2["LonEstacion"]>AristaInferior(ReverseFault2["LatEstacion"])]
# ReverseFault2=ReverseFault2[ReverseFault2["LonEstacion"]>AristaIzquierda(ReverseFault2["LatEstacion"])]
# ReverseFault2=ReverseFault2[ReverseFault2["LonSismo"]<AristaDerecha(ReverseFault2["LatSismo"])]
# ReverseFault2=ReverseFault2[ReverseFault2["LonSismo"]>AristaIzquierda(ReverseFault2["LatSismo"])]



# (X,NS,f,Y,VelMuest)=fourier(ReverseFault2["Base"][110],2)
# plt.plot(X,NS,linewidth=0.3)



# ReverseFault[ReverseFault["Base"].isin(ReverseFault2["Base"])]
fig3 = go.Figure(go.Scattermapbox(
    mode = "markers+lines",
    marker = {'size': 10}))
fig3.update_layout(mapbox_style="open-street-map")
fig3.update_layout(margin={"r":0,"t":0,"l":0,"b":0})



fig3.add_trace(go.Scattermapbox(
    mode = "markers+lines",
    lon = np.array((CU[1],-98.0725)),
    lat = np.array((CU[0],16.4477)),
    marker = {'size': 10}))


fig3.add_trace(go.Scattermapbox(
    mode = "markers+lines",
    lon = np.array((CU[1],-102.9249)),
    lat = np.array((CU[0],18.3125)),
    marker = {'size': 10}))

fig3.add_trace(go.Scattermapbox(
    mode = "markers+lines",
    lon = np.array((-98.0725,-102.9249)),
    lat = np.array((16.4477,18.3125)),
    marker = {'size': 10}))

plot(fig3)



plot(fig3)


CU=[19.4,-99.184008]
CHIL=[17.502333,-99.509542]
CHIL=[17.0,-99.509542]
# CHIL=[15.22333,-99.509542]
OMET=[16.688475,-98.379861]
PINO=[16.338578,-98.031422]
TE85=[18.0,-104.0]
FM1=[17.4678,-102.8287]
FM2=[15.5169,-97.7521]















































