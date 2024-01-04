#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 16:22:18 2023

@author: isaias
"""


###Correr Funciones, SeleccionFinal
#Graficas utiles en Copula, Descriptivo, EligeSismos y Resumen

###############Parametros iniciales
MW=ReverseFault["GMW"]
Deep=ReverseFault["Profundidad_x"]
R=ReverseFault["R"]
ET=ReverseFault["ET"]
fmax=15 #Si se cambia fmax hay que recalcular todo lo que involucre marginal de frecuencia
Ventana=4 #Si se cambia el tiempo de ventaneo hay que recalcular todo
tmax=1000
VelocidadMuestreo=0.01
VelocidadMuestreoFrec=0.1

sns.color_palette("YlOrRd",as_cmap=True)
Colores=sns.color_palette("YlOrRd", len(ReverseFault))


############################################### Se calculan las marginales
flag=0
for idsismo in ReverseFault.index:
    
    f,t,Z=SePSDF(idsismo,"H1",Ventana)            

    ###########Antes no hacia esto
    # Z=Z[f<fmax,:]    
    # f=f[f<fmax]
    ###########
    Marginalf=np.ones(len(f))
    
    for i in np.arange(len(f)):
        Marginalf[i]=integrate.cumtrapz(Z[i,:],t , initial=0)[-1]
        
    Marginalf=Marginalf/integrate.cumtrapz(Marginalf,f,initial=0)[-1]
    
    Marginalt=np.ones(len(t))

    for i in np.arange(len(t)):
        Marginalt[i]=integrate.cumtrapz(Z[:,i],f , initial=0)[-1]
        

    Marginalt=Marginalt/integrate.cumtrapz(Marginalt,t,initial=0)[-1]



    if flag==0:
        Directorio=[np.vstack((f,Marginalf))]
        Directorio2=[np.vstack((t,Marginalt))]
        flag=1
    else :
        Directorio.append(np.vstack((f,Marginalf)))
        Directorio2.append(np.vstack((t,Marginalt)))
############################################### Saber cuanta energia perdemos con el truncamiento de fmax

# plt.plot(Directorio[0][0],Directorio[0][1])

#############No correr a menos que se cambie el modelo en Marginales o la base de datos
##### Cuando se corra esto no correr la siguiente linea en los directorios
#####     Z=Z[f<fmax,:]    ;f=f[f<fmax]
# ReverseFault["ETFmax"]=0.

# for i in range(len(ReverseFault)):
#     ReverseFault["ETFmax"][ReverseFault.index[i]]=integrate.cumtrapz(Directorio[i][1][Directorio[i][0]<fmax],Directorio[i][0][Directorio[i][0]<fmax],initial=0)[-1]

# plt.hist(ReverseFault["ETFmax"])

############################################### Se ajustan las distribuciones

#############No correr a menos que se cambie el modelo en Marginales o la base de datos

# import multiprocessing as mp
# p = mp.get_context("fork").Pool(12)
# p.map(AjustaMarginalFrec,range(len(ReverseFault)))  
# p.close()    


# p = mp.get_context("fork").Pool(12)
# p.map(AjustaMarginalTiempo,range(len(ReverseFault)))  
# p.close()    



######################################################################################################## Graficas de ajustes

#################### Frecuencia


# for i in range(len(ReverseFault)):
#     # AjustaMarginalFrec(i)
#     Corrida=np.genfromtxt(r'./Corrida/CorridaFrecuencia'+ (ReverseFault["Base"][ReverseFault.index[i]]).replace(".","")  +".csv", delimiter=',')
#     plt.figure()
#     # Efectiva=(Corrida[10000:,:3])[0::500]
#     Efectiva=(Corrida[10000:,:2])[0::500]
#     # plt.plot(-(Corrida[10000:,-1]))
    
#     plt.plot(Directorio[i][0][Directorio[i][0]<fmax],Directorio[i][1][Directorio[i][0]<fmax],alpha=0.5)
#     for Efe in range(len(Efectiva)):
#         # a,b,c=Efectiva[Efe,:]
#         a,b=Efectiva[Efe,:]
#         y=gamma.pdf(Directorio[i][0][Directorio[i][0]<fmax], a,scale=1/b)
#         const=gamma.cdf(10,a,scale=1/b)

#         # y=laplace_asymmetric.pdf(Directorio[i][0][Directorio[i][0]<fmax],loc=a,scale=b,kappa=c)
#         plt.plot(Directorio[i][0][Directorio[i][0]<fmax],y/const,alpha=0.1,color="red")
#     plt.plot(Directorio[i][0][Directorio[i][0]<fmax],Directorio[i][1][Directorio[i][0]<fmax],alpha=0.5)
#     plt.title((ReverseFault["Base"][ReverseFault.index[i]])+" R="+str(np.round(ReverseFault["R"][ReverseFault.index[i]],2)))
#     plt.show()

#################### Tiempo



# i=0
# for i in range(len(ReverseFault)):
#     Corrida=np.genfromtxt(r'./Corrida/CorridaTiempo'+ (ReverseFault["Base"][ReverseFault.index[i]]).replace(".","")  +".csv", delimiter=',')
#     plt.figure()
#     Efectiva=(Corrida[10000:,:3])[0::500]
#     plt.plot(Directorio2[i][0][Directorio2[i][0]<tmax],Directorio2[i][1][Directorio2[i][0]<tmax],alpha=0.5)
#     for Efe in range(len(Efectiva)):
#         a,b,c=Efectiva[Efe,:]
#         y=gamma.pdf(Directorio2[i][0][Directorio2[i][0]<tmax], a,scale=1/b)
#         plt.plot(Directorio2[i][0][Directorio2[i][0]<tmax],y,alpha=0.1,color="red")
#     plt.plot(Directorio2[i][0][Directorio2[i][0]<fmax],Directorio2[i][1][Directorio2[i][0]<fmax],alpha=0.2)
#     plt.title((ReverseFault["Base"][ReverseFault.index[i]])+" R="+str(np.round(ReverseFault["R"][ReverseFault.index[i]],2)))
#     plt.ylabel(r"$\frac{cm^2}{s^4 }$")    
#     plt.xlabel("s")    
#     # plt.savefig('/Users/isaias/Desktop/5Jun/Directorio2/'+(ReverseFault["Base"][ReverseFault.index[i]]).replace(".","")+'.png')
#     plt.show()
    



############################################### Se calculan los momentos de las distribuciones


for i in range(len(ReverseFault)):
    print(i)
    Corrida=np.genfromtxt(r'./Corrida/CorridaFrecuencia'+ (ReverseFault["Base"][ReverseFault.index[i]]).replace(".","")  +".csv", delimiter=',')
    Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]
    # a,b,s2=Theta
    a,b,s2=Theta
    # MuestraGamma=gamma.rvs(size=10000, a=a,scale=1/b)
    
    mu=a/b*gamma.cdf(fmax,a=a+1,scale=1/b)/gamma.cdf(fmax,a=a,scale=1/b)
    m2=(a+1)*a/b**2*gamma.cdf(fmax,a=a+2,scale=1/b)/gamma.cdf(fmax,a=a,scale=1/b)    
    var=m2-mu**2



    if i==0:
        # Momentos=(np.mean(MuestraGamma[MuestraGamma<fmax]),np.var(MuestraGamma[MuestraGamma<fmax]))
        Momentos=(mu,var)
    else :
        Momentos=np.vstack((Momentos,(mu,var)))

#############

for i in range(len(ReverseFault)):
    print(i)
    Corrida=np.genfromtxt(r'./Corrida/CorridaTiempo'+ (ReverseFault["Base"][ReverseFault.index[i]]).replace(".","")  +".csv", delimiter=',')
    Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]
    # a,b,s2=Theta
    a,b,s2=Theta
    if i==0:
        MomentosTiempo=gamma.stats(a=a,scale=1/b,moments='mv')
    else :
        MomentosTiempo=np.vstack((MomentosTiempo,gamma.stats(a=a,scale=1/b,moments='mv')))
        

############# Graficas

plt.scatter(Momentos[:,0],Momentos[:,1])
plt.scatter(Momentos[:,0],np.sqrt(Momentos[:,1]))
plt.scatter(MomentosTiempo[:,0],np.log(MomentosTiempo[:,1]))


############################################### Se calcula el parametro rho de la copula gaussiana con su ePSDF

ReverseFault["rho"]=0.0
from scipy.optimize import minimize_scalar
for idsismo in ReverseFault.index:
    f,t,Z=SePSDF(idsismo,"H1",4)#, Muestra=ReverseFault)         
    Z=Z[f<fmax,:][:-1,:-1]
    f=f[f<fmax][:-1]
    t=t[:-1]    
    
    ########Esto es nuevo
    Marginalf=np.ones(len(f))

    for i in np.arange(len(f)):
        Marginalf[i]=integrate.cumtrapz(Z[i,:],t , initial=0)[-1]
        
    Const=integrate.cumtrapz(Marginalf,f,initial=0)[-1]    

    Z=Z/Const
    ########

    
    #GrafCopula.shape,Z.shape
    Graft,Graff,GrafCopula= GaussianaAjustada(idsismo,-0.2)
    
    def Optimiza(b):
        Graft,Graff,GrafCopula= GaussianaAjustada(idsismo,b)
        return np.sum((Z-GrafCopula)**2)
    Optimizado=minimize_scalar(Optimiza, bounds=(-1, 0), method='bounded')
    print("Optimo","\t",Optimizado["x"],"\n",ReverseFault[["Base","R","GMW"]].loc[idsismo])
    ReverseFault["rho"][idsismo]=Optimizado["x"]

####################################################Falta ley de atenuacion para las rho
# plt.hist(ReverseFault["rho"])
# plt.scatter(ReverseFault["GMW"],ReverseFault["rho"],c=ReverseFault["R"],cmap="YlOrRd" )
# plt.scatter((ReverseFault["R"]),ReverseFault["rho"],c=ReverseFault["GMW"],cmap="YlOrRd" )


# plt.scatter(1/R**(2)*np.exp(0.01*R)*MW,np.log(-ReverseFault["rho"]),c=ReverseFault["GMW"],cmap="YlOrRd" )

# plt.scatter(1/(ReverseFault["R"]),ReverseFault["rho"],c=ReverseFault["GMW"],cmap="YlOrRd" )


# import seaborn as sns
# sns.lmplot('R', 'rho', data=ReverseFault, hue='IDSismo', fit_reg=False)
# sns.lmplot('GMW', 'rho', data=ReverseFault, hue='IDSismo', fit_reg=False)
# plt.show()


######################################################################################################## Prediccion ET
####################################################Falta ley de atenuacion para las rho


def PrediceCopula():
    return np.random.choice(ReverseFault["rho"])

########################################################################################################
########################################################################################################
########################################################################################################









########################################################################################################










########################################################################################################
######################################################################################################## Prediccion ET



# theta=np.array((1,2,3,4))
# def Veros(theta):
#     a,b,c,s2=theta
#     return np.sum(((a*MW*np.exp(-R**c)+b)-np.log10(ReverseFault["ET"]))**2)/(2*s2)+len(MW)/2*np.log(s2)

# def ExpSupp(theta):
#     Valor= all(np.abs(theta)<10)
#     Valor= Valor and (np.abs(theta[2])<0.2)
#     Valor= Valor and (theta[-1]>0)
#     return Valor

# import pytwalk


# MCMC= pytwalk.pytwalk( n=len(theta), U=Veros, Supp=ExpSupp)
# MCMC.Run( T=100000, x0=np.array((1.,1.,0.01,1.)), xp0=np.array((1.,1.,0.01,1.))+0.001)

# Corrida=MCMC.Output

# plt.plot(-Corrida[:,-1])
# # plt.plot(-Corrida[10000:,-1])
# Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]
# E1,E2,E3,Es2=np.copy(Theta)
# def PrediceEnergia(MW,R):
#     return E1*MW*np.exp(-R**E3)+E2



# plt.scatter(E1*MW*np.exp(-R**E3)+E2,np.log10(ReverseFault["ET"]),c=R,cmap="YlOrRd" )
# plt.plot( [np.min(E1*MW*np.exp(-R**E3)+E2),np.max(E1*MW*np.exp(-R**E3)+E2)],[np.min(E1*MW*np.exp(-R**E3)+E2),np.max(E1*MW*np.exp(-R**E3)+E2)] )
# plt.xlabel("E1*MW*np.exp(-R**E3)+E2")
# plt.ylabel("log(ET)")


# ######################################################


# ######################################################################################################## Prediccion Frecuencia



# #################### Media
# np.min(R)
# theta=np.array((1,2,3,5))
# def Veros(theta):
#     f1,f2,f3,fs2=theta
#     return np.sum(((f1*R**f2*np.exp(-R*f3))/(100**f2*np.exp(-100*f3))-Momentos[:,0])**2)/(2*fs2)+len(MW)/2*np.log(fs2)

# def ExpSupp(theta):
#     Valor= all(np.abs(theta)<100)
#     Valor= Valor and (theta[1]<0)
#     Valor= Valor and (theta[1]>-1)
#     Valor= Valor and (theta[2]<0.01)
#     Valor= Valor and (theta[2]>0)
#     Valor= Valor and (theta[-1]>0)
#     return Valor

# import pytwalk


# MCMC= pytwalk.pytwalk( n=len(theta), U=Veros, Supp=ExpSupp)
# MCMC.Run( T=100000, x0=np.array((1.,-.5,0.001,1.)), xp0=np.array((1.,-.5,0.001,1.))+0.001)


# Corrida=MCMC.Output
# Corrida0=np.copy(MCMC.Output)
# plt.plot(-Corrida[:,-1])


# Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]
# f1,f2,f3,fs2=np.copy(Theta)



# def PrediceFrecMedia(R):
#     return f1*R**f2*np.exp(-R*f3)/(100**f2*np.exp(-100*f3))



# Predicho=PrediceFrecMedia(R)
# plt.scatter(Predicho,  Momentos[:,0] , c=R,cmap="YlOrRd" )
# plt.plot([np.min(Predicho),np.max(Predicho)],[np.min(Predicho),np.max(Predicho)])
# plt.xlabel("$f_1*R^f_2*e^{-R*f_3}$")
# plt.ylabel("$\mu_t$")


# plt.scatter(R,Predicho-Momentos[:,0])


# #################### Varianza

# theta=np.array((1,2,3,5))
# def Veros(theta):
#     f1,f2,f3,fs2=theta
#     return np.sum(((f1*R**f2*np.exp(-R*f3))/(100**f2*np.exp(-100*f3))-Momentos[:,1])**2)/(2*fs2)+len(MW)/2*np.log(fs2)

# def ExpSupp(theta):
#     Valor= all(np.abs(theta)<100)
#     Valor= Valor and (theta[1]<0)
#     Valor= Valor and (theta[1]>-1)
#     Valor= Valor and (theta[2]<0.01)
#     Valor= Valor and (theta[2]>0)
#     Valor= Valor and (theta[-1]>0)
#     return Valor

# import pytwalk


# MCMC= pytwalk.pytwalk( n=len(theta), U=Veros, Supp=ExpSupp)
# MCMC.Run( T=100000, x0=np.array((1.,-.5,0.001,1.)), xp0=np.array((1.,-.5,0.001,1.))+0.001)


# Corrida=MCMC.Output
# Corrida1=np.copy(MCMC.Output)
# plt.plot(-Corrida[:,-1])


# Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]
# f21,f22,f23,fs22=np.copy(Theta)

# def PrediceFrecVar(R):
#     return     f21*R**f22*np.exp(-R*f23)/(100**f22*np.exp(-100*f23))


# Predicho=PrediceFrecVar(R)
# plt.scatter(Predicho,  Momentos[:,1] , c=R,cmap="YlOrRd" )
# plt.plot([np.min(Predicho),np.max(Predicho)],[np.min(Predicho),np.max(Predicho)])
# plt.xlabel("$f_1*R^f_2*e^{-R*f_3}$")
# plt.ylabel("$\sigma$")


# plt.scatter(-R,Predicho-Momentos[:,1]) 

# ######################################################################################################## Prediccion Tiempo





# #################### Media

# def Veros(theta):
#     a1,a2,s2=theta
#     return np.sum((a1+a2*np.log(R)-MomentosTiempo[:,0])**2)/(2*s2)+len(MomentosTiempo[:,0])/2*np.log(s2)

# def ExpSupp(theta):
#     Valor= all(np.abs(theta)<100)
#     Valor= Valor and theta[-1]<3
#     # Valor= Valor and (np.abs(theta[2])<0.2)
#     # Valor= Valor and (theta[-1]>0)
#     return Valor

# import pytwalk


# theta=np.array((1,2,0.01))
# MCMC= pytwalk.pytwalk( n=len(theta), U=Veros, Supp=ExpSupp)
# MCMC.Run( T=100000, x0=theta, xp0=theta+0.001)

# Corrida=MCMC.Output
# Corrida2=np.copy(MCMC.Output)

# plt.plot(-Corrida[:,-1])
# # plt.plot(-Corrida[10000:,-1])
# Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]
# a1,a2,as2=np.copy(Theta)

# # def PrediceAmpMedia(MW,R):
# #     return      a1+a6*(a2*MW+(a3+a4*MW)*np.log(np.sqrt(R**2+a5**2)))



# def PrediceAmpMedia(R):
#     return      a1+a2*np.log(R)



# Ajuste= PrediceAmpMedia(R)
# plt.scatter(Ajuste,MomentosTiempo[:,0],alpha=0.3 )
# plt.plot( [np.min(Ajuste),np.max(Ajuste)],[np.min(Ajuste),np.max(Ajuste)] )
# # plt.scatter(Ajuste,MomentosTiempo[:,0],alpha=0.3 , c=R,cmap="YlOrRd" )
# plt.xlabel("$a_1+a_6*(a_2*M_W+(a_3+a_4*M_W)*\log(\sqrt{R^2+a_5^2})$ ")
# plt.ylabel("$\mu_t$ ")


# # plt.scatter(MW,Ajuste-MomentosTiempo[:,0])

# # ReverseFault["Predicho"]=Predicho
# # ReverseFault["Observado"]=Momentos[:,0] 
# # sns.lmplot('Observado', 'Predicho',data=ReverseFault, hue='IDSismo', fit_reg=False)





# #################### Varianza



# def Veros(theta):
#     a1,a2,s2=theta
#     return np.sum((a1+a2*np.log(R)-np.sqrt(MomentosTiempo[:,1]))**2)/(2*s2)+len(MomentosTiempo[:,1])/2*np.log(s2)

# def ExpSupp(theta):
#     Valor= all(np.abs(theta)<100)
#     Valor= Valor and theta[-1]<3
#     # Valor= Valor and (np.abs(theta[2])<0.2)
#     # Valor= Valor and (theta[-1]>0)
#     return Valor

# import pytwalk


# theta=np.array((1,2,0.01))
# MCMC= pytwalk.pytwalk( n=len(theta), U=Veros, Supp=ExpSupp)
# MCMC.Run( T=100000, x0=theta, xp0=theta+0.001)

# Corrida=MCMC.Output
# Corrida3=np.copy(MCMC.Output)

# plt.plot(-Corrida[:,-1])
# # plt.plot(-Corrida[10000:,-1])
# Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]
# a21,a22,as22=np.copy(Theta)

# # def PrediceAmpMedia(MW,R):
# #     return      a1+a6*(a2*MW+(a3+a4*MW)*np.log(np.sqrt(R**2+a5**2)))




# def PrediceAmpSD(R):
#     return      a21+a22*np.log(R)



# Ajuste= PrediceAmpSD(R)
# plt.scatter(Ajuste,np.sqrt(MomentosTiempo[:,1]),alpha=0.3 )
# plt.plot( [np.min(Ajuste),np.max(Ajuste)],[np.min(Ajuste),np.max(Ajuste)] )
# # plt.scatter(Ajuste,MomentosTiempo[:,0],alpha=0.3 , c=R,cmap="YlOrRd" )
# plt.xlabel("$a_1+a_6*(a_2*M_W+(a_3+a_4*M_W)*\log(\sqrt{R^2+a_5^2})$ ")
# plt.ylabel("$\sigma_t$ ")

# plt.scatter(R,Ajuste-MomentosTiempo[:,0])
# plt.scatter(MW,Ajuste-MomentosTiempo[:,0])
# ########################################################################################################




# def Predicciones(r):

#     global CopulaDensidad
#     idsismo=ReverseFault.index[r]
#     (X,NS,f,Y,VelMuest)=fourier(ReverseFault["Base"][idsismo], CanalH1)   
#     NS=NS[X>ReverseFault["Inicio"][idsismo]]
#     X=X[X>ReverseFault["Inicio"][idsismo]]
    
    
#     Mag,Dis,EnTot= ReverseFault[["GMW","R","ET"]].loc[idsismo]
#     Dominiof=np.arange(0,fmax,VelocidadMuestreoFrec)
#     Dominio=np.arange(np.min(Directorio2[r][0,:]),np.max(Directorio2[r][0,:]),VelocidadMuestreo)

    
#     PredET=10**PrediceEnergia(Mag,Dis)
#     # EnTot
    
    
#     PredRho=PrediceCopula() #PrediceCopula(EnTot)
    
#     afp,bfp=abFrecPredicho(Dis)
#     atp,btp=abTimePredicho(Dis)
    
    
    
#     #######################################################
    
    
#     DensMarginalf=gamma.pdf(Dominiof, afp,scale=1/bfp)/gamma.cdf(fmax, afp,scale=1/bfp)###############################Marginalf
#     DensMarginalt=gamma.pdf(Dominio, atp,scale=1/btp)###############################Marginalt
    
    
    
#     DistMarginalt=integrate.cumtrapz(DensMarginalt,Dominio)
#     DistMarginalf=integrate.cumtrapz(DensMarginalf,Dominiof)
    
    
#     ##########################3
    
#     b=PredRho ###############################RhoCopula
    
#     Sigma=np.array([[1,b],[b,1]])
#     # Sigma=np.array([[a,b],[b,c]])
#     SigmaInv=np.linalg.inv(Sigma)
#     SigmaDet=np.linalg.det(Sigma)
    
    
    
#     InvNormalf=sp.stats.norm.ppf(DistMarginalf) 
#     InvNormalt=sp.stats.norm.ppf(DistMarginalt) 
    
    
    
#     def CopulaDensidad(X):
#         i=X[0]
#         j=X[1]
#         Vec=np.array([InvNormalf[i],InvNormalt[j]])
#         # Vec=np.array([DistMarginalf[i],DistMarginalt[j]])
#         Vec.shape=(2,1)
#         Valor=0
        
#         if  all(~np.isinf(Vec)) :
#             Valor= np.exp(-Vec.T@(SigmaInv-np.identity(2))@Vec/2)*SigmaDet**(-1/2)* (DensMarginalt[j] * DensMarginalf[i]) 
#         if np.isinf(Valor) or np.isnan(Valor) :
#             Valor=0            
#         return Valor         
    
#     from scipy.ndimage import gaussian_filter
    
#     from multiprocessing import Pool
#     from multiprocessing import get_context
#     import multiprocessing 
    
#     x=np.arange(len(Dominiof)-1)
#     y=np.arange(len(Dominio)-1)    
#     Malla=pd.DataFrame({'x':np.repeat(x,y.shape[0]),
#                   'y':np.tile(y,x.shape[0])})
#     Malla=Malla.to_numpy()
    
    
#     Hilos=multiprocessing.cpu_count()
#     p = get_context("fork").Pool(Hilos)
#     Valores=p.map(CopulaDensidad,Malla)  
#     p.close()    
    
    
#     Conjunta=np.zeros((len(Dominiof)-1,len(Dominio)-1))
    
#     for i in range(Malla.shape[0]):    
#         loc=Malla[i,:]
#         Conjunta[loc[0],loc[1]]=Valores[i]
    
    
#     Sintetico=NS
#     Registro=np.vstack((X,NS))
#     Registro=np.vstack((Registro,Sintetico))
#     Registro=np.vstack((Registro,Sintetico))
#     Registro.shape
    
#     RegistroChristen=Registro.T
#     np.savetxt(r'./RegistroChristen'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv',RegistroChristen,delimiter=',', fmt=('%s'))
#     tst = Acelerogram( r'./RegistroChristen'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv', Dt=-1, skiprows=0)    
#     Real=np.array(Respuesta(tst))
#     os.remove(r'./RegistroChristen'+(ReverseFault["Base"][idsismo]).replace(".","")+ '.csv')
    
    
    
    
    
#     Trayectorias=20
    
#     Const=PredET###############################Et
#     # Const=EnTot
#     PSDFEst=Conjunta*Const
    
    
    
    
#     Proceso,Simulado=RespuestaSimuladaNueva(PSDFEst,Dominio[:-1],Dominiof[:-1],idsismo)
#     for i in range(Trayectorias):
#         print(i+1)
#         np.random.seed(i+1)
#         Proceso,Nuevo=RespuestaSimuladaNueva(PSDFEst,Dominio[:-1],Dominiof[:-1],idsismo)
#         Simulado=np.vstack((Simulado,Nuevo))    
    
#     Simulado=np.array(Simulado)
    
#     #######################################3
    
    
#     f,t,Z=SePSDF(idsismo,"H1",Ventana)       




#     fig, axs = plt.subplots(3, 2,figsize=(10,6))                      
#     axs[0,0].pcolormesh(t,f[f<fmax],Z[f<fmax,:], shading='gouraud', rasterized=True)
#     axs[0,0].set_xlabel("s")
#     axs[0,0].set_ylabel("Hz")
#     # axs[0,0].set_colorbar()
#     # Suave=gaussian_filter(np.log10(Conjunta+np.finfo(float).eps), 1)
#     # Suave=gaussian_filter((Conjunta+np.finfo(float).eps), 1)
#     # plt.contour(Dominio[:-1], Dominiof[:-1], Suave, 10, colors='white',origin="lower",alpha=0.2)
#     axs[0,0].contour(Dominio[:-1], Dominiof[:-1], Conjunta, 10, colors='white',origin="lower",alpha=0.2)
    
    
#     axs[2,0].plot(Real[0],Real[1])
#     for i in range(Trayectorias+1):
#         axs[2,0].plot(Real[0],Simulado[i],color="red",alpha=0.2)
#     axs[2,0].set_xlabel("T")
#     axs[2,0].set_ylabel(r"$\frac{cm}{s^2}$")
    
    
#     y=gamma.pdf(Directorio[r][0], afp,scale=1/bfp)
#     const=gamma.cdf(fmax, afp,scale=1/bfp)
#     axs[0,1].plot(Directorio[r][0],Directorio[r][1])
#     axs[0,1].plot(Directorio[r][0],y/const,color="red")
#     axs[0,1].axvline(x = fmax, color = 'b', label = 'axvline - full height')
#     axs[0,1].set_xlabel("Hz")
#     axs[0,1].set_ylabel("Density")    
    
#     y2=gamma.pdf(Directorio2[r][0], atp,scale=1/btp)
#     axs[1,0].plot(Directorio2[r][0],Directorio2[r][1])
#     axs[1,0].plot(Directorio2[r][0],y2,color="red")
#     axs[1,0].set_xlabel("s")
#     axs[1,0].set_ylabel("Density")

#     axs[1,1].plot(X,NS)
#     axs[1,1].plot(Directorio2[r][0]+ReverseFault["Inicio"].loc[idsismo],np.sqrt(y2*PredET),color="red")
#     axs[1,1].plot(Directorio2[r][0]+ReverseFault["Inicio"].loc[idsismo],-np.sqrt(y2*PredET),color="red")
#     axs[1,1].set_xlabel("s")
#     axs[1,1].set_ylabel(r"$\frac{cm}{s^2}$")
    
#     axs[2,1].plot(X,NS)
#     axs[2,1].plot(Dominio[:-1]+ReverseFault["Inicio"].loc[idsismo],Proceso,alpha=0.5,color="red")
#     axs[2,1].set_xlabel("s")
#     axs[2,1].set_ylabel(r"$\frac{cm}{s^2}$")
    
    
    
#     fig.suptitle(ReverseFault["Base"][idsismo]+' Mw='+str(ReverseFault["GMW"][idsismo])+" R="+str(round(ReverseFault["R"][idsismo],2)), fontsize=16)
#     fig.tight_layout()


#     fig.savefig(r"/Users/isaias/Desktop/Atenuacion/Prediccion/"+(ReverseFault["Base"][idsismo]).replace(".","")+'.png',dpi=200, bbox_inches = "tight")
    
# ###########################



        
## for r in range(len(ReverseFault)):         
##     Predicciones(r)








































































































































































































        