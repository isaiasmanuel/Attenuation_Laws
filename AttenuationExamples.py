#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 14:27:23 2023

@author: isaias
"""
#Correr Funciones, SeleccionFinal, PrediccionFinal antes

# ReverseFault

# ReverseFault.to_csv('./ReverseFaultCompleta2.csv')
# np.savetxt('./MomentosCompleta2.csv', Momentos, delimiter=",")
# np.savetxt('./MomentosTiempoCompleta2.csv', MomentosTiempo, delimiter=",")

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
######################


ReverseFault=pd.read_csv('./ReverseFaultCompleta2.csv', index_col=0)
Momentos=np.loadtxt('./MomentosCompleta2.csv',delimiter=",")
MomentosTiempo=np.loadtxt('./MomentosTiempoCompleta2.csv',delimiter=",")

Trayectorias=0


# plt.scatter(-np.log10(Momentos[:,0]),MomentosTiempo[:,0])
# plt.scatter(-np.log10(Momentos[:,0]),MomentosTiempo[:,0])
# plt.scatter(-np.log10(Momentos[:,0]),np.sqrt(MomentosTiempo[:,1]))

Obs=(np.vstack((Momentos[:,0],np.sqrt(Momentos[:,1])))).T
ObsTiempo=(np.vstack((MomentosTiempo[:,0],np.sqrt(MomentosTiempo[:,1])))).T
Obs=np.delete(Obs,73,0)  #Quite un dato atipico revisar despues
ObsTiempo=np.delete(ObsTiempo,73,0)  #Quite un dato atipico revisar despues
ObsEnergia=np.log10(np.delete(np.array(ReverseFault["ET"]),73,0))  #Quite un dato atipico revisar despues

#######
R2=ReverseFault[ReverseFault.index!=ReverseFault.index[73]]["R"]
MW2=ReverseFault[ReverseFault.index!=ReverseFault.index[73]]["GMW"]
#######
Obs.shape


# plt.scatter(Obs[:,0],Obs[:,1])
# Texto=3
# for i in range(len(ReverseFault)):
#     plt.annotate(str(ReverseFault.index[i])+","+str(i), (Obs[i,0], Obs[i,1]))

#############################
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
        
##################        


# ReverseFault.loc[464]

plt.scatter(Obs[:,0],Obs[:,1])
plt.xlabel("$\mu_2$")
plt.ylabel("$\sigma_2$")


plt.scatter(ObsTiempo[:,0],ObsTiempo[:,1])
plt.xlabel("$\mu_1$")
plt.ylabel("$\sigma_1$")


# plt.scatter(-np.log(Obs[:,0]),ObsTiempo[:,0])

# plt.scatter(ObsTiempo[:,0],np.sqrt(ObsTiempo[:,1]))

# plt.hist(ReverseFault["rho"])


# theta=(f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,s11,s12,s21,s22)
# (f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,s11,s12,s21,s22)=theta


def Predicemf(f11,f12,f13,f14,R):
    return f11*R**(f12)*np.exp(-f13*R)+f14
    
def Predicesf(f21,f22,mf):
    return f21*mf+f22


def Predicemt(t11,t12,R):
    return t11*np.log(R)+t12

# def Predicemt(t11,t12,mf):
#     return t11*np.log(mf)+t12

def Predicest(t21,t22,mt):
    return t21*mt+t22

def PrediceEnergia(E1,E2,E3,MW,R):
    return E1*MW*np.exp(-R**E3)+E2

def PrediceCopula():
    return np.random.choice(ReverseFault["rho"])













def Veros(theta):
    (f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,E1,E2,E3,s11,s12,s21,s22,s00)=theta    
    t1=np.sum(-np.log(sp.stats.norm.cdf((Obs[:,0])/s12)-sp.stats.norm.cdf(-Predicesf(f21,f22,Obs[:,0])/s12))) #sf 
    t2=np.sum(-np.log(sp.stats.norm.cdf((ObsTiempo[:,0])/s22)-sp.stats.norm.cdf(-Predicest(t21,t22,ObsTiempo[:,0])/s22))) #st
    c1=np.sum((Predicemf(f11,f12,f13,f14,R2)-Obs[:,0])**2)/(2*s11)+n/2*np.log(s11)  ####muf
    # c1+=np.sum((Predicesf(f21,f22,Obs[:,0])-Obs[:,1])**2)/(2*s12)+n/2*np.log(s12)   ####sigmaf

    c1+=np.sum((Predicesf(f21,f22,Obs[:,0])-Obs[:,1])**2)/(2*s12)+n/2*np.log(s12)+t1    ####sigmaf

    c1+=np.sum((Predicemt(t11,t12,R2)-ObsTiempo[:,0])**2)/(2*s21)+n/2*np.log(s21) ###mut
    # c1+=np.sum((Predicest(t21,t22,ObsTiempo[:,0])-ObsTiempo[:,1])**2)/(2*s22)+n/2*np.log(s22) ###sigmat
    c1+=np.sum((Predicest(t21,t22,ObsTiempo[:,0])-ObsTiempo[:,1])**2)/(2*s22)+n/2*np.log(s22)+t2 ###sigmat
    c1+=np.sum( (PrediceEnergia(E1,E2,E3,MW2,R2)-ObsEnergia)**2)/(2*s00)+n/2*np.log(s00)  ###ET
    return c1

def ExpSupp(theta):
    (f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,E1,E2,E3,s11,s12,s21,s22,s00)=theta
    Valor = True
    Valor = Valor and all(np.array((s11,s12,s21,s22,s00))>0)
    # Valor = Valor and all(np.array((s11,s12,s21,s22))<20)    

    Valor= Valor and all(np.abs(theta)<100)

    Valor= Valor and (np.abs(E3)<0.2)
    Valor= Valor and (f12<0)
    Valor= Valor and (f12>-1)
    Valor= Valor and (f13<0.01)
    Valor= Valor and (f13>0)
    Valor= Valor and (f21>0)
    Valor= Valor and (t11>0)
    # Valor= Valor and (t11<0)
    Valor= Valor and (t21>0)

   # Valor= Valor and all(theta[-2:]<5)
    return Valor



(f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,E1,E2,E3,s11,s12,s21,s22,s00)=np.ones(18)
n=len(Obs)
f12=-0.5
f13=0.005
E3=0.1
# t11=-1
theta=(f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,E1,E2,E3,s11,s12,s21,s22,s00)
theta=np.array(theta)

np.random.seed(1)
MCMC= pytwalk.pytwalk( n=len(theta), U=Veros, Supp=ExpSupp)
MCMC.Run( T=700000, x0=theta, xp0=theta+0.00001)


MCMC.Ana()
Corrida=MCMC.Output
Corrida.shape


Corrida=np.loadtxt(r'./'+ "Posterior" +".csv",delimiter=',')



from statsmodels.graphics.tsaplots import plot_acf
plt.plot(-Corrida[:,-1])

len(Corrida)
# plot_acf(Corrida[1000:,0], maxlags = 1000)

burnin=10000 #10000
lag=30000 #30000
plot_acf(Corrida[burnin:,0][0::lag,], lags = 10)

CorridaEfectiva=Corrida[burnin::lag,]

CorridaEfectiva.shape
plt.plot(-CorridaEfectiva[:,-1])

# No correr np.savetxt(r'./'+ "Posterior" +".csv",Corrida,delimiter=',', fmt=('%s'))
Corrida.shape
CorridaEfectiva.shape
plt.plot(-Corrida[:,-1])
plt.plot(-CorridaEfectiva[:,-1])


Theta=Corrida[np.argsort(Corrida[:,-1])[0],:-1]

(f11o,f12o,f13o,f14o,f21o,f22o,t11o,t12o,t21o,t22o,E1o,E2o,E3o,s11o,s12o,s21o,s22o,s00o)=Theta

Predichomf=Predicemf(f11o,f12o,f13o,f14o,R2)
Predichosf=Predicesf(f21o,f22o,Obs[:,0])
# Predichomt=Predicemt(t11,t12,Obs[:,0])
Predichomt=Predicemt(t11o,t12o,R2)
Predichost=Predicest(t21o,t22o,Predichomt)
Predichost=Predicest(t21o,t22o,ObsTiempo[:,0])
PredichoET=PrediceEnergia(E1o,E2o,E3o,MW2,R2)


plt.scatter(PredichoET,ObsEnergia)
plt.plot([np.min(PredichoET),np.max(PredichoET)],[np.min(PredichoET),np.max(PredichoET)], color="green",label="f(x)=x")
plt.xlabel(r"$E_1 M_{w,i} e^{-R_i^{E_3}}+E_2$")
plt.ylabel(r"$\log_{10}(E_{T,i})$")
plt.grid(color='gray', linestyle='dashed')


plt.scatter(Predichomf,Obs[:,0])
plt.plot([np.min(Predichomf),np.max(Predichomf)],[np.min(Predichomf),np.max(Predichomf)], color="green")
plt.xlabel(r"$f_{11}R_{i}^{f_{12}}e^{-f_{13}R_{i}}+f_{14}$")
plt.ylabel("$\mu_{1,i}$")
plt.grid(color='gray', linestyle='dashed')

# plt.scatter(Obs[:,0],Predichomf)

plt.scatter(Predichosf,Obs[:,1])
plt.plot([np.min(Predichosf),np.max(Predichosf)],[np.min(Predichosf),np.max(Predichosf)], color="green")
plt.xlabel(r"$f_{21}\hat{\mu}_{1,i}+f_{22}$")
plt.ylabel("$\sigma_{1,i}$")
plt.grid(color='gray', linestyle='dashed')


plt.scatter(Predichomt,ObsTiempo[:,0])
plt.plot([np.min(Predichomt),np.max(Predichomt)],[np.min(Predichomt),np.max(Predichomt)], color="green")
plt.xlabel(r"$t_{11}\log(R_i)+t_{12}$")
plt.ylabel("$\mu_{2,i}$")
plt.grid(color='gray', linestyle='dashed')


plt.scatter(Predichost,ObsTiempo[:,1])
plt.plot([np.min(Predichost),np.max(Predichost)],[np.min(Predichost),np.max(Predichost)], color="green")
plt.xlabel(r"$t_{21}\hat{\mu}_{2,i}+t_{22}$")
plt.ylabel("$\sigma_{2,i}$")
plt.grid(color='gray', linestyle='dashed')  


plt.hist(ReverseFault["rho"])
plt.xlabel(r"$\rho$")


plt.scatter(np.log(R2), ObsTiempo[:,0])





#################################################################################################

pru=20
pre=5
gamma.cdf(fmax, a=pru+1,scale=1/pre)/gamma.cdf(fmax, a=pru,scale=1/pre)
gamma.cdf(fmax, a=pru+2,scale=1/pre)/gamma.cdf(fmax, a=pru,scale=1/pre)



#################################################################################################

# def abFrecPredicho2(X):
#     muF,s2F=X
    
#     def Objetivo(X):
#         alpha,beta=X
#         mu=(alpha/beta)*gamma.cdf(fmax, a=alpha+1,scale=1/beta)/gamma.cdf(fmax, a=alpha,scale=1/beta)
#         mu2=((alpha+1)*alpha/beta**2)*gamma.cdf(fmax, a=alpha+2,scale=1/beta)/gamma.cdf(fmax, a=alpha,scale=1/beta)
#         sigma2=mu2-mu**2
#         return (muF-mu)**2+(s2F-sigma2)**2
    
    
#     Minimiza=minimize(Objetivo, [10,10], method='Nelder-Mead', tol=1e-6, bounds=((1, None), (0, None)))
    
#     alpha,beta=Minimiza["x"]
#     return alpha,beta



def abTimePredicho2(X):
    mu,s2=X
  
    alpha=mu**2/s2
    beta=mu/s2  
    return alpha,beta



# def abTimePredicho2(X):
#     mu,s2=X
#     def Objetivo(X):
#         alpha,beta=X
#         muE=alpha/beta
#         s2E=alpha/beta**2
#         return (muE-mu)**2+(s2E-s2)**2
#     Minimiza=minimize(Objetivo, [2,2], method='Nelder-Mead', tol=1e-6, bounds=((1, None), (0, None)))
#     alpha,beta=Minimiza["x"]
#     return alpha,beta
    

# plt.plot(CorridaEfectiva[:,3])
# CorridaEfectiva=Corrida[10000::10000,]
# CorridaEfectiva.shape

np.random.seed(1)

for i in range(len(R2)):
    plt.plot()
    for j in range(len(CorridaEfectiva)):
        print(j/len(CorridaEfectiva))
        # f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,s11,s12,s21,s22=CorridaEfectiva[j,:-1]
        (f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,E1,E2,E3,s11,s12,s21,s22,s00)=CorridaEfectiva[j,:-1]
        #####################
        # muf=Predicemf(f11,f12,f13,f14,R2[R2.index[i]])+np.random.normal(size=1,loc=0, scale=np.sqrt(s11))
        

        # sf=Predicesf(f21,f22,muf)+np.random.normal(size=1,loc=0, scale=np.sqrt(s12))
        loc=Predicemf(f11,f12,f13,f14,R2[R2.index[i]])
        scale=np.sqrt(s11)
        myclip_a=0
        myclip_b=1000
        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
        muf=sp.stats.truncnorm.rvs(loc=loc,scale=scale,a=a,b=b)

        
        
        
        
        # sf=Predicesf(f21,f22,muf)+np.random.normal(size=1,loc=0, scale=np.sqrt(s12))
        loc=Predicesf(f21,f22,muf)
        scale=np.sqrt(s12)
        myclip_a=0
        myclip_b=np.sqrt(muf)
        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
        sf=sp.stats.truncnorm.rvs(loc=loc,scale=scale,a=a,b=b)

        #####################
        # mut=Predicemt(t11,t12,R2[R2.index[i]])+np.random.normal(size=1,loc=0, scale=np.sqrt(s21))
        
        loc=Predicemt(t11,t12,R2[R2.index[i]])
        scale=np.sqrt(s21)
        myclip_a=0
        myclip_b=1000
        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
        mut=sp.stats.truncnorm.rvs(loc=loc,scale=scale,a=a,b=b)


        
        
        
        
        # st=Predicest(t21,t22,mut)+np.random.normal(size=1,loc=0, scale=np.sqrt(s22))
        loc=Predicest(t21,t22,mut)
        scale=np.sqrt(s22)
        myclip_a=0
        myclip_b=np.sqrt(mut)
        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
        st=sp.stats.truncnorm.rvs(loc=loc,scale=scale,a=a,b=b)
        #####################
        
        
        # EstET=10**PrediceEnergia(E1,E2,E3,MW2[R2.index[i]],R2[R2.index[i]])
        loc=PrediceEnergia(E1,E2,E3,MW2[R2.index[i]],R2[R2.index[i]])
        scale=np.sqrt(s00)
        myclip_a=0
        myclip_b=10**20
        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
        EstET=10**sp.stats.truncnorm.rvs(loc=loc,scale=scale,a=a,b=b)

        
        EstC=PrediceCopula()
        


        
        #####################
        # EstFrec=abFrecPredicho2([muf,sf**2])
        EstFrec=abTimePredicho2([muf,sf**2])
        EstTiem=abTimePredicho2([mut,st**2])
        print(EstFrec)
        print(EstTiem)


        # Sigma=[[s11,0],[0,s22]]
        # g1=f11*R[ReverseFault.index[i]]**f12*np.exp(-R[ReverseFault.index[i]]*f13)+f14
        # g2=f21*g1+f22
        # Media=(np.vstack((g1,g2))).T
        if j==0:
            SimuladosFrec=EstFrec
            SimuladosTiem=EstTiem
            SimuladosET=EstET
            SimuladosC=EstC
            
        else :

            SimuladosFrec=np.vstack((SimuladosFrec,EstFrec))
            SimuladosTiem=np.vstack((SimuladosTiem,EstTiem))
            SimuladosET=np.vstack((SimuladosET,EstET))
            SimuladosC=np.vstack((SimuladosC,EstC))
    
    if i==0:
        DirectorioSimulaciones=[ np.hstack((SimuladosFrec,SimuladosTiem,SimuladosET,SimuladosC))]
        
    else :
        DirectorioSimulaciones.append([ np.hstack((SimuladosFrec,SimuladosTiem,SimuladosET,SimuladosC))])

# np.savetxt(r'./'+ "DirectorioSimulaciones" +".csv",DirectorioSimulaciones,delimiter=',', fmt=('%s'))
        

DirectorioSimulaciones[0].shape
DirectorioSimulaciones2[0].shape
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

###############Zona de pruebas
# np.save('DirectorioSimulaciones2.npy', np.array(DirectorioSimulaciones, dtype=object), allow_pickle=True)
# DirectorioSimulaciones2 = np.load('./DirectorioSimulaciones2.npy', allow_pickle=True)    

##############Bueb gecgi
# np.save('DirectorioSimulaciones.npy', np.array(DirectorioSimulaciones, dtype=object), allow_pickle=True)
DirectorioSimulaciones2 = np.load('./DirectorioSimulaciones.npy', allow_pickle=True)    
DirectorioSimulaciones2.shape

i=0
np.random.seed(1)
for DS in DirectorioSimulaciones2:    
    print(DS)
    if i==0:
        SimuladosFrec=DS[:,0:2]
        SimuladosTiem=DS[:,2:4]
        SimuladosET=DS[:,4]
        SimuladosC=DS[:,5]
    if i!=0:
        SimuladosFrec=DS[0][:,0:2]
        SimuladosTiem=DS[0][:,2:4]
        SimuladosET=DS[0][:,4]
        SimuladosC=DS[0][:,5]
        
    
    plt.plot()
    # try:
    #     sns.kdeplot(SimuladosFrec[:,0], SimuladosFrec[:,1],fill=True,alpha=0.3,color="red")
    # except: 
    #     print("Error")

    # p1=np.median(SimuladosFrec[:,0])
    # p2=np.median(SimuladosFrec[:,1])
    # plt.scatter(p1,p2,label="Estimado: "+ (str(   np.round(p1,2) ) )+ ","+(str(   np.round(p2,2) ) )    )
    # plt.scatter(Obs[i,0],Obs[i,1],label="Observado: "+ (str(   np.round(Obs[i,0],2) ) )+ ","+(str(   np.round(Obs[i,1],2) ) )    )
    # plt.title(ReverseFault["Base"][R2.index[i]]+", R="+str(np.round(ReverseFault["R"][R2.index[i]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[i]],2) ))
    # plt.legend()
    # plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"Frec"+'.png')
    # plt.show()
    
    
    
    # plt.plot()
    # sns.kdeplot(SimuladosTiem[:,0], SimuladosTiem[:,1],fill=True,alpha=0.3,color="red")
    # p1=np.median(SimuladosTiem[:,0])
    # p2=np.median(SimuladosTiem[:,1])
    # plt.scatter(p1,p2,label="Estimado: "+ (str(   np.round(p1,2) ) )+ ","+(str(   np.round(p2,2) ) )    )
    # plt.scatter(ObsTiempo[i,0],ObsTiempo[i,1],label="Observado: "+ (str(   np.round(ObsTiempo[i,0],2) ) )+ ","+(str(   np.round(ObsTiempo[i,1],2) ) )    )
    # plt.title(ReverseFault["Base"][R2.index[i]]+", R="+str(np.round(ReverseFault["R"][R2.index[i]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[i]],2) ))
    # plt.legend()
    # plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"Tiempo"+'.png')
    # plt.show()
    
    plt.plot()
    idsismo=R2.index[i]
    (X,NS,f,Y,VelMuest)=fourier(ReverseFault["Base"][idsismo], CanalH1)   
    NS=NS[X>ReverseFault["Inicio"][idsismo]]
    X=X[X>ReverseFault["Inicio"][idsismo]]
    Mag,Dis,EnTot= ReverseFault[["GMW","R","ET"]].loc[idsismo]
    Dominiot=Directorio2[i][0]
    
    flag=0
    for r in np.random.choice(len(SimuladosFrec),size=100): 
        # afp,bfp=abFrecPredicho2(Simulados[r,:])
        afp,bfp=SimuladosTiem[r,:]
        if afp>0 and bfp>0 :
            # print(afp,bfp)
            DensMarginalf=gamma.pdf(Dominiot, afp,scale=1/bfp)###############################Marginalt
            if  np.max(DensMarginalf)<5:
                if flag==0:
                    plt.plot(Dominiot,DensMarginalf,color="orange",alpha=0.1,label="Simulated")
                    flag=1
                else:
                    plt.plot(Dominiot,DensMarginalf,color="orange",alpha=0.1)

                
        else :
            print("Problema con negativo")
    plt.plot(Directorio2[i][0],Directorio2[i][1],label="Observado")   
    mufo=Predicemt(t11o,t12o,R2[R2.index[i]])
    sfo=Predicest(t21o,t22o,mufo)
    afp,bfp=abTimePredicho2(((mufo,sfo**2)))    
    DensMarginalf=gamma.pdf(Dominiot, afp,scale=1/bfp)###############################Marginalf
    plt.plot(Dominiot,DensMarginalf,color="red",label="MAP")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.title(ReverseFault["Base"][R2.index[i]]+", R="+str(np.round(ReverseFault["R"][R2.index[i]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[i]],2) ))
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"Curvas"+"Tiempo"+'.png')
    plt.show()
    


        
    plt.plot()
    idsismo=R2.index[i]
    (X,NS,f,Y,VelMuest)=fourier(ReverseFault["Base"][idsismo], CanalH1)   
    NS=NS[X>ReverseFault["Inicio"][idsismo]]
    X=X[X>ReverseFault["Inicio"][idsismo]]
    Mag,Dis,EnTot= ReverseFault[["GMW","R","ET"]].loc[idsismo]
    Dominiof=np.arange(0,fmax,VelocidadMuestreoFrec)
    
    flag=0
    for r in np.random.choice(len(SimuladosFrec),size=100): 
        # afp,bfp=abFrecPredicho2(Simulados[r,:])
        afp,bfp=SimuladosFrec[r,:]
        if afp>0 and bfp>0 :
            # print(afp,bfp)
            DensMarginalf=gamma.pdf(Dominiof, afp,scale=1/bfp)/gamma.cdf(fmax, afp,scale=1/bfp)###############################Marginalf
            if  np.max(DensMarginalf)<5:
                if flag==0:
                    plt.plot(Dominiof,DensMarginalf,color="orange",alpha=0.2,label="Simulated")
                    flag=1
                else: 
                    plt.plot(Dominiof,DensMarginalf,color="orange",alpha=0.2)
        else :
            print("Problema con negativo")
    plt.plot(Directorio[i][0][Directorio[i][0]<fmax],Directorio[i][1][Directorio[i][0]<fmax],label="Observado")        
    muto=Predicemf(f11o,f12o,f13o,f14o,R2[R2.index[i]])
    sto=Predicesf(f21o,f22o,muto)
    # atp,btp=abFrecPredicho2(((muto,sto**2)))    
    atp,btp=abTimePredicho2(((muto,sto**2)))    

    DensMarginalf=gamma.pdf(Dominiof, atp,scale=1/btp)/gamma.cdf(fmax, atp,scale=1/btp)###############################Marginalf
    plt.plot(Dominiof,DensMarginalf,color="red",label="MAP")
    plt.xlabel("Frequency [Hz]")
    plt.legend()
    plt.title(ReverseFault["Base"][R2.index[i]]+", R="+str(np.round(ReverseFault["R"][R2.index[i]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[i]],2) ))
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"Curvas"+'.png')
    plt.show()
    
    i+=1
    
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

#### T_list, definida en funciones
Trayectorias=0
i=0

for DS in DirectorioSimulaciones2:    
    print("Numero de registro",i/len(ReverseFault))
    if i==0:
        SimuladosFrec=DS[:,0:2]
        SimuladosTiem=DS[:,2:4]
        SimuladosET=DS[:,4]
        SimuladosC=DS[:,5]
    if i!=0:
        SimuladosFrec=DS[0][:,0:2]
        SimuladosTiem=DS[0][:,2:4]
        SimuladosET=DS[0][:,4]
        SimuladosC=DS[0][:,5]
         
    for j in range(SimuladosFrec.shape[0]):
    # for j in range(50):
        print("Progreso registro", j/SimuladosFrec.shape[0])
        [Proceso,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]=Predicciones2(i,SimuladosFrec[j,0],SimuladosFrec[j,1],SimuladosTiem[j,0],SimuladosTiem[j,1],SimuladosET[j],SimuladosC[j])
        
        if j==0:
            RSSimulado=Simulado
        if j!=0:
            RSSimulado=np.vstack((RSSimulado,Simulado))
    
    if i==0:
        RSSimuladoDirectorio=[RSSimulado]
    if i!=0:
        RSSimuladoDirectorio.append(RSSimulado)
        
    i+=1    
    

# np.save('DirectorioRS.npy', np.array(RSSimuladoDirectorio, dtype=object), allow_pickle=True)
DirectorioRS = np.load('DirectorioRS.npy', allow_pickle=True)    
 
DirectorioRS.shape
    
    
    
######################
######################
######################
######################
# j=2
# DirectorioRS[2]==0
# np.apply_along_axis(all, 1, DirectorioRS[2]==0)
q2=0.2
np.random.seed(1)
for j in range(len(R2)):
    plt.plot()
    Aux= np.apply_along_axis(all, 1, DirectorioRS[j]==0)
    
    for k in range(len(T_list)):
        if k==0:
            inf=np.quantile(DirectorioRS[j][~Aux][:,k],q=q2/2)
            sup=np.quantile(DirectorioRS[j][~Aux][:,k],q=1-q2/2)
            med=np.quantile(DirectorioRS[j][~Aux][:,k],q=0.5)
        else :
            inf=np.vstack((inf, np.quantile(DirectorioRS[j][~Aux][:,k],q=q2/2)))
            sup=np.vstack((sup, np.quantile(DirectorioRS[j][~Aux][:,k],q=1-q2/2)))
            med=np.vstack((med, np.quantile(DirectorioRS[j][~Aux][:,k],q=0.5)))
            
    
    [Proceso,Real,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]=Predicciones(j,1,1,1,1,10**1,PrediceCopula())
    plt.plot(Real[0],Real[1],color="blue", label="$S_a$ observed")
    plt.fill_between(T_list, np.reshape(inf, len(inf)), np.reshape(sup,len(sup)), color='orange', alpha=.35, label="80% Interval")            
    plt.plot(T_list,med, color="red",alpha=1,label="Median")

    # plt.plot(T_list,inf, color="red",alpha=1, label="Quantiles of posterior predictive")
    # plt.plot(T_list,sup, color="red",alpha=1)
    
        
    plt.xlabel("Time [s]")
    plt.ylabel("Acceleration "+r"$\frac{cm}{s^2}$")
    plt.legend()
    plt.title(ReverseFault["Base"][R2.index[j]]+", R="+str(np.round(ReverseFault["R"][R2.index[j]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[j]],2) ))
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[j]]+"RS"+'.png')

    plt.show()
    


# np.apply_along_axis(np.quantile, 1, DirectorioRS[j][~Aux])


for j in range(len(R2)):
    plt.plot()
    Aux= np.apply_along_axis(all, 1, DirectorioRS[j]==0)
        
    for i in range(np.sum(~Aux)):
        plt.scatter(T_list,DirectorioRS[j][~Aux][i], color="red",alpha=0.1)
    
    
    [Proceso,Real,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]=Predicciones(j,afp,bfp,atp,btp,10**1,PrediceCopula())
    plt.plot(Real[0],Real[1],color="blue")
    plt.show()
    
#######################################################

def Predicciones2(r,afp,bfp,atp,btp,PredET,PredRho):

    global CopulaDensidad
    idsismo=R2.index[r]
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

    return [Proceso,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]

    
    
    
    
    




#######################################################
    
    
######################
    
    plt.plot()
    
    
    ###############    
    for j in range(20): # range(len(CorridaEfectiva))
        [Proceso,Real,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]=Predicciones(i,SimuladosFrec[j,0],SimuladosFrec[j,1],SimuladosTiem[j,0],SimuladosTiem[j,1],SimuladosET[j][0],SimuladosC[j][0])
    
        for s in range(Trayectorias+1):
            plt.plot(Real[0],Simulado[s],color="red",alpha=0.1)
        plt.xlabel("T")
        plt.ylabel(r"$\frac{cm}{s^2}$")
        # plt.set_xlabel("T")
        # plt.set_ylabel(r"$\frac{cm}{s^2}$")
        # plt.show(block=False)
        
    ETo=PrediceEnergia(E1o, E2o, E3o, MW2[R2.index[i]], R2[R2.index[i]])

    [Proceso,Real,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]=Predicciones(i,afp,bfp,atp,btp,10**ETo,PrediceCopula())
    plt.plot(Real[0],Real[1],color="blue")
    for s in range(Trayectorias+1):
        plt.plot(Real[0],Simulado[s],color="green",alpha=0.5)
        
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"RS"+'.png')
    plt.show()
    
    




########################################################################################################################################################


np.random.seed(1)
for i in range(len(R2)):
    plt.plot()
    for j in range(len(CorridaEfectiva)):
        print(j/len(CorridaEfectiva))
        # f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,s11,s12,s21,s22=CorridaEfectiva[j,:-1]
        (f11,f12,f13,f14,f21,f22,t11,t12,t21,t22,E1,E2,E3,s11,s12,s21,s22,s00)=CorridaEfectiva[j,:-1]
        #####################
        muf=Predicemf(f11,f12,f13,f14,R2[R2.index[i]])+np.random.normal(size=1,loc=0, scale=np.sqrt(s11))
        # sf=Predicesf(f21,f22,muf)+np.random.normal(size=1,loc=0, scale=np.sqrt(s12))
        loc=Predicesf(f21,f22,muf)
        scale=np.sqrt(s12)
        myclip_a=0
        myclip_b=1000
        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
        sf=sp.stats.truncnorm.rvs(loc=loc,scale=scale,a=a,b=b)

        #####################
        mut=Predicemt(t11,t12,R2[R2.index[i]])+np.random.normal(size=1,loc=0, scale=np.sqrt(s21))
        # st=Predicest(t21,t22,mut)+np.random.normal(size=1,loc=0, scale=np.sqrt(s22))
        loc=Predicest(t21,t22,mut)
        scale=np.sqrt(s22)
        myclip_a=0
        myclip_b=1000
        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
        st=sp.stats.truncnorm.rvs(loc=loc,scale=scale,a=a,b=b)
        #####################
        
        
        # EstET=10**PrediceEnergia(E1,E2,E3,MW2[R2.index[i]],R2[R2.index[i]])
        loc=PrediceEnergia(E1,E2,E3,MW2[R2.index[i]],R2[R2.index[i]])
        scale=np.sqrt(s00)
        myclip_a=0
        myclip_b=10**20
        a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
        EstET=10**sp.stats.truncnorm.rvs(loc=loc,scale=scale,a=a,b=b)

        
        EstC=PrediceCopula()
        


        
        #####################
        EstFrec=abFrecPredicho2([muf,sf**2])
        EstTiem=abTimePredicho2([mut,st**2])


        # Sigma=[[s11,0],[0,s22]]
        # g1=f11*R[ReverseFault.index[i]]**f12*np.exp(-R[ReverseFault.index[i]]*f13)+f14
        # g2=f21*g1+f22
        # Media=(np.vstack((g1,g2))).T
        if j==0:
            SimuladosFrec=EstFrec
            SimuladosTiem=EstTiem
            SimuladosET=EstET
            SimuladosC=EstC
            
        else :
            
            SimuladosFrec=np.vstack((SimuladosFrec,EstFrec))
            SimuladosTiem=np.vstack((SimuladosTiem,EstTiem))
            SimuladosET=np.vstack((SimuladosET,EstET))
            SimuladosC=np.vstack((SimuladosC,EstC))

        
    
        # plt.plot()
        # plt.scatter(Obs[i,0],Obs[i,1])
        # plt.scatter(Simulados[:,0],Simulados[:,1])
        # plt.show()
    # sns.kdeplot(Simulados[:,0], Simulados[:,1],fill=True,alpha=0.3,color="red")  
    try:
        sns.kdeplot(SimuladosFrec[:,0], SimuladosFrec[:,1],fill=True,alpha=0.3,color="red")
    except: 
        print("Error")

    p1=np.median(SimuladosFrec[:,0])
    p2=np.median(SimuladosFrec[:,1])
    plt.scatter(p1,p2,label="Estimado: "+ (str(   np.round(p1,2) ) )+ ","+(str(   np.round(p2,2) ) )    )
    plt.scatter(Obs[i,0],Obs[i,1],label="Observado: "+ (str(   np.round(Obs[i,0],2) ) )+ ","+(str(   np.round(Obs[i,1],2) ) )    )
    plt.title(ReverseFault["Base"][R2.index[i]]+", R="+str(np.round(ReverseFault["R"][R2.index[i]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[i]],2) ))
    plt.legend()
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"Frec"+'.png')
    plt.show()
    
    plt.plot()
    sns.kdeplot(SimuladosTiem[:,0], SimuladosTiem[:,1],fill=True,alpha=0.3,color="red")
    p1=np.median(SimuladosTiem[:,0])
    p2=np.median(SimuladosTiem[:,1])
    plt.scatter(p1,p2,label="Estimado: "+ (str(   np.round(p1,2) ) )+ ","+(str(   np.round(p2,2) ) )    )
    plt.scatter(ObsTiempo[i,0],ObsTiempo[i,1],label="Observado: "+ (str(   np.round(ObsTiempo[i,0],2) ) )+ ","+(str(   np.round(ObsTiempo[i,1],2) ) )    )
    plt.title(ReverseFault["Base"][R2.index[i]]+", R="+str(np.round(ReverseFault["R"][R2.index[i]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[i]],2) ))
    plt.legend()
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"Tiempo"+'.png')

    
    
    plt.show()
    
    plt.plot()
    idsismo=R2.index[i]
    (X,NS,f,Y,VelMuest)=fourier(ReverseFault["Base"][idsismo], CanalH1)   
    NS=NS[X>ReverseFault["Inicio"][idsismo]]
    X=X[X>ReverseFault["Inicio"][idsismo]]
    Mag,Dis,EnTot= ReverseFault[["GMW","R","ET"]].loc[idsismo]
    Dominiot=Directorio2[i][0]

    for r in np.random.choice(len(SimuladosFrec),size=100): 
        # afp,bfp=abFrecPredicho2(Simulados[r,:])
        afp,bfp=SimuladosTiem[r,:]
        if afp>0 and bfp>0 :
            # print(afp,bfp)
            DensMarginalf=gamma.pdf(Dominiot, afp,scale=1/bfp)###############################Marginalt
            if  np.max(DensMarginalf)<5:
                plt.plot(Dominiot,DensMarginalf,color="red",alpha=0.2)
        else :
            print("Problema con negativo")
    plt.plot(Directorio2[i][0],Directorio2[i][1],label="Observado")        
    mufo=Predicemt(t11o,t12o,R2[R2.index[i]])
    sfo=Predicest(t21o,t22o,mut)
    afp,bfp=abTimePredicho2(((mufo,sfo**2)))    
    DensMarginalf=gamma.pdf(Dominiot, afp,scale=1/bfp)###############################Marginalf
    plt.plot(Dominiot,DensMarginalf,color="green",label="MAP")
    plt.legend()
    plt.title(ReverseFault["Base"][R2.index[i]]+", R="+str(np.round(ReverseFault["R"][R2.index[i]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[i]],2) ))
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"Curvas"+"Tiempo"+'.png')
    plt.show()
    


        
    plt.plot()
    idsismo=R2.index[i]
    (X,NS,f,Y,VelMuest)=fourier(ReverseFault["Base"][idsismo], CanalH1)   
    NS=NS[X>ReverseFault["Inicio"][idsismo]]
    X=X[X>ReverseFault["Inicio"][idsismo]]
    Mag,Dis,EnTot= ReverseFault[["GMW","R","ET"]].loc[idsismo]
    Dominiof=np.arange(0,fmax,VelocidadMuestreoFrec)

    for r in np.random.choice(len(SimuladosFrec),size=100): 
        # afp,bfp=abFrecPredicho2(Simulados[r,:])
        afp,bfp=SimuladosFrec[r,:]
        if afp>0 and bfp>0 :
            # print(afp,bfp)
            DensMarginalf=gamma.pdf(Dominiof, afp,scale=1/bfp)/gamma.cdf(fmax, afp,scale=1/bfp)###############################Marginalf
            if  np.max(DensMarginalf)<5:
                plt.plot(Dominiof,DensMarginalf,color="red",alpha=0.2)
        else :
            print("Problema con negativo")
    plt.plot(Directorio[i][0],Directorio[i][1],label="Observado")        
    muto=Predicemf(f11o,f12o,f13o,f14o,R2[R2.index[i]])
    sto=Predicesf(f21o,f22o,muto)
    atp,btp=abFrecPredicho2(((muto,sto**2)))    
    DensMarginalf=gamma.pdf(Dominiof, atp,scale=1/btp)/gamma.cdf(fmax, atp,scale=1/btp)###############################Marginalf
    plt.plot(Dominiof,DensMarginalf,color="green",label="MAP")
    plt.legend()
    plt.title(ReverseFault["Base"][R2.index[i]]+", R="+str(np.round(ReverseFault["R"][R2.index[i]],2))+", M=" +str(np.round(ReverseFault["GMW"][R2.index[i]],2) ))
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"Curvas"+'.png')
    plt.show()
    
    
    
    plt.plot()
    
    
    ###############    
    for j in range(20): # range(len(CorridaEfectiva))
        [Proceso,Real,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]=Predicciones(i,SimuladosFrec[j,0],SimuladosFrec[j,1],SimuladosTiem[j,0],SimuladosTiem[j,1],SimuladosET[j][0],SimuladosC[j][0])
    
        for s in range(Trayectorias+1):
            plt.plot(Real[0],Simulado[s],color="red",alpha=0.1)
        plt.xlabel("T")
        plt.ylabel(r"$\frac{cm}{s^2}$")
        # plt.set_xlabel("T")
        # plt.set_ylabel(r"$\frac{cm}{s^2}$")
        # plt.show(block=False)
        
    ETo=PrediceEnergia(E1o, E2o, E3o, MW2[R2.index[i]], R2[R2.index[i]])

    [Proceso,Real,Simulado,t,f,Z,Dominio,Dominiof, Conjunta]=Predicciones(i,afp,bfp,atp,btp,10**ETo,PrediceCopula())
    plt.plot(Real[0],Real[1],color="blue")
    for s in range(Trayectorias+1):
        plt.plot(Real[0],Simulado[s],color="green",alpha=0.5)
        
    plt.savefig("./Regiones/"+ReverseFault["Base"][R2.index[i]]+"RS"+'.png')
    plt.show()
        
    # fig, axs = plt.subplots(2, 1,figsize=(10,6))                      
    # axs[0,0].pcolormesh(t,f[f<fmax],Z[f<fmax,:], shading='gouraud', rasterized=True)
    # axs[0,0].set_xlabel("s")
    # axs[0,0].set_ylabel("Hz")
    # # axs[0,0].set_colorbar()
    # # Suave=gaussian_filter(np.log10(Conjunta+np.finfo(float).eps), 1)
    # # Suave=gaussian_filter((Conjunta+np.finfo(float).eps), 1)
    # # plt.contour(Dominio[:-1], Dominiof[:-1], Suave, 10, colors='white',origin="lower",alpha=0.2)
    # axs[0,0].contour(Dominio[:-1], Dominiof[:-1], Conjunta, 10, colors='white',origin="lower",alpha=0.2)
    
    

        
    # axs[0].plot(X,NS)
    # axs[0].plot(Dominio[:-1]+ReverseFault["Inicio"].loc[idsismo],Proceso,alpha=0.5,color="red")
    # axs[0].set_xlabel("s")
    # axs[0].set_ylabel(r"$\frac{cm}{s^2}$")
    
    
    
    # fig.suptitle(ReverseFault["Base"][idsismo]+' Mw='+str(ReverseFault["GMW"][idsismo])+" R="+str(round(ReverseFault["R"][idsismo],2)), fontsize=16)
    # fig.tight_layout()
    
    
    # fig.savefig(r"/Users/isaias/Desktop/Atenuacion/Prediccion/"+(ReverseFault["Base"][idsismo]).replace(".","")+'.png',dpi=200, bbox_inches = "tight")
    
    
    
    
    
    
    
################################################################    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    





















