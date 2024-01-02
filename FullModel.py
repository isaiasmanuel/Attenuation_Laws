#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 17:55:12 2022

@author: isaias
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 12:47:23 2022

@author: isaias
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 13:00:25 2022

@author: isaias

"""
#PruebaRecortadaConReal,Errores,Correr,GraficasArticulo,Funciones

print("Entramos")




def EvaluaZhat(TiempoReducido):
    global CopulaDensidad
    tmax=max(X)
    # plt.plot(DensMarginalt)
    DensMarginalf1=scipy.stats.lognorm.pdf(f2, s=a1,scale=np.exp(b1))/scipy.stats.lognorm.cdf(fmax, s=a1,scale=np.exp(b1))*(f2<fmax)
    DistMarginalf1=scipy.stats.lognorm.cdf(f2, s=a1,scale=np.exp(b1))/scipy.stats.lognorm.cdf(fmax, s=a1,scale=np.exp(b1))*(f2<fmax)
    
    DensMarginalf21=scipy.stats.lognorm.pdf(f2, s=a21,scale=np.exp(b21))/scipy.stats.lognorm.cdf(fmax, s=a21,scale=np.exp(b21))*(f2<fmax)
    DensMarginalf22=scipy.stats.lognorm.pdf(f2, s=a22,scale=np.exp(b22))/scipy.stats.lognorm.cdf(fmax, s=a22,scale=np.exp(b22))*(f2<fmax)
    
    
    DensMarginalf2=m2*DensMarginalf21+(1-m2)*DensMarginalf22
    DistMarginalf2=integrate.cumtrapz(DensMarginalf2,f2 , initial=0)
    DensMarginalt1=scipy.stats.lognorm.pdf(TiempoReducido, s=am1, scale=np.exp(bm1),loc=dm1)/(scipy.stats.lognorm.cdf(TiempoReducido[-1], s=am1, scale=np.exp(bm1),loc=dm1))
    DensMarginalt2=scipy.stats.lognorm.pdf(TiempoReducido, s=am2, scale=np.exp(bm2),loc=dm1+alpha)/(scipy.stats.lognorm.cdf(TiempoReducido[-1], s=am2, scale=np.exp(bm2),loc=dm1+alpha))
    
    
    DensMarginalTotal=DensMarginalt1*m1+DensMarginalt2*(1-m1)
    
    DensMarginalt1=DensMarginalTotal*(TiempoReducido<dm1+alpha)
    DensMarginalt2=DensMarginalTotal*(TiempoReducido>dm1+alpha)
    DistMarginalt1=integrate.cumtrapz(DensMarginalt1,TiempoReducido , initial=0)
    DistMarginalt2=integrate.cumtrapz(DensMarginalt2,TiempoReducido , initial=0)

    C1=DistMarginalt1[-1]
    C2=DistMarginalt2[-1]
    DensMarginalt1=DensMarginalt1/DistMarginalt1[-1]
    DensMarginalt2=DensMarginalt2/DistMarginalt2[-1]
    DistMarginalt1=DistMarginalt1/DistMarginalt1[-1]
    DistMarginalt2=DistMarginalt2/DistMarginalt2[-1]
    
    
    InvNormalf1=sp.stats.norm.ppf(DistMarginalf1) 
    
    
    
    InvNormalf2=sp.stats.norm.ppf(DistMarginalf2) 
    
    
    
    Sigma=np.array([[1,b],[b,1]])
    # Sigma=np.array([[a,b],[b,c]])
    SigmaInv=np.linalg.inv(Sigma)
    SigmaDet=np.linalg.det(Sigma)
    
    
    Const=sp.special.gamma((nu+2)/2)*sp.special.gamma((nu)/2)/(sp.special.gamma((nu+1)/2)**2*np.sqrt(SigmaDet))
    
    
    

    def CopulaDensidad(X):
        i=X[0]
        j=X[1]
        Vec=np.array([InvNormalf[i],InvNormalt[j]])
        # Vec=np.array([DistMarginalf[i],DistMarginalt[j]])
        Vec.shape=(2,1)
        Valor=0
        
        if  all(~np.isinf(Vec)) :
            Valor= np.exp(-Vec.T@(SigmaInv-np.identity(2))@Vec/2)*SigmaDet**(-1/2)* (DensMarginalt[j] * DensMarginalf[i]) 
        if np.isinf(Valor):
            Valor=0            
        return Valor         

    if len(TiempoReducido)==Reducido.shape[1]:        
        Zhat=np.copy(Reducido)
    if len(TiempoReducido)==t2.shape[0]:        
        Zhat=np.copy(PSDFSuavizada)
 
    Zhat1=np.copy(Zhat-Zhat)
    Zhat2=np.copy(Zhat1)
    
    
    x = np.array(np.arange(np.sum(f2<fmax)))
    y1=np.arange(np.sum(TiempoReducido<dm1),np.min((np.sum(DistMarginalt1<cuantil),np.sum(TiempoReducido<tmax))))    
    
    y1=np.arange(np.sum(TiempoReducido<dm1),np.min((np.sum(TiempoReducido<dm1+alpha),np.sum(DistMarginalt1<cuantil),np.sum(TiempoReducido<tmax))))        
    y2=np.arange(np.sum(TiempoReducido<dm1+alpha),np.min((np.sum(DistMarginalt2<cuantil),np.sum(TiempoReducido<tmax))))
    
    y=y1    
    Malla=pd.DataFrame({'x':np.repeat(x,y.shape[0]),
                  'y':np.tile(y,x.shape[0])})
    Malla=Malla.to_numpy()
    
    DensMarginalt=DensMarginalt1
    DistMarginalt=DistMarginalt1
    DensMarginalf=DensMarginalf1
    DistMarginalf=DistMarginalf1
    InvNormalf=InvNormalf1
    
    InvNormalt=sp.stats.norm.ppf(DistMarginalt1) 
    
    
    p = get_context("fork").Pool(Hilos)
    Valores=p.map(CopulaDensidad,Malla)  
    p.close()    
    
    for i in range(Malla.shape[0]):    
        loc=Malla[i,:]
        Zhat1[loc[0],loc[1]]=Valores[i]
    
    y=y2
    Malla=pd.DataFrame({'x':np.repeat(x,y.shape[0]),
                  'y':np.tile(y,x.shape[0])})
    Malla=Malla.to_numpy()
        
    DensMarginalt=DensMarginalt2
    DistMarginalt=DistMarginalt2
    DensMarginalf=DensMarginalf2
    DistMarginalf=DistMarginalf2
    InvNormalf=InvNormalf2


    InvNormalt=sp.stats.norm.ppf(DistMarginalt2) 
    p = get_context("fork").Pool(Hilos)
    Valores=p.map(CopulaDensidad,Malla)  
    p.close()    
    
    for i in range(Malla.shape[0]):    
        loc=Malla[i,:]
        Zhat2[loc[0],loc[1]]=Valores[i]
        
    Zhat=C1*Zhat1+(C2)*Zhat2
    Zhat[np.isnan(Zhat)]=0
    
    return Zhat




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
    
    
    
    
    Sintetico=Proceso
    Registro=np.vstack((t2,Sintetico))
    Registro=np.vstack((Registro,Sintetico))
    Registro=np.vstack((Registro,Sintetico))
    Registro.shape
    ############### Sismos registrados lejanos
    
    # Muestra.sort_values(by=["R"])[-10:]
    # STFT(59,10,90) #Canditato
    # STFT(52,10,90) #Canditato
    # STFT(116,10,90) #Canditato
    # STFT(168,10,90) #Canditato
    # Muestra.sort_values(by=["R"])[-30:-20]
    
    
    RegistroChristen=Registro.T
    np.savetxt(r'./RegistroChristen2'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv',RegistroChristen,delimiter=',', fmt=('%s'))
    # np.savetxt(r'./RegistroChristen2.csv',RegistroChristen,delimiter=',', fmt=('%s'))
    tst = Acelerogram( r'./RegistroChristen2'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv', Dt=-1, skiprows=1)
    
    return Proceso,np.array(Respuesta(tst))[1]


import pycopula
import os
import matplotlib.pyplot as plt
from datetime import datetime

try:
    os.chdir('/home/isaias.ramirez/')
except:
    os.chdir('/Users/isaias/Desktop/Servidor/')



try:
    os.mkdir("./Corrida")
    print("Directorio Corrida nuevo creado")
except:
    print("Directorio Corrida ya existe")

try:
    os.mkdir("./Graficas")
    print("Directorio Graficas nuevo creado")
except:
    print("Directorio Graficas ya existe")


from Funciones import *

TiempoInicial=pd.read_csv("TiempoInicial.csv", delimiter=",")

Muestra

Muestra=pd.merge(Muestra, TiempoInicial, on="Base")

# Muestra.loc[Muestra["Base"]=="SUCH8509.191"]


# NS=np.genfromtxt("./ImperialCM.csv",delimiter=',')
# NS=NS[1:]
# NS=np.reshape(NS,(NS.shape[0]*NS.shape[1]))
# NS=NS*980.665
# X=np.arange(NS.shape[0])*0.005
# plt.plot(X,NS)
# for i in np.hstack((np.array([0,1,2,143,238,247]),Muestra.index.values[3:7])):
    
# for i in np.hstack((np.array([0,1,2,143,238,247]),Muestra.index.values[3:20])):
# i=169
# (169, 184, 82, 75, 135, 57, 82)

# for i in Muestra.index.values[:]:    
# for i in (Muestra[Muestra["IDEstacion"].str.contains("CUP")]).index.values[:] :
# for i in np.arange(11,25):



Intervalo=0#IntervaloFinal

# for i in np.arange(Intervalo*len(Muestra.index.values[:])//10,(Intervalo+1)*len(Muestra.index.values[:])//10):
for i in np.arange(len(Muestra.index.values[:])):
    
    if not os.path.isfile(r'./Corrida/Corrida'+ (Muestra["Base"][i]).replace(".","")  +".csv"):
    
        try:    
            print(Muestra["Base"][i])
            idsismo=np.copy(i)
            CanalV=np.int(np.float(RecuperaOrientacion(Muestra["Base"][i])[7]))
            CanalH1=np.int(np.float(RecuperaOrientacion(Muestra["Base"][i])[8]))
            CanalH2=np.int(np.float(RecuperaOrientacion(Muestra["Base"][i])[9]))
            
            (X,NS,f,Y,VelMuest)=fourier(Muestra["Base"][i], CanalH1) #Usar registrocrudo cuando no quiero normalizado y centrado
            
            
            LenMCMC=10000
            Hilos=46 ###Elegir 12 si  es local y 48 en el servidor
            fmax=15
            Trayectorias=20
            tmax=max(X)
            tmin=00
            nu=30        #Grados de libertad
            cuantil=1 #La muestra se corta en este cuantil de la marginal en tiempo
            Reduccion=int(1*VelMuest )#Numero de muestras que se juntan, no hay interseccion
            
            Dominio=X>Muestra["Inicio"][idsismo]
            X=X[Dominio]
            NS=NS[Dominio]
            
            # plt.plot(X,NS,linewidth=.5)
            # print(VelMuest)
            
            ###############Creamos registro
            Sintetico=NS
            Registro=np.vstack((X,NS))
            Registro=np.vstack((Registro,Sintetico))
            Registro=np.vstack((Registro,Sintetico))
            Registro.shape
            
            RegistroChristen=Registro.T
            
            np.savetxt(r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv',RegistroChristen,delimiter=',', fmt=('%s'))
            
            
            ####### Calculamos STFT
            NS=signal.detrend(NS)
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
            batch=np.sum(X<4) #50 es el bueno
            ####### Calculamos STFT
            
            
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
            
            # fmax=np.max(f)
            
            
            ########Mantenemos energia
            
            # print(np.max(f),f[1]-f[0],"\n",np.max(t),t[1]-t[0])
            
            Marginalf=np.ones(len(f))
            for i in np.arange(len(f)):
                Marginalf[i]=integrate.cumtrapz(Z[i,:],t , initial=0)[-1]
            
            Normalizacion=integrate.cumtrapz(Marginalf,f , initial=0)[-1]/integrate.cumtrapz(NS**2,X , initial=0)[-1]
            
            
            PSDF=np.copy(Z)*Normalizacion**(-1)
            
            
            PSDF[:,0]=0
            
            ######### Interpolamos
            from scipy.interpolate import griddata
            import matplotlib.pyplot as plt
            
            Puntos=np.ones((len(f)*len(t),2))
            Valores=np.ones((len(f)*len(t)))
            Suavizada=scipy.interpolate.RectBivariateSpline(f,t,PSDF, kx=4, ky=4)
            
            contador=0
            
            for i in f:
                for j in t:        
                    Puntos[contador,:]=[i,j]    
                    Valores[contador]=Suavizada(i,j)
                    if i==0:
                        Valores[contador]=0
                    contador+=1
            
            
            grid_x, grid_y = np.meshgrid(f, X, indexing='ij')
            
            
            f2=grid_x[:,0]
            t2=grid_y[0,:]
            
            grid_z0 = griddata(Puntos, Valores, (grid_x, grid_y), method='linear')
            
            PSDFSuavizada=np.copy(grid_z0)
            
            # plt.pcolormesh(t2, f2, PSDFSuavizada, shading='gouraud')
            
            # plt.pcolormesh(t, f, PSDF, shading='gouraud')
            
            
            
            #####################
            
            
            #################Normalizado y calculando marginales
            
            
            Marginalf=np.ones(len(f2))
            for i in np.arange(len(f2)):
                Marginalf[i]=integrate.cumtrapz(PSDFSuavizada[i,:],t2 , initial=0)[-1]
            
            Normalizacion=integrate.cumtrapz(Marginalf,f2 , initial=0)[-1]
            
            # print(Normalizacion/integrate.cumtrapz(NS**2,X , initial=0)[-1])
            
            
            #f,t,Z=STFT(idsismo,desv=10,over=1/(X[1]-X[0])*.9)
            
            Marginalf=np.ones(len(f2))
            for i in np.arange(len(f2)):
                Marginalf[i]=integrate.cumtrapz(PSDFSuavizada[i,:]/Normalizacion,t2 , initial=0)[-1]
                
            
            Marginalt=np.ones(len(t2))
            for i in np.arange(len(t2)):
                Marginalt[i]=integrate.cumtrapz(PSDFSuavizada[:,i]/Normalizacion,f2 , initial=0)[-1]
            
            
            # integrate.cumtrapz(Marginalt,t2 , initial=0)[-1]
            # integrate.cumtrapz(Marginalf,f2 , initial=0)[-1]
            
            
            PSDFNorm=np.copy(PSDFSuavizada)/Normalizacion
            
            # plt.plot(t2,Marginalt)    
            # plt.plot(f2,Marginalf)    
            
            
            
            # plt.scatter(T_list,Real[1])
            # for i in range(Trayectorias):
            #     plt.scatter(T_list ,RSReal[i][1],color="orange",alpha=0.2)
            
    
            
            
            ########################### Reducida
            
            
            import torch
            import torch.nn
            
            
            PSDFNorm2=PSDFNorm#[f2<fmax,:][:,t2<tmax]
            Tensor=torch.from_numpy(PSDFNorm2)
            Tensor=torch.reshape(Tensor, (1,PSDFNorm2.shape[0], PSDFNorm2.shape[1]))
            
            # plt.pcolormesh(PSDFNorm2, shading='gouraud')
            # plt.pcolormesh(Tensor[0,:,:], shading='gouraud')
            
            
            m=torch.nn.AvgPool2d((1,Reduccion))
            Reducido=m(Tensor)
            Reducido=Reducido[0,:,:].numpy()
            TiempoReducido=t2[np.arange(0,Reducido.shape[-1])*Reduccion]
            TiempoReducido=TiempoReducido+(TiempoReducido[1]-TiempoReducido[0])/2
            
            
            
            ###########################
            
            
            import seaborn as sn
            
            locs=np.reshape(f2[f2<fmax],(np.sum(f2<fmax),1))
            K = skMatern(nu=1/2, length_scale=2)
            # Mat=np.diag(np.ones(len(f2)))
            Mat1=K(locs)#@(Pol/2)
            # sn.heatmap(Mat)
            #sn.heatmap(Mat)
            #sn.heatmap(K)
            #Mat=np.diag(np.ones(len(locs)))
            V=np.linalg.inv(Mat1)
            
            
            
            locs=np.reshape(TiempoReducido,(len(TiempoReducido),1))
            K = skMatern(nu=1/2, length_scale=8)
            # Mat=np.diag(np.ones(len(f2)))
            Mat2=K(locs)#@(Pol/2)
            # sn.heatmap(Mat)
            
            #Mat=np.diag(np.ones(len(locs)))
            U=np.linalg.inv(Mat2)
            
            # sn.heatmap(Mat1)
            # sn.heatmap(Mat2)
            
            # U.shape
            # V.shape
            # A.shape
            
            
            
            ######################
            
            
            Dist=Muestra["R"][idsismo]
            
            TiempoEsp=(Dist/3.5-Dist/6.1)
            
            
            
            
            
            # Theta=np.array([1.1,1,0.01,#Gamma1
            #                 1,1,18,
            #                 0.2,#Parametro mezclante
            #                 0.97,0.96,#Lognormal
            #                 0.8,#Clayton
            #                 # 1,0.8,1,#Clayton
            #                 0.001#Error              
            #                 ])    
            
            # Theta=np.array([0.5,2,0.01,#Gamma1
            #                 0.4,2.4,13,
            #                 0.05,#Parametro mezclante
            #                 0.5,0.8,#Lognormal
            #                 0.1,#Clayton
            #                 0.1#Error              
            #                 ])   
            
            
            
            ###########################
            
            MarginalfReduc=np.ones(len(f2))
            for i in np.arange(len(f2)):
                MarginalfReduc[i]=integrate.cumtrapz(Reducido[i,:],TiempoReducido , initial=0)[-1]
                
            MarginaltReduc=np.ones(len(TiempoReducido))
            for i in np.arange(len(TiempoReducido)):
                MarginaltReduc[i]=integrate.cumtrapz(Reducido[:,i],f2 , initial=0)[-1]
            
            
            
            ###########################
            
            
            def EnSuppT(Theta):
                #(am1,bm1,dm1,am2,bm2,alpha,m1)
                #  0   1    2  3   4    5   6
                # a1,b1
                # 7, 8
                # b
                # 9  
                # sigma
                # 10
                #Falta condicion de picos
                am1,bm1,dm1,am2,bm2,alpha,m1,sigma1,sigma2=Theta
            
                
                DensMarginalt1=scipy.stats.lognorm.pdf(t2, s=am1, scale=np.exp(bm1),loc=dm1)
                DensMarginalt2=scipy.stats.lognorm.pdf(t2, s=am2, scale=np.exp(bm2),loc=dm1+alpha)
                
            
                Moda1= t2[np.argmax(DensMarginalt1)] 
                Moda2= t2[np.argmax(DensMarginalt2)]     
                
                Var1=(np.exp(am1**2)-1)*np.exp(2*bm1+am1**2)
                Var2=(np.exp(am2**2)-1)*np.exp(2*bm2+am2**2)
                # Media1=np.exp((bm1)+am1**2/2)+dm1
                # Media2=np.exp((bm2)+am2**2/2)+dm1+alpha
                Pico1=np.max(DensMarginalt1*m1)
                Pico2=np.max(DensMarginalt2*(1-m1))
                # Pico1=m1*scipy.stats.lognorm.pdf(Moda1, s=am1, scale=np.exp(bm1),loc=dm1)
                # Pico2=(1-m1)*scipy.stats.lognorm.pdf(Moda2, s=am2, scale=np.exp(bm2),loc=(dm1+alpha))
            
                
                Verificacion= all(Theta>0) and all (Theta<100)
                # Verificacion=Verificacion and  (Media1<Media2)
                Verificacion= Verificacion and (dm1+alpha>Moda1)
                Verificacion= Verificacion and (TiempoEsp*0.5<Moda2-Moda1)
                Verificacion= Verificacion and (TiempoEsp*0.5<alpha)
                Verificacion= Verificacion and (Moda2-Moda1<TiempoEsp*1.5)        
                Verificacion= Verificacion and (alpha<TiempoEsp*1.5)        
                # Verificacion= Verificacion and ((Dist/3.5-Dist/6.1)*0.3<Media2-Media1)
                # Verificacion= Verificacion and (Media2-Media1<(Dist/3.5-Dist/6.1)*2)        
                Verificacion= Verificacion and (m1<0.5)
                # Verificacion= Verificacion and (bm1<np.log(20))
                # Verificacion= Verificacion and (bm2<np.log(20))
                Verificacion= Verificacion and (Pico1<Pico2)
                Verificacion= Verificacion and (Var1<50**2)
                Verificacion= Verificacion and (Var2<50**2)
                Verificacion= Verificacion and (sigma1<sigma2)
                Verificacion= Verificacion and dm1<2
                Verificacion= Verificacion and np.sum(DensMarginalt1>0)
    
                return Verificacion    
            
            
            def EnergiaT(Theta):
                
                # plt.plot(DensMarginalt)
                
                # Sesgo1=(np.exp(am1**2)+2)*np.sqrt(np.exp(am1**2)-1)
                # Sesgo2=(np.exp(am2**2)+2)*np.sqrt(np.exp(am2**2)-1)
                am1,bm1,dm1,am2,bm2,alpha,m1,sigma1,sigma2=Theta
                DensMarginalt1=scipy.stats.lognorm.pdf(TiempoReducido, s=am1, scale=np.exp(bm1),loc=dm1)
                DensMarginalt2=scipy.stats.lognorm.pdf(TiempoReducido, s=am2, scale=np.exp(bm2),loc=dm1+alpha)
                
                Var1=(np.exp(am1**2)-1)*np.exp(2*bm1+am1**2)
                Var2=(np.exp(am2**2)-1)*np.exp(2*bm2+am2**2)
                
                # Salida=np.sum((Marginalt[t2<am1+alpha]-(m1*DensMarginalt1+(1-m1)*DensMarginalt2)[t2<am1+alpha])**2)/(sigma1**2)+np.sum([t2<am1+alpha])*np.log(sigma1)+np.sum((Marginalt[t2>=am1+alpha]-(m1*DensMarginalt1+(1-m1)*DensMarginalt2)[t2>=am1+alpha])**2)/(sigma2**2)+np.sum([t2>=am1+alpha])*np.log(sigma2)
                # Salida=np.sum((MarginaltReduc[TiempoReducido<dm1+alpha]-(m1*DensMarginalt1+(1-m1)*DensMarginalt2)[TiempoReducido<dm1+alpha])**2)/(sigma1**2)+np.sum([TiempoReducido<dm1+alpha])*np.log(sigma1)+np.sum((MarginaltReduc[TiempoReducido>=dm1+alpha]-(m1*DensMarginalt1+(1-m1)*DensMarginalt2)[TiempoReducido>=dm1+alpha])**2)/(sigma2**2)+np.sum([TiempoReducido>=dm1+alpha])*np.log(sigma2)    
                    
                Diferencia1=MarginaltReduc[(TiempoReducido<dm1+alpha)]-(m1*DensMarginalt1+(1-m1)*DensMarginalt2)[(TiempoReducido<dm1+alpha)]
                Diferencia2=MarginaltReduc[(TiempoReducido>=dm1+alpha)* (TiempoReducido<=tmax)]-(m1*DensMarginalt1+(1-m1)*DensMarginalt2)[(TiempoReducido>=dm1+alpha)* (TiempoReducido<=tmax)]
                # Salida=np.sum(Diferencia1**2)/(sigma1**2)+np.sum([(TiempoReducido<dm1+alpha)])*np.log(sigma1)+np.sum(Diferencia2**2)/(sigma2**2)+np.sum([(TiempoReducido>=dm1+alpha)* (TiempoReducido<=tmax)])*np.log(sigma2)    
                
                
                Sobreajuste=1
                Salida= np.sum(Sobreajuste*(Diferencia1>=0)*Diferencia1**2+1*(Diferencia1<0)*Diferencia1**2)/(sigma1**2)+np.sum([(TiempoReducido<dm1+alpha)])*np.log(sigma1)+np.sum(Sobreajuste*(Diferencia2>=0)*Diferencia2**2+1*(Diferencia2<0)*Diferencia2**2)/(sigma2**2)+np.sum([(TiempoReducido>=dm1+alpha)* (TiempoReducido<=tmax)])*np.log(sigma2)#+Var1+Var2
                
                
                
                # Salida=np.trace((A.T)@A)/(2*sigma**2)+(len(TiempoReducido)*len(f2))*np.log(sigma)/2+am1**2/2+am2**2/2+a1**2/2 
                # Salida=np.trace((A1.T)@A1)/(sigma1**2)+np.trace((A2.T)@A2)/(sigma2**2)+((np.sum(TiempoReducido<dm1+alpha))*len(f2))*np.log(sigma1)+((np.sum(TiempoReducido>=dm1+alpha))*len(f2))*np.log(sigma2)
                return Salida 
            
            
            
            def EnSuppF(Theta):
                a,b,sigma=Theta
                Var=(np.exp(a**2)-1)*np.exp(2*b+a**2)
                Verificacion= all(Theta>0) and all (Theta<100)
                Verificacion= Verificacion and (Var>4 )
                # Verificacion= Verificacion and (Var<10**2)        
                
                return Verificacion    
            
            
            def EnergiaF(Theta):
                
                a1,b1,sigma=Theta
                
    
                
                DensMarginalf=scipy.stats.lognorm.pdf(f2, s=a1,scale=np.exp(b1))/scipy.stats.lognorm.cdf(fmax, s=a1,scale=np.exp(b1))*(f2<fmax)
                    
                Salida=np.sum((Marginalf[f2<fmax]-DensMarginalf[f2<fmax])**2)/(sigma**2)+np.sum([f2<fmax])*np.log(sigma)    
            
                return Salida 
            
            # Theta=np.array([5.00000000e-01, 1.00000000e+00, 1.00000000e-02, 1.40000000e+00,
            #        2.00000000e+00, 1.09490391e+01, 5.00000000e-02, 5.00000000e-01,
            #        1.00000000e+00, 3.00000000e-01, 1.00000000e-03, 1.00000000e-03])
            
            
            
            def EnSupp(Theta):
                
                #(am1,bm1,dm1,am2,bm2,alpha,m1)
                #  0   1    2  3   4    5   6
                # a1,b1
                # 7, 8
                # b
                # 9  
                # sigma
                # 10
                #Falta condicion de picos
            
                am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a21,b21,a22,b22,m2,b,sigma1,sigma2=Theta
                
                # scipy.stats.lognorm.stats(s=am1, scale=np.exp(bm1),loc=dm1, moments='mvsk')
                
                DensMarginalt1=scipy.stats.lognorm.pdf(TiempoReducido, s=am1, scale=np.exp(bm1),loc=dm1)
                DensMarginalt2=scipy.stats.lognorm.pdf(TiempoReducido, s=am2, scale=np.exp(bm2),loc=dm1+alpha)
                
                
                DensMarginalf21=scipy.stats.lognorm.pdf(f2, s=a21,scale=np.exp(b21))/scipy.stats.lognorm.cdf(fmax, s=a21,scale=np.exp(b21))*(f2<fmax)
                DensMarginalf22=scipy.stats.lognorm.pdf(f2, s=a22,scale=np.exp(b22))/scipy.stats.lognorm.cdf(fmax, s=a22,scale=np.exp(b22))*(f2<fmax)
                # dm1+alpha>Moda1
                Moda1= TiempoReducido[np.argmax(DensMarginalt1)] 
                Moda2= TiempoReducido[np.argmax(DensMarginalt2)]    
                Moda21=f2[np.argmax(DensMarginalf21)] 
                Moda22=f2[np.argmax(DensMarginalf22)] 
                
                Var1=(np.exp(am1**2)-1)*np.exp(2*bm1+am1**2)
                Var2=(np.exp(am2**2)-1)*np.exp(2*bm2+am2**2)
                
                Var21=(np.exp(a21**2)-1)*np.exp(2*b21+a21**2)
                Var22=(np.exp(a22**2)-1)*np.exp(2*b22+a22**2)
                
                # Media1=np.exp((bm1)+am1**2/2)+dm1
                # Media2=np.exp((bm2)+am2**2/2)+dm1+alpha
                Pico1=np.max(DensMarginalt1*m1)
                Pico2=np.max(DensMarginalt2*(1-m1))
                # Pico1=m1*scipy.stats.lognorm.pdf(Moda1, s=am1, scale=np.exp(bm1),loc=dm1)
                # Pico2=(1-m1)*scipy.stats.lognorm.pdf(Moda2, s=am2, scale=np.exp(bm2),loc=(dm1+alpha))
                
                
                Verificacion= all(Theta>0) and all (Theta<100)
                # Verificacion= all(Theta[np.arange(len(Theta))!=9]) and all (Theta<100)
                
                # Verificacion=Verificacion and  (Media1<Media2)
                # Verificacion= Verificacion and (dm1+alpha>Moda1)
                Verificacion= Verificacion and (TiempoEsp*0.4<Moda2-Moda1)
                Verificacion= Verificacion and (TiempoEsp*0.4<alpha)
                Verificacion= Verificacion and (Moda2-Moda1<TiempoEsp*1.5)        
                Verificacion= Verificacion and (alpha<TiempoEsp*1.5)        
                # Verificacion= Verificacion and ((Dist/3.5-Dist/6.1)*0.3<Media2-Media1)
                # Verificacion= Verificacion and (Media2-Media1<(Dist/3.5-Dist/6.1)*2)        
                Verificacion= Verificacion and (m1<1)
                Verificacion= Verificacion and (m2<1)
                # Verificacion= Verificacion and (bm1<np.log(20))
                # Verificacion= Verificacion and (bm2<np.log(20))
                Verificacion= Verificacion and (b<1)
                Verificacion= Verificacion and (Pico1<Pico2)
                Verificacion= Verificacion and (sigma1<sigma2)
                Verificacion= Verificacion and dm1<2            
                Verificacion= Verificacion and (Var21>4 and Var22<2)
                
                # Verificacion= Verificacion and (Sesgo1<Sesgo2)
                # Verificacion= Verificacion and (Var1<Var2)
                return Verificacion    
            
            
            
            def Energia(Theta):
        
                global CopulaDensidad
                am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a21,b21,a22,b22,m2,b,sigma1,sigma2=Theta
                b=-b    
                # plt.plot(DensMarginalt)
            
                        
                DensMarginalf1=scipy.stats.lognorm.pdf(f2, s=a1,scale=np.exp(b1))/scipy.stats.lognorm.cdf(fmax, s=a1,scale=np.exp(b1))*(f2<fmax)
                DistMarginalf1=scipy.stats.lognorm.cdf(f2, s=a1,scale=np.exp(b1))/scipy.stats.lognorm.cdf(fmax, s=a1,scale=np.exp(b1))*(f2<fmax)
                
                DensMarginalf21=scipy.stats.lognorm.pdf(f2, s=a21,scale=np.exp(b21))/scipy.stats.lognorm.cdf(fmax, s=a21,scale=np.exp(b21))*(f2<fmax)
                DensMarginalf22=scipy.stats.lognorm.pdf(f2, s=a22,scale=np.exp(b22))/scipy.stats.lognorm.cdf(fmax, s=a22,scale=np.exp(b22))*(f2<fmax)
                
                
                DensMarginalf2=m2*DensMarginalf21+(1-m2)*DensMarginalf22
                DistMarginalf2=integrate.cumtrapz(DensMarginalf2,f2 , initial=0)
                DensMarginalt1=scipy.stats.lognorm.pdf(TiempoReducido, s=am1, scale=np.exp(bm1),loc=dm1)/(scipy.stats.lognorm.cdf(TiempoReducido[-1], s=am1, scale=np.exp(bm1),loc=dm1))
                DensMarginalt2=scipy.stats.lognorm.pdf(TiempoReducido, s=am2, scale=np.exp(bm2),loc=dm1+alpha)/(scipy.stats.lognorm.cdf(TiempoReducido[-1], s=am2, scale=np.exp(bm2),loc=dm1+alpha))
                
                
                DensMarginalTotal=DensMarginalt1*m1+DensMarginalt2*(1-m1)
                DensMarginalt1=DensMarginalTotal*(TiempoReducido<dm1+alpha)
                DensMarginalt2=DensMarginalTotal*(TiempoReducido>dm1+alpha)
                DistMarginalt1=integrate.cumtrapz(DensMarginalt1,TiempoReducido , initial=0)
                DistMarginalt2=integrate.cumtrapz(DensMarginalt2,TiempoReducido , initial=0)
            
                C1=DistMarginalt1[-1]
                C2=DistMarginalt2[-1]
                DensMarginalt1=DensMarginalt1/DistMarginalt1[-1]
                DensMarginalt2=DensMarginalt2/DistMarginalt2[-1]
                DistMarginalt1=DistMarginalt1/DistMarginalt1[-1]
                DistMarginalt2=DistMarginalt2/DistMarginalt2[-1]
                
                
                InvNormalf1=sp.stats.norm.ppf(DistMarginalf1) 
    
                
                
                InvNormalf2=sp.stats.norm.ppf(DistMarginalf2) 
    
                
                
                Sigma=np.array([[1,b],[b,1]])
                # Sigma=np.array([[a,b],[b,c]])
                SigmaInv=np.linalg.inv(Sigma)
                SigmaDet=np.linalg.det(Sigma)
                
                
                # Const=sp.special.gamma((nu+2)/2)*sp.special.gamma((nu)/2)/(sp.special.gamma((nu+1)/2)**2*np.sqrt(SigmaDet))
                
                
                def CopulaDensidad(X):
                    i=X[0]
                    j=X[1]
                    Vec=np.array([InvNormalf[i],InvNormalt[j]])
                    # Vec=np.array([DistMarginalf[i],DistMarginalt[j]])
                    Vec.shape=(2,1)
                    Valor=0
                    
                    if  all(~np.isinf(Vec)) :
                        Valor= np.exp(-Vec.T@(SigmaInv-np.identity(2))@Vec/2)*SigmaDet**(-1/2)* (DensMarginalt[j] * DensMarginalf[i]) 
                    if np.isinf(Valor):
                        Valor=0            
                    return Valor            
                    
                Zhat=np.copy(Reducido)
                Zhat1=np.copy(Zhat-Zhat)
                Zhat2=np.copy(Zhat1)
                
                
                x = np.array(np.arange(np.sum(f2<fmax)))
                y1=np.arange(np.sum(TiempoReducido<dm1),np.min((np.sum(DistMarginalt1<cuantil),np.sum(TiempoReducido<tmax))))    
                
                y1=np.arange(np.sum(TiempoReducido<dm1),np.min((np.sum(TiempoReducido<=dm1+alpha),np.sum(DistMarginalt1<cuantil),np.sum(TiempoReducido<tmax))))        
                y2=np.arange(np.sum(TiempoReducido<=dm1+alpha),np.min((np.sum(DistMarginalt2<cuantil),np.sum(TiempoReducido<tmax))))
                
                y=y1    
                Malla=pd.DataFrame({'x':np.repeat(x,y.shape[0]),
                              'y':np.tile(y,x.shape[0])})
                Malla=Malla.to_numpy()
                
                DensMarginalt=DensMarginalt1
                DistMarginalt=DistMarginalt1
                DensMarginalf=DensMarginalf1
                DistMarginalf=DistMarginalf1
                InvNormalf=InvNormalf1
    
                
                InvNormalt=sp.stats.norm.ppf(DistMarginalt1) 
    
                
                
                p = get_context("fork").Pool(Hilos)
                Valores=p.map(CopulaDensidad,Malla)  
                p.close()    
                
                for i in range(Malla.shape[0]):    
                    loc=Malla[i,:]
                    Zhat1[loc[0],loc[1]]=Valores[i]
                
                y=y2
                Malla=pd.DataFrame({'x':np.repeat(x,y.shape[0]),
                              'y':np.tile(y,x.shape[0])})
                Malla=Malla.to_numpy()
                    
                DensMarginalt=DensMarginalt2
                DistMarginalt=DistMarginalt2
                DensMarginalf=DensMarginalf2
                DistMarginalf=DistMarginalf2
                InvNormalf=InvNormalf2
    
            
                InvNormalt=sp.stats.norm.ppf(DistMarginalt2) 
                p = get_context("fork").Pool(Hilos)
                Valores=p.map(CopulaDensidad,Malla)  
                p.close()    
                
                for i in range(Malla.shape[0]):    
                    loc=Malla[i,:]
                    Zhat2[loc[0],loc[1]]=Valores[i]
                    
                Zhat=C1*Zhat1+(C2)*Zhat2
                Zhat[np.isnan(Zhat)]=0
                A=Reducido-Zhat
                A=A[:,TiempoReducido<tmax][f2<fmax,:]
                # Salida=np.trace(U@(A.T)@V@A)/(2*sigma1*sigma2)+U.shape[0]**2/2*np.log(sigma1)+V.shape[0]**2/2*np.log(sigma2)+Sesgo1**2/2+Sesgo2**2/2+a1**2/2
                A.shape
                A1=A[:,TiempoReducido[TiempoReducido<tmax]<=dm1+alpha]
                A2=A[:,TiempoReducido[TiempoReducido<tmax]>=dm1+alpha]
                
                # plt.pcolormesh(A2, shading='gouraud')
                # plt.pcolormesh(Zhat, shading='gouraud')
                # plt.pcolormesh(Reducido, shading='gouraud')
            
                # Salida=np.trace((A1.T)@A1)/(sigma1**2)+np.trace((A2.T)@A2)/(sigma2**2)+((np.sum(TiempoReducido<dm1+alpha))*len(f2))*np.log(sigma1)+((np.sum(TiempoReducido>dm1+alpha))*len(f2))*np.log(sigma2)+Sesgo1**2+Sesgo2**2+a1**2+bm1**2+bm2**2
               
                
                Salida=np.trace((A1.T)@A1)/(sigma1**2)+np.trace((A2.T)@A2)/(sigma2**2)+((np.sum(TiempoReducido<=dm1+alpha))*len(f2))*np.log(sigma1)+((np.sum(TiempoReducido>=dm1+alpha))*len(f2))*np.log(sigma2)
                
                return Salida 
            
            
            
            tmax=max(X)
            
            # tmax=60
            
            Theta=np.array([0.5,1,0.01,#Gamma1
                            1.4,2,TiempoEsp*0.6,
                            0.05,#Parametro mezclante
                            0.001,0.002#Error              
                            ])    
            
            
            
            # Theta=np.array([0.5,1,0.01,#Gamma1
            #                 1.4,2,TiempoEsp*0.9,
            #                 0.05,#Parametro mezclante
            #                 0.001,0.001#Error              
            #                 ])    
            
            
            am1,bm1,dm1,am2,bm2,alpha,m1,sigma1,sigma2=Theta
            
            
            np.random.seed(1)
            TwEnergia = pytwalk.pytwalk( n=len(Theta), U=EnergiaT, Supp=EnSuppT)
            TwEnergia.Run( T=100000, x0=Theta+0.0001, xp0=Theta)
            Corrida=TwEnergia.Output
            
            # plt.plot(-Corrida[:,-1])
            FinalT= Corrida[-1,:-1]
            am1,bm1,dm1,am2,bm2,alpha,m1,sigma1,sigma2=FinalT
            
            DensMarginalt1=scipy.stats.lognorm.pdf(t2, s=am1, scale=np.exp(bm1),loc=dm1)
            DensMarginalt2=scipy.stats.lognorm.pdf(t2, s=am2, scale=np.exp(bm2),loc=dm1+alpha)
            
            # plt.plot(t2,Marginalt)
            # plt.plot(t2,m1*DensMarginalt1+(1-m1)*DensMarginalt2)
            
            # plt.plot(TiempoReducido,MarginaltReduc)
            
            
            # plt.plot(t2[t2<dm1+alpha],Marginalt[t2<dm1+alpha])
            # plt.plot(t2[t2<dm1+alpha],m1*DensMarginalt1[t2<dm1+alpha])
            
            
            # plt.plot(t2[t2>dm1+alpha],Marginalt[t2>=dm1+alpha])
            # plt.plot(t2[t2>dm1+alpha],(1-m1)*DensMarginalt2[t2>=dm1+alpha])
            
            #############################################
            
                
            Marginalf1Reduc=np.ones(len(f2<dm1+alpha))
            for i in np.arange(len(f2)):
                Marginalf1Reduc[i]=integrate.cumtrapz(Reducido[i,TiempoReducido<dm1+alpha],TiempoReducido[TiempoReducido<dm1+alpha] , initial=0)[-1]
            
            Marginalf2Reduc=np.ones(len(f2>dm1+alpha))
            for i in np.arange(len(f2)):
                Marginalf2Reduc[i]=integrate.cumtrapz(Reducido[i,TiempoReducido>dm1+alpha],TiempoReducido[TiempoReducido>dm1+alpha] , initial=0)[-1]
            
            
            # plt.plot(f2,MarginalfReduc)
            # plt.plot(f2,Marginalf2Reduc)
            
            # plt.plot(f2,Marginalf2Reduc/Marginalf2Reduc[-1])
            # plt.plot(f2,Marginalf1Reduc/Marginalf1Reduc[-1])
            
            
            # Theta=np.array([9.64750012e-01, 3.84607773e+00, 1.94015608e-02, 9.93997099e-01,
            #  1.14207148e+00, 2.55522188e+01, 5.43143683e-01, 6.33069153e-01,
            #  2.90704382e+00, 9.39894005e-01, 7.99796340e-01, 4.08040326e+00,
            #  8.83946853e-02])
            
            
            
            Theta=np.array([1,1,0.001])    
            
            np.random.seed(1)
            TwEnergia = pytwalk.pytwalk( n=len(Theta), U=EnergiaF, Supp=EnSuppF)
            TwEnergia.Run( T=50000, x0=Theta+0.0001, xp0=Theta)
            Corrida=TwEnergia.Output
            
            # plt.plot(-Corrida[:,-1])
            FinalF= Corrida[-1,:-1]
            a1,b1,sigma=FinalF
            
            DensMarginalf=scipy.stats.lognorm.pdf(f2, s=a1,scale=np.exp(b1))/scipy.stats.lognorm.cdf(fmax, s=a1,scale=np.exp(b1))*(f2<fmax)
            
            # plt.plot(f2,Marginalf)
            # plt.plot(f2,DensMarginalf)
            
            
            #############################################
            
            
            
            
            am1,bm1,dm1,am2,bm2,alpha,m1,sigma1,sigma2=FinalT
            a1,b1,sigma=FinalF
            # Theta=np.array([am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a2,b2,a2,b2,0.4,0.7,sigma1,sigma2])
            
            
            # Theta=np.array([am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a1,b1,0.1,2,0.4,0.7,sigma1,sigma2])    
            Theta=np.array([am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a1,b1,0.1,1,0.4,0.7,sigma1,sigma2])    
            am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a21,b21,a22,b22,m2,b,sigma1,sigma2=Theta
        
            
            np.random.seed(1)
            TwEnergia = pytwalk.pytwalk( n=len(Theta), U=Energia, Supp=EnSupp)
            TwEnergia.Run( T=LenMCMC, x0=Theta+0.0001, xp0=Theta)
            Corrida=TwEnergia.Output
            
            # plt.plot(-Corrida[:,-1])
            Final= Corrida[-1,:-1]
            # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a2,b2,b,sigma1,sigma2=Final
            
            am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a21,b21,a22,b22,m2,b,sigma1,sigma2=Final
            
            DensMarginalt1=scipy.stats.lognorm.pdf(t2, s=am1, scale=np.exp(bm1),loc=dm1)
            DensMarginalt2=scipy.stats.lognorm.pdf(t2, s=am2, scale=np.exp(bm2),loc=dm1+alpha)
            # plt.plot(t2,Marginalt)
            # plt.plot(TiempoReducido,MarginaltReduc)
            # plt.plot(t2,m1*DensMarginalt1+(1-m1)*DensMarginalt2)
            
            
        
            DensMarginalf1=scipy.stats.lognorm.pdf(f2, s=a1,scale=np.exp(b1))/scipy.stats.lognorm.cdf(fmax, s=a1,scale=np.exp(b1))*(f2<fmax)
            DensMarginalf21=scipy.stats.lognorm.pdf(f2, s=a21,scale=np.exp(b21))/scipy.stats.lognorm.cdf(fmax, s=a21,scale=np.exp(b21))*(f2<fmax)
            DensMarginalf22=scipy.stats.lognorm.pdf(f2, s=a22,scale=np.exp(b22))/scipy.stats.lognorm.cdf(fmax, s=a22,scale=np.exp(b22))*(f2<fmax)
            # plt.plot(f2,Marginalf)
            # plt.plot(f2,DensMarginalf21*m2+DensMarginalf22*(1-m2))
            # plt.plot(f2,DensMarginalf21)
            # plt.plot(f2,DensMarginalf22)
        
            # plt.plot(f2,MarginalfReduc)
            # plt.plot(f2,DensMarginalf1)
            
            
            ##############
            
            
            # now = datetime.now()
            # current_time = now.strftime("%H:%M:%S")
            np.savetxt(r'./Corrida/Corrida'+ (Muestra["Base"][idsismo]).replace(".","")  +".csv",Corrida,delimiter=',', fmt=('%s'))
            os.remove(r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')
        except:
            print("No se pudo "+ Muestra["Base"][idsismo])
            os.remove(r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')            
            file1 = open("Problemas.txt", "a")  # append mode
            file1.write((Muestra["Base"][idsismo])+"\n")
            file1.close()
            
    # try:
        
    #     Corrida=genfromtxt(r'./Corrida/Corrida'+ (Muestra["Base"][idsismo]).replace(".","")  +".csv", delimiter=',')
    #     # Corrida=genfromtxt(r"./Imperial"+".csv", delimiter=',')
    #     MAP=np.argmin(Corrida[:,-1])
    #     # Final= Corrida[-1,:-1]
    #     Final= Corrida[MAP,:-1]
    #     # plt.plot(-Corrida[:,-1])
    
    #     # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,b,sigma=Final
    #     # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,b,sigma1,sigma2=Final
    #     # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a2,b2,b,sigma1,sigma2=Final
    #     am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a21,b21,a22,b22,m2,b,sigma1,sigma2=Final
    #     # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,b,sigma11,sigma12,sigma21,sigma22=Final
    #     b=-b
        
    #     # plt.plot(-Corrida[:,-1])

    #     Zhat=EvaluaZhat(t2)
        
    #     MarginalfZhat=np.ones(len(f2))
    #     for i in np.arange(len(f2)):
    #         MarginalfZhat[i]=integrate.cumtrapz(Zhat[i,:],t2 , initial=0)[-1]
       
    #     NormalizacionZhat=integrate.cumtrapz(MarginalfZhat,f2 , initial=0)[-1]
       
    #     Zhat=Zhat*Normalizacion/NormalizacionZhat
       
    #     PSDFEst=np.copy(Zhat)
    
       
    #     MarginalfZhat=np.ones(len(f2))
    #     for i in np.arange(len(f2)):
    #         MarginalfZhat[i]=integrate.cumtrapz(Zhat[i,:],t2 , initial=0)[-1]
       
       
    #     MarginaltZhat=np.ones(len(t2))
    #     for i in np.arange(len(t2)):
    #         MarginaltZhat[i]=integrate.cumtrapz(Zhat[:,i],f2 , initial=0)[-1]
       
                
    #     np.random.seed(1)
    #     Proceso,Simulado=RespuestaSimulada(PSDFEst)
        
        
        ## tst = Acelerogram( "./RegistroChristen.csv", Dt=-1, skiprows=0)
        # tst = Acelerogram( r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv', Dt=-1, skiprows=0)
            
    #     Real=np.array(Respuesta(tst))
        
            
    #     for i in range(Trayectorias):
    #         print(i+1)
    #         Proceso,Nuevo=RespuestaSimulada(PSDFEst)
    #         Simulado=np.vstack((Simulado,Nuevo))    
    #     Simulado=np.array(Simulado)
        
    #     ProcesoReal,SimuladoReal=RespuestaSimulada(PSDFSuavizada)
        
    #     for i in range(Trayectorias):
    #         print(i+1)
    #         ProcesoReal,Nuevo=RespuestaSimulada(PSDFSuavizada)
    #         SimuladoReal=np.vstack((SimuladoReal,Nuevo))    
    #     SimuladoReal=np.array(SimuladoReal)
        
        
                
        
    #     #############################
    #     import matplotlib.pyplot as plt
    #     import numpy as np
          
        
    #     plt.rcParams.update({'font.size': 26})
        
    #     fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(36, 20))
        
    #     ax[0,0].plot(X,NS,alpha=0.5,label="Real")
    #     ax[0,0].plot(t2,Proceso,alpha=0.5,label="Fitted")
    #     ax[0,0].set_title('Record '+Muestra["Base"][idsismo]+" $,M_w$ "+str(Muestra["GMW"][idsismo])+", Distance " +str(round(Muestra["R"][idsismo],2)) +"km ")
    #     ax[0,0].set_xlabel("Time [sec]")
    #     ax[0,0].set_ylabel("Acceleration [$\\frac{cm}{s^2}$]")
    #     ax[0,0].yaxis.grid(color='gray', linestyle='dashed')
    #     ax[0,0].xaxis.grid(color='gray', linestyle='dashed')
    
    #     ax[0,0].legend()
        
        
    #     ax[0,1].pcolormesh(t[t<tmax], f[f<fmax], PSDF[f<fmax,:][:,t<tmax], shading='gouraud')
    #     ax[0,1].set_title('Evolutionary Power Spectral Density Function (ePSDF)')
    #     ax[0,1].set_ylabel('Frequency [Hz]')
    #     ax[0,1].set_xlabel('Time [sec]')
        
        
    #     ax[0,2].pcolormesh(t2[t2<tmax], f2[f2<fmax], PSDFSuavizada[f2<fmax,:][:,t2<tmax], shading='gouraud')
    #     ax[0,2].contour(t2[t2<tmax], f2[f2<fmax], gaussian_filter(Zhat[f2<fmax,:][:,t2<tmax]+np.finfo(float).eps, 1), 10, colors='w',origin="lower")
    #     ax[0,2].set_title('Fitted model over ePSDF')
    #     ax[0,2].set_ylabel('Frequency [Hz]')
    #     ax[0,2].set_xlabel('Time [sec]')
        
    #     # ax[1,0].plot(f2[f2<fmax],Marginalf[f2<fmax]*Normalizacion)
    #     # ax[1,0].plot(f2[f2<fmax],MarginalfZhat[f2<fmax])
    #     ax[1,0].plot(f2,Marginalf*Normalizacion)
    #     ax[1,0].plot(f2,MarginalfZhat)
    #     ax[1,0].set_axisbelow(True)
    #     ax[1,0].yaxis.grid(color='gray', linestyle='dashed')
    #     ax[1,0].xaxis.grid(color='gray', linestyle='dashed')
    #     # ax[1,0].plot(f2,scipy.stats.lognorm.pdf(f2, s=a1,scale=np.exp(b1)  ))
    #     # ax[0].set_title('Frequency Marginal density')
    #     ax[1,0].set_xlabel('Frequency [Hz]')
    #     ax[1,0].set_title('Frequency Marginal density')
    #     ax[1,0].set_xlabel('Frequency [Hz]')
        
        
    #     ax[1,1].plot(t2,Marginalt*Normalizacion)
    #     ax[1,1].plot(t2,MarginaltZhat)
    #     ax[1,1].yaxis.grid(color='gray', linestyle='dashed')
    #     ax[1,1].xaxis.grid(color='gray', linestyle='dashed')
    #     ax[1,1].set_title('Time Marginal Density')
    #     ax[1,1].set_xlabel('Time [sec]')
        
        
    #     for i in range(len(Simulado)):
    #             if i==0:
    #                 plt.scatter(Real[0],Simulado.T[:,i], alpha=0.3,label="Simulated Model",color="red")
    #                 plt.scatter(Real[0],SimuladoReal.T[:,i], alpha=0.3,label="Simulated",color="tab:purple")
    #             else :
    #                 plt.scatter(Real[0],Simulado.T[:,i], alpha=0.3,color="red")
    #                 plt.scatter(Real[0],SimuladoReal.T[:,i], alpha=0.3,color="tab:purple")
                    
                    
    #     ax[1,2].scatter(Real[0],Real[1],label="Real")
    #     ax[1,2].scatter(Real[0],np.median(Simulado,0),label="Median of simulations", color="black")
    #     ax[1,2].scatter(Real[0],np.median(SimuladoReal,0),label="Median of simulations", color="tab:cyan")        
    #     ax[1,2].legend()
    #     ax[1,2].set_title("Response spectra")
    #     ax[1,2].set_xlabel("Time [sec]")
    #     ax[1,2].set_ylabel("Acceleration [$\\frac{cm}{s^2}$]")
    #     ax[1,2].yaxis.grid(color='gray', linestyle='dashed')
    #     ax[1,2].xaxis.grid(color='gray', linestyle='dashed')
        
    #     fig.tight_layout()
    
    #     plt.show()  # doctest: +SKIP
        
        
    #     fig.savefig(r"./Graficas/"+(Muestra["Base"][idsismo]).replace(".","")+'.png')
    #     # fig.savefig("Imperial"+'.png')
        
    
      
       
        
    #     fig1, ax1 = plt.subplots()
    #     ax1.scatter(np.log10(Real[0]),np.log10(Real[1]/980.665),label="Real",s=100)

    #     for i in range(Trayectorias):
    #             if i==0:
    #                 plt.scatter(np.log10(Real[0]),np.log10(Simulado.T[:,i]/980.665), alpha=0.3,label="Simulated",color="red")
    #                 plt.scatter(np.log10(Real[0]),np.log10(SimuladoReal.T[:,i]/980.665), alpha=0.3,label="Simulated",color="tab:purple")

    #             if i!=0:
    #                 plt.scatter(np.log10(Real[0]),np.log10(Simulado.T[:,i]/980.665), alpha=0.3,color="red")
    #                 plt.scatter(np.log10(Real[0]),np.log10(SimuladoReal.T[:,i]/980.665), alpha=0.3,color="tab:purple")
                    
                    
    #     ax1.scatter(np.log10(Real[0]),np.log10(np.median(Simulado,0)/980.665),alpha=0.7,label="Median Model", color="black")
    #     ax1.scatter(np.log10(Real[0]),np.log10(np.median(SimuladoReal,0)/980.665),alpha=0.7,label="Median Real", color="tab:cyan")
    #     ax1.legend(prop={'size': 6})
    #     ax1.set_title("Response spectra")
    #     ax1.set_xlabel("Time [sec]")
    #     ax1.set_ylabel("Acceleration [$g$]")
    #     ax1.yaxis.grid(color='gray', linestyle='dashed')
    #     ax1.xaxis.grid(color='gray', linestyle='dashed')
    #     fig.tight_layout()
    #     plt.show()
        
    #     fig1.savefig(r"./Graficas/"+(Muestra["Base"][idsismo]).replace(".","")+"loglog"+'.png')
    #     #os.remove(r'./RegistroChristen.csv')
    #     os.remove(r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')
    #     os.remove(r'./RegistroChristen2'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')
    #     # fig1.savefig("ImperialLogLog"+'.png')
    # except:
    #     print("No se pudo graficar")
        # os.remove(r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')
        # os.remove(r'./RegistroChristen2'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')


        
        
        









