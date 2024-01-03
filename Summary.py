###Al final vienen separadas

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
TiempoOptimo=pd.read_csv("TiempoOptimo.csv", delimiter=",",header=None)
TiempoOptimo.set_axis({'Base', 'TiempoOptimo'}, axis=1, inplace=True)


Muestra

Muestra=pd.merge(Muestra, TiempoInicial, on="Base")
# Muestra=pd.merge(Muestra, TiempoOptimo, on="Base")

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



Intervalo=0 






Muestra[Muestra["Base"]=="UNIO0301.221"]
Muestra[Muestra["Base"]=="ATYC8509.191"]

s=65

s=140


# np.arange(Intervalo*len(Muestra.index.values[:])//10,(Intervalo+1)*len(Muestra.index.values[:])//10)
Muestra.index.values[5]
# for i in np.arange(Intervalo*len(Muestra.index.values[:])//10,(Intervalo+1)*len(Muestra.index.values[:])//10):
for s in Muestra.index.values[:]:
    print(Muestra["Base"][s])
    idsismo=np.copy(s)
    if os.path.isfile(r'./Corrida/Corrida'+ (Muestra["Base"][idsismo]).replace(".","")  +".csv"): # and (not os.path.isfile(r'./RSSimulado/'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')):
        try:    
            # print(Muestra["Base"][i])
            # idsismo=np.copy(i)
            CanalV=np.int(np.float(RecuperaOrientacion(Muestra["Base"][idsismo])[7]))
            CanalH1=np.int(np.float(RecuperaOrientacion(Muestra["Base"][idsismo])[8]))
            CanalH2=np.int(np.float(RecuperaOrientacion(Muestra["Base"][idsismo])[9]))
            
            (X,NS,f,Y,VelMuest)=fourier(Muestra["Base"][idsismo], CanalH1) #Usar registrocrudo cuando no quiero normalizado y centrado
            
            
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
            
            
    
        except:
            print("No se pudo "+ Muestra["Base"][idsismo])
            
            
    try:
        
        Corrida=genfromtxt(r'./Corrida/Corrida'+ (Muestra["Base"][idsismo]).replace(".","")  +".csv", delimiter=',')
        # Corrida=genfromtxt(r"./Imperial"+".csv", delimiter=',')
        MAP=np.argmin(Corrida[:,-1])
        # Final= Corrida[-1,:-1]
        Final= Corrida[MAP,:-1]
        # plt.plot(-Corrida[:,-1])
    
        # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,b,sigma=Final
        # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,b,sigma1,sigma2=Final
        # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a2,b2,b,sigma1,sigma2=Final
        am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,a21,b21,a22,b22,m2,b,sigma1,sigma2=Final
        # am1,bm1,dm1,am2,bm2,alpha,m1,a1,b1,b,sigma11,sigma12,sigma21,sigma22=Final
        b=-b
        
        # plt.plot(-Corrida[:,-1])

        Zhat=EvaluaZhat(t2)
        
        MarginalfZhat=np.ones(len(f2))
        for i in np.arange(len(f2)):
            MarginalfZhat[i]=integrate.cumtrapz(Zhat[i,:],t2 , initial=0)[-1]
       
        NormalizacionZhat=integrate.cumtrapz(MarginalfZhat,f2 , initial=0)[-1]
       
        Zhat=Zhat*Normalizacion/NormalizacionZhat
       
        PSDFEst=np.copy(Zhat)
    
       
        MarginalfZhat=np.ones(len(f2))
        for i in np.arange(len(f2)):
            MarginalfZhat[i]=integrate.cumtrapz(Zhat[i,:],t2 , initial=0)[-1]
       
       
        MarginaltZhat=np.ones(len(t2))
        for i in np.arange(len(t2)):
            MarginaltZhat[i]=integrate.cumtrapz(Zhat[:,i],f2 , initial=0)[-1]
       
                
        np.random.seed(1)
        Proceso,Simulado=RespuestaSimulada(PSDFEst)
        
        
        # tst = Acelerogram( "./RegistroChristen.csv", Dt=-1, skiprows=0)
        tst = Acelerogram( r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv', Dt=-1, skiprows=0)
            
        Real=np.array(Respuesta(tst))
        
            
        for i in range(Trayectorias):
            print(i+1)
            np.random.seed(i+1)
            Proceso,Nuevo=RespuestaSimulada(PSDFEst)
            Simulado=np.vstack((Simulado,Nuevo))    
        Simulado=np.array(Simulado)
        
        np.random.seed(1)
        ProcesoReal,SimuladoReal=RespuestaSimulada(PSDFSuavizada)
        
        for i in range(Trayectorias):
            print(i+1)
            np.random.seed(i+1)
            ProcesoReal,Nuevo=RespuestaSimulada(PSDFSuavizada)
            SimuladoReal=np.vstack((SimuladoReal,Nuevo))    
        SimuladoReal=np.array(SimuladoReal)
        
        np.savetxt(r'./RSSimulado/'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv',Simulado,delimiter=',', fmt=('%s'))
        np.savetxt(r'./RSReal/'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv',SimuladoReal,delimiter=',', fmt=('%s'))


        # Simulado=np.loadtxt(r'./RSSimulado/'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv',delimiter=",")
        # SimuladoReal=np.loadtxt(r'./RSreal/'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv',delimiter=",")
        #############################
        import matplotlib.pyplot as plt
        import numpy as np
          
        
        plt.rcParams.update({'font.size': 26})
        
        fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(36, 20))
        
        ax[0,0].plot(X,NS,alpha=0.5,label="Real")
        ax[0,0].plot(t2,Proceso,alpha=0.5,label="Simulated")
        ax[0,0].set_title('Record '+Muestra["Base"][idsismo]+" $,M_w$ "+str(Muestra["GMW"][idsismo])+", Distance " +str(round(Muestra["R"][idsismo],2)) +"km ")
        ax[0,0].set_xlabel("Time [s]")
        ax[0,0].set_ylabel("Acceleration [$\\frac{cm}{s^2}$]")
        ax[0,0].yaxis.grid(color='gray', linestyle='dashed')
        ax[0,0].xaxis.grid(color='gray', linestyle='dashed')
    
        ax[0,0].legend()
        
        
        ax[0,1].pcolormesh(t[t<tmax], f[f<fmax], PSDF[f<fmax,:][:,t<tmax], shading='gouraud')
        ax[0,1].set_title('Evolutionary Power Spectral Density Function (ePSDF)')
        ax[0,1].set_ylabel('Frequency [Hz]')
        ax[0,1].set_xlabel('Time [s]')
        
        
        ax[0,2].pcolormesh(t2[t2<tmax], f2[f2<fmax], PSDFSuavizada[f2<fmax,:][:,t2<tmax], shading='gouraud')
        ax[0,2].contour(t2[t2<tmax], f2[f2<fmax], gaussian_filter(Zhat[f2<fmax,:][:,t2<tmax]+np.finfo(float).eps, 1), 10, colors='w',origin="lower")
        ax[0,2].set_title('Fitted model over ePSDF')
        ax[0,2].set_ylabel('Frequency [Hz]')
        ax[0,2].set_xlabel('Time [s]')
        
        # ax[1,0].plot(f2[f2<fmax],Marginalf[f2<fmax]*Normalizacion)
        # ax[1,0].plot(f2[f2<fmax],MarginalfZhat[f2<fmax])
        ax[1,0].plot(f2,np.sqrt(Marginalf*Normalizacion))
        ax[1,0].plot(f2,np.sqrt(MarginalfZhat))
        ax[1,0].set_axisbelow(True)
        ax[1,0].yaxis.grid(color='gray', linestyle='dashed')
        ax[1,0].xaxis.grid(color='gray', linestyle='dashed')
        # ax[1,0].plot(f2,scipy.stats.lognorm.pdf(f2, s=a1,scale=np.exp(b1)  ))
        # ax[0].set_title('Frequency Marginal density')
        ax[1,0].set_xlabel('Frequency [Hz]')
        ax[1,0].set_title('Frequency Marginal density')
        ax[1,0].set_xlabel('Frequency [Hz]')
        
        
        ax[1,1].plot(t2,np.sqrt(Marginalt*Normalizacion))
        ax[1,1].plot(t2,np.sqrt(MarginaltZhat))
        ax[1,1].yaxis.grid(color='gray', linestyle='dashed')
        ax[1,1].xaxis.grid(color='gray', linestyle='dashed')
        ax[1,1].set_title('Time Marginal Density')
        ax[1,1].set_xlabel('Time [s]')
        
        
        for i in range(len(Simulado)):
                if i==0:
                    plt.scatter(Real[0],Simulado.T[:,i], alpha=0.3,label="Simulated Model",color="red")
                    # plt.scatter(Real[0],SimuladoReal.T[:,i], alpha=0.3,label="Simulated",color="tab:purple")
                else :
                    plt.scatter(Real[0],Simulado.T[:,i], alpha=0.3,color="red")
                    # plt.scatter(Real[0],SimuladoReal.T[:,i], alpha=0.3,color="tab:purple")
                    
                    
        ax[1,2].scatter(Real[0],Real[1],label="Real")
        ax[1,2].scatter(Real[0],np.median(Simulado,0),label="Median of simulations", color="black")
        # ax[1,2].scatter(Real[0],np.median(SimuladoReal,0),label="Median of simulations", color="tab:cyan")        
        ax[1,2].legend()
        ax[1,2].set_title("Response spectra")
        ax[1,2].set_xlabel("Time [s]")
        ax[1,2].set_ylabel("Acceleration [$\\frac{cm}{s^2}$]")
        ax[1,2].yaxis.grid(color='gray', linestyle='dashed')
        ax[1,2].xaxis.grid(color='gray', linestyle='dashed')
        
        fig.tight_layout()
    
        plt.show()  # doctest: +SKIP
        
        
        fig.savefig(r"./Graficas/"+(Muestra["Base"][idsismo]).replace(".","")+'.png')
        # fig.savefig("Imperial"+'.png')
        
    
      
       
        
        fig1, ax1 = plt.subplots()
        ax1.scatter(np.log10(Real[0]),np.log10(Real[1]/980.665),label="Real",s=100)

        for i in range(Trayectorias):
                if i==0:
                    plt.scatter(np.log10(Real[0]),np.log10(Simulado.T[:,i]/980.665), alpha=0.3,label="Simulated",color="red")
                    plt.scatter(np.log10(Real[0]),np.log10(SimuladoReal.T[:,i]/980.665), alpha=0.3,label="Simulated",color="tab:purple")

                if i!=0:
                    plt.scatter(np.log10(Real[0]),np.log10(Simulado.T[:,i]/980.665), alpha=0.3,color="red")
                    plt.scatter(np.log10(Real[0]),np.log10(SimuladoReal.T[:,i]/980.665), alpha=0.3,color="tab:purple")
                    
                    
        ax1.scatter(np.log10(Real[0]),np.log10(np.median(Simulado,0)/980.665),alpha=0.7,label="Median Model", color="black")
        ax1.scatter(np.log10(Real[0]),np.log10(np.median(SimuladoReal,0)/980.665),alpha=0.7,label="Median Real", color="tab:cyan")
        ax1.legend(prop={'size': 6})
        ax1.set_title("Response spectra")
        ax1.set_xlabel("Time [s]")
        ax1.set_ylabel("Acceleration [$g$]")
        ax1.yaxis.grid(color='gray', linestyle='dashed')
        ax1.xaxis.grid(color='gray', linestyle='dashed')
        fig.tight_layout()
        plt.show()
        
        fig1.savefig(r"./Graficas/"+(Muestra["Base"][idsismo]).replace(".","")+"loglog"+'.png')
        #os.remove(r'./RegistroChristen.csv')
        os.remove(r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')
        os.remove(r'./RegistroChristen2'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')
        # fig1.savefig("ImperialLogLog"+'.png')
    except:
        
        try:
            os.remove(r'./RegistroChristen'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')
            os.remove(r'./RegistroChristen2'+(Muestra["Base"][idsismo]).replace(".","")+ '.csv')
            print("No se pudo graficar")
        except:
            print("No se pudo graficar")








################### Separadas

Fuente=20

plt.rcParams.update({'font.size': Fuente})



plt.figure()
plt.plot(X,NS,alpha=0.5,label="Observed")
plt.plot(t2,Proceso,alpha=0.5,label="Simulated")
plt.title(Muestra["Base"][idsismo], fontsize=Fuente)
plt.xlabel("Time [s]")
plt.ylabel("Acceleration [$\\frac{cm}{s^2}$]")
plt.grid(color='gray', linestyle='dashed')
plt.legend( prop={'size': 15})
#plt.tight_layout()
plt.savefig("/Users/isaias/Desktop/Articulo/"+Muestra["Base"][idsismo][:4]+"6"+"1"+".pdf",dpi=300, bbox_inches = "tight")



#######################
plt.figure()
plt.pcolormesh(t[t<tmax], f[f<fmax], Z[f<fmax,:][:,t<tmax], shading='gouraud', rasterized=True)
plt.colorbar()
# plt.title('ePSDF', fontsize=Fuente)
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [s]')
#plt.tight_layout()
plt.savefig("/Users/isaias/Desktop/"+Muestra["Base"][idsismo][:4]+"6"+"2"+".pdf",dpi=300, bbox_inches = "tight")
######################
plt.figure()
plt.pcolormesh(t2[t2<tmax], f2[f2<fmax], PSDFSuavizada[f2<fmax,:][:,t2<tmax], shading='gouraud', rasterized=True)
plt.contour(t2[t2<tmax], f2[f2<fmax], gaussian_filter(Zhat[f2<fmax,:][:,t2<tmax]+np.finfo(float).eps, 1), 10, colors='w',origin="lower")
plt.title('Fitted model over ePSDF', fontsize=Fuente)
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [s]')
#plt.tight_layout()
plt.savefig("/Users/isaias/Desktop/Articulo/"+Muestra["Base"][idsismo][:4]+"6"+"3"+".pdf",dpi=300, bbox_inches = "tight")


####################

# ax[1,0].plot(f2[f2<fmax],Marginalf[f2<fmax]*Normalizacion)
# ax[1,0].plot(f2[f2<fmax],MarginalfZhat[f2<fmax])
plt.figure()
plt.plot(f2[f2<25],np.sqrt(Marginalf[f2<25]*Normalizacion))
plt.plot(f2[f2<25],np.sqrt(MarginalfZhat[f2<25]))
plt.grid(color='gray', linestyle='dashed')

plt.grid(color='gray', linestyle='dashed')
plt.grid(color='gray', linestyle='dashed')
# ax[1,0].plot(f2,scipy.stats.lognorm.pdf(f2, s=a1,scale=np.exp(b1)  ))
# ax[0].set_title('Frequency Marginal density')
plt.xlabel('Frequency [Hz]')
plt.title('Frequency Marginal density', fontsize=Fuente)
plt.xlabel('Frequency [Hz]')
#plt.tight_layout()
plt.savefig("/Users/isaias/Desktop/Articulo/"+Muestra["Base"][idsismo][:4]+"6"+"4"+".pdf",dpi=300, bbox_inches = "tight")


#################
plt.figure()
plt.plot(t2,np.sqrt(Marginalt*Normalizacion))
plt.plot(t2,np.sqrt(MarginaltZhat))
plt.grid(color='gray', linestyle='dashed')
plt.title('Time Marginal Density', fontsize=Fuente)
plt.xlabel('Time [s]')
#plt.tight_layout()
plt.savefig("/Users/isaias/Desktop/Articulo/"+Muestra["Base"][idsismo][:4]+"6"+"5"+".pdf",dpi=300, bbox_inches = "tight")

##################
plt.figure()
for i in range(len(Simulado)):
        if i==0:
            plt.scatter(Real[0],Simulado.T[:,i], alpha=0.3,label="Simulated Model",color="red")
            # plt.scatter(Real[0],SimuladoReal.T[:,i], alpha=0.3,label="Simulated",color="tab:purple")
        else :
            plt.scatter(Real[0],Simulado.T[:,i], alpha=0.3,color="red")
            # plt.scatter(Real[0],SimuladoReal.T[:,i], alpha=0.3,color="tab:purple")
            
            
plt.scatter(Real[0],Real[1],label="Observed",s=70,alpha=0.7)
plt.scatter(Real[0],np.median(Simulado,0),label="Median of simulations", color="black",s=70,alpha=0.7)
plt.xticks(Real[0])
plt.yticks(np.round(np.arange(6)/5*np.ceil(np.max((np.max(Real[1]),np.max(Simulado))))))
# ax[1,2].scatter(Real[0],np.median(SimuladoReal,0),label="Median of simulations", color="tab:cyan")        
plt.legend( prop={'size': 15})
plt.title("Response spectra", fontsize=Fuente)
plt.xlabel("Time [s]")
plt.ylabel("Acceleration [$\\frac{cm}{s^2}$]")
plt.grid(color='gray', linestyle='dashed')
#plt.tight_layout()
plt.savefig("/Users/isaias/Desktop/Articulo/"+Muestra["Base"][idsismo][:4]+"6"+"6"+".pdf",dpi=300, bbox_inches = "tight")

# plt.tight_layout()

# plt.show()  # doctest: +SKIP





#####################################################################


plt.figure()
plt.plot(X,NS,alpha=0.5,label="Real")
plt.xlabel("Time [s]")
plt.ylabel("Acceleration [$\\frac{cm}{s^2}$]")
plt.grid(color='gray', linestyle='dashed')
#plt.legend()
#plt.tight_layout()
plt.savefig("/Users/isaias/Desktop/Articulo/"+"Acelerograma"+Muestra["Base"][idsismo][:4]+".pdf",dpi=300, bbox_inches = "tight")



plt.figure()
cmap=plt.pcolormesh(t[t<tmax], f[f<fmax], Z[f<fmax,:][:,t<tmax], shading='gouraud', rasterized=True)
plt.colorbar(cmap)
plt.xlabel("Time [s]")
plt.ylabel("Frequency [$Hz$]")
plt.savefig("/Users/isaias/Desktop/Articulo/"+"STFTBarra"+(Muestra["Base"][idsismo]).replace(".","")+'.pdf',dpi=300, bbox_inches = "tight")

# plt.tight_layout()





######################################################################



np.random.seed(1)
Pro1, Res=RespuestaSimulada(PSDFSuavizada)
np.random.seed(2)
Pro2, Res=RespuestaSimulada(PSDFSuavizada)
np.random.seed(3)
Pro3, Res=RespuestaSimulada(PSDFSuavizada)



plt.figure()
plt.plot(X,NS,linewidth=.5, label="Simulated X(t)", color="tab:blue")
plt.grid(color='gray', linestyle='dashed')
plt.ylabel("Acceleration"+r"$[\frac{cm}{s^2}]$")
plt.savefig("/Users/isaias/Desktop/Articulo/"+"Simulado0"+(Muestra["Base"][idsismo])[:4]+'.pdf',dpi=300, bbox_inches = "tight")


plt.figure()
plt.plot(X,Pro1,linewidth=.5, label="Simulated X(t)", color="tab:orange")
plt.grid(color='gray', linestyle='dashed')
plt.ylabel("Acceleration"+r"$[\frac{cm}{s^2}]$")
plt.savefig("/Users/isaias/Desktop/Articulo/"+"Simulado1"+(Muestra["Base"][idsismo])[:4]+'.pdf',dpi=300, bbox_inches = "tight")

plt.figure()
plt.plot(X,Pro2,linewidth=.5, label="Simulated X(t)", color="tab:orange")
plt.ylabel("Acceleration"+r"$[\frac{cm}{s^2}]$")
plt.grid(color='gray', linestyle='dashed')
plt.xlabel("Time [$s$]")
plt.savefig("/Users/isaias/Desktop/Articulo/"+"Simulado2"+(Muestra["Base"][idsismo])[:4]+'.pdf',dpi=300, bbox_inches = "tight")



plt.figure()
plt.plot(X,Pro3,linewidth=.5, label="Simulated X(t)", color="tab:orange")
plt.grid(color='gray', linestyle='dashed')
plt.xlabel("Time [$s$]")
plt.ylabel("Acceleration"+r"$[\frac{cm}{s^2}]$")
plt.savefig("/Users/isaias/Desktop/Articulo/"+"Simulado3"+(Muestra["Base"][idsismo])[:4]+'.pdf',dpi=300, bbox_inches = "tight")






























