import numpy as np
from geopack import geopack
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from matplotlib.patches import Patch
from matplotlib.patches import Wedge, Circle
import time
from utilidades import zonaOvalo
from utilidades import leerOMNI
from utilidades import proyeccionOrtografica

#%%

# Configuración inicial

# modelo T89: parametros = kp
# modelo T96: parametros = np.array([Pdyn, Dst, ByIMF, BzIMF, 5, 6, 7, 8, 9, 10])
# modelo T01: parametros = np.array([Pdyn, Dst, ByIMF, BzIMF, G1, G2, 7, 8, 9, 10])
# modelo T03: parametros = np.array([Pdyn, Dst, ByIMF, BzIMF, W1, W2, W3, W4, W5, W6])

def contornos(RS, RC, step):
    
    # 1. Definir contornos ecuatoriales: son las curvas alrededor de la tierra, por las cuales pasan las lineas
    # de campo magnetico que se considerar que desvian a las particulas que dan lugar a las auroras
    # la función devuelve 3 listas: 
    #   R_interno = distancias al centro de la tierra de cada punto del ovalo interno en Radios terrestres(Re)
    #   R_externo = distancias al centro de la tierra de cada punto del ovalo externo en Radios terrestres(Re
    #   phi       = longitudes GSM en radianes
    
    RS, RC = radios(RS, RC)
    phi = np.linspace(0, 2*np.pi, step)  # Longitud GSM en radianes
    RS1, RS2 = RS[0], RS[1]  # Radios subsolares (Re)
    RC1, RC2 = RC[0], RC[1]   # Radios en la cola (Re)
    R_interno = RS1 + (RC1 - RS1) * np.sin(phi/2)**2
    R_externo = RS2 + (RC2 - RS2) * np.sin(phi/2)**2
    
    return R_interno, R_externo, phi

def radios(RS=[10, 10.5], RC =[6, 12.0]):
    
      # función que permite modificar los radios que se utilizan en los contornos (no muy util por separado de contornos) 
  
      RS = [RS[0], RS[1]]  # Radios subsolares (Re)
      RC = [RC[0], RC[1]]   # Radios en la cola (Re)
      return RS, RC

def seguirLineasGSM(modelo, polo, ovalo, parametros, ut, RS, RC, step=60):
    
    # 2. Función para trazar líneas de campo hasta la ionosfera:
    
    # Llamo a la función contornos, obtengo los puntos de los ovalos
    R_interno, R_externo, phi = contornos(RS, RC, step)
    
    # Selecciono cual contorno utilizar para la simulación
    if ovalo == 'ext':
        x = R_externo*np.cos(phi)
        y = R_externo*np.sin(phi)
    elif ovalo == 'int':
        x = R_interno*np.cos(phi)
        y = R_interno*np.sin(phi)
    else:
        return print('Parametro de ovalo no valido')
    
    z = 0 # Los contornos se definen el ecuador z=0
    
    # Seleccióno la dirección del siguimiento de las lineas de campo
    if polo == 'sur':
        dir = 1
    elif polo == 'norte':
        dir = -1
    else:
        dir =-1            # default: norte
    
    
    # polo norte dipolar (lat, lon en grados)
    # parametros físicos para la simulación
    R0 = 1.02    # distancia del centro de la tierra, donde el seguimiento de linea tiene que terminar
    RL = 100     # radio limite de la simulación, si al seguir una linea se supuera, el codigo termina
    XF = []
    YF = []
    ZF = []
    
    # para cada valor de angulo en phi, se llama a la función trace
    # las coordenadas iniciales son la de el punto del contorno correspondiente
    # a ese valor de angulo, los coordenadas finales de cada seguimiento de linea 
    # se guardan en 3 listas
    for i in range(len(phi)):
        xf, yf, zf, xx, yy, zz = geopack.trace(
        x[i], y[i], z, #coordenadas iniciales en gsm
        dir=dir, 
        parmod = parametros, 
        r0     = R0, 
        exname = modelo, 
        inname = 'igrf', 
        rlim   = RL
        )

    #Elimina puntos que quedaron fuera del limite de 10 radios terrestres  
        R = np.sqrt(xf**2 + yf**2 + zf**2)
        if R < 5: 
            XF.append(xf)
            YF.append(yf)
            ZF.append(zf)
        else:
            continue
        
    return np.array(XF), np.array(YF), np.array(ZF)

def coord(XF,YF,ZF):
    # Cambio de coordenadas a geocentricas
    Xgeo, Ygeo, Zgeo = geopack.geogsm(XF, YF, ZF, -1)
    r = np.sqrt(Xgeo**2 + Ygeo**2 + Zgeo**2)
    lat = np.arcsin(Zgeo / r)
    lon = np.arctan2(Ygeo, Xgeo)
    
    return  np.degrees(lat), np.degrees(lon), r

def seguirLineasGEO(modelo, polo, ovalo, kp, fecha, RS, RC, step):

    
    # una función que combina las dos anteriores simplemente por comodidad
    
    XF,YF,ZF = seguirLineasGSM(modelo, polo, ovalo, kp, fecha, RS, RC, step)
    
    return coord(XF, YF, ZF)

def dual_half_circle(center=(0, 0), radius=1, angle=90, ax=None, colors=('#FDD922', '#0F1E84', 'white'),
                     **kwargs):
    
    # crea una figura con un circulo central para el grafico lateral
    # es una modificación de un codigo dentro de un ejemplo del modulo geopack
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    # w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], **kwargs)
    # w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], **kwargs)

    w1 = Wedge(center, radius, theta1, theta2, fc=colors[1], **kwargs, zorder=3, label="Lado noche")
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[0], **kwargs, zorder=3, label="Lado día")

    cr = Circle(center, radius, fc=colors[2], fill=False, **kwargs, zorder = 3)
    for wedge in [w1, w2, cr]:
        ax.add_artist(wedge)
    return [w1, w2, cr]

def setup_fig(xlim=(15, -20), ylim=(15, -15), xlabel='$X_{GSM} [Re]$', ylabel='$Y_{GSM}[Re]$'):
    
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Times New Roman'
    plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
    plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)
    # ax.axvline(0, ls=':', color='k', zorder=0)
    # ax.axhline(0, ls=':', color='k', zorder=0)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_aspect('equal')
    w1, w2, cr = dual_half_circle(ax=ax)
    
    #Grosor de los bordes
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    #Grosor de los ticks
    plt.tick_params(width=1.5, length=5,labelsize=15, top=True, right=True, direction="in",color="w") 
    ax.set_facecolor('k')
    plt.tight_layout()
    return ax
      
def graficaLateral():
       

    ax = setup_fig(xlim=(15, -30), ylim=(-15, 15), xlabel='$X_{GSM} [Re]$', ylabel='$Z_{GSM} [Re]$')


    i = 15
    fecha, ut, parametros, kps = leerOMNI("2025_01_01.lst", i)
    j = 0
    _, _, parametrosComp, _ = leerOMNI("2025_01_01.lst", j)
                                                  
    # Parametros del campo
    dip = geopack.recalc(ut)
    kp = parametros
    CampoExterno = 't96'
    CampoInterno = 'igrf'
    RS = [8.5, 9.5]
    RC = [5, 30]


    #%%
    # para graficar correctamente la magnetosfera se hace por separado las zonas
    # de lineas de campo cerradas y abiertas.

    #Lineas de ovalo lado dia
    for i in np.linspace(RS[0],RS[1], 10):
        
            xfOS1, yfOS1, zfOS1, xOS1, yOS1, zOS1 = geopack.trace(
                i, 0, 0, dir=1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=100)
            xfOS2, yfOS2, zfOS2, xOS2, yOS2, zOS2 = geopack.trace(
                i, 0, 0, dir=-1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=100)
            
            if i == RS[0]:
                
                ax.plot(xOS1, zOS1, color="#00F964", label="Zona Auroral")
                ax.plot(xOS2, zOS2, color="#00F964")
            else:
                ax.plot(xOS1, zOS1, color="#00F964")
                ax.plot(xOS2, zOS2, color="#00F964")
            
    # Lineas de campo cerradas dia 
    for i in np.linspace(1, 9.6, 10):   
        
        xf1, yf1, zf1, xx1, yy1, zz1 = geopack.trace(
            i, 0, 0, dir=1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=50)
        xf2, yf2, zf2, xx2, yy2, zz2 = geopack.trace(
            i, 0, 0, dir=-1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=50) 
        
        if RS[0]< i < RS[1]:
            
            ax.plot(xx1, zz1, color="#00F964")
            ax.plot(xx2, zz2, color="#00F964")
            
        else:
            ax.plot(xx1, zz1, color="#2B76C1")
            ax.plot(xx2, zz2, color="#2B76C1")   

    # lineas de campo cerradas noche 
    for i in np.linspace(1.02, 5, 15):

        xf3, yf3, zf3, xx3, yy3, zz3 = geopack.trace(
            -i, 0, 0, dir=1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=50)
        xf4, yf4, zf4, xx4, yy4, zz4 = geopack.trace(
            -i, 0, 0, dir=-1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=50)
        
        if RC[0] < i <RC[1]:
             ax.plot(xx3, zz3, color="#00F964")
             ax.plot(xx4, zz4, color="#00F964")
             
        else:
            ax.plot(xx3, zz3, color="#2B76C1")
            ax.plot(xx4, zz4, color="#2B76C1")

    for i in np.linspace(10, 30, 15):

        xf3, yf3, zf3, xx3, yy3, zz3 = geopack.trace(
            -i, 0, 0, dir=1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=50)
        xf4, yf4, zf4, xx4, yy4, zz4 = geopack.trace(
            -i, 0, 0, dir=-1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=50)
        
        if RC[0] < i <RC[1]:
             ax.plot(xx3, zz3, color="#00F964")
             ax.plot(xx4, zz4, color="#00F964")
             
        else:
            ax.plot(xx3, zz3, color="#2B76C1")
            ax.plot(xx4, zz4, color="#2B76C1")
            
     
    # lineas de campo abiertas desde z=0 hasta z=20       
    for i in np.linspace(-1, 25, 30):
        
        xf5, yf5, zf5, xx5, yy5, zz5 = geopack.trace(
            -30, 0, i, dir=-1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=50)
        
        ax.plot(xx5, zz5, color="#2B76C1")
        
    #lineas de campo abiertas desde z=-5 hasta -20
    for i in np.linspace(1, 25, 30):
        
        xf6, yf6, zf6, xx6, yy6, zz6 = geopack.trace(
            -30, 0, -i, dir=1.0, parmod=kp, exname=CampoExterno, inname='igrf', rlim=50)
        
        ax.plot(xx6, zz6, color="#2B76C1")
        
    dt = datetime.fromtimestamp(ut, tz=timezone.utc)
    horaUT = dt.strftime("%H:%M:%S")
    fechaTor = fecha.strftime("%Y/%m/%d")

    ax.text(
        x=0.025, y=0.25,          # posición relativa en axes (0 a 1)
        s=f"$UT:{horaUT}$ \n $P_d={kp[0]}\,nPa$ \n $Dst={kp[1]}\,nT$ \n $By={kp[2]}\,nT$ \n $Bz={kp[3]}\,nT$",
        transform=ax.transAxes,   # coordenadas relativas al axes
        fontsize=25,
        verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",edgecolor="grey", alpha=1)
    )
    ax.text(
        x=0.50, y=0.11,          # posición relativa en axes (0 a 1)
        s=f"$R_S: \; int={RS[0]} \; ; \; ext={RS[1]} $ \n $R_C: \; int={RC[0]} \; ; \; ext={RC[1]} $",
        transform=ax.transAxes,   # coordenadas relativas al axes
        fontsize=25,
        verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",edgecolor="grey", alpha=1)
    )
    ax.text(
        x=0.025, y=0.97,          # posición relativa en axes (0 a 1)
        s=f"$Simulación \;de \;la \;magnetosfera \;para \;la \; tormenta \; del \;{fechaTor}\; con \;los \;modelos \;T96 \;e \;IGRF$",
        transform=ax.transAxes,   # coordenadas relativas al axes
        fontsize=25,
        verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",edgecolor="grey", alpha=1)
    )
    ax.legend(
        fontsize = 20,           # tamaño del texto del label
        facecolor="white",       # color de fondo
        edgecolor="grey",     # borde de la leyenda
        labelcolor="black",     # color del texto de los labels
        loc="lower right",     # ubicación
        framealpha=1         # transparencia del fondo
    )
    plt.show(ax)

    end = time.time()        
   

def main():
        
    for polo in ["norte", "sur"]:
        
        # llamar a leerOMNI para obtener parametros de tormenta
        hora = 15
        fecha, ut, parametros, kps = leerOMNI("2025_01_01.lst", hora)
        
        # llamo a recalc para actualizar los parametros de la simulación de campo
        dipAngle = geopack.recalc(ut)  
        dt = datetime.fromtimestamp(ut, tz=timezone.utc)
        horaUT = dt.strftime("%H:%M:%S")
        fechaTor = fecha.strftime("%Y/%m/%d")
        
        # seleccionar radios de contorno
        RS = [8.5, 9.5]
        RC = [5, 30]
        
        # llamar a seguirlineas para seguir las lineas del campo hasta la atmosfera
        puntosExt = seguirLineasGEO("t96", polo, "ext", parametros, ut, RS, RC, 100)
        puntosInt = seguirLineasGEO("t96", polo, "int", parametros, ut, RS, RC, 100)
        
        # Segunda simulación

        hora2 = 0
        _, _, parametrosComp, kps2 = leerOMNI("2025_01_01.lst", hora2)
        
        RS = [8, 9]
        RC = [6, 15]      
        
        puntosExt2 = seguirLineasGEO("t96", polo, "ext", parametrosComp, ut, RS, RC, 100)
        puntosInt2 = seguirLineasGEO("t96", polo, "int", parametrosComp, ut, RS, RC, 100)
        
        # Crear figura/axes con tamaño fijo
        fig = plt.figure(figsize=(7, 8), dpi=150)
    
        # Llamar a proyeccionOrtografica para obtener axes con proyección
        ax = proyeccionOrtografica(fig, polo, [43, 90], fecha, kps[1])
        
        if polo == "norte":
            zonaOvalo(ax, puntosExt, puntosInt)
            zonaOvalo(ax, puntosExt2, puntosInt2, color="orange")
        elif polo == "sur":    
            zonaOvalo(ax, puntosInt, puntosExt)
            zonaOvalo(ax, puntosInt2, puntosExt2, color="orange")
        else:
            raise NameError
            
        # figure config
        ax.text(
             x=0.025, y=0.98,          # posición relativa en axes (0 a 1)
             s=f"$UT={horaUT}$",
             transform=ax.transAxes,   # coordenadas relativas al axes
             fontsize=20,
             verticalalignment='top',
             bbox=dict(boxstyle="round,pad=0.2", facecolor="white",edgecolor="k", alpha=1)
         )
    
        fech = fecha.strftime("%Y/%m/%d")
        ax.text(
             x=0.2, y=0.1,          # posición relativa en axes (0 a 1)
             s=f"$ Simulación \; con \; el\;modelo \;de \;Tsyganenko\;T96$ \n $para \;la \; Tormenta \;del\; {fech}$",
             transform=ax.transAxes,   # coordenadas relativas al axes
             fontsize=16,
             weight=700,
             zorder=4,
             verticalalignment='top',
             bbox=dict(boxstyle="round,pad=0.2", facecolor="white",edgecolor="k", alpha=1)
         )
        ovalo1 = Patch(facecolor='#19DA40', edgecolor='#6F1EC0', label=f'Kp = {kps[1]}')
        ovalo2 = Patch(facecolor='orange', edgecolor='#6F1EC0', label=f'Kp = {kps2[1]}')


        handles, labels = ax.get_legend_handles_labels()
        handles.append(ovalo2)
        labels.append(f'$Kp = {kps2[1]}$')
        handles.append(ovalo1)
        labels.append(f'$Kp = {kps[1]}$')
        
        # Mostrar leyenda combinada
        ax.legend(handles=handles, 
                  labels=labels, 
                  loc='upper right',
                  fontsize=15,
                  framealpha=1,
                  edgecolor="k",
                  alignment="right"
                  )
        
        plt.show()
        

    return 0


if __name__ == '__main__':
    
    start = time.time()
    
    main()
    
    end = time.time()
    print(f"Tiempo de ejecución: {end - start:.3f} segundos")