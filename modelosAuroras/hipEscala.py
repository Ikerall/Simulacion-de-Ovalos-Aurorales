import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from matplotlib.patches import Patch
from utilidades import proyeccionOrtografica as po
from utilidades import graficarOvalo as go
from utilidades import leerOMNI
from starkov import ovalos
from utilidades import zonaOvalo
#%%
# Funcion que calcula el cambio de latitudes por la hipotesis de escala en función de la presion
def escalar(lat1, P1, P2):
    
    # Las latitudes de entrada son en grados, se las pasa a radianes

    
    lat1 = np.deg2rad(lat1)
    
    # esta expresion devuelve el cos de la latitud nueva,z=cos(lat2)
    z = np.cos(lat1) * (P2/P1)**( 1/(12) )
    
    lat2 = np.rad2deg( np.arccos(z) )


    return lat2


def main():
    
    for polo in ["norte","sur"]:
        
        # Valores del ovalo antes del escalado
        i = 5
        fecha, ut, parametrosEsc, kps = leerOMNI("2025_01_01.lst", i)
        
        dt = datetime.fromtimestamp(ut, tz=timezone.utc)
        horaUT = dt.strftime("%H:%M:%S")
        fechaTor = fecha.strftime("%Y/%m/%d")
        
        kpf = 2 # kpf no real solo para comparar
        ptsPo = ovalos(0, kpf, fecha, polo)
        ptsEq = ovalos(1, kpf, fecha, polo)


        
        # a la función de escalar se le pasa los pontos de los limites del 
        # ovalo del otro modelo, la presion para el kp aproximado del otro modelo
        # y la presión con la que se va a escalar el ovalo
        
        P1 = 0.50 # presion solo para comprar aprox a kp anterior
        P2 = parametrosEsc[0]
        
        latPoEsc = escalar(ptsPo[0], P1, P2)    
        latEqEsc = escalar(ptsEq[0], P1, P2)    
                
           
        # Definir las latitudes y longitudes de los ovalos escalados
        if polo == "sur":
            ptsPoEsc = [-latPoEsc, ptsPo[1]]
            ptsEqEsc = [-latEqEsc, ptsEq[1]]
        else:          
            ptsPoEsc = [latPoEsc, ptsPo[1]]
            ptsEqEsc = [latEqEsc, ptsEq[1]]
        
        df = pd.DataFrame({
            "latPolar"   : ptsPo[0],
            "latEscPolar": ptsPoEsc[0],
            "latEscuador": ptsEq[0],
            "latEscEcuad": ptsEqEsc[0]
            })
        
        print(df)
        # graficarOvalo(ax, puntos, tipo, color="red", linestyle="-",label=None):

        fig1 = plt.figure(figsize=(7, 8), dpi=150)
        ax = po(fig1, polo, [45, 90], fecha, kps[1])
        
        zonaOvalo(ax, ptsPoEsc, ptsEqEsc, color="red")
        zonaOvalo(ax, ptsPo, ptsEq, color = "#0F2BE2")      
 
        
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
            x=0.25, y=0.09,          # posición relativa en axes (0 a 1)
            s=f"$ Simulación \; con \; la\;hipotesis \;de \;escala$\n $para \;la\;Tormenta \;del\; {fech}$",
            transform=ax.transAxes,   # coordenadas relativas al axes
            fontsize=16,
            weight=700,
            verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",edgecolor="k", alpha=1)
        )
        
        #Label del ovalo azul 
        ovalo2 = Patch(facecolor='#0F2BE2', edgecolor='#6F1EC0', label=f'Kp = {kpf}, Pd = {P1}')
        
        #label del ovalo rojo
        ovalo1 = Patch(facecolor='red', edgecolor='#6F1EC0', label=f'Kp = {kps[1]}, Pd ={P2}')


        handles, labels = ax.get_legend_handles_labels()
        handles.append(ovalo2)
        labels.append(f'$Kp = {kpf},\; P_d = {P1}\;nPa$')
        handles.append(ovalo1)
        labels.append(f'$Kp = {kps[1]},\; P_d ={P2}\;nPa$')
        
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


if __name__ == "__main__":
    
    main()