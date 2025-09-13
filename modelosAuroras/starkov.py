import matplotlib.pyplot as plt
from datetime import datetime, timezone
from matplotlib.patches import Patch
from utilidades import proyeccionOrtografica as po
from utilidades import leerOMNI
from utilidades import zonaOvalo
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_sun, ITRS
import astropy.units as u
#%%

#  Calcular el punto subsolar 
def subsolar_geodetic(dt_utc):
    t = Time(dt_utc, scale='utc')
    sun = get_sun(t)
    sun_itrs = sun.transform_to(ITRS(obstime=t))
    x, y, z = [c.to(u.m).value for c in sun_itrs.cartesian.xyz]
    r = np.sqrt(x*x + y*y + z*z)
    lat = np.degrees(np.arcsin(z / r))
    lon = np.degrees(np.arctan2(y, x))
    lon = ((lon + 180.0) % 360.0) - 180.0
    return lat, lon    # degrees (geo), lon east-positive

def _make_R_for_pole(pole='norte'):
    if pole.lower() == 'norte':
        lon0 = np.deg2rad(-82.86)
        lat0 = np.deg2rad(82.41)
    elif pole.lower() == 'sur':
        lon0 = np.deg2rad(107.0)
        lat0 = np.deg2rad(-80.65)
    else:
        raise ValueError("polo debe ser 'norte' o 'sur'")
    lmda = np.pi/2 - lat0
    R = np.array([
        [np.cos(lon0)*np.cos(lmda), -np.sin(lon0), np.cos(lon0)*np.sin(lmda)],
        [np.sin(lon0)*np.cos(lmda),  np.cos(lon0), np.sin(lon0)*np.sin(lmda)],
        [-np.sin(lmda),              0.0,           np.cos(lmda)]
    ])
    return R

# devuelve Δφ en (coordenadas geomagneticas en grados) 
def delta_phi_magnetic(dt_utc, pole='norte'):
    """
    Devuelve (delta_deg, lonmag_subsolar_deg)
    - delta_deg: Δφ en grados magnéticos ((-180,180])
    - lonmag_subsolar_deg: longitud magnética del subsolar (deg, (-180,180])
    Uso: lonMag_deg = 360*(t/24) + delta_deg  (t en horas MLT)
    """
    # 1) subsolar en GEO (deg)
    ss_lat_deg, ss_lon_deg = subsolar_geodetic(dt_utc)

    # 2) vector unitario subsolar en GEO (radianes internamente)
    lat_r = np.deg2rad(ss_lat_deg)
    lon_r = np.deg2rad(ss_lon_deg)
    v_geo = np.array([np.cos(lat_r)*np.cos(lon_r),
                      np.cos(lat_r)*np.sin(lon_r),
                      np.sin(lat_r)])   # shape (3,)

    # 3) construir R (la misma convención que en tu CGCtoGEO)
    R = _make_R_for_pole(pole)

    # 4) expresar v_geo en coordenadas MAG: v_mag = R^T * v_geo
    v_mag = R.T @ v_geo

    # 5) longitud magnética del subsolar (rad -> deg)
    lonmag_subsolar = np.arctan2(v_mag[1], v_mag[0])   # radians in (-pi,pi]
    lonmag_subsolar_deg = np.rad2deg(lonmag_subsolar)

    # 6) Δφ en rad: queremos Δφ tal que phi_subsolar = pi + Δφ -> Δφ = phi_subsolar - pi
    delta_rad = lonmag_subsolar - np.pi
    # normalizar a (-pi, pi]
    delta_rad = (delta_rad + np.pi) % (2*np.pi) - np.pi
    delta_deg = np.rad2deg(delta_rad)

    return float(delta_deg), float(lonmag_subsolar_deg)

# Entrega los coeficientes para calcular las amplitudes y las fases del modelo
#de Starkov, en funcion del parametro m, polo=0, ecuador=1, difuso=2
def coeficientes(m):
    
    if m == 0:
        b0 = np.fromstring(" -0.07 -10.06 -4.44 -3.77 -6.61   6.37 -4.48", sep=" ")
        b1 = np.fromstring(" 24.54  19.83  7.47  7.90 10.17  -1.10 10.16", sep=" ")
        b2 = np.fromstring("-12.53 -9.33  -3.01 -4.73 -5.80   0.34 -5.87", sep=" ")
        b3 = np.fromstring("  2.15  1.24   0.25  0.91  1.19  -0.38  0.98", sep=" ")
        
    elif m == 1:
        b0 = np.fromstring("  1.61 -9.59 -12.07 -6.56 -2.22 -23.98 -20.07", sep=" ")
        b1 = np.fromstring(" 23.21 17.78  17.49 11.44  1.50  42.79  36.67", sep=" ")
        b2 = np.fromstring("-10.97 -7.20  -7.96 -6.73 -0.58 -26.96 -24.20", sep=" ")
        b3 = np.fromstring("  2.03  0.96   1.15  1.31  0.08   5.56   5.11", sep=" ") 
        
    elif m == 2:
        b0 = np.fromstring("  3.44 -2.41 -0.74 -2.12 -1.68   8.69  8.61", sep=" ")
        b1 = np.fromstring(" 29.77  7.89  3.94  3.24 -2.48 -20.73 -5.34", sep=" ")
        b2 = np.fromstring("-16.38 -4.32 -3.09 -1.67  1.58  13.03 -1.36", sep=" ")
        b3 = np.fromstring("  3.35  0.87  0.72  0.31 -0.28  -2.14  0.76", sep=" ")
        
    else:
        return "Valor no permitido, m = 0,1,o 2"
    
    B = np.array([b0,b1,b2,b3])
    
    return np.transpose(B)

# Esta función calcula el indice Al para el metodo de Starkov, recibe el indice Kp=[0-9]
def calcularAL(kp):
    
    AL = 18.0 - 12.3*kp +  27.2*kp**2 - 2.0*kp**3
    
    return AL

# Esta función calcula las amplitudes y las fases usadas en la expresión de la
# colatitud de los ovalos del modelo de Starkov, recibe el indice Al, y el parametro
# m, para definir cual limite del ovalo se calcula: m=0 polar, m=1 limite ecuatorial
# m=2 limite difuso.
def amplitudes(AL, m):
    
    B = coeficientes(m) #Entrega una lista de listas de los coeficientes para un m
    A = np.zeros(4) 
    a = np.zeros(3)
    
    for i in range(7):
        
        b = B[i] # Toma solo una lista
        
        if i < 4:                  
            A[i] =  b[0] + b[1]*np.log10(abs(AL)) + b[2]*(np.log10(abs(AL)))**2 + b[3]*(np.log10(abs(AL)))**3
            
        elif i >= 4:
            a[i-4] =  b[0] + b[1]*np.log10(abs(AL)) + b[2]*(np.log10(abs(AL)))**2 + b[3]*(np.log10(abs(AL)))**3
            
    return A, a

# Devuele la latidud y la longitud en coordenadas centradas en el dipolo magnetico
# el calculo se hace en radianes y el resultado se transforma a coordenadas 
# geofracicas en grados usando la funcion CGCtoGEO
def ovalos(m, kp, fecha, polo):
    
         AL = calcularAL(kp) # caculo Al con el kp
         
         A, a = amplitudes(AL, m)   # saco los coeficientes para este caso
         
         t = np.arange(0,24, 0.5)   # una lista para la hora local magnetica
         
         # correcciones = deltaphi(fecha)   #Correción con la hora local magnetica
         # delta2 = correcciones["delta_phi_deg"]
         # MLT = correcciones["mlt_apexpy_h"]
         colatMg = np.zeros(len(t)) #Lista vacia para la colatitud 
         lonMg = np.zeros(len(t))   #Lista vacia para la longitud
         
         # if polo == "norte":
         #     delta = delta_phi_geom(fecha, -82.86)
         # elif polo == "sur":
         #     delta = delta_phi_geom(fecha, 107.0)
             
         correccion = delta_phi_magnetic(fecha, polo)
         delta = correccion[0]
             
         #Calculo las colatitudes 
           
         colatMg = A[0] + A[1]*np.cos(np.deg2rad(15*(t + a[0])))+ A[2]*np.cos( np.deg2rad(15*(2*t + a[1]))) + A[3]*np.cos(np.deg2rad(15*(3*t + a[2])))

         lonMg   = 2*np.pi*(t)/24 + np.deg2rad(delta)
            
         # Convierto la lista de colatitudes a latitudes    
         latMg = np.pi/2 - np.deg2rad(colatMg)
                    
         lat, lon = CGCtoGEO(latMg, lonMg, polo)

         return lat, lon

# cambiar de coordenadas geomaneticas a geocentricas
def CGCtoGEO(latMg, lonMg, polo):
    
    # #Convertir a radianes
    # latMg = np.deg2rad(latMg)
    # lonMg = np.deg2rad(lonMg)

    # Calculo las coordenadas cartesianas magneticas
    xm = np.cos(latMg) * np.cos(lonMg)
    ym = np.cos(latMg) * np.sin(lonMg)
    zm = np.sin(latMg)
    
    if polo == "norte":
        #lon0 es la longitud del polo norte magnetico en coordenadas geograficas
        lon0 = np.deg2rad(-82.86) 
    
        #lat0 es la latitud del polo norte magnetico en coordenadas geograficas
        lat0 = np.deg2rad(82.41) 
    
        #lmda es la diferencia longitudinal entre el polo norte geografico y magnetico
        lmda = np.pi/2 - lat0
    elif polo == "sur":
        
        #lon0 es la longitud del polo sur magnetico en coordenadas geograficas
        lon0 = np.deg2rad(107.0) 
     
        #lat0 es la latitud del polo sur magnetico en coordenadas geograficas
        lat0 = np.deg2rad(-80.65) 
    
        #lmda es la diferencia longitudinal entre el polo sur geografico y magnetico
        lmda = np.pi/2 - lat0
        
    # Defino la matriz de rotación
    R = np.array([[np.cos(lon0)*np.cos(lmda)  , -np.sin(lon0)  , np.cos(lon0)*np.sin(lmda)] ,
                  [np.sin(lon0)*np.cos(lmda)  , np.cos(lon0)   , np.sin(lon0)*np.sin(lmda)] , 
                  [-np.sin(lmda)              , 0              , np.cos(lmda)               ]])
    
    # defino dos listas vacias, donde se pondran la latidud y la longitud en coordenadas GEO
    latGeo = np.zeros(len(xm))
    lonGeo   = np.zeros(len(xm)) 
    
    # En este loop hago el producto de la matriz de rotación por un vector con las
    # componentes magneticas de cada punto. Las coordenadas obtenidas en GEO, se guaradan
    # en tres variables x,y,z y luego se calculan las coordendas esfericas de ese punto.
    
    for i in range(len(xm)):
        
        n = np.matmul(R, np.array([xm[i], ym[i], zm[i]]))
        
        x = n[0]
        y = n[1]
        z = n[2]
    
        latGeo[i] = np.arcsin(z)
        lonGeo[i] = np.arctan2(y, x)
    
    return np.rad2deg(latGeo), np.rad2deg(lonGeo)


def grafStarkov():
    
    for polo in ["norte","sur"]:
        i = 23
        fecha, ut, parametros, kps = leerOMNI("2025_01_01.lst", i)
    
        dt = datetime.fromtimestamp(ut, tz=timezone.utc)
        horaUT = dt.strftime("%H:%M:%S")
        fechaTor = fecha.strftime("%Y/%m/%d")
        
        ptsPo   = ovalos(0, kps[1], fecha, polo)
        ptsEq   = ovalos(1, kps[1], fecha, polo)
        kpf = 2
        ptsPo2   = ovalos(0, kpf, fecha, polo)
        ptsEq2   = ovalos(1, kpf, fecha, polo)

        # ptsDiff = ovalos(2, kps[1], fecha)
        
        fig1 = plt.figure(figsize=(7, 8), dpi=150)
        ax = po(fig1, polo, [45, 90], fecha, kps[1])
            
        zonaOvalo(ax, ptsPo, ptsEq)
        zonaOvalo(ax, ptsPo2, ptsEq2 , color="#0F2BE2", alpha=0.6)

        
        ax.text(
            x=0.025, y=0.98,          # posición relativa en axes (0 a 1)
            s=f"$UT={horaUT}$",
            transform=ax.transAxes,   # coordenadas relativas al axes
            fontsize=20,
            verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",edgecolor="k", alpha=1)
        )
        # ax.set_title(
        #     "Simulación con el modelo de Starkov para la Tormenta del " + fecha.strftime("%Y/%m/%d"), 
        #     fontname='serif', 
        #     size="large", 
        #     weight=700
        #     )
        fech = fecha.strftime("%Y/%m/%d")
        ax.text(
            x=0.25, y=0.09,          # posición relativa en axes (0 a 1)
            s=f"$ Simulación \; con \; el\;modelo \;de \;Starkov$\n $para \;la \; Tormenta \;del\; {fech}$",
            transform=ax.transAxes,   # coordenadas relativas al axes
            fontsize=16,
            weight=700,
            verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",edgecolor="k", alpha=1)
        )

        ovalo1 = Patch(facecolor='#19DA40', edgecolor='#6F1EC0', label=f'Kp = {kps[1]}')
        ovalo2 = Patch(facecolor='#0F2BE2', edgecolor='#6F1EC0', label=f'Kp = {kpf}')


        handles, labels = ax.get_legend_handles_labels()
        handles.append(ovalo2)
        labels.append(f'$Kp = {kpf}$')
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


#%%  
def main():
    
    grafStarkov()

    
    return 0



if __name__ == "__main__":
    main()
    
