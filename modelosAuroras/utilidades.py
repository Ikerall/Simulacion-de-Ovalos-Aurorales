import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.feature.nightshade import Nightshade
import pandas as pd
from datetime import datetime, timedelta, timezone
from shapely.geometry import Polygon
import warnings


# Esta funcion plotea los puntos sobre la proyección, recibe el axes,
# los puntos en coordenadas geo, el tipo de grafica, la fecha y el color.
def graficarOvalo(ax, puntos, tipo, color="red", linestyle="-",label=None):
    
    
    # Sacamos las latitudes y longitudes por separado
    lats = puntos[0]
    lons = puntos[1]
    
    if tipo == "puntos":
            
        ax.scatter(
            lons, lats,
            transform=ccrs.Geodetic(),
            color=color,
            linewidth=0.5,
            linestyle=linestyle,
            label=label,
            zorder=0,
            s = 6,
            marker = "o",
        )
        
    elif tipo == "lineas":        
        ax.plot(
            lons, lats,
            transform=ccrs.Geodetic(),
            color=color,
            linewidth=2,
            linestyle=linestyle,
            label=label,
            zorder=0,
        )
        
    return ax

# genera una figura con shapely para la zona entre los ovalos
def zonaOvalo(ax, pInterno,pExterno, color="#19DA40",edgecolor="#6F1EC0", alpha=0.5):
    
    warnings.filterwarnings("ignore", message="Approximating coordinate system.*")
                            
    latExt = pExterno[0]
    lonExt = pExterno[1]
    latInt = pInterno[0]
    lonInt = pInterno[1]
    

    poly = Polygon(shell = list(zip(lonExt, latExt)), holes = [list(zip(lonInt[::-1], latInt[::-1]))])
    feature = cfeature.ShapelyFeature([poly], crs=ccrs.Geodetic(), 
                                      facecolor=color,
                                      linewidth=1.8, 
                                      edgecolor=edgecolor, 
                                      alpha=alpha, 
                                      zorder=3)
    ax.add_feature(feature)

# función para leer datos desde OMNI web
def leerOMNI(archivo, indice):
    
    # Leer datos desde el archivo bajado en OMNI web
    data = pd.read_csv(archivo,
                       sep="\s+", 
                       names=["Año","Día","Hora","ByIMF(GSM)","BzIMF(GSM)","Pdyn(nPa)","kp*10","Dst(nT)"])

    #Saco un fila de parametros correspondientes a una hora
    fila = np.array(data.iloc[indice].tolist())

    fecha = fila[:3]                                    
    parametros = [fila[5], fila[7], fila[3], fila[4]]
    kpMedido = int(fila[6]/10)
    
    # restricción del kp para el modelo de t86, sujeto a modificación para los otros modelos
    if 0 <= kpMedido < 1:
        kp = 1
    elif 1 <= kpMedido < 2:
        kp = 2
    elif 2 <= kpMedido < 3:
        kp = 3
    elif 3 <= kpMedido < 4:
        kp = 4
    elif 4 <= kpMedido < 5:
        kp = 5
    elif 5 <= kpMedido < 6:
        kp = 6
    elif 6 <= kpMedido:
        kp = 7
    
    kps =[kp, kpMedido]
    
    año, dia, hora = map(int, fecha)  
    fecha = datetime(año, 1, 1, tzinfo=timezone.utc) + timedelta(days=dia-1, hours=hora)
    ut = int(fecha.timestamp())  
    
    return fecha, ut, parametros, kps


# Esta función crea y devuelve el axes de matplot con la proyección Ortografica de 
# cartopy. Los argumentos que recibe son el polo de la grafica, "sur" o "norte".
# fecha es un objeto tipo datetime con la fecha de los datos. Rango es un
# array con el rango inferior y superior de latitudes [latmin , latmax] 
def proyeccionOrtografica(fig, polo, rango, fecha, info=None):
    
    # segun el parametro sea sur o norte la proyección se centra en un polo o otro.
    if polo == "sur":
        lat = -90
        lon = 180
        la1 = -rango[1]
        la2 = -rango[0]
        latitudes = np.arange(-90, -30 + 1, 15)  # latitudes visibles
        poloMG =[-80.65, 107.0]
    else:
        lat = 90
        lon = 0
        la1  = rango[0]
        la2  = rango[1]
        latitudes = np.arange(30, 90 + 1, 15)     # latitudes visibles
        poloMG =[82.41, -82.86]
        
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Times New Roman'
    plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
    plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'
    
    # creo el objeto axes de matplotlib, en este caso es una projección de cartopy, de tipo ortografica.
    ax = fig.add_subplot(projection=ccrs.Orthographic(lon, lat))

    # detalles para la proyección
    ax.stock_img()  # fondo satelital (opcional)
    ax.coastlines(resolution='110m', linewidth=0.3 ,zorder=2)  # líneas de costa
    ax.add_feature(cfeature.BORDERS.with_scale('110m'), linewidth=0.3, zorder=3) # límites de países
    ax.add_feature(Nightshade(fecha), alpha=0.5,zorder=2) # Sombra de la noche en función de la hora
    # ax.set_title("Tormenta del " + fecha.strftime("%Y/%m/%d"), fontname='serif', size="large", weight=700)

    #modificar el rango de latitudes de la imagen
    ax.set_extent([-180, 180, la1, la2], crs=ccrs.PlateCarree())
    
    #grilla de puntos
    gl = ax.gridlines(linewidth=0.2, linestyle=(5,(10, 3)), color="k" ,draw_labels=False, zorder=1)
    lon_step = 45
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, lon_step))
    lat_step = 15
    gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, lat_step))
    
    # Etiquetas de latitud sobre los paralelos
    lon_label = 45  # meridiano donde se colocan las etiquetas
    for lat_val in latitudes:
        ax.text(
            lon_label, lat_val, f"{abs(lat_val)}°",
            transform=ccrs.Geodetic(),
            horizontalalignment='right',
            verticalalignment='bottom' if lat_val >= 0 else 'top',
            fontsize=20,
            fontname='Times New Roman',
            color='black',
            zorder=4
        )

    #Grosor de los bordes
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    #Grosor de los ticks
    plt.tick_params(width=1.5, length=5,labelsize=15, top=True, right=True, direction="in",color="w") 
    plt.tight_layout()
    
    ax.scatter(poloMG[1], poloMG[0], transform=ccrs.PlateCarree(), color="red", marker="D",s=100, alpha=1, zorder=3,label=f"$Polo\;{polo}\;magnético$")
    return ax