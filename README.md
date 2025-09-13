# Simulacion de Ovalos Aurorales
Author: Iker Algañaraz, mail: alganaraziker@gmail.com

Source code usado para la simulación del ovalos aurorales utilizando los modelos de Ley de Escala, Starkov y Tsyganenko T96. 
Los codigos son funciones sencillas para el caso de la Ley de Escala y el modelo de Starkov. Para el modelo de Tsyganenko se usa principalmente el modulo geopack: https://github.com/tsssss/geopack/tree/master , el cual incluye hasta el modelo t04 de Tsyganenko, y los coeficientes para IGRF.
# Modo de uso
Hay 4 archivos de codigo, cada modelo tiene el suyo y el archvio de utilidades contiene ciertas funciones que sirven para las graficas o para leer datos y son utilizadas en todos los modelos. Tambien hay 3 archivos ".lst"con datos de tormentas y que se usan a modo de ejemplo en los codigos.

## Starkov
Un codigo sencillo basado completamente en la metodologia descrita en el paper de Fred Sigernes para al app de AuroraForecast, se puede descargar la app y acceder al paper desde: http://aurora.unis.no/Forecast3D.html. La unica parte que no se detalla en el paper es el calculo de la diferencia longitudinal entre el punto subsolar y los polos magneticos. Esto se realizo utilizando astropy y se documenta en el codigo.

## T96
La metodologia de este codigo esta basada en el paper de Tsyganenko de 2019 "Tsyganenko, N. A., Secular drift of the auroral ovals: How fast do they actually move?, Geophysical Research Letters, 46, 3017-3023, 2019.". Leyendo de ahi y siguiendo el codigo es claro el uso de cada funcion. Es importante la elección de los limites de contorno a la hora del funcionamiento, ya que puede dar lugar a ovalos muy poco definidos si se toman muchas lineas que no son cerradas. Para encontrar los ovalos se utiliza principalmente la función "trace" del modulo de geopack, que permite seguir las lineas de campo, usando eso el codigo es sencillo y puede ser adaptado a cualquiera de los otros modelos de campo externo disponibles en geopack.

## Ley de Escala
Es solo una función que realiza el escalado con la presion dinamica, necesita obtener las latitudes y longitudes de otro modelo primero, puede ser el de Starkov o T96.

## Utilidades
La función mas importante es "leerOMNI" ya que es la que permite leer datos de tormentas, en su definicion esta claro los parametros puestos y su orden, la misma se puede modificar facilmente para otros parametros que se deseen. Las otras funciones sirven para graficar, en particular "proyecciónOrtografica" automatiza el tener que generar proyecciónes sobre los polos, en este caso lo hace usando una transformación Ortografica, para otras funciones se uso de transformación de coordenadas las "Geodetic" o las "PlateCarree". "graficarOvalo" sencillamente grafica una curva de puntos unidas por lineas para los limites de las auroras, "zonaOvalo" por su parte plotea una zona pintada entre los dos limites para denotar el ovalo, lo hace usando shapely. en algunos casos puede dar error sobre todo al graficar a la vez en el polo norte y en el sur, una solucion momentanea fue cambiar el orden de limite interno y externo que se le entragaba a la funcion.
