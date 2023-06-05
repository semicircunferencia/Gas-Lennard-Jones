import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "DejaVu Sans"
#fig = plt.figure(figsize=(11.,8.))      # figura (ancho, largo)
ax = plt.subplot()      # subfigura

# Datos
data = np.loadtxt('EnergiasTyP.dat')
tiempo = data[0:1000,0] # cojo el tiempo de la primera columna
energiakin = data[0:1000,1]
energiapot = data[0:1000,2]
energia = data[0:1000,3]

# configurar ejes
ax.set_ylabel('Energía', fontname='DejaVu Sans', fontsize='12')
ax.set_xlabel('Tiempo', fontname='DejaVu Sans', fontsize='12')

#Cambiar ticks
#for label in ax.get_xticklabels():
    #label.set_fontproperties('Times New Roman')
#plt.xticks(fontsize='15')
#for label in ax.get_yticklabels():
    #label.set_fontproperties('Times New Roman')
#plt.yticks(fontsize='15')

# Creación de la gráfica
ax.plot(tiempo, energiakin, linestyle='-', marker='', markersize=4, color='#B4045F', label="E. cinética", linewidth=1.0)  #marker=puntos
ax.plot(tiempo, energiapot, linestyle='-', marker='', markersize=4, color='#5FB404', label="E. potencial", linewidth=1.0)
ax.plot(tiempo, energia, linestyle='-', marker='', markersize=4, color='#045FB4', label="E. total", linewidth=1.0)

plt.legend(loc="center right")


# Guardar la gráfica
plt.savefig('grafica_energia.pdf',dpi=300, bbox_inches = "tight")

# Mostrarla en pantalla
#plt.show()