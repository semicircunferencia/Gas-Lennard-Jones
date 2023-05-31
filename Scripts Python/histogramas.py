import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "DejaVu Sans"
#fig = plt.figure(figsize=(11.,8.))      # figura (ancho, largo)
ax = plt.subplot()      # subfigura

# Calculo la temperatura
'''archivotemp = np.loadtxt('EnergiasTyP.dat')
temp1=archivotemp[10:300,4]
temp2=archivotemp[350:600,4]
temp3=archivotemp[650:900,4]
temp4=archivotemp[950:1200,4]
temp5=archivotemp[1250:1500,4]
temp6=archivotemp[1150:1740,4]

T1=np.mean(temp1)
T2=np.mean(temp2)
T3=np.mean(temp3)
T4=np.mean(temp4)
T5=np.mean(temp5)
T6=np.mean(temp6)

print(T1, " ")
print(T2, " ")
print(T3, " ")
print(T4, " ")
print(T5, " ")
print(T6, " ")'''


# Datos
archivodatos = np.loadtxt('Velocidades.dat')
datos = archivodatos[:,0] # cojo el modulo de la la velocidad

# x = np.linspace(min(datos), max(datos), 100)
# y= x
# y=(2*3.141593*T)**(-0.5)*np.exp(-(x-0)**2/(2*T))

# configurar ejes
ax.set_ylabel('Densidad de partículas $g(r)$', fontname='DejaVu Sans', fontsize='12')
ax.set_xlabel('Distancia $r$', fontname='DejaVu Sans', fontsize='12')

#Cambiar ticks
#for label in ax.get_xticklabels():
    #label.set_fontproperties('Times New Roman')
#plt.xticks(fontsize='15')
#for label in ax.get_yticklabels():
    #label.set_fontproperties('Times New Roman')
#plt.yticks(fontsize='15')

# Creación de la gráfica
plt.hist(datos, bins="auto", density=True, color='#CBC3E3')
# ax.plot(x, y, linestyle='-', marker='', markersize=4, color='#B4045F')  #marker=puntos


# Guardar la gráfica
plt.savefig('Histograma_y.pdf',dpi=300, bbox_inches = "tight")

# Mostrarla en pantalla
#plt.show()