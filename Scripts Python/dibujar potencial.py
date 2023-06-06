import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

plt.rcParams["font.family"] = "DejaVu Sans"
# Create the figure and axis
fig, ax = plt.subplots()

x=np.linspace(0.5,4.0,200)
y=4*(1.0/np.power(x,12)-1.0/np.power(x,6))
ax.plot(x, y, linestyle='-', marker='', markersize=4, color='#B4045F')

# Set limits
ax.set_ylim(-1.5,2)
ax.set_xlim(0,4)

ax.set_ylabel('$V(r)$', fontname='DejaVu Sans', fontsize='12')
ax.set_xlabel('$r$', fontname='DejaVu Sans', fontsize='12')

# Set labels and title

# Show the plot
plt.savefig('lenn.pdf',dpi=300, bbox_inches = "tight")
