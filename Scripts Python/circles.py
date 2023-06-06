import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

plt.rcParams["font.family"] = "DejaVu Sans"
# Create the figure and axis
fig, ax = plt.subplots()

# Plot the circles
circle1 = patches.Circle((5, 5), radius=5, edgecolor='#B4045F', facecolor='none')
circle2 = patches.Circle((5, 5), radius=6, edgecolor='#5FB404', facecolor='none')
ax.add_patch(circle1)
ax.add_patch(circle2)

ax.plot([5], [5], linestyle='', marker='.', markersize=8, color='#3A405A', linewidth=1.0)

# Define the region to color


# Create a polygonal patch for the colored region

# Set limits
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

# Set labels and title
ax.set_aspect("equal")

# Show the plot
plt.savefig('circulos.pdf',dpi=300, bbox_inches = "tight")
