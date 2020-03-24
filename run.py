from ctypes import *
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Set simulation parameters
size = 100      # Grid size
n = 1000        # Number of time steps
h = 0.02        # Grid spacing
U = 1.0         # Top boundary velocity
dt = 0.01       # Time step
nu = 0.003      # Kinematic viscosity

# Set initial fluid volume
phi = np.zeros((size, size), dtype=np.float64)
for i in range(size):
    for j in range(size):
        phi[i,j] = size / 2 * h - j * h

# Compile c code
os.system('cc -Ofast -fPIC -shared -o vsfluid.so vsfluid.c')

# Set up interface with c library
c_library = CDLL('./vsfluid.so')
c_library.solve.restype = POINTER(c_double * (size * size * n))
c_library.solve.argtypes = [POINTER(c_double * (size * size)), c_int, c_int, c_double, c_double, c_double, c_double]

# Call c code
output = c_library.solve(phi.ctypes.data_as(POINTER(c_double * (size * size))), size, n, h, dt, nu, U)

data = np.reshape(output.contents, (n, size, size))

# Define animation loop
def update(i):
    field = data[i, :, :]
    img.set_data(field)
    #img.set_clim(vmin=field.min(), vmax=field.max())
    img.set_clim(vmin=-0.02, vmax=0.02)
    return img,

# Display animation
fig = plt.figure()
img = plt.imshow(data[0, :, :], animated=True)
ani = animation.FuncAnimation(fig, update, frames=range(n), interval=1000/60, blit=True)
plt.axis('off')
plt.show()

# Save animated gif (must have imagemagick installed)
#fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
#ani.save('render.gif', writer='imagemagick')
