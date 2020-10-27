from ctypes import Structure, c_int, c_double, POINTER, CDLL
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

nx = 200           # Grid width
ny = 200           # Grid height
nt = 2000          # Number of time steps
h = 0.02           # Grid spacing
dt = 0.01          # Time step
nu = 0.003         # Kinematic viscosity
epsilon = 2 * h    # Half interface thickness
sigma = 0          # Surface tension (TODO)
g = 1              # Gravity (TODO)

# Set initial fluid volume
phi = np.zeros((ny, nx), dtype=np.float64)

for i in range(ny):
    for j in range(nx):
        # Circle in center distance function
        #dist = np.sqrt((i-100)*h*(i-100)*h+(j-100)*h*(j-100)*h) - 50*h

        # Half half distance function
        dist = (nx / 2 - j) * h
        
        level = 2 * epsilon / np.pi

        if dist > epsilon:
            phi[i,j] = -level
        elif dist < -epsilon:
            phi[i,j] = level
        else:
            phi[i,j] = -level * np.sin(dist / level)

class Params(Structure):
    _fields_ = [('nx', c_int),
        ('ny', c_int),
        ('nt', c_int),
        ('h', c_double),
        ('dt', c_double),
        ('nu', c_double),
        ('epsilon', c_double),
        ('sigma', c_double),
        ('g', c_double),
        ('phi0', POINTER(c_double * (nx * ny)))]

params = Params(nx,
    ny,
    nt,
    h,
    dt,
    nu,
    epsilon,
    sigma,
    g,
    phi.ctypes.data_as(POINTER(c_double * (nx * ny))))

# Set up interface with c library
c_library = CDLL('./fluid.so')
c_library.solve.restype = POINTER(c_double * (nx * ny * nt))
c_library.solve.argtypes = [Params]

# Call C code
output = c_library.solve(params)

data = np.reshape(output.contents, (nt, nx, ny))

# Define animation loop
def update(i):
    field = data[i, :, :]
    img.set_data(field)
    img.set_clim(vmin=field.min(), vmax=field.max())
    return img,

# Display animation
fig = plt.figure()
img = plt.imshow(data[0, :, :], animated=True)
anim = ani.FuncAnimation(fig, update, frames=range(nt), interval=1000/60, blit=True)
plt.axis('off')
plt.show()

# Save animated gif (must have imagemagick installed)
#fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
#ani.save('render.gif', writer='imagemagick')
