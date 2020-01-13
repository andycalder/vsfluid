import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

data = np.fromfile('./output', dtype=np.float64)
n = int(data.size / (100 * 100))
data = np.reshape(data, (n, 100, 100))

fig = plt.figure()
img = plt.imshow(data[0, :, :], animated=True)

def update(i):
    field = data[i, :, :]
    img.set_data(field)
    #img.set_clim(vmin=field.min(), vmax=field.max())
    img.set_clim(vmin=-0.06, vmax=0.06)
    return img,

ani = animation.FuncAnimation(fig, update, frames=range(n), interval=1000/60, blit=True)

plt.axis('off')
plt.show()
