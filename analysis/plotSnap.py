import numpy as np
import matplotlib.pyplot as plt
import h5py


# parameters
N       = 64     # resolution
beta    = 0.06   # gravity
Lbox    = 1.0/beta

b = 1
G = 1

snap = 20

simDir = '../output/'

# read snapshot
snapFile = simDir + 'N' + str(N) + 'box' + '%0.4f' % (Lbox) + 'B' + str(beta) + 'b' + str(b) + 'G' + str(G) + '/snap_' + '%03d' % (snap) + '.h5' 

f = h5py.File(snapFile, 'r')

rho = f.get('/psiRe')[()]**2 + f.get('/psiIm')[()]**2



# plot snapshot

plt.imshow(np.log10(np.mean(rho,2)), cmap="plasma")
plt.axis('off')
plt.savefig("snap" + '%03d' % (snap) + ".png", bbox_inches='tight')
plt.show()
