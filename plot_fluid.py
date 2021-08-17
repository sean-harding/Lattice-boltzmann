import matplotlib.pyplot as py
from matplotlib.animation import FuncAnimation
import scipy as sp
import numpy as np
#Need to add a header to the velocity text file to say how many sites there are
data = np.loadtxt('velocities.txt',delimiter='\n',dtype='str')
size = data[0].split(",")
size = [int(s) for s in size]
N = size[0]*size[1]
data_parsed = sp.recarray((int(len(data)/N),N),dtype=[('vx','float'),('vy','float')]) #Shape is nFrames,nSites
#Don't plot boundaries
for i,d in enumerate(data):
    if i == 0:
        pass
    else:
        data_parsed[(i-1)//N,(i-1)%N] = tuple(d.split(','))

fig,ax = py.subplots()
ax.set(xlim=(0,size[0]),ylim=(0,size[1]))
x = sp.linspace(0,size[0]-1,size[0])
y = sp.linspace(0,size[1]-1,size[1])
X,Y = sp.meshgrid(x,y)

qax = ax.quiver(X,Y,data_parsed[1]['vx'],data_parsed[1]['vy'])
def for_animate(i):
    qax.set_UVC(data_parsed[i]['vx'],data_parsed[i]['vy'])
anim = FuncAnimation(fig,for_animate,frames=len(data_parsed)-1,repeat=False)
py.draw()
py.show()