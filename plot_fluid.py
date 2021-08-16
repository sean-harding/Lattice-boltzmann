import matplotlib.pyplot as py
from matplotlib.animation import FuncAnimation
import scipy as sp
import numpy as np
#Need to add a header to the velocity text file to say how many sites there are
N = 25
data = np.loadtxt('velocities.txt',delimiter='\n',dtype='str')

data_parsed = sp.recarray((int(len(data)/25),25),dtype=[('vx','float'),('vy','float')]) #Shape is nFrames,nSites

for i,d in enumerate(data):
    data_parsed[i//25,i%25] = tuple(d.split(','))

fig,ax = py.subplots()
ax.set(xlim=(0,5),ylim=(0,5))
x = sp.linspace(0,4,5)
y = sp.linspace(0,4,5)
X,Y = sp.meshgrid(x,y)

qax = ax.quiver(X,Y,data_parsed[0]['vx'],data_parsed[0]['vy'])
def for_animate(i):
    qax.set_UVC(data_parsed[i]['vx'],data_parsed[i]['vy'])
anim = FuncAnimation(fig,for_animate,frames=len(data_parsed)-1)
py.draw()
py.show()