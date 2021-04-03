#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

D = 0.25 #m
z = 0.0 #m
g = 9.8 #m/s/s
T = 820+273.15 #K
R = 8.314 #J/mol/K
A = 0.25*np.pi*D**2
Ko = np.array([4.65E13,8.75E8,3.85E11,9.81E8,5.87E4,1.03E12,7.08E13])
E = np.array([65.2,32.7,65.25,36.92,7.04,41.26,60.43])
KT = np.multiply(Ko,np.exp(-E/(R*T)))
FT0 = 100 #mol/s
fT0 = FT0/A
Po = 11*101325 #Pa
vo = (4*FT0*R*T/(np.pi*D**2*Po))
Mi = np.array([44.1E-3,30.07E-3,28.05E-3,2.01E-3,16.04E-3,42.08E-3,26.04E-3,54.09E-3,18.02E-3])
yo = np.array([0.0,0.594,0.006,0.0,0.0,0.0,0.0,0.0,0.4])
Go = np.array([0.0,Po,vo,sum(np.multiply(yo,Mi))])
Co = Po*yo/(R*T)
y = np.concatenate((Go,Co))


def func(y):
    Af = np.identity(np.size(y))
    b = np.zeros(np.size(y))
    b[0] = 0
    b[1] = 0
    b[2] = 0
    b[3] = 0
    b[4] = -0.5*KT[2]*y[5]
    b[5] = -(KT[1]*y[6]*y[7]-KT[0]*y[5]-KT[2]*y[5])
    b[6] = -(KT[0]*y[5]-KT[1]*y[6]*y[7]-KT[5]*y[10]*y[6])
    b[7] = -(KT[0]*y[5]-KT[1]*y[6]*y[7])
    b[8] = -(0.5*KT[2]*y[5]+KT[3]*y[9]-KT[4]*y[10]*y[8]+KT[6]*y[5]*y[6])
    b[9] = -(KT[4]*y[8]*y[10]-KT[3]*y[9]+KT[6]*y[5]*y[6])
    b[10] = -(KT[3]*y[9]-KT[4]*y[8]*y[10]-KT[5]*y[6]*y[10])
    b[11] = -(KT[5]*y[6]*y[10])
    b[12] = 0.0
    Af[1][1] = 2*R*T+y[3]*g*z+y[3]*y[2]**2 #dp
    Af[1][3] = 2*g*z*y[1]+y[1]*y[2]**2
    Af[1][2] = 2*y[1]*y[2]*y[3]
    Af[2][1] = y[2]*y[3] #dv
    Af[2][2] = y[3]*y[1]
    Af[2][3] = y[2]*y[1]
    Af[3][3] = fT0
    for i in range(0,np.size(Mi)):
        Af[3][i+4] = -Mi[i]   
    f = np.linalg.solve(Af, b)
    return f;

xi = 0.0
xf = 95.0
N = 1000
h = (xf-xi)/N

for i in range(0,N):
    if i==0:
        fi = y + h*func(y)
        y = np.vstack((y,fi))
    else:
        fi = y[i,:] + h*func(y[i,:])
        y = np.vstack((y,fi))

# Plot the data
#plt.plot(y[:,0], y[:,0], label='L')
#plt.plot(y[:,0], y[:,1], label='P')
#plt.plot(y[:,0], y[:,2], label='v')
#plt.plot(y[:,0], y[:,3], label='Mp')
plt.plot(y[:,0], y[:,4], label='P')
plt.plot(y[:,0], y[:,5], label='E')
plt.plot(y[:,0], y[:,6], label='Y')
plt.plot(y[:,0], y[:,7], label='H')
plt.plot(y[:,0], y[:,8], label='M')
plt.plot(y[:,0], y[:,9], label='L')
plt.plot(y[:,0], y[:,10], label='A')
plt.plot(y[:,0], y[:,11], label='B')
plt.plot(y[:,0], y[:,12], label='W')

# Show the plot
plt.xlabel('Lenght (m)')
plt.ylabel('Concentracion')
plt.legend()
plt.show()