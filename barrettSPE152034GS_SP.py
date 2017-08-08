#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 16:49:21 2017
Modified on Fri Aug 4 13:26:30 2017
@author: SP,GS
"""

import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

def TyP(w,Tb,T0,P0,zInic,zFinal):
	# Defición del sistema de ODEs
	def func(y, z, a1, a2, Cj, b1, b2, b3):
		P, T = y
	        f = 1 # factor Fanning: hay que calcularlo
                Z = 1 # factor de compresibilidad del gas: hay que calcularlo
	        Cj = 0.0 # Ojo: se debe calcular. Depende de P y T.
		H = a1*f*Z*T/P + (a2/Z)*P/T
		dydz = [H, Cj*H + b1*T + b2*z + b3]
		return dydz

	# parámetros y condiciones de borde
	rto = (2*0.0254+(7/8.0)*0.0254)/2.0
	rti = rto - 0.3
        A = np.pi*rti**2
	Ma = 28.97 #peso molecular del aire. De Barrett SPE152034
        R = 8.31445 # constante gases
        Gammag = 0.5537 # gravedad específica metano

        a1 = -(w**2)*R/((A**2)*rti*Ma*Gammag)
	#a1 = -.001

	g = 9.8
	a2 = -Ma*Gammag*g/R
	#a2 = -0.01


	Ut = 1/((rto/2)*np.log(200.0/rto))
	cp = 2260.0 # asumimos Cp constante, pero depende de la temperatura T. Es el Cp del metano
	Lr = (2*np.pi*rto*Ut)/(cp*w)

	b1 = -Lr

	gG = 1/27.5
        
	b2 = -Lr*gG

	Teb = Tb
	
	b3 = Lr*Teb - g/cp

	y0 = [P0, T0]

	z = np.linspace(zInic,zFinal,int((zFinal-zInic)/0.001)) # malla espacial


	# solver ODEs
	sol = [odeint(func,y0,z,args=(a1,a2,Cj,b1,b2,b3)), z]
	return sol

# Temperatura de la geoterma
def Tzeta(z,Tb,gG):
	T = np.zeros(len(z));
	for i in range(0,len(z)):
		T[i] = Tb - gG*z[i]
	return T

gG = 1/27.5
finalPozo = 1100.0
tempSup = 20.0

Tb = gG*finalPozo + tempSup

#perforaciones
zPerf = [100.0,150.0,200.0,300.0]		# posiciones de las perforaciones en metros, empezando por la última
#fPerf = [0.2,0.4,0.1,0.3]			# fracciones de caudal másico, deben sumar 1
fPerf = [0.2,0.4,0.3,0.1]			# fracciones de caudal másico, deben sumar 1
nPerf = len(zPerf)

wTot = 0.656*100000/86400.0			# caudal másico total: 100 mil m³ por día, en kg/s, a presión 1 atm

wPerf = np.zeros(len(fPerf));
for i in range(0,nPerf):
	wPerf[i] = wTot*fPerf[i]

w = 0

# suponemos que los gases entran a la temperatura de la geoterma
TPerf = Tzeta(zPerf,Tb,gG)

for i in range(0,nPerf):
	print 'paso ',i
	#P0=1
	if (i==0):
		T0 = TPerf[0]
		P0 = 1.e6 #presión primera perforación 
	else:
		T0 = (TPerf[i]*wPerf[i] + TSol[len(TSol)-1]*w)/(wPerf[i]+w) 
		P0 = PSol[len(PSol)-1]
	zInic = zPerf[i]
	if (i==nPerf-1):
		zFinal = finalPozo	# si es el último tramo, calcula hasta arriba (z=finalPozo)
	else:
		zFinal = zPerf[i+1]	# si no, calcula hasta la próxima perforación
	w = w+wPerf[i]

	TyPSoli,zi = TyP(w,Tb,T0,P0,zInic,zFinal)
	if (i==0):
		PSol = TyPSoli[:,0]
		TSol = TyPSoli[:,1]
		z = zi
	else:
		PSol = np.concatenate((PSol,TyPSoli[:,0]))
		TSol = np.concatenate((TSol,TyPSoli[:,1]))
		z = np.concatenate((z,zi))

#z = np.linspace(zInic,zFinal,10000)
# Graficos
#plt.plot(z, PSol, 'b', label=u'Presión')
plt.plot(finalPozo-z, TSol, 'g', label=u'Temperatura')
geoterma = np.array([ [finalPozo,0], [Tb,tempSup] ])
plt.plot( geoterma[0], geoterma[1], label=u'Geoterma' )
plt.legend(loc='best')
plt.xlabel('Profundidad [m]')
plt.ylabel('Temperatura [C]')
plt.xlim([0, finalPozo])
plt.ylim([10, Tb+10])
plt.xticks(np.arange(0, finalPozo, 100))
plt.yticks(np.arange(10, Tb+10, 5))
plt.grid()
plt.show()
