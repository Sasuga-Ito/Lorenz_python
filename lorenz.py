"""
Lorenz63 model simulation


"""

__author__ = 'balshark(Twitter: @balsharkPhD)'
__version__ = '1.0.0'
__date__ = '03/15/2022'
__status__ = 'Development'

from mpl_toolkits.mplot3d import Axes3D
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

class Lorenz63():

	def __init__(self):
		# time step
		self.t0 = 0
		self.tn = 30
		self.dt = 0.01
		self.tt = np.arange(self.t0,self.tn,self.dt)

		# initial condition
		self.init_x = 1.5
		self.init_y =-1.5
		self.init_z = 25.0

		# lorenz63 paramter
		self.s = 10.0
		self.r = 28.0
		self.b = 8.0/3.0

	def main(self):
		print('Lorenz63 simulation')

		# insert initial condition
		xk = [self.init_x]
		yk = [self.init_y]
		zk = [self.init_z]
		xyz = [self.init_x,self.init_y,self.init_z]
		dt = self.dt

		# time integration using 4th order Runge-Kutta
		for t in self.tt:
			k1 = self.cal_lnz(xyz)
			k2 = self.cal_lnz(xyz+dt*k1*0.5)
			k3 = self.cal_lnz(xyz+dt*k1*0.5)
			k4 = self.cal_lnz(xyz+dt*k3)
			xyz = xyz + dt*(k1+k2+k3+k4)/6
			xk.append(xyz[0])
			yk.append(xyz[1])
			zk.append(xyz[2])

		self.output(xk,yk,zk)

		print('done')

	def cal_lnz(self,xyz):
		# insert input parameter
		x,y,z = xyz
		s,r,b = self.s,self.r,self.b

		# calculate difference term
		dxdt =-s*(x-y)
		dydt = x*(r-z) - y
		dzdt = x*y - b*z

		return np.array([dxdt,dydt,dzdt])

	def output(self,xk,yk,zk):
		# output plot and figure
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')
		ax.plot(xk,yk,zk,color='green')
		plt.show()
		fig.savefig('lorenz63.png')


if __name__ == '__main__':
	proc = Lorenz63()
	proc.main()
