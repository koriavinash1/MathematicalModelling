import numpy as np
import matplotlib.pyplot as plt


###################################################################
#               Initialize Space and Time steps                   #
#                                                                 #
###################################################################
k = 0.005
h = 0.05
Tmax = 2.0

Xmin   = -0.55
Xmax   = 0.5
deltaX = h  
X0     = np.linspace(Xmin, Xmax, int((Xmax - Xmin)/deltaX))

Tmin = 0
Tmax = Tmax
deltaT = k
T      = np.linspace(Tmin, Tmax, int((Tmax - Tmin)/deltaT))



###################################################################
#            Function Defination for Burger's equation            #
#                                                                 #
###################################################################
# U_t + F(U(X,t))_x = 0
# U(X, T) is constant

# initilize F
rhoMax = 1.0
F = lambda U: U**2/2
FD = lambda U: U
FStarSolve = lambda: 0

LB, RB = 1.0, 0.0  
U2rho = lambda U: (1.0 - U)*rhoMax/2.0
rho2U = lambda rho: (1.0 - rho*2.0/rhoMax)


###################################################################
#                 Initial Condition Defination                    #
#				     Conditions on Rho not U                      #
###################################################################

def initial_conditions_red(x):
	if x <= 0.0:
		return 1.0
	elif x > 0.0:
		return -1.0	


initial_conditions = initial_conditions_red
# boundary contitions 
BL, AL = 1, 0 # before and after light
# initial U ...
U0 = rho2U(0.55)*np.ones(len(X0))
# U0 = np.array([initial_conditions(x) for x in X0])


###################################################################
#       Plot Characteristics of Burgers Equation (Static)         #
###################################################################

def characteristic_solution(x0, t):
	x = []
	for i in range(len(x0)):
		x.append(rho2U(initial_conditions(x0[i])*t + x0[i]))
	return np.array(x)


def plot_characteristic(x, t):
	for i in range(x.shape[0]):
		plt.plot(x[i], t, 'g')
	plt.show()


# Rankine-Hugoniot condition
def RHCondition(F, U):
	s = []
	for i in range(len(U) - 1):
		s.append((F(U[i+1]) - F(U[i]))/ (U[i+1] - U[i]))
	return s

# plot_characteristic(characteristic_solution(X0, T), T)



###################################################################
#           Numerical Schemes for solving Burger's Eqn            #
###################################################################

def Upwind_Method(F, FD, U, k, h):
	temp = np.zeros_like(U)
	for i in range(1, len(U) - 1):
		if FD(U[i]) >= 0:
			temp[i] = U[i] - (k/h)*(F(U[i]) - F(U[i-1]))
		else:
			temp[i] = U[i] - (k/h)*(F(U[i +1]) - F(U[i]))
	return temp


def Lax_Friedrichs_scheme(F, FD, U, k, h):
	temp = np.zeros_like(U)
	for i in range(1, len(U) - 1):
		temp[i] =0.5*(U[i+1] + U[i-1]) - (k/(2.*h))*(F(U[i +1]) - F(U[i]))
	return temp


def Mac_Cornack_scheme(F, FD, U, k, h):
	Ustar = lambda u, up1: u - (k/h)*(F(up1) - F(u))
	temp = np.zeros_like(U)
	for i in range(1, len(U) - 1):
		temp[i] =0.5*(U[i] + Ustar(U[i], U[i+1])) -\
		 (k/h)*(F(Ustar(U[i], U[i+1])) - F(Ustar(U[i-1], U[i])))
	return temp


def Richtmyer_two_step_Lax_Wendroff_scheme(F, FD, U, k, h):
	Ustar = lambda u, up1: 0.5*(u + up1) - (k/h)*(F(up1) - F(u))
	temp = np.zeros_like(U)
	for i in range(1, len(U) - 1):
		temp[i] = U[i] - (k/h)*(F(Ustar(U[i], U[i+1])) - F(Ustar(U[i-1], U[i])))
	return temp


def Gudonov_Method(F, FD, FStarSolve, U, k, h):
	def Speed(u, up1):
		return float(F(up1) - F(u))/(up1 - u)

	def Ustar(u, up1):
		if (FD(u) >= 0) and (FD(up1) >= 0):
			return u
		elif (FD(u) < 0) and (FD(up1) < 0):
			return up1
		elif (FD(u) >= 0) and (FD(up1) < 0):
			if Speed(u, up1) >= 0:
				return u
			else:
				return up1
		elif (FD(up1) >= 0) and (FD(u) < 0):
			return FStarSolve()

	temp = np.zeros_like(U)
	for i in range(1, len(U) - 1):
		temp[i] = U[i] -(k/h)*(F(Ustar(U[i], U[i+1])) - F(Ustar(U[i-1], U[i])))
	return temp


###################################################################
#                     Solve Burgers Equation                      #
###################################################################


def find_solution(U0, T, nsteps, k, h, method='Gudonov_Method'):
	tsteps = int(T/k) + 1
	U = np.zeros((len(U0), tsteps))
	U[:, 0] = U0
	U[0, :] = U[0, 0]  
	idx = len(U0)//2 if len(U0) %2 == 0 else len(U0)//2 + 1
	print "idx: {}, len(X): {}".format(idx, len( U0[U0==1.0]))
	# solver
	for tt in range(tsteps - 1):
		# boundary conditions
		# first half time red and next half time green light condition
		U[idx, :tsteps//2]     = rho2U(BL)
		U[idx + 1, :tsteps//2] = rho2U(AL)

		U[0, :]  = U[1, :]
		U[-1, :] = U[-2, :]
		
		for xx in range(nsteps):
			if method == 'Upwind_Method':
				U[:, tt + 1] = Upwind_Method(F, FD, U[:, tt], k, h)
			elif method == 'Mac_Cornack_scheme':
				U[:, tt + 1] = Mac_Cornack_scheme(F, FD, U[:, tt], k, h)
			elif method == 'Lax_Friedrichs_scheme':
				U[:, tt + 1] = Lax_Friedrichs_scheme(F, FD, U[:, tt], k, h)
			elif method == 'Richtmyer_two_step_Lax_Wendroff_scheme':
				U[:, tt + 1] = Gudonov_Method(F, FD, FStarSolve, U[:, tt], k, h)
			elif method == 'Gudonov_Method':
				U[:, tt + 1] = Gudonov_Method(F, FD, FStarSolve, U[:, tt], k, h)

		if tt % 200 == 0:
			print "[INFO] tt: {}: Utt: {}".format(tt*k,list(U2rho(U[1:-1, tt])))


	return U



###################################################################
#                     Solve Burgers Equation                      #
###################################################################

legend_array = []
nsteps = int(len(X0)/h)
U = find_solution(U0, Tmax, nsteps, k, h)

plt.ion()
for tt in range(int(Tmax/k)):
	# if tt % 20 == 0:

	plt.clf()
	if tt > int(Tmax/k)//2:
		plt.plot(X0[1:-1], U2rho(U[1:-1, tt]), 'g')
		plt.title("GREEN Light: {}/{}".format(tt,int(Tmax/k)))
	else:
		plt.plot(X0[1:-1], U2rho(U[1:-1, tt]), 'r')
		plt.title("RED Light: {}/{}".format(tt,int(Tmax/k)))
        
	plt.ylim([-0.5, 1.5])
        plt.xlabel('x')
        plt.ylabel('density')
	plt.pause(0.0001)
	plt.savefig('../imgs/tl-gif/'+str(tt)+'.png')
        if tt % 100 == 99: 
            plt.savefig('../imgs/tl-'+str(tt)+'.png')
        	# legend_array.append('t = {}'.format(tt*k))
		# plt.legend(legend_array)	
plt.show()

import imageio
import os
png_dir = '../imgs/tl-gif/'
images = []
for tt in range(int(Tmax/k)):
    # if file_name.endswith('.png'):
    file_path = os.path.join(png_dir, str(tt) + '.png')
    images.append(imageio.imread(file_path))
imageio.mimsave('../imgs/tl-movie.gif', images, fps=50)
