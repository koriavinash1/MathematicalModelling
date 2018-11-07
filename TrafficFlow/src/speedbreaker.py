import numpy as np
import matplotlib.pyplot as plt


###################################################################
#               Initialize Space and Time steps                   #
#                                                                 #
###################################################################
k = 0.005
h = 0.05
Tmax = 1.0

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
vmax   = 1.0
F = lambda rho: rho*vmax*(1.0 - rho/rhoMax)
FD = lambda rho: vmax - 2*rho/rhoMax
FStarSolve = lambda: vmax*rhoMax/2.0
 
v2rho = lambda v: (-1*v/vmax + 1.0)*rhoMax/2.0 
rho2v = lambda rho: vmax*(1.0 - 2*rho/rhoMax) 

U2rho = lambda U: (1.0 - U)*rhoMax/2.0
rho2U = lambda rho: (1.0 - rho*2.0/rhoMax)
###################################################################
#                 Initial Condition Defination                    #
#	                 Conditions on Rho not U                      #
###################################################################
def initial_conditions(x):
	if x <= 0.0:
		return 1.0
	elif x > 0.0:
		return 0.0
BSB, ASB = 1.0, 0.1	 # rho before and after speedBreaker

# initial rho ...
# rho0 = 0.55*np.ones(len(X0))
rho0 = np.array([initial_conditions(x) for x in X0])


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
	return U


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


def find_solution(rho0, T, nsteps, k, h, method='Gudonov_Method'):
	tsteps = int(T/k) + 1
	rho = np.zeros((len(rho0), tsteps))
	rho[:, 0] = rho0
	rho[0, :] = rho[0, 0]  
	idx = len(rho0)//2 if len(rho0) %2 == 0 else len(rho0)//2 + 1
	print "idx: {}, len(X): {}".format(idx, len( rho0[rho0==1.0]))
	# solver
	for tt in range(tsteps - 1):
		# boundary conditions
		rho[idx, :] = v2rho(BSB)
		rho[idx + 1, :] = v2rho(ASB)

		rho[0, :]  = rho[1, :]
		rho[-1, :] = rho[-2, :]
		
		for xx in range(nsteps):
			if method == 'Upwind_Method':
				rho[:, tt + 1] = Upwind_Method(F, FD, rho[:, tt], k, h)
			elif method == 'Mac_Cornack_scheme':
				rho[:, tt + 1] = Mac_Cornack_scheme(F, FD, rho[:, tt], k, h)
			elif method == 'Lax_Friedrichs_scheme':
				rho[:, tt + 1] = Lax_Friedrichs_scheme(F, FD, rho[:, tt], k, h)
			elif method == 'Richtmyer_two_step_Lax_Wendroff_scheme':
				U[:, tt + 1] = Gudonov_Method(F, FD, FStarSolve, U[:, tt], k, h)
			elif method == 'Gudonov_Method':
				rho[:, tt + 1] = Gudonov_Method(F, FD, FStarSolve, rho[:, tt], k, h)

		if tt % 20 == 0:
			print "[INFO] tt: {}: Utt: {}".format(tt*k,list(rho[1:-1, tt]))


	return rho


# def find_solution(U0, T, nsteps, k, h, method='Gudonov_Method'):
# 	tsteps = int(T/k) + 1
# 	U = np.zeros((len(U0), tsteps))
# 	U[:, 0] = U0
# 	U[0, :] = U[0, 0]  
# 	idx = len(U0)//2 if len(U0) %2 == 0 else len(U0)//2 + 1
# 	print "idx: {}, len(X): {}".format(idx, len( U0[U0==1.0]))
# 	# solver
# 	for tt in range(tsteps - 1):
# 		# boundary conditions
# 		U[idx, :]     = rho2U(v2rho(BSB))
# 		U[idx + 1, :] = rho2U(v2rho(ASB))

# 		U[0, :]  = U[1, :]
# 		U[-1, :] = U[-2, :]
		
# 		for xx in range(nsteps):
# 			if method == 'Upwind_Method':
# 				U[:, tt + 1] = Upwind_Method(F, FD, U[:, tt], k, h)
# 			elif method == 'Mac_Cornack_scheme':
# 				U[:, tt + 1] = Mac_Cornack_scheme(F, FD, U[:, tt], k, h)
# 			elif method == 'Lax_Friedrichs_scheme':
# 				U[:, tt + 1] = Lax_Friedrichs_scheme(F, FD, U[:, tt], k, h)
# 			elif method == 'Richtmyer_two_step_Lax_Wendroff_scheme':
# 				U[:, tt + 1] = Gudonov_Method(F, FD, FStarSolve, U[:, tt], k, h)
# 			elif method == 'Gudonov_Method':
# 				U[:, tt + 1] = Gudonov_Method(F, FD, FStarSolve, U[:, tt], k, h)

# 		if tt % 200 == 0:
# 			print "[INFO] tt: {}: Utt: {}".format(tt*k,list(U2rho(U[1:-1, tt])))


# 	return U

###################################################################
#                     Solve Burgers Equation                      #
###################################################################
legend_array = []
nsteps = int(len(X0)/h)
rho = find_solution(rho0, Tmax, nsteps, k, h)


plt.ion()
for tt in range(int(Tmax/k)):
	plt.clf()
	plt.plot(X0[1:-1], rho[1:-1, tt])
	plt.title("Speed Breaker: {}/{}".format(tt,int(Tmax/k)))
	plt.ylim([-0.5, 1.5])
	plt.pause(0.001)
        plt.xlabel('x')
        plt.ylabel('density')
        plt.savefig('../imgs/sb-gif/'+str(tt)+'.png')
        if tt % 100 == 99: 
            plt.savefig('../imgs/sb-'+str(tt)+'.png')
	# legend_array.append('t = {}'.format(tt*k))
	# plt.legend(legend_array)	
# plt.show()


import imageio
import os
png_dir = '../imgs/sb-gif/'
images = []
for tt in range(int(Tmax/k)):
    # if file_name.endswith('.png'):
    file_path = os.path.join(png_dir, str(tt) + '.png')
    images.append(imageio.imread(file_path))
imageio.mimsave('../imgs/sb-movie.gif', images, fps=50)