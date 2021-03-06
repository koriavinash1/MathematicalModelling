import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

Xmin   = -0.2
Xmax   = 1.6
deltaX = 0.05  
X0      = np.linspace(Xmin, Xmax, int((Xmax - Xmin)/deltaX))

Tmin = 0
Tmax = 2
deltaT = 0.15
T      = np.linspace(Tmin, Tmax, int((Tmax - Tmin)/deltaT))

# U_t + F(U(X,t))_x = 0
# U(X, T) is constant

# initilize F
F = lambda U: U**2/2
FD = lambda U: U
FStarSolve = lambda: 0

def initial_conditions(x):
	if x <= 0.0:
		return 1.0
	elif 0.0< x <=1.0:
		return 1.0 - x
	elif x > 0.0:
		return 0.0

def characteristic_solution(x0, t):
	x = []
	for i in range(len(x0)):
		x.append(initial_conditions(x0[i])*t + x0[i])
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
plot_characteristic(characteristic_solution(X0, T), T)

def Upwind_Method(F, FD, U, k, h):
	for i in range(1, len(U) - 1):
		if FD(U[i]) >= 0:
			U[i] = U[i] - (k/h)*(F(U[i]) - F(U[i-1]))
		else:
			U[i] = U[i] - (k/h)*(F(U[i +1]) - F(U[i]))
	return U

def Lax_Friedrichs_scheme(F, FD, U, k, h):
	for i in range(1, len(U) - 1):
		U[i] =0.5*(U[i+1] + U[i-1]) - (k/(2.*h))*(F(U[i +1]) - F(U[i]))
	return U

def Mac_Cornack_scheme(F, FD, U, k, h):
	Ustar = lambda u, up1: u - (k/h)*(F(up1) - F(u))
	for i in range(1, len(U) - 1):
		U[i] =0.5*(U[i] + Ustar(U[i], U[i+1])) -\
		 (k/h)*(F(Ustar(U[i], U[i+1])) - F(Ustar(U[i-1], U[i])))
	return U

def Richtmyer_two_step_Lax_Wendroff_scheme(F, FD, U, k, h):
	Ustar = lambda u, up1: 0.5*(u + up1) - (k/h)*(F(up1) - F(u))
	for i in range(1, len(U) - 1):
		U[i] = U[i] - (k/h)*(F(Ustar(U[i], U[i+1])) - F(Ustar(U[i-1], U[i])))
	return U

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

	for i in range(1, len(U) - 1):
		U[i] = U[i] -(k/h)*(F(Ustar(U[i], U[i+1])) - F(Ustar(U[i-1], U[i])))
	return U


def find_solution(U0, T, nsteps, k, h, method='Richtmyer_two_step_Lax_Wendroff_scheme'):
	tsteps = int(T/k) + 1
	U = np.zeros((len(U0), tsteps))
	U[:, 0] = U0

	# print ("[INFO] tsteps: {}, nsteps: {}".format(tsteps, nsteps))
	for tt in tqdm(range(tsteps - 1)):
		for xx in range(nsteps):
			# print "[Secondary INFO] tt: {}: Utt: {}".format(tt, U[:, tt])
			if method == 'upwind':
				U[:, tt + 1] = Upwind_Method(F, FD, U[:, tt], k, h)
			elif method == 'Mac_Cornack_scheme':
				U[:, tt + 1] = Mac_Cornack_scheme(F, FD, U[:, tt], k, h)
			elif method == 'Lax_Friedrichs_scheme':
				U[:, tt + 1] = Lax_Friedrichs_scheme(F, FD, U[:, tt], k, h)
			elif method == 'Richtmyer_two_step_Lax_Wendroff_scheme':
				U[:, tt + 1] = Gudonov_Method(F, FD, FStarSolve, U[:, tt], k, h)
			elif method == 'Gudonov_Method':
				U[:, tt + 1] = Gudonov_Method(F, FD, FStarSolve, U[:, tt], k, h)

	return U




U0 = np.array([initial_conditions(x) for x in X0])
legend_array = []
for T in [0, 0.005, 0.5, 1.0]:
	U = find_solution(U0, T, 2, 0.005, 0.14)
	plt.plot(X0, U[:, -1])
	legend_array.append('t = {}'.format(T))
plt.legend(legend_array)	
plt.show()
"""
legend_array = []
k = 0.005
h = 0.05
T = 1.0
U0 = np.array([initial_conditions(x) for x in X0])
U = find_solution(U0, T, 25, k, h)


# plt.ion()
for tt in range(int(T/k)):
	if tt % 20 == 0:
		# plt.clf()
		plt.plot(X0[1:-1], (1.0 - U[1:-1, tt])/2.0)
		# plt.pause(0.2)
		legend_array.append('t = {}'.format(tt*k))
		plt.legend(legend_array)	
plt.show()
"""
