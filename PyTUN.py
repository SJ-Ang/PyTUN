#!/usr/bin/python

# Comments are incredibly appreciated (e-mail: a0072211@u.nus.edu)
#######################################################################
#                              PyTUN.py                               #
#  Evaluation of Tunneling Rates and Associated Effective Barrier     #
#        heights from Gaussian 09 using Eckart Potential.             #
#######################################################################
####################  Written by: Shi Jun Ang #########################
################  Last modified:  Dec 04, 2017 ########################
#######################################################################

import sys, math, time
from array import array

# PHYSICAL CONSTANTS
GAS_CONSTANT = 8.3144621
PLANCK_CONSTANT = 6.62606957e-34
BOLTZMANN_CONSTANT = 1.3806488e-23
SPEED_OF_LIGHT = 2.99792458e10
AVOGADRO_CONSTANT = 6.0221415e23
AMU_to_KG = 1.66053886E-27
autokcal = 627.509541
kjtokcal = 4.184
atmos = 101.325
PI = 3.14159265359
k = 3.1668114E-6 #Boltzmann Constant in atomic units 

#10-point Gauss-Legendre Quadrature abscissa and weight (exact solution for up to 21st order polynomial)
x = array('d',[-0.9739065285,-0.8650633667,-0.6794095683,-0.4333953941,-0.1488743390,0.1488743390,0.4333953941,0.6794095683,0.8650633667,0.9739065285])
w = array('d',[0.0666713443,0.1494513492,0.2190863625,0.2692667193,0.2955242247,0.2955242247,0.2692667193,0.2190863625,0.1494513492,0.0666713443])


#Read gaussian output for the level of theory and basis set used
def level_of_theory(file):
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()
   level = "none"
   bs = "none"
   for i in range(0,len(inlines)):
      if inlines[i].strip().find('\\Freq\\') > -1:
          if len(inlines[i].strip().split("\\")) > 5:
              level = (inlines[i].strip().split("\\")[4])
              bs = (inlines[i].strip().split("\\")[5])
   return level+"/"+bs

#Read gaussian output for the Final SCF Energy
def final_scf_energy(file):
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()
   scf = "none"
   for i in range(0,len(inlines)):
      if inlines[i].strip().startswith('SCF Done:'):
              scf = (inlines[i].strip().split()[4])
   return scf

#Read gaussian output for Zero Point Energy
def zero_point_energy(file):
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()
   zpe = "none"
   for i in range(0,len(inlines)):
      if inlines[i].strip().startswith('Zero-point correction'):
              zpe = (inlines[i].strip().split()[2])
   return zpe

#Read gaussian output for Thermal Correction to Gibbs Free Energy
def g_corr(file):
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()
   gcorr = "none"
   for i in range(0,len(inlines)):
      if inlines[i].strip().startswith('Thermal correction to Gibbs Free Energy='):
              gcorr = (inlines[i].strip().split()[6])
   return gcorr

#Read gaussian output for Temperature
def temperature(file):
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()
   temperature = "none"
   for i in range(0,len(inlines)):
      if inlines[i].strip().startswith('Temperature'):
              temperature = (inlines[i].strip().split()[1])
   return temperature

#Read gaussian output for reduced mass of normal mode of transition
def reduced_mass(file):
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()
   mu = "none"
   for i in range(0,len(inlines)):
      if inlines[i].strip().startswith('Red. masses'):
	if mu is "none":
          mu = (inlines[i].strip().split()[3])
   return mu

#Read gaussian output for force constant of normal mode of transition
def force_constant(file):
   g09_output = open(file, 'r')
   inlines = g09_output.readlines()
   fc = "none"
   for i in range(0,len(inlines)):
      if inlines[i].strip().startswith('Frc consts'):
	if fc is "none":
          fc = (inlines[i].strip().split()[3])
   return fc

#Parameters B, ALPHA, a, b, d of Eckart Potential
def Bee(V_max,V_r, V_p):
	bee = (V_max ** 0.5 + (V_max - (V_p - V_r)) ** 0.5) ** 2
	return bee

def ALPHA(B,F_s,V_max,V_r, V_p):
	alpha = (B * F_s / (2 * V_max * (V_max - (V_p - V_r)))) ** 0.5
	return alpha

def A(E,mu,ALPHA):
	a =  2 * PI * (2 * mu * E)**0.5 / ALPHA
	return a

def B(E,mu,V_p,V_r,ALPHA):
	b =  2 * PI * (2 * mu * (E - (V_p - V_r)))**0.5 / ALPHA
	return b

def D(bee,mu,ALPHA):
	d =  2 * PI * abs((2 * mu * bee - (ALPHA/2)**2))**0.5 / ALPHA
	return d

#Calculation of Transmission Probabilty of Eckart Potential
def T(a,b,d):
	T = (math.cosh(a+b) - math.cosh(a-b))/(math.cosh(a+b) + math.cosh(d))
	return T

#Calculation of SINH function of Kappa
def S(V_max,E):
	S = math.sinh(((V_max-E)) / (TEMPERATURE*k))
	return S
TEMPERATURE = float(temperature('ts.out'))
#TEMPERATURE = 200
V_r = float(final_scf_energy('reactant_sp.out')) + float(zero_point_energy('reactant.out'))
V_p = float(final_scf_energy('product_sp.out')) + float(zero_point_energy('product.out'))
V_max = float(final_scf_energy('ts_sp.out')) + float(zero_point_energy('ts.out'))
G_r = float(final_scf_energy('reactant_sp.out')) + float(g_corr('reactant.out'))
G_p = float(final_scf_energy('product_sp.out')) + float(g_corr('product.out'))
G_ts = float(final_scf_energy('ts_sp.out')) + float(g_corr('ts.out'))

if V_r > V_p:
	E_o = V_r
else:
	E_o = V_p

#Scaling of Energies(define V_r == 0)
V_max = V_max - V_r
V_p = V_p - V_r
E_o = E_o - V_r
V_r = V_r - V_r
#V_p = 0.000/autokcal
#V_r = 0.00
#V_max = 9.80/autokcal
#E_o = 0.00/autokcal
y = (V_max - E_o)/2.0 
z = (V_max + E_o)/2.0 

# Specifing Parameters for the Eckart Potential
mu = float(reduced_mass('ts.out'))*1836
F_s = float(force_constant('ts.out'))/15.569141
bee = Bee(V_max,V_r,V_p)
alpha = ALPHA(bee,F_s,V_max,V_r,V_p)
d = D(bee,mu,alpha)

#Calculation of Uncorrected Gibbs Free Energy Barrier in kJ/mol and Rate
delta_g_dagg = (G_ts - G_r)*2625.5
print "Uncorrected delta G dagger is %f" % delta_g_dagg
ln_uncorr_rate = math.log(BOLTZMANN_CONSTANT * TEMPERATURE / (PLANCK_CONSTANT))-delta_g_dagg*1000/(GAS_CONSTANT*TEMPERATURE)  
uncorr_rate = "{:.10E}".format(math.exp(ln_uncorr_rate))
print "Uncorrected rate is %s" % uncorr_rate

#Calculation of Wigner tunneling correction
kappa_w = 1+(( (F_s/mu)**0.5 / (k*TEMPERATURE) )**2) / 24
print "The Wigner tunneling correction is %f" % kappa_w

#Calculation of Skodje tunneling correction
alpha_s = (2*PI/((F_s/mu)**0.5))
beta_s = (1/ (k*TEMPERATURE))
if V_p < V_r:
	Vee = 0
else:
	Vee = V_p-V_r
if alpha_s > beta_s:
	kappa_s = (beta_s*PI/alpha_s)/(math.sin(beta_s*PI/alpha_s))+(beta_s/(alpha_s-beta_s)*math.exp((beta_s-alpha_s)*(V_max-Vee)))
else:
	kappa_s = (beta_s/(beta_s-alpha_s)*(math.exp((beta_s-alpha_s)*(V_max-Vee))-1))
print "The Skodje tunneling correction is %f" % kappa_s


#Calculation of Eckart tunneling correction using 10-point Gauss-Legendre Quadrature
kappa = 1
for i in range(0,10):
	a = A((x[i] * y + z),mu,alpha)
	b = B((x[i] * y + z),mu,V_p,V_r,alpha)
	kappa = (2 * y  / (TEMPERATURE*k) * w[i] * S((V_max),(x[i] * y + z)) * T(a,b,d)) + kappa
print "The Eckart tunneling correction is %f" % kappa

#Calculation of Wigner Apparent Gibbs Free Energy Barrier in kJ/mol and Rate
print "For Wigner Tunneling:"
corr_rate_w = "{:.10E}".format(kappa_w*math.exp(ln_uncorr_rate))
print "Corrected rate is %s" % corr_rate_w
delta_g_dagg_app_w = GAS_CONSTANT * TEMPERATURE * (math.log(BOLTZMANN_CONSTANT * TEMPERATURE / (PLANCK_CONSTANT)) - math.log(kappa_w) - ln_uncorr_rate) /1000
print "Apparent delta G dagger is %f" % delta_g_dagg_app_w

#Calculation of Skodje Apparent Gibbs Free Energy Barrier in kJ/mol and Rate
print "For Skodje Tunneling:"
corr_rate_s = "{:.10E}".format(kappa_s*math.exp(ln_uncorr_rate))
print "Corrected rate is %s" % corr_rate_s
delta_g_dagg_app_s = GAS_CONSTANT * TEMPERATURE * (math.log(BOLTZMANN_CONSTANT * TEMPERATURE / (PLANCK_CONSTANT)) - math.log(kappa_s) - ln_uncorr_rate) /1000
print "Apparent delta G dagger is %f" % delta_g_dagg_app_s

#Calculation of Eckart Apparent Gibbs Free Energy Barrier in kJ/mol and Rate
print "For Eckart Tunneling:"
corr_rate = "{:.10E}".format(kappa*math.exp(ln_uncorr_rate))
print "Corrected rate is %s" % corr_rate
delta_g_dagg_app = GAS_CONSTANT * TEMPERATURE * (math.log(BOLTZMANN_CONSTANT * TEMPERATURE / (PLANCK_CONSTANT)) - math.log(kappa) - ln_uncorr_rate) /1000
print "Apparent delta G dagger is %f" % delta_g_dagg_app
