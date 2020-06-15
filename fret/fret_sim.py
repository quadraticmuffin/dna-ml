import sys
import itertools as it
import numpy as np
import scipy.special, scipy.integrate, scipy.interpolate

PI = np.pi
AVOGADRO = 6.022e23 # molecules / mole
LIGHT_SPEED = 3e8 # m/s

def random_gaussian(num_elements):
  mean = np.random.uniform(0, num_elements)
  std = np.random.uniform(0,num_elements/4.)
  x = np.arange(num_elements)
  y = 1./np.sqrt(2*PI*std**2) * np.exp(-(x - mean)**2 / (2*std**2))
  return y / y.sum()

def load_spectrum(filepath, min_wavelength = 400, max_wavelength = 800, step_size = 1.0):
  wavelengths = np.arange(min_wavelength, max_wavelength+step_size/2, step_size)

  file = open(filepath)
  
  spectrum = []
  for line in file:
    if not line.startswith('#'):
      spectrum.append([max(0, float(v)) for v in line.split()])

  spectrum_interp = scipy.interpolate.interp1d(*zip(*spectrum), fill_value = 0, bounds_error=False)(wavelengths)
  return spectrum_interp

def load_emission_spectra(filepaths, min_wavelength = 400e-9, max_wavelength = 800e-9, step_size = 1.0e-9):
  # Assumes wavelengths in data files are in nanometers
  spectra = [load_spectrum(f, 10**9*min_wavelength, 10**9*max_wavelength, 10**9*step_size) for f in filepaths]

  # Normalize spectra
  spectra_norm = [spectrum/spectrum.sum()/step_size for spectrum in spectra]
  return spectra_norm

def load_absorption_spectra(filepaths, min_wavelength = 400e-9, max_wavelength = 800e-9, step_size = 1.0e-9):
  # Assumes wavelengths in data files are in nanometers
  spectra = [load_spectrum(f, 10**9*min_wavelength, 10**9*max_wavelength, 10**9*step_size) for f in filepaths]

  # Convert units of spectra from cm^-1 * M^-1 to m^2/mol
  return [0.1*spectrum for spectrum in spectra]
  

# Forster radius
def calc_forster_radius(donor_emission = None, acceptor_absorbance = None, wavelength_low = None, wavelength_high = None, kappa2 = 2./3, quant_yield = 1.0, refractive_index = 1.4):
  # donor_emission: 1D numpy array, units of inverse distance
  # acceptor_absorbance: 1D numpy array, in terms of wavelength, units of inverse concentration * inverse distance
  # wavelength_low, wavelength_high: float, units of distance
  # kappa2 (default 2/3): float, unitless
  # quant_yield (default 1.0): float, unitless
  # refractive_index: float, unitless

  assert(donor_emission is not None)
  assert(acceptor_absorbance is not None)
  assert(len(donor_emission) == len(acceptor_absorbance))

  coeff = 9. * np.log(10) * kappa2 * quant_yield / (128. * PI**5 * AVOGADRO * refractive_index**4)
  step_size = (wavelength_high - wavelength_low) / (len(donor_emission)-1)
  integral = scipy.integrate.simps(donor_emission * acceptor_absorbance * np.arange(wavelength_low, wavelength_high+step_size/2, step_size)**4) * step_size
  #print coeff, integral
  return (coeff * integral)**(1./6)

def calc_fret_rate_poisson3D(acceptor_conc, donor_lifetime, forster_radius = None, **kwargs):
  # acceptor_conc: float, units of concentration
  # donor_lifetime: float, units of time
  # if forster_radius not given, it is calculated with calc_forster_radius() using **kwargs

  # calculates using the formula:
  #   E = k / (k + 1/tau)

  E = calc_fret_efficiency_poisson3D(acceptor_conc, forster_radius, **kwargs)
  k = E/(1-E) / donor_lifetime
  return k


# Naive estimate of FRET transfer
def calc_fret_efficiency_poisson3D(acceptor_conc, forster_radius = None, **kwargs):
  # acceptor_conc: float, units of concentration
  # forster_radius (if given): float, units of distance
  #   if not given: additional keyword arguments are given to calc_forster_radius()

  if forster_radius is None:  forster_radius = forster_radius(**kwargs)

  c_0 = 3. / (2. * PI**(3./2) * AVOGADRO ** forster_radius**3)
  x = acceptor_conc / c_0
  return PI ** (.5) * x * np.exp(x**2) * scipy.special.erfc(x)
def calc_fret_efficiency_poisson3D_multiacceptor(acceptor_concs, forster_radii):
  # acceptor_concs: list of float, units of molar concentration
  # forster_radii: list of float, units of distance

  x = (2./3) * PI**(1.5) * np.array(forster_radii)**3 * np.array(acceptor_concs) * AVOGADRO
  #print x
  return PI**(.5) * x * np.exp(x.sum()**2) * scipy.special.erfc(x.sum())



def calc_transfer_efficiencies(fluorophore_concs, fluorophore_emissions, fluorophore_absorbances, quantum_yields, **kwargs):
  # fluorophore_concs: list of concentrations of each fluorophore
  # fluorophore_emissions: list of emission spectra for each fluorophore
  # fluorophore_absorbances: list of absorbances spectra for each fluorophore
  # **kwargs: keyword args to be given to calc_fret_efficiency_poisson3D

  num_fluorophores = len(fluorophore_concs)

  # Calculate matrix of transfer efficiencies
  # Columns are donors, rows acceptors
  transfer_matrix = np.empty(shape = (num_fluorophores, num_fluorophores))
  for i in range(num_fluorophores):
    args = kwargs.copy()
    args['quant_yield'] = quantum_yields[i]
    forster_radii = [calc_forster_radius(fluorophore_emissions[i],fluorophore_absorbances[j],**args) for j in range(num_fluorophores)]
    #print i, forster_radii
    transfer_matrix[:,i] = calc_fret_efficiency_poisson3D_multiacceptor(fluorophore_concs, forster_radii)
#  transfer_matrix /= transfer_matrix.sum(axis = 0) + donor_lifetimes

  P_adj = transfer_matrix - np.identity(num_fluorophores)
  B = np.zeros((num_fluorophores, num_fluorophores))
  for i in range(num_fluorophores):
    B[i,i] = transfer_matrix[:,i].sum() - 1

  Q = np.matmul(B, np.linalg.inv(P_adj))

  return transfer_matrix,Q
  

def calc_transfer_efficiencies_MC(iters, fluorophore_concs, fluorophore_emissions, fluorophore_absorbances, quantum_yields, **kwargs):

  num_fluorophores = len(fluorophore_concs)

  # Calculate forster radii between every pair of fluorophores
  # Column = donor, row = acceptor
  forster_radii = np.empty((num_fluorophores, num_fluorophores))
  for i,j in it.product(range(num_fluorophores), range(num_fluorophores)):
    args = kwargs.copy()
    args['quant_yield'] = quantum_yields[i]
    forster_radii[j][i] = calc_forster_radius(fluorophore_emissions[i], fluorophore_absorbances[j], **args)

  transfer_matrix_all = np.empty((num_fluorophores, num_fluorophores, iters))
  Q_all = np.empty((num_fluorophores, num_fluorophores, iters))
  # Perform <iters> simulations.
  # Each simulation involves:
  #  1. For each ordered pair of fluorophores (D,A):
  #    a. Choose a reaction volume (cube of side length 4*forster radius)
  #    b. Sample the number N_A of A molecules from Poisson distribution with mean (A_conc*(4r)**3)
  #    c. Randomly sample N_A positions for the A molecules uniformly in the cube.
  #    d. Calculate the FRET rate from D to each A molecule (assume D is at position (0,0,0))
  #      - Assume a donor lifetime of 1 since this gets normalized out anyway
  #    e. Store this efficiency in rates_matrix[D][A]
  #  2. Turn rates_matrix into a matrix of transfer_efficiencies:
  #    a. For each column (donor) calculate the sum of the rates in that column (total propensity)
  #    b. Divide each entry by 1+sum of that entry's column, and store in transfer_matrix
  #  3. Calculate the Q matrix (transfer_matrix**infty)
  #    a. Q = (P - I)^(-1)
  for iter in range(iters):
    rates_matrix = np.empty((num_fluorophores, num_fluorophores))
    for D,A in it.product(range(num_fluorophores), range(num_fluorophores)):
      forster_radius = forster_radii[A][D]
      box_sidelength = 4*forster_radius
      N_A_avg = AVOGADRO * fluorophore_concs[A]*box_sidelength**3
      N_A = np.random.poisson(lam = N_A_avg)
#      print D,A,N_A_avg,N_A
      A_pos = np.random.uniform(low = -box_sidelength/2, high = +box_sidelength/2, size = (N_A, 3))
      A_dist = ((A_pos**2).sum(axis=1))**(.5)
      rates_matrix[A][D] = ((forster_radius / A_dist)**6).sum()

    transfer_matrix = np.empty((num_fluorophores, num_fluorophores))
    for D in range(num_fluorophores):
      transfer_matrix[:,D] = rates_matrix[:,D] / (rates_matrix[:,D].sum() + 1)

    #print rates_matrix
    #print transfer_matrix - np.identity(num_fluorophores)
    transfer_matrix_adj = np.empty((num_fluorophores, num_fluorophores))
    for D,A in it.product(range(num_fluorophores),range(num_fluorophores)):
      if D != A:
        transfer_matrix_adj[A,D] = rates_matrix[A,D] / (rates_matrix[:,D].sum() + 1)
      else:
        transfer_matrix_adj[A,D] = (rates_matrix[A,D] - rates_matrix[:,D].sum() - 1) / (rates_matrix[:,D].sum() + 1)
    B = np.zeros((num_fluorophores, num_fluorophores))
    for i in range(num_fluorophores):  B[i,i] = 1./ (rates_matrix[:,i].sum() + 1)
    Q = -np.matmul(B, np.linalg.inv(transfer_matrix_adj))

    transfer_matrix_all[:,:,iter] = transfer_matrix
    Q_all[:,:,iter] = Q

    print "\r  {}%".format(100*iter/iters),
    sys.stdout.flush()
  
  return transfer_matrix_all.sum(axis=2)/iters, Q_all.sum(axis=2)/iters, transfer_matrix_all, Q_all

def calc_barcode(emissions, absorptions, wavelengths, excitation_wavelengths, detection_wavelengths, transfer_matrix, quantum_yields, concentrations, path_length = None, detection_width = 10e-9):
  emissions_interp = [scipy.interpolate.interp1d(wavelengths, e, kind='linear', bounds_error=False, fill_value=0) for e in emissions]
  absorptions_interp = [scipy.interpolate.interp1d(wavelengths, a, kind='linear', bounds_error=False, fill_value=0) for a in absorptions]

  step_size = wavelengths[1] - wavelengths[0]

  if path_length is None:
    # Average path length through a sphere of radius R is 4R/3
    # Radius of sphere is 0.5um
    path_length = 4.*(.5e-6)/3

  absorption_by_fluorophore = np.empty((len(detection_wavelengths), len(emissions)))
  emission_by_fluorophore = np.empty((len(detection_wavelengths), len(emissions)))
  for i, (excitation, detection) in enumerate(zip(excitation_wavelengths, detection_wavelengths)):
    absorbances = concentrations * np.array([interp(excitation) for interp in absorptions_interp]) * path_length
    net_absorbance = sum(absorbances)
    if net_absorbance > 0:
      fluorophore_absorption = np.reshape(absorbances / net_absorbance * (1 - 10**-net_absorbance), (len(absorbances), 1))
    else:
      fluorophore_absorption = np.zeros((len(absorbances), 1))
    fluorophore_emission = np.reshape(np.matmul(transfer_matrix, fluorophore_absorption), len(absorbances))
    #print transfer_matrix, fraction_absorbed, fluorophore_emission
    absorption_by_fluorophore[i,:] = np.reshape(fluorophore_absorption, len(absorbances))
    emission_by_fluorophore[i,:] = quantum_yields * np.array([scipy.integrate.simps(interp(np.arange(detection-detection_width/2, detection+detection_width/2, step_size))) for interp in emissions_interp]) * step_size * fluorophore_emission

  net_emission = emission_by_fluorophore.sum(axis = 1)

  return net_emission, absorption_by_fluorophore, emission_by_fluorophore

def plot_transfer_matrix(mat, fluorophores, xlabel='donor', ylabel='acceptor', title = None):
  plt.figure()
  plt.imshow(mat, interpolation='nearest', vmin = 0, vmax = 1)
  cbar = plt.colorbar()
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if title is not None:  plt.title(title)
  plt.xticks(np.arange(len(fluorophores)), fluorophores)
  plt.yticks(np.arange(len(fluorophores)), fluorophores)
  cbar.set_label('probability')


## TEST DATA
NUM_FLUOROPHORES = 3
FLUOROPHORES = ['fluorescein','rhodamine','texasred']

#CONCENTRATIONS = np.array([1e3, 1e3, 1e3]) # moles per m^3
LAMBDA_MIN = 400e-9 # meters
LAMBDA_MAX = 800e-9 # meters
LAMBDA = np.arange(LAMBDA_MIN, LAMBDA_MAX+5e-10, 1e-9)
EMISSIONS = load_emission_spectra(['{}_emission.txt'.format(f) for f in FLUOROPHORES], LAMBDA_MIN, LAMBDA_MAX)
ABSORPTIONS = load_absorption_spectra(['{}_absorption.txt'.format(f) for f in FLUOROPHORES], LAMBDA_MIN, LAMBDA_MAX)
QUANTUM_YIELDS = np.array([0.79, 0.7, 0.9]) #[np.random.uniform(0.2,0.8) for _ in range(NUM_FLUOROPHORES)]
KAPPA2 = 2./3
REFRACTIVE_INDEX = 1.4


import matplotlib.pyplot as plt
from matplotlib.colors import Colormap as cm
plt.ion()

plt.figure()
for i in range(NUM_FLUOROPHORES):
  plt.subplot(3,1,i+1)
  plt.plot(LAMBDA, EMISSIONS[i]/EMISSIONS[i].max())
  plt.plot(LAMBDA, ABSORPTIONS[i]/ABSORPTIONS[i].max())
  plt.xlabel(FLUOROPHORES[i])
  plt.ylabel('a.u.')

#transfer_matrix_analytic,Q_analytic = calc_transfer_efficiencies(CONCENTRATIONS, EMISSIONS, ABSORPTIONS, QUANTUM_YIELDS, wavelength_low=LAMBDA_MIN, wavelength_high=LAMBDA_MAX, kappa2=KAPPA2, refractive_index=REFRACTIVE_INDEX)
#plot_transfer_matrix(transfer_matrix_analytic, ['F','R','T'], xlabel = 'donor', ylabel = 'acceptor', title = 'Single Transfer (analytic)')
#plot_transfer_matrix(Q_analytic, ['F','R','T'], xlabel = 'first donor', ylabel = 'final acceptor', title = 'Overall Transfer (analytic)')

#passes = 13
#transfer_matrix_MC,Q_MC,_,_ = calc_transfer_efficiencies_MC(1, CONCENTRATIONS, EMISSIONS, ABSORPTIONS, QUANTUM_YIELDS, wavelength_low=LAMBDA_MIN, wavelength_high=LAMBDA_MAX, kappa2=KAPPA2, refractive_index=REFRACTIVE_INDEX)
#for p in range(passes):
#  transfer_matrix_MC_temp,Q_MC_temp,_,_ = calc_transfer_efficiencies_MC(2**p, CONCENTRATIONS, EMISSIONS, ABSORPTIONS, QUANTUM_YIELDS, wavelength_low=LAMBDA_MIN, wavelength_high=LAMBDA_MAX, kappa2=KAPPA2, refractive_index=REFRACTIVE_INDEX)
#  transfer_matrix_MC_new = (transfer_matrix_MC + transfer_matrix_MC_temp)/2
#  Q_MC = (Q_MC + Q_MC_temp)/2
#  print ((transfer_matrix_MC_new - transfer_matrix_MC)**2).sum()
#  transfer_matrix_MC = transfer_matrix_MC_new
#
#plot_transfer_matrix(transfer_matrix_MC, ['F','R','T'], xlabel = 'donor', ylabel = 'acceptor', title = 'Single Transfer (MC)')
#
#plot_transfer_matrix(Q_MC, ['F','R','T'], xlabel = 'first donor', ylabel = 'final acceptor', title = 'Overall Transfer (MC)')



total_conc = 1000.0 # moles/m^3

def label_to_concs(frt):
  fw = 1. if 'f' in frt else 5. if 'F' in frt else 0
  rw = 1. if 'r' in frt else 5. if 'R' in frt else 0
  tw = 1. if 't' in frt else 5. if 'T' in frt else 0
  return [fw/(fw+rw+tw), rw/(fw+rw+tw), tw/(fw+rw+tw)]
  

#labels = [f+r+t for f,r,t in it.product(['','F'],['','R'],['','T']) if f+r+t!='']
labels = [f+r+t for f,r,t in it.product(['','f','F'],['','r','R'],['','t','T']) if f+r+t!='']
concentrations = total_conc * np.array([label_to_concs(l) for l in labels])
#concentrations = total_conc * np.array([[f/(f+r+t), r/(f+r+t), t/(f+r+t)] for f,r,t in it.product([0,1.,5.],[0,1.,5.],[0,1.,5.]) if (f,r,t)!=(0,0,0)])
#concentrations = total_conc * np.array([[f/(f+r+t), r/(f+r+t), t/(f+r+t)] for f,r,t in it.product([0,5.],[0,5.],[0,5.]) if (f,r,t)!=(0,0,0)])
MC_iters = 2**13

excitation_wavelengths = np.array([488, 532, 561]) * 1e-9
detection_wavelengths = np.array([513, 576, 622]) * 1e-9
barcodes = np.empty((len(concentrations), len(excitation_wavelengths)))

for i in range(concentrations.shape[0]):
  concs = concentrations[i,:]

  print "{}: Testing concentrations {}".format(i, concs)

  transfer_matrix_analytic,Q_analytic = calc_transfer_efficiencies(concs, EMISSIONS, ABSORPTIONS, QUANTUM_YIELDS, wavelength_low=LAMBDA_MIN, wavelength_high=LAMBDA_MAX, kappa2=KAPPA2, refractive_index=REFRACTIVE_INDEX)
  transfer_matrix_MC,Q_MC,_,_ = calc_transfer_efficiencies_MC(MC_iters, concs, EMISSIONS, ABSORPTIONS, QUANTUM_YIELDS, wavelength_low=LAMBDA_MIN, wavelength_high=LAMBDA_MAX, kappa2=KAPPA2, refractive_index=REFRACTIVE_INDEX)

  print "  Quality check: analytic transfer matches MC with error {}".format(((transfer_matrix_analytic - transfer_matrix_MC)**2).sum())

  barcode, fluorophore_absorptions, fluorophore_emissions = calc_barcode(
    emissions = EMISSIONS,
    absorptions = ABSORPTIONS,
    wavelengths = LAMBDA,
    excitation_wavelengths = excitation_wavelengths,
    detection_wavelengths = detection_wavelengths,
    transfer_matrix = Q_MC,
    quantum_yields = QUANTUM_YIELDS,
    concentrations = concs
  )
  barcodes[i,:] = barcode

  print "  Barcode: {}".format(barcode)

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xs,ys,zs = zip(*barcodes)
for (f,r,t),(x,y,z) in zip(concentrations, barcodes):
  ax.scatter(x,y,zs=z, color=(f/total_conc, r/total_conc, t/total_conc))
ax.legend(labels)
ax.set_xlabel('{}nm'.format(detection_wavelengths[0]*1e9))
ax.set_ylabel('{}nm'.format(detection_wavelengths[1]*1e9))
ax.set_zlabel('{}nm'.format(detection_wavelengths[2]*1e9))
#ax.scatter(xs,ys,zs=zs)
  

