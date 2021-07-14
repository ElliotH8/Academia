#!/usr/bin/env python
# coding: utf-8

# # RAYS FOR EXCELLENCE 2021

# # Fitting Komrad data with Band-function

# In[952]:


#Necessary things to import
from threeML import *
from threeML.io.package_data import get_path_of_data_file

from astromodels import *
from astromodels.functions.function import Function1D, FunctionMeta, ModelAssertionViolation
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from threeML.plugins.DispersionSpectrumLike import DispersionSpectrumLike
from threeML.utils.OGIP.response import OGIPResponse
import random

import warnings
warnings.simplefilter('ignore')

import os
os.environ["MKL_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'
os.environ["NUMEXPR_NUM_THREADS"] = '1'

silence_logs()


# ## Reading the Komrad spectra

# ### Function that reads each spectral file

# In[953]:


####################################
## READ TEXT FILE INTO PARAMETERS ##
####################################
def read_file(file_name):

    f = open(file_name)
    header_line = True
    skip_next_line = False
    i = 0
    for line in f.readlines():
        line_s = line.split()
        if skip_next_line == True:
            skip_next_line = False

        elif header_line == True:
            if len(line_s) == 5:
                if str(line_s[1]) == "ENERGIES" and str(line_s[3]) == "SPECTRUM":
                    header_line = False
                    skip_next_line = True
            elif "nBins" == str(line_s[0]):
                nBins = int(line_s[2])
                energies = np.zeros(nBins)
                spectra  = np.zeros(nBins)
            elif "tau" == str(line_s[0]):
                tau_init = float(line_s[2])
            elif "theta_RMS" == str(line_s[0]):
                theta_R = float(line_s[2])
            elif "R_theta" == str(line_s[0]):
                R = float(line_s[2])
            elif "y_RMS" == str(line_s[0]):
                yr = float(line_s[2])
            elif "<epsilon>" == str(line_s[0]):
                mean_epsilon = float(line_s[6])
            elif "Spectral" == str(line_s[0]):
                if "Nnu" == str(line_s[3]):
                    spectralInput = 0
                elif "Fnu" == str(line_s[3]):
                    spectralInput = 1
                elif "nuFnu" == str(line_s[3]):
                    spectralInput = 2
        else:
            energies[i] = float(line_s[0])
            spectra[i]  = float(line_s[1])
            i = i+1
    f.close

    return([nBins, tau_init, theta_R, R, yr,
                mean_epsilon, spectralInput, energies, spectra])


# ### Function that reads the parent file

# In[954]:


def read_parent_file(parent_file, save_name, model_name, renormalize):

    f = open(parent_file)
    nSpectra = 0
    #For loop to get the number of spectra
    for line in f.readlines():
        nSpectra = nSpectra+1
    f.close

    tauThetas = np.zeros(nSpectra)
    Rs        = np.zeros(nSpectra)
    yrs       = np.zeros(nSpectra)
    gridSize  = int(round(nSpectra**(1./3.)))

    #Read all spectra
    i = 0
    f = open(parent_file)
    for line in f.readlines():
        [nBins, tau_init, theta_R, R, yr,
                    mean_epsilon, spectralInput, energies, spectra] = read_file(line[:-1])
        #Create Nss that holds all imported spectra
        if i == 0:
            Nss = np.zeros((nSpectra, nBins))
        #Adds the parameters for each spectra
        tauThetas[i] = tau_init*theta_R
        Rs[i]        = R
        yrs[i]       = yr
        #Adds the spectra to Nss
        Nss[i]       = spectra
        i = i+1
    f.close
    
    if renormalize:
        index_energy_1 = np.where(energies == 1)
        
    #The software needs to create a "TemplateModelFactory" to work with.
    temp_factory = TemplateModelFactory(save_name, model_name, energies, ['tauTheta', 'R', 'yr'])
    tauTheta_grid  = np.logspace(np.log10(min(tauThetas)), np.log10(max(tauThetas)), gridSize)
    R_grid         = np.logspace(np.log10(min(Rs)),        np.log10(max(Rs)),        gridSize)
    yr_grid        = np.logspace(np.log10(min(yrs)),       np.log10(max(yrs)),       gridSize)
    temp_factory.define_parameter_grid('tauTheta', tauTheta_grid)
    temp_factory.define_parameter_grid('R', R_grid)
    temp_factory.define_parameter_grid('yr', yr_grid)
    
    #Below adds all spectra from Komrad into the TemplateModel
    n = 0
    fig, ax = plt.subplots(1,figsize=(8, 5))
    Nss_normalized = []
    for tau in tauTheta_grid:
        for R in R_grid:
            for yr in yr_grid:
                if renormalize:
                    Nss_normalized.append(Nss[n]/Nss[n][index_energy_1])
                    temp_factory.add_interpolation_data(Nss[n]/Nss[n][index_energy_1], tauTheta=tau, R=R, yr=yr)
                    #Makes a plot of the spectra
                    ax.plot(energies, Nss[n]/Nss[n][index_energy_1], color = (n*1.0/(1.3*nBins), 0.1, 0.1))
                else:
                    temp_factory.add_interpolation_data(Nss[n], tauTheta=tau, R=R, yr=yr)
                    #Makes a plot of the spectra
                    ax.plot(energies, Nss[n], color = (n*1.0/(1.3*nBins), 0.1, 0.1))
                n = n+1
    temp_factory.save_data(overwrite=True)
    #Below only relevant for the plot of the spectra
    ax.set_xlim([1e-1, 1e2])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([1e-4, 2])
    plt.show()
    
    if renormalize:
        return([nBins , nSpectra, energies, Nss_normalized, tauThetas, Rs, yrs])
    else:
        return([nBins , nSpectra, energies, Nss, tauThetas, Rs, yrs])


# ### Inporting the Komrad spectra

# In[955]:


#parent_file should contain the path to where the file is located at your computer. 
#If it is in the same folder as the script the full path is not needed, parent_file is simply"Parentfile.txt"
parent_file = "Parentfile.txt"
[nBins, nSpectra, energies, Nss, tauTheta_grid, R_grid, yr_grid] = read_parent_file(parent_file, save_name="RaysModel", model_name="RMS_Model", renormalize=True)

RMS_Model = TemplateModel('RaysModel')


# # Generating synthetic spectra

# In[956]:


#Import response function of GBM
response = OGIPResponse("glg_cspec_n4_bn180427442_v01.rsp2")

#Set the background function to be a powerlaw in energy
initial_K = 10^3
background_function = Powerlaw(K=initial_K, index=-1.5, piv=100.0)
SNR = 100000
RMS_Model.K = initial_K*SNR
RMS_Model.scale = 10.0
RMS_K = RMS_Model.K.value

#The arrays below will contain the paramters for each randomly generated spectra
nr_fake_spectra = 100
random_tauThetas = np.zeros(nr_fake_spectra)
random_Rs        = np.zeros(nr_fake_spectra)
random_yrs       = np.zeros(nr_fake_spectra)
scales           = np.zeros(nr_fake_spectra)
fake_spectra = []

#Randomizes the spectra
for i in range(nr_fake_spectra):
    #Picks out random values (in log-space) between the min and the max
    #random_tauThetas[i]  = np.exp(np.log(min(tauTheta_grid)) + (np.log(max(tauTheta_grid)) - np.log(min(tauTheta_grid)))*random.random())
    random_tauThetas[i] = 5
    #random_Rs[i]         = np.exp(np.log(min(R_grid))        + (np.log(max(R_grid))        - np.log(min(R_grid)))*random.random())
    random_Rs[i] = np.array([30+2.7*i])
    #random_yrs[i]        = np.exp(np.log(min(yr_grid))       + (np.log(max(yr_grid))       - np.log(min(yr_grid)))*random.random())
    random_yrs[i] = 1.5
    scales[i] = RMS_Model.scale.value
    
    #Say these are the values for the current spectrum
    RMS_Model.tauTheta = random_tauThetas[i]
    RMS_Model.R        = random_Rs[i]
    RMS_Model.yr       = random_yrs[i]
    
    #Print for us to see
    #print("random_tauTheta = %4.4e"%(random_tauThetas[i]),"random_R = %4.4e"%(random_Rs[i]),"random_yr = %4.4e"%(random_yrs[i]))
    
    #Generate a spectrum with the current parameter values, response function, and background
    spectrum_generator = DispersionSpectrumLike.from_function('fake',
                                                source_function=RMS_Model,
                                                response=response, 
                                                background_function = background_function)   
    #Save the spectrum to fake_spectra
    fake_spectra.append(spectrum_generator)
    #View the spectrum
    #spectrum_generator.view_count_spectrum()
    #plt.show()


# # Fitting the synthetic spectra

# ## With the Band function

# In[957]:


#Add a Band function model with specific values of the parameters
band = Band()
band.alpha = -1.
band.beta  = -3.
band.xp    = 300.
band.K     = 1e3
band_model = Model(PointSource("fake", 0.0, 0.0, spectral_shape = band))

#Fit the fake spectra using the Band function
saved_fits = []
for i in range(nr_fake_spectra):
    #Print the "real" parameter values for the spectra
    #print("random_tauTheta = %4.4e"%(random_tauThetas[i]),"random_R = %4.4e"%(random_Rs[i]),"random_yr = %4.4e"%(random_yrs[i]))
    
    #Fit the spectra using the band model and save it
    jl = JointLikelihood(band_model, DataList(fake_spectra[i]))
    saved_fits.append(jl)
    
    #Show the best fit parameters
    _ = jl.fit()
    #_=display_spectrum_model_counts(jl, show_legend=True, min_rate=1, figsize = (10,7))
    #plt.show()


# ## Access your best fit parameters

# In[958]:


for i in range(nr_fake_spectra):
    alpha_i = saved_fits[i].results._values[1]
    beta_i  = saved_fits[i].results._values[3]
    xp_i    = saved_fits[i].results._values[2]
    #print("Best fit values for fake_spectra %d are"%i)
    #print("alpha = %3.3e, beta = %3.3e, xp = %3.3e\n"%(alpha_i, beta_i, xp_i))


# In[959]:


#making arrays out of the obtained values of alpha, beta, and xp
alphaarray = []
betaarray = []
xparray = []
band_Ks = []
for i in range(nr_fake_spectra):
    alpha_i = saved_fits[i].results._values[1]
    beta_i  = saved_fits[i].results._values[3]
    xp_i    = saved_fits[i].results._values[2]
    band_Ki = saved_fits[i].results._values[0]
    arrays = np.array([(alpha_i), (beta_i), (xp_i), (band_Ki)])
    alphaarray.append([arrays[0]])
    betaarray.append([arrays[1]])
    xparray.append([arrays[2]])
    band_Ks.append([arrays[3]])
alphaarray = np.array(alphaarray)
betaarray = np.array(betaarray)
xparray = np.array(xparray)
band_Ks = np.array(band_Ks)
#print(alphaarray)
#print(betaarray)
#print(xparray)


# In[960]:


plt.plot(random_Rs, alphaarray)
plt.title('R vs alpha')
plt.grid()
plt.xlabel('R (ThetaR/ThetaU)')
plt.ylabel('alpha parameter of Band function')
plt.savefig('RvA.png')


# In[961]:


#quadratic regression
random_Rs = np.squeeze(random_Rs)
alphaarray = np.squeeze(alphaarray)

quadraticmodel = np.poly1d(np.polyfit(random_Rs, alphaarray, 2))

polyline = np.linspace(30, 300, 100)
plt.scatter(random_Rs, alphaarray)
plt.plot(polyline, quadraticmodel(polyline))
plt.grid()
plt.title('Relationship between R and Alpha Parameters')
plt.xlabel('R')
plt.ylabel('Alpha')
plt.savefig("regression.png")

print("Quadratic model:")
print(quadraticmodel)

#correlation
corr_matrix = np.corrcoef(alphaarray, quadraticmodel(polyline))
corr = corr_matrix[0,1]
R_sq = corr**2

print("R squared:")
print(R_sq)


# In[962]:


alphaarrayadded = np.array(np.asarray(alphaarray)+ 0.6)
#exponential fit
plt.scatter(random_Rs, alphaarrayadded)
#fit and plot the model
fit = np.polyfit(random_Rs, np.log(alphaarrayadded), 1)
#fit = np.array(np.exp(np.asarray(fit)))
y_temp = np.exp(fit[1])*np.exp(fit[0]*random_Rs)
print(fit)
plt.plot(random_Rs, y_temp)
plt.title('Relationship between R and Alpha Parameters')
plt.grid()
plt.xlabel('R')
plt.ylabel('Alpha')
plt.show()

#Calculating R Squared
corr_matrix = np.corrcoef(alphaarray, y_temp)
corr = corr_matrix[0,1]
R_sq = corr**2

print("R squared:")
print(R_sq)


# In[963]:


#double linear fit
x1 = np.arange(30, 100.2, 2.7)
x2 = np.arange(100.2, 300, 2.7)
alphaarray1 = alphaarray[0:26]
alphaarray2 = alphaarray[25:99]

linmodel1 = np.poly1d(np.polyfit(x1, alphaarray1, 2))
linmodel2 = np.poly1d(np.polyfit(x2, alphaarray2, 1))


#Calculating R Squared
corr_matrix1 = np.corrcoef(alphaarray1, linmodel1(x1))
corr1 = corr_matrix1[0,1]
R_sq1 = corr1**2

corr_matrix2 = np.corrcoef(alphaarray2, linmodel2(x2))
corr2 = corr_matrix2[0,1]
R_sq2 = corr2**2

fig, ax = plt.subplots(1,figsize=(8, 5))
plt.scatter(random_Rs, alphaarray)
plt.plot(x1, linmodel1(x1), 'r', label='R squared: {}'.format(R_sq1))
plt.plot(x2, linmodel2(x2), 'b', label='R squared: {}'.format(R_sq2))
plt.grid()
plt.title('Relationship between R and α Parameters')
plt.xlabel('R', size=15)
plt.ylabel('α', size=20)
ax.set_ylim([-0.6, -0.4])
plt.legend()
plt.show()

print(linmodel1)
print(linmodel2)


# In[964]:


#finding ranges for alpha, beta and xp
alpharange = np.ptp(alphaarray)
betarange = np.ptp(betaarray)
xprange = np.ptp(xparray)
print(alpharange, betarange, xprange)


# In[965]:


#Choose what spectra to plot
what_spectras = [0, 49, 99]
#Change these to edit the axis of the plot
y_lim_lower = 1e2
y_lim_upper = 3e7
x_lim_lower = 1
x_lim_upper = 1e5


fig, ax = plt.subplots(1,figsize=(8, 5))
i = -1
for what_spectra in what_spectras:
    i = i+1
    #Extract the closest input spectra
    tmp = abs(tauTheta_grid/random_tauThetas[what_spectra] - 1.)
    tauTheta_idx = np.where(tmp == min(tmp))[0][0]
    tmp = abs(R_grid/random_Rs[what_spectra] - 1.)
    R_idx = np.where(tmp == min(tmp))[0][0]
    tmp = abs(yr_grid/random_yrs[what_spectra] - 1.)
    yr_idx = np.where(tmp == min(tmp))[0][0]
    idx = tauTheta_idx + R_idx + yr_idx

    #Easier with shorter names
    alpha  = alphaarray[what_spectra]
    beta   = betaarray[what_spectra]
    band_K = band_Ks[what_spectra]
    scale  = scales[what_spectra]
    E_0    = xparray[what_spectra]/(alpha+2)
    #if (100 <= (alpha-beta)*E_0):
    A = band_K/np.exp(-100./E_0)
    #    print("I go here!")
    #else:
    #    A = band_K/(((alpha-beta)*E_0/100.)**(alpha - beta) * np.exp(beta-alpha))
    print("tauTheta = %4.4e"%(tauTheta_grid[tauTheta_idx]), "R = %4.4e"%(R_grid[R_idx]),"yr = %4.4e"%(yr_grid[yr_idx]),"scale = %4.4e"%(scale))
    Es = scale*energies
    xs = Es
    ys_rms = (Es**2.)*Nss[idx]*RMS_K/(scale**1.)
    ys_band = (Es**2.)*band_K*( (Es <= (alpha-beta)*E_0) *((Es/100.)**(alpha)) * np.exp(-Es/E_0) +                 (Es > (alpha-beta)*E_0) * (((alpha-beta)*E_0/100.)**(alpha - beta)) * ((Es/100.)**beta) * np.exp(beta-alpha) )
    ys_bkg = initial_K*100.**2.*(energies/100.)**(1./2.)

    #NOTE! I have no good reason for doing below. I must figure out how K_band and K_RMS works.
    #ys_rms = ys_rms/max(ys_rms)*max(ys_band)
    #if beta >= -2:
    #    ys_rms = ys_rms/max(ys_rms)*max((Es**2.)*A*(Es <= (alpha-beta)*E_0) *((Es/100.)**(alpha)) * np.exp(-Es/E_0))    
    #else:
    #    ys_rms = ys_rms/max(ys_rms)*max(ys_band)
    
    
    if i == 0:
        ax.plot(xs, ys_rms, 'r-', xs, ys_band, 'r:', xs, ys_bkg, 'k-.')
    elif i == 1:
        ax.plot(xs, ys_rms, 'b-', xs, ys_band, 'b:')
    elif i == 2:
        ax.plot(xs, ys_rms, 'g-', xs, ys_band, 'g:')
    elif i == 3:
        ax.plot(xs, ys_rms, 'c-', xs, ys_band, 'c:')
    elif i == 4:
        ax.plot(xs, ys_rms, 'm-', xs, ys_band, 'm:')
    else:
        print("NOTE: I can only make 5 pretty colors")
        ax.plot(xs, ys_rms, 'm-', xs, ys_band, 'm:')
    #ax.plot(xs, ys_band)
    

#ax.legend([r'$L_1$'])
GBM_low_lim = [8, 8]
GBM_upp_lim = [1e3, 1e3]
ax.plot(GBM_low_lim, [y_lim_lower, y_lim_upper], 'k--', GBM_upp_lim, [y_lim_lower, y_lim_upper], 'k--')


ax.set_ylim([y_lim_lower, y_lim_upper])
ax.set_xlim([x_lim_lower, x_lim_upper])
ax.set_xscale("log")
ax.set_yscale("log")
plt.xlabel('keV', fontsize=16)
plt.ylabel(r'$E^2 *$ counts s$^{-1}$ keV$^{-1}$', fontsize=16)
#plt.tick_params(
#    axis='x',          # changes apply to the x-axis
#    which='both',      # both major and minor ticks are affected
#    bottom=False,      # ticks along the bottom edge are off
#    top=False)         # ticks along the top edge are off
#ax.set_xticklabels([])
plt.tight_layout()
plt.show()


# In[ ]:




