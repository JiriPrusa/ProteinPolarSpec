from dipole_lib import importDipoleData
from MEM_burg import burg_AR
from MEM_burg import cic
from spectrum import arma2psd
import numpy as np
from numpy.fft import fft
from scipy.constants import k, epsilon_0, c
import matplotlib.pyplot as plt

###################MEM SPECTRA CALCULATION #############################
########################################################################
def MEM():
    timestep = 60 * 1e-15  # in sec
    nfft = 4096  # number of points to calculate FT for
    p = 75  # order of Burg
    alpha = 0.5  # damping factor for coefficients
    s_num = 0
    t_num = 50
    ###############
    plot_style = 1  # 0 == Frequency, 1 == wavelength
    # Create a for loop to read in each csv file and append the results to a list
    data_list = []
    tic = time.time()
    # import multiple CSV files
    for i in range(s_num, t_num):
        print(i)
        dipFileName = "/home/jirka/WORK/UFE/Jessica/lsozyme_2021/dip_data/SUMDIP_%i.csv" % i
        # import dipole file
        dip_data, V = importDipoleData(dipFileName)
        dip = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
    #    dip_array = np.transpose(dip[:,:])
        dip_array = np.transpose(dip[:,1])
        data_list.append(dip_array)

    dip_data = np.vstack(data_list)
    print("Available number of time series: %i with %i frames" % np.shape(dip_data))
    # uncomment next line to shrink the input data (Don't forget to change timestep if you do so)
    # dip_data = np.array(dip_data)[0, :]
    print("You calculate spectra from %i time series of %i length." % np.shape(dip_data))
    print("With timestep of %.4e sec between frames." % timestep)
    # Create an array of derivatives
    array_of_derivatives = np.diff(dip_data)
    # array_of_derivatives = np.atleast_2d(array_of_derivatives)
    # feed array into burg_AR class
    coeffs, RC, sigma = burg_AR(p, array_of_derivatives)
    rc = np.zeros(p + 1)
    for k in range(0, p):
        rc[k + 1] = RC[k]

    res = sigma * np.cumprod(1 - rc ** 2)
    rc[0] = 1
    cic_ar, pe_est = cic(res, 166667)
    # AR, rho, ref = spectrum.arburg(array_of_derivatives, 100)
    toc = time.time()
    print("time: " + str(toc - tic))
    print("order selection due to cic: " + str(np.argmin(cic_ar)))
    # p = spectrum.pburg(array_of_derivatives, 8, scale_by_freq=False);
    # p.plot()
    # print(coeffs)
    ###############################################################
    # Damp the coefficients to make it smooth for Fourier Transform
    T = len(coeffs)
    coeffs_dmp = np.zeros(T)
    for i in range(T):
        damp = np.exp(-i / (alpha * T))
        coeffs_dmp[i] = coeffs[i] * damp
    # calculate power spectral density as FT of Burg coeffitients utilizing spectra package
    psd = arma2psd(A=coeffs_dmp, rho=sigma, T=1 / timestep, NFFT=nfft, sides='default', norm=False)
    freq = np.fft.fftfreq(len(psd), d=timestep)
    # strip the left (negative) side
    psd = psd[1:len(psd) // 2]
    freq = freq[1:len(freq) // 2]
    ##################################################################
    if plot_style == 0:
        x_ticks = freq
    else:
        x_ticks = freq / c / 100

    fig, axs = plt.subplots(2)
    axs[0].plot(x_ticks, psd)
    axs[1].plot(coeffs)
    plt.show()
    np.savetxt("/home/jirka/WORK/UFE/alignByQueen/TEST/coefs_openmm.csv", coeffs)
    np.savetxt("/home/jirka/WORK/UFE/alignByQueen/TEST/spec_openmm.csv", np.vstack((freq, psd)).T, header='Freq  Intensity', delimiter="\t")

########################################################################################################################
########################################################################################################################

MEM()
