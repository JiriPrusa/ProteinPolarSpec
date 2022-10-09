# Reads in dipole data from non-polarizable FF and performs ACF
import pandas as pd
import io
import numpy as np


def importDipoleData(dipFileName):
    with open(dipFileName, 'r') as f:
        file_content = f.read()  # Read whole file in the file_content string
    f = io.StringIO(file_content)
    vol_comment = f.readline()
    V = float(vol_comment.strip("#").split(" ")[2]) * 1e-30  # convert A^3 to m^3
    df = pd.read_csv(dipFileName, skiprows=0, header=1, delimiter="\s+", escapechar='#', index_col=False)
    return df, V

def exportDipoleData(dipData, Volume, oFileName):
    f = open(oFileName, "w")
    V = Volume / 1e-30 # convert m^3 to A^3
    f.write("# BoxVolume: %.4f        A^3\n# " % V)
    a = pd.DataFrame(dipData, columns=['time', 'dip_x', 'dip_y', 'dip_z', '|dip|'])
    a.to_csv(f, sep='\t', encoding='utf-8', index=False)
    f.close()

def importMtot(xvgFileName):
    data = np.loadtxt(xvgFileName, comments=["@", "#"])
    frames = data[:, 0]
    dip = data[:, 1:-1]
    return frames, dip

def monoExp(t, tau):
    return np.exp(-t/tau)

def biExp(t, A1, tau1, tau2):
    return A1 * np.exp(-t / tau1) + (1-A1) * np.exp(-t / tau2)

def numFLap(funct, time, freq):
    # Returns the image of $funct in $freq area.Where $funct is function of $time( in time domain).
    # ang_freq = freq. * 2 * pi;
    ang_freq = 2*np.pi*freq
    dir_phi_n = funct
    FreqImage = np.zeros(len(ang_freq), dtype = 'complex_')
    for omega in range(len(ang_freq)):
        integral = 0 + 0j
        for clock in range(len(time) - 2):
            # preforms the Fourier - Laplace transform
            integral += np.exp(1j * ang_freq[omega] * time[clock]) * dir_phi_n[clock]
        FreqImage[omega] = integral
    return FreqImage
