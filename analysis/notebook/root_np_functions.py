import ROOT
import numpy as np
from ROOT import TGraphErrors
from ROOT import TVector
import math

def TVT_to_numpy(file,TVectorT_name):
    TVectorT = file.Get(TVectorT_name)
    N = TVectorT.GetNoElements()
    np_bins = np.zeros(N)
    for i in range(N):
        np_bins[i] = TVectorT[i]
    return np_bins,N-1 #return number of bin centers

def TH1_to_numpy_wErrors(file, TH1_name, normalize=False, zero_suppress=False, scale_factor=1.):    
    th1 = file.Get(TH1_name)
    #th1.Scale(scale_factor)
    #scale_factor = 1.
    Nx = th1.GetNbinsX()
    th1_array = np.zeros(Nx)
    th1_errors = np.zeros(Nx)
    
    for i in range(Nx):
        val = th1.GetBinContent(i+1)
        error = th1.GetBinError(i+1) * math.sqrt(scale_factor)
        if zero_suppress and val == 0:
            val = np.nan
        if not(zero_suppress) and val!=0 and error/val >= 0.5:
            val=np.nan #skip converged fits with poor statistics
        th1_array[i] = val
        th1_errors[i] = error
    if normalize:
        N = np.sum(th1_array)*scale_factor
        th1_array = th1_array/N # Equivalent to (1/sqrt(scale_factor))
        th1_errors = th1_errors/N
    return th1_array, th1_errors

def get_th1_binning_np(file, TH1_name, dirname=""):
    if dirname:
        TH1_name = dirname+"/"+TH1_name
    th1 = file.Get(TH1_name)
    N = th1.GetNbinsX()
    min = th1.GetXaxis().GetXmin()
    max = th1.GetXaxis().GetXmax()
    bins = np.linspace(min,max,N+1)
    centers = (bins[1:] + bins[:-1]) /2
    widths = [(j-i)/2 for i, j in zip(bins[:-1], bins[1:])]
    if (len(centers)!=N):
        print("get_th1_binning_np: something went wrong")
    return bins,centers,widths


def get_binning_from_edges(bin_edges):
    centers = (bin_edges[1:] + bin_edges[:-1]) /2
    widths = [(j-i)/2 for i, j in zip(bin_edges[:-1], bin_edges[1:])]
    return centers,widths

def TH1_to_numpy_wErrors_old(file, TH1_name, cut_low_stat=True, Normalize = False, dirname=""):    
    if dirname:
        TH1_name = dirname+"/"+TH1_name
    th1 = file.Get(TH1_name)
    N = th1.GetNbinsX()
    th1_array = np.zeros(N)
    th1_errors = np.zeros(N)
    for i in range(N):
        val = th1.GetBinContent(i+1)
        error = th1.GetBinError(i+1)
        if cut_low_stat:
            if val == 0: #ROOT for some reason puts unconverged fits from poor stats to 0. Should be NaN
                val = np.nan
            if error/val >= 0.3: 
                val=np.nan #skip converged fits with poor statistics
        th1_array[i] = val
        th1_errors[i] = error
    return th1_array, th1_errors