from scipy import signal
import numpy as np
import math

def find_peaks(counts,hist_binwidth):
  # call with histogram counts as array
  # returns index of detected peaks
  no_bins = math.ceil(0.02/hist_binwidth)
  peaks_idx = signal.find_peaks_cwt(counts, np.arange(1,no_bins,1)) # searched between 1 and 10 binwidths
  return (peaks_idx+1)
