import numpy as np
from peakdetect import peakdetect
from scipy.signal import butter, lfilter, lfilter_zi
from biosppy import storage
from biosppy.signals import ecg
from scipy.signal import correlate
from scipy.signal import hilbert, chirp
from matplotlib import pylab as plt

def invertation(ppg): # type  of ppg is np.array
	new_ppg = []
	moda = sum(ppg)/len(ppg)
	for i in -ppg:
		new_ppg.append(i + 2*moda)
	return new_ppg

def ppg_sistolic_peaks(ppg_signal):
    peaks = peakdetect(ppg_signal,lookahead=5)
    peaks_max = []
    peaks_max_ind = []
    for i in peaks[0]:
        peaks_max.append(i[1])
    for j in peaks[0]:
        peaks_max_ind.append(j[0])
    return peaks_max_ind, peaks_max

def ppg_min_peaks(ppg_signal):
    peaks = peakdetect(ppg_signal,lookahead=5)
    peaks_min = []
    for i in peaks[1]:
        peaks_min.append(i[1])
    peaks_min_ind = []
    for j in peaks[1]:
        peaks_min_ind.append(j[0])
    return peaks_min_ind, peaks_min


def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    zi = lfilter_zi(b, a)
    y, z = lfilter(b, a, data, zi=zi*data[0])
    return y

def read_file(file_name):
    ch1 = []
    ch2 = []
    with open(file_name, 'r') as f:
        for row in f:
            row_str = row[:-1]
            row_splt = row_str.split('\t')
            ch1.append(float(row_splt[0]))
            ch2.append(float(row_splt[1]))
    offset = int(abs(len(ch1)-3000)/2)
    offset = 5
    return ch1[offset:-offset-1], ch2[offset:-offset-1]


def top_peaks(ppg_peaks, ppg_peaks_ind):
	tops = []
	tops_ind = []
	i = 0
	try:
	    for ind, val in enumerate(ppg_peaks):
		    if val > ppg_peaks[ind + 1] and val > ppg_peaks[ind -1]:
			    tops.append(val)
			    tops_ind.append(ppg_peaks_ind[i])
			    i += 1
		    else:
			    i += 1
	except IndexError:
		pass
	return tops, tops_ind

def top_min_peaks(ppg_peaks, ppg_peaks_ind):
	tops = []
	tops_ind = []
	i = 0
	try:
	    for ind, val in enumerate(ppg_peaks):
		    if val < ppg_peaks[ind + 1] and val < ppg_peaks[ind -1]:
			    tops.append(val)
			    tops_ind.append(ppg_peaks_ind[i])
			    i += 1
		    else:
			    i += 1
	except IndexError:
		pass
	return tops, tops_ind

def get_indices_of_R_peaks(ecg_signal):
	rp_ind = ecg.engzee_segmenter(signal=ecg_signal, sampling_rate=300.0, threshold=0.48)
	rp_peaks = []
	for i in list(rp_ind[0]):
		rp_peaks.append(ecg_signal[i])
	return list(rp_ind[0]), rp_peaks


def get_shift_by_peaks(ecg_p, ppg_p, nyqyist_freq):
    if ecg_p[0] > ppg_p[0]:
	    ppg_p = ppg_p[1:len(ecg_p)]
    else:
	    ppg_p = ppg_p[0:len(ecg_p)]
    ecg_p = ecg_p[0:len(ppg_p)]
    ecg_p = np.asarray(ecg_p)
    ppg_p = np.asarray(ppg_p)
    diff = np.mean(ecg_p - ppg_p) / nyqyist_freq
    # print("true shift by peaks ", round(diff, 3))

    shift = round(diff, 3)


    return shift


def BPM_calculation(ppg_signal):
	ppg_x, ppg_y = ppg_sistolic_peaks(ppg_signal)
	ppg_peaks = top_peaks(ppg_y, ppg_x)[0]
	k = 60/(len(ppg_signal)/300)
	beats = k * len(ppg_peaks)
	return round(beats)

def get_attitude(ppg_signal):
    ppg_signal = butter_bandpass_filter(ppg_signal, 0.1, 5, 300, 3)

    ppg_x, ppg_y = ppg_sistolic_peaks(ppg_signal)
    p_high, p_high_ind = top_peaks(ppg_y, ppg_x)

    ppg_x_min, ppg_y_min = ppg_min_peaks(ppg_signal)
    p_low, p_low_ind = top_min_peaks(ppg_y_min, ppg_x_min)

    y_max = np.array([])
    for i in range(len(p_high)):
        if i == len(p_high) - 1:
            break
        apendix = np.linspace(p_high[i], p_high[i+1], p_high_ind[i+1] - p_high_ind[i] + 1)
        if i != len(p_high) - 2:
            apendix = np.delete(apendix, i + 1)
        y_max = np.concatenate([y_max, apendix])

    y_min = np.array([])
    for i in range(len(p_low)):
        if i == len(p_low) - 1:
            break
        apendix = np.linspace(p_low[i], p_low[i+1], p_low_ind[i+1] - p_low_ind[i] + 1)
        if i != len(p_high) - 2:
            apendix = np.delete(apendix, i + 1)
        y_min = np.concatenate([y_min, apendix])

    y_max_ind = np.linspace(p_high_ind[0], p_high_ind[len(p_high)-1], p_high_ind[len(p_high)-1] - p_high_ind[0] + 1  )
    y_min_ind = np.linspace(p_low_ind[0], p_low_ind[len(p_low)-1], p_low_ind[len(p_low)-1] - p_low_ind[0] + 1  )

    L = min(len(y_max), len(y_min))
    y_max = y_max[0:L]
    y_min = y_min[0:L]

    L_ind = min(len(y_max_ind), len(y_min_ind))
    y_max_ind = y_max_ind[0:L_ind]
    y_min_ind = y_min_ind[0:L_ind]

    attitude = y_max/(y_max - y_min)
    attitude = np.mean(attitude)
    attitude = round(attitude, 3)
        

    x = np.arange(0, len(ppg_signal))
    y = np.array(ppg_signal[:len(ppg_signal)])
    plt.plot(x, y)
    plt.plot(y_max_ind, y_max)
    plt.plot(y_min_ind, y_min)
    plt.show()

    return attitude





# def invertation(ppg): # type  of ppg is np.array
# 	new_ppg = []
# 	moda = sum(ppg)/len(ppg)
# 	for i in -ppg:
# 		new_ppg.append(i + 2*moda)
# 	return new_ppg

# def write_file(file_name, shift, BPM, attitude):
# 	sig1, sig2 = read_file(file_name)
#     with open(file_name, 'w') as f:
#         for i in range(10):
#             row = str(shift) + '\t' + str(BPM) + '\t' + str(attitude) + '\r'
#             f.write(row)
#     return




