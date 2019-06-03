import numpy as np
from peakdetect import peakdetect
from scipy.signal import butter, lfilter, lfilter_zi
from biosppy import storage
from biosppy.signals import ecg
from scipy.signal import correlate

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

def get_indices_of_R_peaks(ecg_signal):
	rp_ind = ecg.engzee_segmenter(signal=ecg_signal, sampling_rate=300.0, threshold=0.48)
	rp_peaks = []
	for i in list(rp_ind[0]):
		rp_peaks.append(ecg_signal[i])
	return list(rp_ind[0]), rp_peaks

def get_phase_shift(sig1, sig2):
    cross_corr = correlate(sig2, sig1, 'full', 'direct')
    # plt.plot(range(len(cross_corr)), cross_corr)
    # plt.show()
    shift = np.argmax(cross_corr) - len(sig1)
    DESCRETISATION_PERIOD = 1/300
    return shift*DESCRETISATION_PERIOD

def get_shift_by_peaks(ecg_p, ppg_p, nyqyist_freq):
	if ecg_p[0] > ppg_p[0]:
		ppg_p = ppg_p[1:len(ecg_p)]
	else:
		ppg_p = ppg_p[0:len(ecg_p)]
	ecg_p = ecg_p[0:len(ppg_p)]
	ecg_p = np.asarray(ecg_p)
	ppg_p = np.asarray(ppg_p)
	diff = np.mean(ecg_p - ppg_p) / nyqyist_freq

	print("true shift by peaks ", diff)


	return nyqyist_freq
