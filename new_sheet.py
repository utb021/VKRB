from matan import *
from matplotlib import pylab as plt
import os

sig1, sig2 = read_file('Temirlan1.txt')
new_ppg = invertation(np.array(sig1[:len(sig1)]))
# new_ppg = butter_bandpass_filter(new_ppg, 0.1, 5, 300, 3)
ppg_x, ppg_y = ppg_sistolic_peaks(new_ppg)
ppg_sist_peaks, ppg_sist_ind = top_peaks(ppg_y, ppg_x)
ppg_x_min, ppg_y_min = ppg_min_peaks(new_ppg)
ppg_diast_peaks, ppg_diast_peaks_ind = top_min_peaks(ppg_y_min, ppg_x_min)

ecg_peaks_ind, ecg_peaks = get_indices_of_R_peaks(sig2)

shift = get_shift_by_peaks(ecg_peaks, ppg_sist_peaks, 300)
beats = BPM_calculation(new_ppg)
attitude = get_attitude(new_ppg)
print(shift, beats, attitude)

# x = np.arange(0, len(new_ppg))
# y = np.array(new_ppg[:len(new_ppg)])
# a = np.arange(0, len(sig2))
# b = np.array(sig2[:len(sig2)])
# plt.plot(a, b)
# plt.plot(ecg_peaks_ind, ecg_peaks, 'x')
# plt.plot(x, y)
# plt.plot(ppg_sist_ind, ppg_sist_peaks, 'x')
# plt.plot(ppg_diast_peaks_ind, ppg_diast_peaks, 'x')
# plt.show()


