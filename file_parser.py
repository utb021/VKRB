from matan import *
import os
folders = os.listdir()
k = 0
for i in range(7):
	with open (folders[i].upper, 'w') as f:

        for folder in folders:
        	if k == 10:
        		continue
	        for file in os.listdir(folder):
	            sig1, sig2 = read_file(file)
	            new_ppg = invertation(np.array(sig1[:len(sig1)]))
	            ppg_x, ppg_y = ppg_sistolic_peaks(new_ppg)
	            ppg_sist_peaks, ppg_sist_ind = top_peaks(ppg_y, ppg_x)
	            ecg_peaks_ind, ecg_peaks = get_indices_of_R_peaks(sig2)
	            shift = get_shift_by_peaks(ecg_peaks, ppg_sist_peaks, 300)
                beats = BPM_calculation(new_ppg)
                attitude = get_attitude(new_ppg)
                row = str(shift) + '\t' + str(BPM) + '\t' + str(attitude) + '\r'
                f.write(row)
                k += 1


         


