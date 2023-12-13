# COMCAL-cirrus-cloud-retrieval-scheme
COMCAL cirrus cloud retrieval scheme
The ComCAL is a traditional Raman lidar. It has been introduced here:
https://epic.awi.de/id/eprint/15547/1/Imm2006l.pdf
Since 2018 it is located at Palau (7.3°N, 134.5°E). Its most important technical parameters are:
Laser: 20Hz with 120 mJ / 180 mJ and 65 mJ for the three colors of 355 nm, 532 nm and 1064 nm
Telescope: 40cm with 0.83mrad fov
Detection: 355 and 532 in “p” and “s” polarization. Further a 387nm channel is installed (not used here)
The data is in netcdf format with 10min / 60m resolution. The variables and their units are described there. Data is stored up to that altitude in which the SNR of the lidar signal drops below 3 
Due to legal constraints, daytime measurements between 8 a.m. and 5 p.m. local time are not permitted.  
The evaluation has been done according to Klett 1985: http://www.lalinet.org/viiiwlmla/Courses/Lecture4/Klett_Applied%20optics_1985.pdf
The 1064nm channel is currently weaker than expected.

The "comcal_daily_analysis.m" is the main program for the analysis of the lidar data. It's used for reading vertical and time resolution as the input for the data interactive process. 

The "dattoaer_iterative2.m" is the basic code for the iteractive process (Klett, 1985) to retrieve the optical properties from COMCAL. The LR is always iteractively inspected so each day has its own 'dattoaer_iterative' code, and it can be provided by Xiaoyu Sun (xiaoyu_sun@iup.physik.uni-bremen.de).

The "detection_of_cirrus_layers.py" is the code for the detection of the cirrus cloud base and top height by the Wavelet Covariance Transform (WCT) method Gamage and Hagelberg (1993). 

More information can be provided by: 
xiaoyu_sun@iup.physik.uni-bremen.de
christoph.ritter@awi.de
