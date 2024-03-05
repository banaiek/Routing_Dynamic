%% This script collects input parameters for LFP burst detection. The input parameters include:

% Fsin - sampling frequency of the input signal.
% Par.fds - frequency at which the input signal should be downsampled.
% Par.lpf - cutoff frequency for the low-pass filter.
% Par.Fband - defines narrow frequency bands for burst detection. For example, [1.5 4; 4 8; 6 12; 8 16; 15 25 ;65 115].
% Par.NC - sets the moothing cycle factor (depending on the frequency band-width between 1 to 2).
% Par.mCs - establishes the minimum number of prominent cycles following detection.
% Par.Csep - denotes the minimum separating cycles.
% Par.MaxCycl - sets the maximum number of cycles on each side of the burst.

% The output is a structured array, which includes:

% out.Phase_bp - phase of the band-pass filtered signal.
% out.Band_Passed_F - band-pass filtered signal.
% out.Logical_P - logical value indicating if a burst was detected.
% out.Wrap - envelope of the signal.
% out.ThPS - signal that passed the threshold.
% out.Os_Cp - center point of the detected burst.
% out.Os_Ct - center time of the burst.
% out.Time - time series of the downsampled signal.
% out.OS_StartEndtimes - start and end times of each detected burst.
% out.OS_StartEndPoints - start and end points of each burst.
% out.NumC - total number of detected bursts.
% out.Os_CM - center of the burst.

%***********************************
% Kianoush Banaie Boroujeni (Kia) - 2023
%***********************************
%example
%parameters

Fsin = 1000;
Par.fds=1000; 
Par.lpf=1000; 
Par.Fband=[1.5 4; 4 8; 6 12; 8 16; 15 25 ;65 115]; 
Par.NC=[1, 1, 1, 1, 1.5, 1]; 
Par.mCs=[1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
Par.Csep=[1, 1, 2, 2, 2, 5]; 
Par.MaxCycl=[5, 5, 10, 10, 10, 10];

out=Osc_burst_detection_V01(X,time,Fsin,Par);

