%% This function is butterworth IIR filter, which we use to make a sequence of filters and with sequential
%% subtraction we obtain band-pass filtered signals


%%----------
%written by Kia. Banaie Boroujeni 2018

function[yF, yrF]= BP_filter_S_01(Fc,data,Fs) 
   NB = 0; NA = 2; Wn = Fc/(Fs); 
   [B,A] = maxflat(NB,NA,Wn);
   yF = filtfilt(B,A,data); % filtered data with sampling frequecy at FS
    yrF=data-yF; %subtraction for the next level of filtering
end