%% Script for calculation of Spike-triggered Multiunit Activity

% Inputs:

%         WB :    single trial of a Wide-Band signal of  Multiunit channel
%         Spks :  single trial vector of Spike locations ('not time')
%         Par :   Parameters for specification of stMUA
%                Par.Fs  : Sampling frequency, Default: 32000 Hz
%                Par.WL  : Window length around spikes (msec.), Default : 75 msec
%                Par.Dur : half duration of the sliding window (msec.), Default: 2.5 msec
%                Par.sw  : step of the sliding window (msec.), Default : 0.5 msec
%                Par.Fband : frequency bands for band-pass filtering, Default : [800 3000]

% Output:
%         MUA : a matrix of MUAs where rows denote each spike triggered MUA
% How to call Example :
%         WB =randn(32000,1);
%         Spks =randi(32000,30);
%         MUA=st_MUA(WB,Spks);

%%%%%%%%%%%%%%%%%-------------------%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Kianoush Banaie Boroujeni 2021 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%-------------------%%%%%%%%%%%%%%%%%%%%

function  [MUA]=st_MUA(WB,Spks,Par)

if exist('Par.Fs')
    Fs=Par.Fs;
else
    Fs=32E3;
end
if exist('Par.WL')
    WL=Par.WL;
else
    WL=75;
end
if exist('Par.Fband')
    Fb=Par.Fband;
else
    Fb=[800 3000];
end
if exist('Par.sw')
    sw=Par.sw;
else
    sw=0.5;
end
if exist('Par.Dur')
    Dur=Par.Dur;
else
    Dur=2.5;
end

WL=75;
Dur=2.5;
Snum=round(sw*Fs/1000);
TW=round(Dur*Fs/1000);
TL=round(WL*Fs/1000);

yrF=bandpass(WB,Fb,Fs);
RyrF=abs(yrF);
L=length(WB);

numSpk=length(Spks);
LineS=1:Snum:2*TL+1;
Npoint=length(LineS);

SPK_Med=nan(numSpk,Npoint);        
        for vr=1:numSpk %finding trough of each spike to match with LFP
            if Spks(vr)>TL+TW+1 && Spks(vr)<L-TL-TW-1
            for j=1:Npoint
                R=RyrF((Spks(vr)-TL+LineS(j)-1)-TW:(Spks(vr)-TL+LineS(j)-1)+TW);
                R1=RyrF((Spks(vr)-TL+LineS(j)-2)-TW:(Spks(vr)-TL+LineS(j)-2)+TW);
                R2=RyrF((Spks(vr)-TL+LineS(j))-TW:(Spks(vr)-TL+LineS(j))+TW);
                Peaks=R(R>R1 & R>R2);

                if ~isempty(Peaks)
                    SPK_Med(vr,j)=median(Peaks);
                end
            end
            end
        end
        MUA=SPK_Med;

end
