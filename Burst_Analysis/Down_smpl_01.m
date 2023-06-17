% This function down-samples the data y with Freq Fs to FS/D

function [yds,Time2]=Down_smpl_01(y,Time,dFs,Fs)
D=round(Fs/dFs);
tf=range(Time);
T1=Time(1);
Tl=Time(end);
L = fix(size(y,1)/D);
yr = y(1:D*L,:); % recorded data
 ym = reshape(yr,D,L,size(y,2));
%   ym = vec2mat(yr,D); %reshaped
yds = squeeze(nanmean(ym,1)); % peak detection
Fs2 = Fs/D; % new sampling frequency
L2=size(yds,2);
Time2=T1:tf/(L2-1):Tl;
end
