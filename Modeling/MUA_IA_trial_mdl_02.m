function  [MUA,Ratio]=MUA_IA_trial_mdl_02(WB,Spks,Fs,sw,Fb)

WL=100;
Dur=2.5;
Snum=round(sw*Fs/1000);
TW=Dur*Fs/1000;
TL=WL*Fs/1000;

yrF=bandpass(WB,Fb,Fs);
RyrF=abs(yrF);
L=length(WB);
Ratio=zeros(1,L);
% rms=sqrt(mean((yrF).^2));
numSpk=length(Spks);
LineS=1:Snum:2*TL+1;
Npoint=length(LineS);
time=linspace(-WL,WL,Npoint);
pre=find(time<0 & time>-20);
post=find(time<20 & time>0);
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
            Ratio(1,Spks(vr)) = nanmean(SPK_Med(vr,post)) - nanmean(SPK_Med(vr,pre));
            end
        end
        MUA=SPK_Med;
        
        

end
