%% Burst Detection Function:
%**************************
%% Kianoush Banaie Boroujeni (Kia)- 2023
%**************************

function [out]=Osc_burst_detection_V01(data,time,Fsin,Par)
%parameters
Fs=Par.fds; %down sample freq
lpf=Par.lpf; % low-pass filter freq
Fband=Par.Fband; %frequency bands



Nb=size(Fband,1); %Number of bands

Timein=time;
yin=data;



[yds,Time]=Down_smpl_01(yin',Timein,Fs,Fsin);
[y,~]=BP_filter_S_01(lpf,yds,Fs);
L1=length(y);


for fr=1:Nb
    clearvars -except fr Phase_bp ThPS Wrap Logical_P Band_Passed_F lpf L1 y yds Time Nb Fband Par Fs out
    NC=Par.NC(fr); %min number of prominent cycles
    mCs=Par.mCs(fr);
    Csep=Par.Csep(fr); %seperating min cycles
    MaxCycl=Par.MaxCycl(fr); %Max cycle on each side
    bp=z_Phase_bp_01(y,Fband(fr,:),Fs);
    hy=hilbert(bp);
    Ph_y = atan2(imag(hy),real(hy));


    Cf=mean(Fband(fr,:)); %Center freq

    %% Rectification and peak detection
    Abp=abs(bp);
    R1=Abp(1:end-2);
    R2=Abp(3:end);
    R=Abp(2:end-1);
    peaks=find(R>R1 & R>R2)+1;

    peak_Values(peaks)=Abp(peaks);


   %% Phase discontinuity
   d_phase = diff(Ph_y);
   dX_unwrapped = unwrap(d_phase);
   threshold = pi;
   Ph_discontinuity = [0,(abs(dX_unwrapped) > threshold)];

    %% Setting the window size
    WL = round(NC*Fs/Cf); % Window length for amp. conv.
    mWL = round(mCs*Fs/Cf); %min window length
    hWL = round(Fs/Cf); % mean window length
    Ana_bp=abs(hy);

    %% RMS thresholding and min prom. cycle detection
    Wrp=Ana_bp;
    abs_bp=abs(Ana_bp); % abs Ana. signal amp.
    
    E_Ana_bp=Ana_bp.^2; % Energy

    Med_abp=sqrt(nanmean(E_Ana_bp));
    labp=abs_bp(abs_bp<Med_abp);
    habp=abs_bp(abs_bp<Med_abp);
    rms1=Med_abp+3.3*nanstd(labp);
    rms2=Med_abp+3.3*nanstd(habp);
    
    A_Pks=findpeaks(E_Ana_bp);
    Pk_tr=nanmean(A_Pks)+nanstd(A_Pks); % Setting a peak threshold of one SD above the mean
    [~,locs] = findpeaks(-E_Ana_bp);
    D_locs=abs(diff(locs));
    W_Thr = min(nanmedian(D_locs),mWL); % adjusted window threshold 
    

    rms = min(rms1,rms2);
    rmst = 2*(nanmedian(peak_Values.^2)); % peak Energy threshold of 2 RMS

    A=nan(size(bp));
    A2=nan(size(bp));
    A2(E_Ana_bp-Pk_tr>0)=1;
    A(Wrp-1*rms>0 & E_Ana_bp-rmst>0 & A2>0 & ~Ph_discontinuity )=1; %logical value of threshold passed signal
    WC=ones(WL,1)/(WL);

    DOs1=(conv(A.*E_Ana_bp,WC,'same'));
    D1=DOs1(1:end-2);
    D2=DOs1(2:end-1);
    D3=DOs1(3:end);

    OS_D=round(Csep*Fs/Cf); %min distance of oscillations
    OsC=find(D2-D1>0 & D2-D3>0)+1;
    OsC(OsC<=hWL | OsC>=L1-hWL)=[];

    %removing invalid oscillations (too close and spourios local peaks)
    invp=find(Wrp(OsC)<Wrp(OsC-round(hWL)) | Wrp(OsC)<Wrp(OsC+round(hWL)) | isnan(A(OsC))  );
    OsC(invp)=[];

    %Removing too close oscillations
    %length(A_Pks)./length(P_V)
    Clp=find(diff(OsC)<OS_D);

    while ~isempty(Clp)
        [~,mWrp]=min([Wrp(OsC(Clp(1))),Wrp(OsC(Clp(1)+1))]);
        OsC(Clp(1)+mWrp-1)=[];
        Clp=find(diff(OsC)<OS_D);
    end

    OsC=reshape(OsC,1,numel(OsC)); % oscillation centers

    % Portions with valid oscillation

    logic_p=zeros(size(bp));





    NOsc=length(OsC);
    OStp=[];


    cutoff = max(rmst,Med_abp.^2);
    Wrp2 = E_Ana_bp-cutoff;
    
    for i=1:NOsc
        t1=max(1,OsC(i)-MaxCycl*hWL); %starting limit of Burst
        t2=min(L1,OsC(i)+MaxCycl*hWL); %end limit of Burst

        OSt1=[];
        OSt2=[];

        [~,OSt2]=find(Wrp2(OsC(i):t2)<0,1,'first');
        [~,OSt1]=find(Wrp2(t1:OsC(i)-1)<0,1,'last');
        if isempty(OSt1)

            OSt1=1;%
        end
        if isempty(OSt2)

            OSt2=min(L1-OsC(i),MaxCycl*hWL)-1;

        end

        OStp(i,1)=t1+OSt1(end);
        OStp(i,2)=OsC(i)+OSt2(1)-1;

        if i<NOsc
            if OStp(i,2)>OsC(i+1)
                [~,ml]=min(Wrp(OsC(i):OsC(i+1)));
                OStp(i,2)=OsC(i)+ml-1;
            end
        end

        if i>1
            if OStp(i,1)<OsC(i-1)
                [~,ml]=min(Wrp(OsC(i-1):OsC(i)));
                OStp(i,1)=OsC(i-1)+ml+1;
            end
        end

        if OStp(i,2)-OStp(i,1)>=W_Thr
            logic_p(OStp(i,1):OStp(i,2))=1;
        end
    end

    L_osc=diff(OStp,[],2);
    inv_osc=find(L_osc<W_Thr);
    OStp(inv_osc,:)=[];
    OsC(inv_osc)=[];


    if ~isempty(OStp)
        OS_times{fr}(:,1)=Time(OStp(:,1));
        OS_times{fr}(:,2)=Time(OStp(:,2));
    else
        OS_times{fr}=[];
    end

    out.Phase_bp(fr,:)=Ph_y;
    out.Band_Passed_F(fr,:)=bp;
    out.Logical_P(fr,:)=logic_p;
    out.Wrap(fr,:)=Wrp;
    out.ThPS(fr,:)=A;
    out.Os_Cp{fr}=OsC;
    out.Os_Ct{fr}=Time(OsC);
    out.Time=Time;
    out.OS_StartEndtimes{fr}=OS_times{fr};
    out.OS_StartEndPoints{fr}=OStp;
    out.NumC(fr)=NOsc;
    out.Os_CM{fr}=Wrp(OsC)/nanmean(Wrp);
end
