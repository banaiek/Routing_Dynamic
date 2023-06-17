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


Band_Passed_F=nan(Nb,L1); %band-passed signal
Logical_P=nan(Nb,L1); %logical burst portion
Wrap=nan(Nb,L1); %warp of bp signal
ThPS=nan(Nb,L1); %primary threshold passed
Phase_bp=nan(Nb,L1); %Phase of band passed signal

for fr=1:Nb
    clearvars -except fr Phase_bp ThPS Wrap Logical_P Band_Passed_F lpf L1 y yds Time Nb Fband Par Fs out
    NC=Par.NC(fr); %min number of prominent cycles
    mCs=Par.mCs(fr);
    Csep=Par.Csep(fr); %seperating min cycles
    MaxCycl=Par.MaxCycl(fr); %Max cycle on each side
    bp=bandpass(y,Fband(fr,:),Fs);
    hy=hilbert(bp);
    Ph_y = atan2(imag(hy),real(hy));


    Cf=mean(Fband(fr,:)); %Center freq

    %% Rectification and peak detection
    Abp=abs(bp);
    R1=Abp(1:end-2);
    R2=Abp(3:end);
    R=Abp(2:end-1);
    peaks=find(R>R1 & R>R2)+1;

    P_V(peaks)=Abp(peaks);

    %% Median Filtering
    WL=round(NC*Fs/Cf);
    mWL=round(mCs*Fs/Cf);
    hWL=round(Fs/Cf);
    M_Abp=abs(hy);

    %% RMS thresholding and min prom. cycle detection
    Wrp=M_Abp;
    abp=abs(M_Abp);
    M_Abp=M_Abp.^2;
    Med_abp=nanmean(M_Abp);
    labp=abp(abp<Med_abp);
    diffSD =((nanstd(labp)/sqrt(1-2/pi))-(nanmean(labp)/(sqrt(2)/pi)))/sqrt(length(P_V)/2);
    habp=abp(abp<Med_abp);
    sd=min(nanstd(labp),nanstd(habp));
    rms1=Med_abp+3.3*nanstd(labp);%+2*nanstd(P_V)/sqrt(length(P_V));
    rms2=Med_abp+3.3*nanstd(habp);
    A_Pks=findpeaks(M_Abp);
    Pk_tr=nanmean(A_Pks)+2*nanstd(A_Pks);
    [~,locs] = findpeaks(-M_Abp);
    D_locs=abs(diff(locs));
    W_Thr=min(nanmedian(D_locs),mWL);
    mWL=W_Thr;

    rms=min(rms1,rms2);
    rmst=2*(nanmean(P_V.^2));
    A=nan(size(bp));
    A2=nan(size(bp));
    A2(M_Abp-Pk_tr>0)=1;
    A(M_Abp-1*rms>0 & M_Abp-rmst>0 )=1; %logical value of threshold passed signal
    WC=ones(WL,1)/(WL);
    WL2=round(WL/2);
    WC2=ones(WL2,1)/(WL2);

    DOs1=(conv(A.*M_Abp,WC,'same'));
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

    A_O=nan(size(bp));
    logic_p=zeros(size(bp));

    for i=1:length(OsC)

        A_O(max(1,OsC(i)-MaxCycl*hWL):min(OsC(i)+MaxCycl*hWL,L1))=1;
    end

    A_OS=A_O.*A; % logical burst portions



    dPh_y=diff(Ph_y);
    [~,h]=findpeaks(abs(dPh_y),'Threshold',pi/2,'MinPeakDistance',1);
    h(h>L1-4 | h<4)=[];
    for inp=1:length(h)
        dPh_y(h(inp)-2:h(inp)+2)=(dPh_y(h(inp)-3)+dPh_y(h(inp)+3))/2;
    end
    fdPh_y=dPh_y;
    bk=find(abs(fdPh_y(2:end-1))>pi/2)+1;

    fdPh_y(bk)=(fdPh_y(bk+1)+fdPh_y(bk-1))/2;
    Pt=(nanstd((fdPh_y)));
    fdPh_y=smooth(fdPh_y,3);
    [~,Phase_disc]=findpeaks(abs(diff(fdPh_y)),'MinPeakProminence',Pt,'MinPeakDistance',1);
    Ph_disc=zeros(size(bp));
    Ph_disc(Phase_disc)=1;

    NOsc=length(OsC);
    OStp=[];


    cutoff=max(rmst.^2,Med_abp);
    Wrp2=M_Abp-cutoff;
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

        if OStp(i,2)-OStp(i,1)>=mWL
            logic_p(OStp(i,1):OStp(i,2))=1;
        end
    end

    L_osc=diff(OStp,[],2);
    inv_osc=find(L_osc<mWL);
    OStp(inv_osc,:)=[];
    OsC(inv_osc)=[];


    if ~isempty(OStp)
        OS_points{fr}=OStp;
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
