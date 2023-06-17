%% Spiking network modeling and dynamic routing

% This script simulates the biophysical interplay of two spiking networks emulating
% the behaviors of LPFC and ACC. The model, which gets network configuration,
% sets parameters for specific neuron-type interactions within the network, and
% generates spike times along with local field potential (LFP) geometrically.
% The main script establishes EI type population sizes, network dimensions, sampling
% frequency, and trial duration. It also sets the overall connectivity weights between
% the E/E of each individual within the two networks.
% All other parameters tied to the neuron-type specific network designs
% are defined in the set_initial_network_V01 script.
% The network employs a function, InputFun_02, to create input pulses
% with adjustable levels of coherency for each network.
% The network also illustrates the lead/lag routing dynamics of the model...

%******************* 
% Kianoush Banaie Boroujeni (Kia)- 2023
%*******************

clear 
close all



Net.NumNets = 2; %number of networks
Net.NumCells = 1000; %number of neurons, an ordinal exponent of 3

  % number of excitatory and inhibitory neurons

Net.EI.prop{1,1}= [.76 .24  ];          %ratio of number of excitatory to inhibitory neurons in network 1
Net.EI.prop{2,1} = [.72 .28];          %ratio of number of excitatory to inhibitory neurons in network 2

Net.EItypes.prop{1,1} = [.85 .1 .05];
Net.EItypes.prop{1,2} = [.30 .09 .35 .26 ]; % prop for PV, CCK, CB, CR,  in LPFC

Net.EItypes.prop{2,1} = [.85 .1 .05];
Net.EItypes.prop{2,2} = [.22 .18 .40 .20 ]; % prop for PV, CCK, CB, CR in ACC


Fs=4E3;
tau=1E3/Fs;
trialtime=4500;
W11 = .25;
W22 = .25;
Nmua=1000;

% proportion of internetwork projection type
inter_E = .95;
inter_I = .05;

WIA = .075; % inter-network connection ratio
W12 = WIA;
W21 = WIA;

w12 = 0.06./WIA;
w21 = 0.06./WIA;

load('WF_mat.mat')




%% Making a cubic structure for each network with recording sites
Net_Label={'Network1','Network2'};

NNs = Net.NumNets; %number of networks
N = Net.NumCells; %number of neurons, an ordinal exponent of 3
Ntot = NNs*N; % total number of neurons in the simulation

Np=nthroot(N,3);
Sp=round(nthroot(N,3)/2);
Ep=(Np-1); %end neuron point
Nd=-Ep:2:Ep;
Ed=2*(-Sp+1:2:Sp-1);
SpN=length(Ed);

[nx,~,~]=meshgrid(Nd);

[ex,~,~]=meshgrid(Ed);

% current nodes
nx=repmat(Nd,Np.^2,1);
nx=reshape(nx,numel(nx),1);
ny=meshgrid(Nd);
ny=repmat(reshape(ny,numel(ny),1),Np,1);
nz=repmat(Nd',Np.^2,1);
Nc=[nx,ny,nz];

% recording nodes
ex=repmat(Ed,SpN.^2,1);
ex=reshape(ex,numel(ex),1);
ey=meshgrid(Ed);
ey=repmat(reshape(ey,numel(ey),1),SpN,1);
ez=repmat(Ed',SpN.^2,1);
Ec=[ex,ey,ez];
NE=length(ex);

%pair wise distance of R nodes to C nodes
dE=zeros(NE,N);

for i=1:NE
dE(i,:)=sqrt(sum((Ec(i,:)-Nc).^2,2));
end
for i=1:N
dNs(i,:)=sqrt(sum((Nc(i,:)-Nc).^2,2));
end
  gC=exp(-(dNs/std(dNs,[],'all')).^2)+diag(-1*ones(1,N));

gC(isinf(gC))=0;
gC=1*gC./max(gC(gC>0),[],'all');

[rn,~,r_ind] = unique(dE);
r_ind = reshape(r_ind,NE,N);




cntr1=0;
cntr2=0;

Level1_name=[]; 
Level2_name=[];
Level3_name=[];

      Net.Level1.Net_ID = [1:NNs]';
      Net.Level1.Net_Size = repmat(N,NNs,1);
      Net.Level1.Net_Ind = [N:N:Ntot]';
      Level1_name=fields(Net.Level1);
      Net.Level1.Net_Spec = table(Net.Level1.Net_ID,Net.Level1.Net_Size,Net.Level1.Net_Ind,'VariableNames',Level1_name);
      Net.Level1.Net_Mat = [Net.Level1.Net_ID,Net.Level1.Net_Size,Net.Level1.Net_Ind];
      Net.Level1.Net_Mat=[zeros(1,size(Net.Level1.Net_Mat,2));Net.Level1.Net_Mat];
      
for i=1:NNs % network #

    for j=1:2 % E/I  type
        cntr1=cntr1+1;
          Net.Level2.Net_ID   (cntr1,1) = i;
          Net.Level2.Net_Size (cntr1,1) = Net.Level1.Net_Size (i);
          Net.Level2.Net_Ind  (cntr1,1) = Net.Level1.Net_Ind (i);
          Net.Level2.EI_ID   (cntr1,1)  = j;
          Net.Level2.EI_Size (cntr1,1)  = Net.EI.prop{i,1}(j) * Net.Level1.Net_Size(i);
          if cntr1>1 
          Net.Level2.EI_Ind  (cntr1,1)  = Net.Level2.EI_Ind(cntr1-1,1) + Net.Level2.EI_Size(cntr1,1);
          else
              Net.Level2.EI_Ind  (cntr1,1)  = Net.Level2.EI_Size(cntr1,1);
          end
                  
        for k=1:length(Net.EItypes.prop{i,j}) % E/I subtypes
            cntr2=cntr2+1;
                  Net.Level3.Net_ID       (cntr2,1) = i;
                  Net.Level3.Net_Size     (cntr2,1) = Net.Level1.Net_Size (i);
                  Net.Level3.Net_Ind      (cntr2,1) = Net.Level1.Net_Ind (i);
                  Net.Level3.EI_ID        (cntr2,1) = j;
                  Net.Level3.EI_Size      (cntr2,1) = Net.Level2.EI_Size (cntr1,1);
                  Net.Level3.EI_Ind       (cntr2,1) = Net.Level2.EI_Ind  (cntr1,1);
                  Net.Level3.EItypes_ID   (cntr2,1) = k;
                  Net.Level3.EItypes_Size (cntr2,1) = round(Net.EItypes.prop{i,j}(k) * Net.Level2.EI_Size (cntr1,1)); 
                 
                  if cntr2>1 
                      Net.Level3.EItypes_Ind  (cntr2,1) = Net.Level3.EItypes_Ind(cntr2-1,1) + Net.Level3.EItypes_Size(cntr2,1);
                  else
                      Net.Level3.EItypes_Ind  (cntr2,1)  = Net.Level3.EItypes_Size(cntr2,1);
                  end
                                              
        end
    end
        Net.rand_Ind {i} = randperm(N); %randomly assigning each neuron to a an index of the geometric structure
        [~,Net.rand_IndS2C{i}]= sort(Net.rand_Ind {i}); %index of the structure element assigned to the Neuron ID
end

Level2_name=fields(Net.Level2);
Net.Level2.EI_Spec = table(Net.Level2.Net_ID,Net.Level2.Net_Size,Net.Level2.Net_Ind,Net.Level2.EI_ID,...
   Net.Level2.EI_Size,Net.Level2.EI_Ind,'VariableNames',Level2_name);
Net.Level2.Net_Mat = [Net.Level2.Net_ID,Net.Level2.Net_Size,Net.Level2.Net_Ind,...
   Net.Level2.EI_ID,Net.Level2.EI_Size,Net.Level2.EI_Ind];
Net.Level2.Net_Mat=[zeros(1,size(Net.Level2.Net_Mat,2));Net.Level2.Net_Mat];
          
Level3_name=fields(Net.Level3);
Net.Level3.EI_Spec = table(Net.Level3.Net_ID,Net.Level3.Net_Size,Net.Level3.Net_Ind,Net.Level3.EI_ID,...
                      Net.Level3.EI_Size,Net.Level3.EI_Ind,Net.Level3.EItypes_ID,Net.Level3.EItypes_Size,Net.Level3.EItypes_Ind,...
                      'VariableNames',Level3_name);
Net.Level3.Net_Mat  = [Net.Level3.Net_ID,Net.Level3.Net_Size,Net.Level3.Net_Ind,Net.Level3.EI_ID,...
                      Net.Level3.EI_Size,Net.Level3.EI_Ind,Net.Level3.EItypes_ID,Net.Level3.EItypes_Size,Net.Level3.EItypes_Ind];
Net.Level3.Net_Mat  = [zeros(1,size(Net.Level3.Net_Mat,2)); Net.Level3.Net_Mat];

Wcon=nan(Ntot);
l_ind=[0;Net.Level1.Net_Ind];
for i=1:NNs
    
   Wcon(l_ind(i)+1:l_ind(i+1),l_ind(i)+1:l_ind(i+1)) = gC(Net.rand_Ind {i},Net.rand_Ind {i});
end
GL = Wcon;
GL (GL==0) = 1;
GL(isnan(GL)) = 0;
Wcon(isnan(Wcon))=0.1; % basedline inter-network connection probability



NWE12 = randperm(Net.Level2.EI_Size(1),round(Net.Level2.EI_Size(1)*(1-W12*inter_E)));
NWE21 = randperm(Net.Level2.EI_Size(3),round(Net.Level2.EI_Size(3)*(1-W21*inter_E)))+1000;

NWI12 = randperm(Net.Level2.EI_Size(2),round(Net.Level2.EI_Size(2)*(1-W12*inter_I)))+Net.Level2.EI_Size(1);
NWI21 = randperm(Net.Level2.EI_Size(4),round(Net.Level2.EI_Size(4)*(1-W21*inter_I)))+1000+Net.Level2.EI_Size(3);

% discunnected neurons
NW12 = [NWE12,NWI12];
NW21 = [NWE21,NWI21];

figure
l_EI_Ind = [0;Net.Level2.EI_Ind];
l_Net_ID = Net.Level2.Net_ID; 
l_EI_ID = Net.Level2.EI_ID;  
l_Net_Ind = [0;Net.Level1.Net_Ind];  
for i=1:NNs
    l_ind = Net.rand_Ind {i};
    subplot(1,NNs,i)
    Einds = find(l_Net_ID==i & l_EI_ID==1);
    Iinds = find(l_Net_ID==i & l_EI_ID==2);
    EInds = (l_EI_Ind(Einds)+1:l_EI_Ind(Einds+1)) - l_Net_Ind(i);
    IInds = (l_EI_Ind(Iinds)+1:l_EI_Ind(Iinds+1)) - l_Net_Ind(i);

scatter3(nx(l_ind(EInds)),ny(l_ind(EInds)),nz(l_ind(EInds)),1,'b','filled')
hold on
scatter3(nx(l_ind(IInds)),ny(l_ind(IInds)),nz(l_ind(IInds)),1,'r','filled')
scatter3(ex,ey,ez,8,'k','filled','s')

axis off
view(-63.4166*((-1).^(i+1)),38.79*((-1).^(i+1)))
title(Net_Label{i})
end
set(gca,'tickdir','out')
set(gcf,'position',[72   612   428   190])

%% main connectivity 
Net=set_initial_network_V01(W11,W22,w12,w21,Net);

%   set_initial_network_08
l_ind = [0;Net.Level3.EItypes_Ind];
l_size = Net.Level3.EItypes_Size; 
l_ID = Net.Level3.EI_ID;  
[I,J] = size(Net.C);

% connectivity mat

for i = 1:I
    for j = 1:J
        PSCSi=((-1).^(l_ID(i)-1));
        PSCSj=((-1).^(l_ID(i)-1));
        Net.Cmat(l_ind(i)+1:l_ind(i+1),l_ind(j)+1:l_ind(j+1)) = Net.C(i,j).*(.9+.1*rand(l_size(i),l_size(j))).*PSCSi;         
    end
end



%% parameter settings 
SPK_time=[];             % spike-times
tspan = 0:tau:trialtime;
Nsample = length(tspan);
SPK_time = zeros(Ntot,Nsample);
SPK_timeZ = SPK_time;
SPK_Ones = ones(Ntot,Nsample);
ext_Input = zeros(Ntot,Nsample);

for ipar=1:9
    for j=1:J
        if ipar<9
       Net.ParMat(l_ind(j)+1:l_ind(j+1),ipar) = Net.par(ipar,j).* ones(l_size(j),1)+Net.par(ipar+9,j).*(.5-rand(l_size(j),1));
        else
        ext_Input(l_ind(j)+1:l_ind(j+1),:) =  Net.par(ipar,j).* ones(l_size(j),Nsample)+.2*Net.par(ipar+9,j).*(.5-rand(l_size(j),Nsample)); 
        end
    end
end

 a = Net.ParMat(:,1);
 b = Net.ParMat(:,2);
 c = Net.ParMat(:,3);
 d = Net.ParMat(:,4);
 tr = Net.ParMat(:,5);
 td = Net.ParMat(:,6);
 Syn_amp = Net.ParMat(:,7);
 v = Net.ParMat(:,8);
  

cn = dsp.ColoredNoise(2,'SamplesPerFrame',Nsample,...
    'NumChannels',Ntot);
 
Bnoise = cn()';
Bnoise = Bnoise./max(abs(Bnoise));
Bnoise = Bnoise(:,1:Nsample);

Syn_cn = dsp.ColoredNoise(2,'SamplesPerFrame',Nsample,...
    'NumChannels',Ntot);
NG_syn = Syn_cn()';
NG_syn=NG_syn./max(abs(NG_syn));
NG_syn=(NG_syn/50);
  
%% Simulation 
C=Net.Cmat.*Wcon;
PC=rand(size(C))<abs(C);
IA1 = 1:1000;
IA2 = 1001:2000;
C(NW12,1001:2000)=0;
C(NW21,1:1000)=0;
IA_C1 = find(sum(C(1:1000,1001:2000),2));
IA_C2 = 1000+find(sum(C(1001:2000,1:1000),2));
nIA_C1 = setdiff(IA1,IA_C1);
nIA_C2 = setdiff(IA2,IA_C2);

PC(IA_C1,nIA_C2) = rand(size(C(IA_C1,nIA_C2)))<abs(C(IA_C1,nIA_C2));
PC(IA_C2,nIA_C1) = rand(size(C(IA_C2,nIA_C1)))<abs(C(IA_C2,nIA_C1));

PC(nIA_C1,IA_C1) = rand(size(C(nIA_C1,IA_C1)))<abs(C(nIA_C1,IA_C1));
PC(nIA_C2,IA_C2) = rand(size(nIA_C2,IA_C2))<abs(C(nIA_C2,IA_C2));

C=sign(C).*PC;

%% External Input 
t1s = linspace(-2*pi,2*pi,Nsample/50);
t2s = linspace(-pi,pi,Fs/2);
x1 = square(t2s,70)+2;
x2 = -square(t2s,70)+2;
 GWin=gausswin(100./tau,5)';

Input_pulse1 = [];
Input_pulse2 = [];

% desynchronized input ratio
i1=[99  5 99];
i2=[99 99 5];

InputL = length(i1);

for i=1:InputL
in1 = InputFun_02 (Ntot/2, (Nsample-1)/InputL, Fs, 5, i1(i));
in2 = InputFun_02 (Ntot/2, (Nsample-1)/InputL, Fs, 5, i2(i));
if i==InputL

in1 = InputFun_02 (Ntot/2, (Nsample-1)/InputL+1, Fs, 5,  i1(i));
in2 = InputFun_02 (Ntot/2, (Nsample-1)/InputL+1, Fs, 5,  i2(i));
end
Input_pulse1 = [Input_pulse1,in1];
Input_pulse2 = [Input_pulse2,in2];
end
    
Input_pulse=[Input_pulse1;Input_pulse2];
Input_pulse=conv2(Input_pulse,GWin,'same')+rand(Ntot, Nsample);
Input_pulse=Input_pulse./max([Input_pulse,ones(Ntot,1)],[],'all');


%%

vv=[];
VVrmv=[];
Nspiked=[];
u=b.*v;
AttenL=25*Fs/1000;
Atten=[1-hamming(2*AttenL+1)]';

%setting delay for long range connections
int_ind=[IA_C1;IA_C2];
v_int=zeros(Ntot,5./tau);
v_Nint=zeros(Ntot,1./tau);

%initialization of synaptic parameters
T=1E6*ones(Ntot,1);
v_syn=-(-(1./tr).*exp((-T-tau)./tr)+(1./td).*exp((-T-tau)./td));
G_syn=zeros(Ntot,1);
Gout = G_syn;
tm=log(tr./td)./(1./td -1./tr);
Mm=-exp(-tm./tr)+exp(-tm./td);

%initialization function of Synaptic voltage gradient
init_V = @(Tr,Td,Ts) -(-(1./Tr).*exp((-Ts+tau)./Tr)+(1./Td).*exp((-Ts+tau)./Td));

for t=1:Nsample
    
    
    I=4*ext_Input(:,t).*Input_pulse(:,t)+(Bnoise(:,t));
    vs=v;
    vs=v_Nint(:,t+1);
    vs(int_ind)=v_int(int_ind,t+1);
    Nspiked=find(v > 30);
    Nspiked2=find(vs > 30);
    SPK_time(Nspiked,t)=1;

    I=I+C'*((Gout+NG_syn(:,t)).*Syn_amp*1.3);
    v = v + tau*(0.04*v.^2+5*v+140-u+I);
    u = u + tau*a.*(b.*v-u);
    v(Nspiked)=c(Nspiked);
    u(Nspiked)= u(Nspiked) + d(Nspiked);

    T(Nspiked2)=tr(Nspiked2).*(td(Nspiked2)./(tr(Nspiked2)+td(Nspiked2))).*(Gout(Nspiked2)/(.63));
    v_syn(Nspiked2)=init_V(tr(Nspiked2),td(Nspiked2),T(Nspiked2));
    v_syn=v_syn+tau*(-v_syn.*((1./td)+(1./tr)) - ((1./(td.*tr)).*G_syn));
    G_syn=G_syn+v_syn.*tau;
    Gout = G_syn./Mm;
    GV(:,t) = Gout;
  
    
    VVrmv =[VVrmv,v];
    v_int = [v_int,v];
    v_Nint = [v_Nint,v];
    
    vv = [vv,v];
    
    mInd=min([AttenL,t-1]);
    MInd=min([AttenL,Nsample-t]);    
    SPK_timeZ(Nspiked,t-mInd:t+MInd)=1;
    SPK_Ones(Nspiked,t-mInd:t+MInd)=SPK_Ones(Nspiked,t-mInd:t+MInd).*Atten(AttenL+1-mInd:AttenL+1+MInd);

     
end
%% plot rasters
Acolor{1,1} = [.1 .1 1; .2 .2 .7; .3 .3 .5];
Acolor{1,2} = [1 .1 .1;.9 .1 .5; .7 .2 .2; .5 .3 .3];

l_ind  = [0;Net.Level3.EItypes_Ind];  
l_EItype_ID = Net.Level3.EItypes_ID;
l_EI_ID = Net.Level3.EI_ID;
l_Net_ID = Net.Level3.Net_ID;  

figure
for i=1:NNs
    iNet=find(l_Net_ID==i);
    subplot(NNs,1,i)
for j=[iNet']
    
    Inds = l_ind(j)+1:l_ind(j+1);
    rownum = Inds- l_Net_Ind(i);
    Est=SPK_time(Inds,:).*rownum';

    Est(Est==0)=nan;
    col=Acolor{1,l_EI_ID(j)}(l_EItype_ID(j),:);
    plot(tspan,Est,'.','markersize',2,'markerfacecolor',col','color',col')

    hold on
    set(gca,'tickdir','out','ylim',[1 1000])
    xlabel 'time(msec.)'
    ylabel 'Neuron#'
    title(Net_Label{i})
end
end

set(gcf,'position',[72          55        1166         738])

%% rasters for nearby neurons
r_ind1=r_ind(:,Net.rand_IndS2C{1});
r_ind2=r_ind(:,Net.rand_IndS2C{2}); 

figure
inn=27;
subplot(2,1,1)
rownum = 1:inn;
[val,ind]=sort(Wcon(500,1:1000));
inds=(ind(1000-inn+1:1000));
Est=SPK_time(inds,:).*rownum';
plot(tspan,Est,'.','color','k')
inds=(ind(1:inn));
hold on
plot(get(gca,'xlim'),[inn+1 inn+1],'k','linewidth',1)
Est=SPK_time(inds,:).*(rownum'+(inn+1));
plot(tspan,Est,'.','color','k')

subplot(2,1,2)
rownum = 1:inn;
[val,ind]=sort(Wcon(1501,1001:2000));
inds=1000+(ind(1000-inn+1:1000));
Est=SPK_time(inds,:).*rownum';
plot(tspan,Est,'.','color','k')
inds=1000+(ind(1:inn));
hold on
Est=SPK_time(inds,:).*(rownum'+(inn+1));
plot(tspan,Est,'.','color','k')
plot(get(gca,'xlim'),[inn+1 inn+1],'k','linewidth',1)

%% adding_Spike triggered average 

figure
Vout0=-(VVrmv+70);
Vout=Vout0;
for i=1:2000
A(i,:)=conv(SPK_time(i,:),WF_mat(i,:),'same');
end
sd=std(Vout0,[],2);
B=70.*A;

PinkN=dsp.ColoredNoise(1,'SamplesPerFrame',Nsample,...
    'NumChannels',Ntot)
CC=PinkN()';
cn = dsp.ColoredNoise(2,'SamplesPerFrame',Nsample,...
    'NumChannels',Ntot);

spkt_N1=[];
spkt_B1=[];
spkt_N2=[];
spkt_B2=[];

while (length(spkt_N1)<2 || length(spkt_B1)<2 || length(spkt_N2)<2 || length(spkt_B2)<2 )
iE=randi(100,1)+10;

[~,Bind_1] = min(r_ind1(iE,1:Net.Level2.EI_Ind(1)));
[~,Nind_1] = min(r_ind1(iE,Net.Level2.EI_Ind(1)+1:1000));

[~,Bind_2] = min(r_ind2(iE,1:Net.Level2.EI_Ind(3)-1000));
[~,Nind_2] = min(r_ind2(iE,Net.Level2.EI_Ind(3)-1000+1:1000));


Nind_1 = Nind_1 + Net.Level2.EI_Ind(1);
Nind_2 = Nind_2 + Net.Level2.EI_Ind(3);
Bind_2 = Bind_2 + Net.Level2.EI_Ind(2);

spkt_N1=tspan(find(SPK_time(Nind_1,:))-1);
spkt_B1=tspan(find(SPK_time(Bind_1,:))-1);
spkt_N2=tspan(find(SPK_time(Nind_2,:))-1);
spkt_B2=tspan(find(SPK_time(Bind_2,:))-1);
end

CC=CC(1:Ntot,1:Nsample);
CC=CC+cn()';
normCC=(CC./range(CC,2));
sdN=nanstd(normCC,[],2);

Spkwin=find(SPK_timeZ);
Vout(Spkwin)=Vout0(Spkwin).*SPK_Ones(Spkwin)+B(Spkwin);

Vout=Vout+10*(normCC./sdN);

%
subplot(7,2,1)
plot(tspan,nanmean(Vout(1:1000,:)),'k')
hold on
plot(tspan,nanmean(GV(1:1000,:)),'k')

plot(tspan,nanmean(Vout(1001:2000,:)),'b')
hold on
plot(tspan,nanmean(GV(1001:2000,:)),'b')

set(gca,'tickdir','out')
% plot(nanmean(-(vv+65)))


subplot(7,2,3)
hold on
%  plot(tspan,(Vout0(Nind_1,:)),'r')
plot(tspan,(Vout(Nind_1,:)),'k')
set(gca,'tickdir','out')
plot([spkt_N1', spkt_N1'],get(gca,'ylim'),'r')

subplot(7,2,5)
hold on
% plot(tspan,(Vout0(100,:)),'r')
plot(tspan,(-GV(Nind_1,:)),'r')
set(gca,'tickdir','out')
plot([spkt_N1', spkt_N1'],get(gca,'ylim'),'k')

subplot(7,2,7)
hold on
% plot(tspan,(Vout0(1100,:)),'r')
plot(tspan,(Vout(Bind_1,:)),'k')
set(gca,'tickdir','out')
plot([spkt_B1', spkt_B1'],get(gca,'ylim'),'b')

subplot(7,2,9)
hold on
% plot(tspan,(Vout0(1100,:)),'r')
plot(tspan,(GV(Bind_1,:)),'b')
set(gca,'tickdir','out')
plot([spkt_B1', spkt_B1'],get(gca,'ylim'),'b')

subplot(7,2,2)
plot(tspan,nanmean(Vout(1000:2000,:)))
hold on
plot(tspan,nanmean(GV(1000:2000,:)))
set(gca,'tickdir','out')
% plot(nanmean(-(vv+65)))

subplot(7,2,4)
hold on
% plot(tspan,(Vout0(100,:)),'r')
plot(tspan,(Vout(Nind_2,:)),'k')
set(gca,'tickdir','out')
plot([spkt_N2', spkt_N2'],get(gca,'ylim'),'r')

subplot(7,2,6)
hold on
% plot(tspan,(Vout0(100,:)),'r')
plot(tspan,(-GV(Nind_2,:)),'r')
set(gca,'tickdir','out')
plot([spkt_N2', spkt_N2'],get(gca,'ylim'),'k')

subplot(7,2,8)
hold on
% plot(tspan,(Vout0(1100,:)),'r')
plot(tspan,(Vout(Bind_2,:)),'k')
set(gca,'tickdir','out')
plot([spkt_B2', spkt_B2'],get(gca,'ylim'),'b')

subplot(7,2,10)
hold on
% plot(tspan,(Vout0(1100,:)),'r')
plot(tspan,(GV(Bind_2,:)),'b')
set(gca,'tickdir','out')
plot([spkt_B2', spkt_B2'],get(gca,'ylim'),'k')







set(gcf,'position',[42 394 1066 555]);

% Ohmic LFP signals on the recording sites

Vr1 = rn(r_ind1)./min(rn(r_ind1),[],2);
Vr2 = rn(r_ind2)./min(rn(r_ind2),[],2);
% Ohmic Monopole
mVEs1 = ((1./Vr1).^1)* (Vout(1:1000,:)+(sign(nansum(C(1:1000,:),2)).*GV(1:1000,:)));
mVEs2 = ((1./Vr2).^1)* (Vout(1001:2000,:)+(sign(nansum(C(1001:2000,:),2)).*GV(1001:2000,:)));
% Ohmic Dipole
dVEs1 = ((1./Vr1).^2)* Vout(1:1000,:)  + ((1./Vr1).^1)*(sign(nansum(C(1:1000,:),2)).*GV(1:1000,:));
dVEs2 = ((1./Vr2).^2)* Vout(1001:2000,:) + ((1./Vr2).^1)*(sign(nansum(C(1001:2000,:),2)).*GV(1001:2000,:));

subplot(7,2,11)
hold on
% plot(tspan,(Vout0(1100,:)),'r')
plot(tspan,(mVEs1(iE,:)),'k')
set(gca,'tickdir','out')
plot([spkt_B1', spkt_B1'],get(gca,'ylim'),'b')
plot([spkt_N1', spkt_N1'],get(gca,'ylim'),'r')

subplot(7,2,12)
hold on
% plot(tspan,(Vout0(1100,:)),'r')
plot(tspan,(mVEs2(iE,:)),'k')
set(gca,'tickdir','out')
plot([spkt_B2', spkt_B2'],get(gca,'ylim'),'b')
plot([spkt_N2', spkt_N2'],get(gca,'ylim'),'r')


subplot(7,2,13)
hold on
% plot(tspan,(Vout0(1100,:)),'r')
plot(tspan,(dVEs1(iE,:)),'k')
set(gca,'tickdir','out')
plot([spkt_B1', spkt_B1'],get(gca,'ylim'),'b')
plot([spkt_N1', spkt_N1'],get(gca,'ylim'),'r')

subplot(7,2,14)
hold on
% plot(tspan,(Vout0(1100,:)),'r')
plot(tspan,(dVEs2(iE,:)),'k')
set(gca,'tickdir','out')
plot([spkt_B2', spkt_B2'],get(gca,'ylim'),'b')
plot([spkt_N2', spkt_N2'],get(gca,'ylim'),'r')


figure
subplot(4,2,1)
dis=prctile(mVEs1,95,'all');
mkdis=(1:NE)*dis;
VEp1=mVEs1+mkdis';
plot(tspan,VEp1(1:10,:),'k')
title 'LPFC Ohmic, Monopole'
% set(gca,'ylim',[0.1 .6],'tickdir','out')
subplot(4,2,3)
dis=prctile(mVEs2,95,'all');
mkdis=(1:NE)*dis;
VEp2=mVEs2+mkdis';
plot(tspan,VEp2(1:10,:),'k')
title 'ACC Ohmic, Monopole'

% set(gca,'ylim',[0.1 .6],'tickdir','out')



subplot(4,2,5)
dis=prctile(mVEs1,95,'all');
mkdis=(1:NE)*dis;
VEp1=dVEs1+mkdis';
plot(tspan,VEp1(1:10,:),'k')
title 'LPFC Ohmic, dipole'

subplot(4,2,7)
dis=prctile(dVEs2,95,'all');
mkdis=(1:NE)*dis;
VEp2=dVEs2+mkdis';
plot(tspan,VEp2(1:10,:),'k')
title 'ACC Ohmic, dipole'

mean(sum(SPK_time(1:1000,50:end),2))
mean(sum(SPK_time(1001:2000,50:end),2))

%       return
%%
%computing st_MUA
%Intra
warning off

 
[~,uints1]=min(rn(r_ind1),[],2);
[~,uints2]=min(rn(r_ind2),[],2);
uints2=uints2+1000;
% Fs=8000;
sw=.5;

Fb=[800 3000];
cntri=0;
irep=1;
RatioL11 = zeros(Nmua,Nsample);
RatioL22 = zeros(Nmua,Nsample);
RatioL12 = zeros(Nmua,Nsample);
RatioL21 = zeros(Nmua,Nsample);
for j=1:irep
for i=1:Nmua
    uN=i;%randi(1000,1);
    cntri=cntri+1;
    Spks1=find(SPK_time(uN,:));
    ch1=randi(125,1);
    WB1=dVEs1(ch1,:);
    [MUA11{cntri,1},RatioL11(i,:)]=MUA_IA_trial_mdl_02(WB1,Spks1,Fs,sw,Fb);
    
    Spks2=find(SPK_time(1000+uN,:));
    ch2=randi(125,1);
    WB2=dVEs2(ch2,:);
    [MUA22{cntri,1},RatioL22(i,:)]=MUA_IA_trial_mdl_02(WB2,Spks2,Fs,sw,Fb);
    

    [MUA12{cntri,2},RatioL12(i,:)]=MUA_IA_trial_mdl_02(WB2,Spks1,Fs,sw,Fb);
    

    [MUA21{cntri,1},RatioL21(i,:)]=MUA_IA_trial_mdl_02(WB1,Spks2,Fs,sw,Fb);
end
end


%%
Nc=40;
timeW=-100:.5:100;
PreTime=find(timeW>-20 & timeW<0);
PostTime=find(timeW>0 & timeW<20);

a=cellfun(@(x) nanmean(x,1), MUA11, 'UniformOutput', false);
p=cell2mat(a);
np11=((((p-nanmean(p,2))./nanstd(p,[],2))));
np11=conv2(np11,ones(1,Nc),'same')/Nc;
S_se11=nanstd(np11)/sqrt(irep*Nmua);
S_m11=nanmean(np11);
figure, subplot(2,2,1)
shadedErrorBar(-100:.5:100,S_m11,S_se11)
set(gca,'tickdir','out','ylim',[-.5 .5],'xlim',[-75 75])
hold on
plot([0 0],get(gca,'ylim'),'color','r','linestyle','--')
plot(get(gca,'xlim'),[0 0],'color','r','linestyle','--')
Nc=40;
a=cellfun(@(x) nanmean(x,1), MUA22, 'UniformOutput', false);
p=cell2mat(a);
np22=((((p-nanmean(p,2))./nanstd(p,[],2))));
np22=conv2(np22,ones(1,Nc),'same')/Nc;
S_se22=nanstd(np22)/sqrt(irep*Nmua);
S_m22=nanmean(np22);
subplot(2,2,4)
shadedErrorBar(-100:.5:100,S_m22,S_se22)
set(gca,'tickdir','out','ylim',[-.5 .5],'xlim',[-75 75])
hold on
plot([0 0],get(gca,'ylim'),'color','r','linestyle','--')
plot(get(gca,'xlim'),[0 0],'color','r','linestyle','--')

a=cellfun(@(x) nanmean(x,1), MUA12, 'UniformOutput', false);
p=cell2mat(a);
np12=((((p-nanmean(p,2))./nanstd(p,[],2))));
np12=conv2(np12,ones(1,Nc),'same')/Nc;
S_se12=nanstd(np12)/sqrt(irep*Nmua);
S_m12=nanmean(np12);
 subplot(2,2,2)
shadedErrorBar(-100:.5:100,S_m12,S_se12)
set(gca,'tickdir','out','ylim',[-.5 .5],'xlim',[-75 75])
hold on
plot([0 0],get(gca,'ylim'),'color','r','linestyle','--')
plot(get(gca,'xlim'),[0 0],'color','r','linestyle','--')

a=cellfun(@(x) nanmean(x,1), MUA21, 'UniformOutput', false);
p=cell2mat(a);
np21=((((p-nanmean(p,2))./nanstd(p,[],2))));
np21=conv2(np21,ones(1,Nc),'same')/Nc;
S_se21=nanstd(np21)/sqrt(irep*Nmua);
S_m21=nanmean(np21);
subplot(2,2,3)
shadedErrorBar(-100:.5:100,S_m21,S_se21)
set(gca,'tickdir','out','ylim',[-.5 .5],'xlim',[-75 75])
hold on
plot([0 0],get(gca,'ylim'),'color','r','linestyle','--')
plot(get(gca,'xlim'),[0 0],'color','r','linestyle','--')
%%
S_MnPre=[nanmean(S_m11(PreTime)),nanmean(S_m12(PreTime));nanmean(S_m21(PreTime)),nanmean(S_m22(PreTime))];
S_MnPost=[nanmean(S_m11(PostTime)),nanmean(S_m12(PostTime));nanmean(S_m21(PostTime)),nanmean(S_m22(PostTime))];
S_SEPre=[nanmean(S_se11(PreTime)),nanmean(S_se12(PreTime));nanmean(S_se21(PreTime)),nanmean(S_se22(PreTime))];
S_SEPost=[nanmean(S_se11(PostTime)),nanmean(S_se12(PostTime));nanmean(S_se21(PostTime)),nanmean(S_se22(PostTime))];

 b=(1:-.05:.5)'.*([0.2 0.2 1]);
 r=(.5:.05:1)'.*([1 0.2 0.2]);
 bk2=(.25:.05:.5)'.*[0.7 0.2 0.2];
 bk1=(.5:-.05:.25)'.*[0.2 0.2 0.8];
 map=[b; bk1; 0.05 0.05 0.05; bk2; r]; 


hylim=1.2*max( [S_MnPost+S_SEPost,S_MnPre+S_SEPre],[],'all');
lylim=1.2*min( [S_MnPost-S_SEPost,S_MnPre-S_SEPre,0*S_MnPre],[],'all');



S_AreaLabel={'LPFC-LPFC','LPFC-ACC','ACC-LPFC','ACC-ACC'};
S_AreaLabelS={'LPFC','ACC'};
figure
subplot(1,3,1)
imagesc(S_MnPost-S_MnPre)
set(gca,'tickdir','out','xtick',[1:2],'xticklabel',S_AreaLabelS,'ytick',[1:2],'yticklabel',S_AreaLabelS,'clim',[-.3 .3])
subplot(1,3,2)
imagesc(S_MnPost)
set(gca,'tickdir','out','xtick',[1:2],'xticklabel',S_AreaLabelS,'ytick',[1:2],'yticklabel',S_AreaLabelS)
subplot(1,3,3)
imagesc(S_MnPre)
set(gca,'tickdir','out','xtick',[1:2],'xticklabel',S_AreaLabelS,'ytick',[1:2],'yticklabel',S_AreaLabelS)
 colormap(map)
 
S_locs=[1:4];

figure
subplot(2,3,4)
scatter(S_locs,reshape(S_MnPost,1,4),12,'filled','r')
hold on
plot([S_locs; S_locs],[reshape(S_MnPost+S_SEPost,1,4); reshape(S_MnPost-S_SEPost,1,4)],'color','r','linestyle','--')
set(gca,'tickdir','out','ylim',[lylim hylim],'xlim',[0 5],'xtick',[1:4],'xticklabel',S_AreaLabel)

subplot(2,3,5)
scatter(S_locs,reshape(S_MnPre,1,4),12,'filled','b')
hold on
plot([S_locs; S_locs],[reshape(S_MnPre+S_SEPre,1,4); reshape(S_MnPre-S_SEPre,1,4)],'color','b','linestyle','--')
set(gca,'tickdir','out','ylim',[lylim hylim],'xlim',[0 5],'xtick',[1:4],'xticklabel',S_AreaLabel)

subplot(2,3,6)
scatter(S_locs,reshape(S_MnPost-S_MnPre,1,4),12,'filled','k')
hold on
plot([S_locs; S_locs],[reshape(S_MnPost-S_MnPre-S_SEPost,1,4); reshape(S_MnPost-S_MnPre+S_SEPost,1,4)],'color','k','linestyle','--')

plot(get(gca,'xlim'),[0 0],'color','r','linestyle','--')
set(gca,'tickdir','out','xlim',[0 5],'xtick',[1:4],'xticklabel',S_AreaLabel)


 
 

%% Power Plot
Nc=40;
nnn=find(isnan(sum(np11,2)));
np11(nnn,:)=[];
nnn=find(isnan(sum(np21,2)));
np21(nnn,:)=[];
nnn=find(isnan(sum(np22,2)));
np22(nnn,:)=[];
nnn=find(isnan(sum(np12,2)));
np12(nnn,:)=[];

figure
ntp=3; [pxx,f]=pmtm(np22',(ntp:-1:1)/sum(1:ntp),'Tapers','sine',[1:100],2000);
  [pxx,f] =pspectrum(np22',timeW/1000,'leakage',1,'frequencylimit',[1 100]);
hold on;
s1=shadedErrorBar(f,nanmean(pxx'.*f'),std(pxx'.*f')/sqrt(1000),{'color',[1 .2 .2]},1)
ntp=3; [pxx,f]=pmtm(np11',(ntp:-1:1)/sum(1:ntp),'Tapers','sine',[1:100],2000);
  [pxx,f] =pspectrum(np11',timeW/1000,'leakage',1,'frequencylimit',[1 100]);
s1=shadedErrorBar(f,nanmean(pxx'.*f'),std(pxx'.*f')/sqrt(1000),{'color',[.2 .2 1]},1)


%  return
%% Time resolved Power
M = 2000;
L = 1900;
g = flattopwin(M);
Ndft = 995;
neven = ~mod(Ndft,2);
Spect1 = [];
Spect2 = [];

figure
for i=1:50
    i
[spect1,sfreqs,stimes] = spectrogram(dVEs1(i,100:end),g,L,[1:100],Fs);
[spect2,sfreqs,stimes] = spectrogram(dVEs2(i,100:end),g,L,[1:100],Fs);
Spect1(i,:,:) = abs(spect1);
Spect2(i,:,:) = abs(spect2);
end

subplot(2,2,1)
imagesc(stimes,sfreqs,squeeze(nanmean(Spect1)).*sfreqs)
  axis xy

colormap(jet)
subplot(2,2,2)
imagesc(stimes,sfreqs,squeeze(nanmean(Spect2)).*sfreqs)
 axis xy

colormap(jet)

RatioS12=conv(nanmean(RatioL12),ones(1,50*Nc),'same');
RatioS21=conv(nanmean(RatioL21),ones(1,50*Nc),'same');

RatioS12(RatioS12==0)=nan;
RatioS21(RatioS21==0)=nan;
subplot(2,2,3)
plot(tspan,RatioS12,'k')
axis tight
subplot(2,2,4)
plot(tspan,RatioS21,'k')
axis tight

