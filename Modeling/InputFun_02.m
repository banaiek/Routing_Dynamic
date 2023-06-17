%%
function [Vout] = InputFun_02 (Ntot, Nsample, fs, f, Var)

p=(f/fs);
Nr=round(Ntot*(Var/100));
Rep=ceil(Ntot/Nr);
v0=round(rand(Nr,Nsample)-(.5-p));
vs=repmat(v0,Rep,1);
vs=vs(1:Ntot,:);


Vout=vs(randperm(Ntot),:);
end
