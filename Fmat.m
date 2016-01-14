function [Fmat_vector,Fmat_tensor]=Fmat(Flist);

Fmat_vector=F_vector(Flist);
Fmat_tensor=F_tensor(Flist);


function Fmat_vector=F_vector(Flist);
%% This function calculates the Fu components for certain F that will be an
%% input.
%% Inputs:  F is the hyperfine quantum number
%% Outputs: Fmat_vector is a 1x3 cell where each element is a matrix corresponding to
%%          the vector component that multiplies the electric field components in cartesian basis
%% Each matrix writen in the |F,mF> basis, where the eigenvectors are arranged such that |F,F>=(1 0 0 ...)'

% The cell Fmat_vector contains the components (u)=(x,y,z) 

Fmat_vector=cell(numel(Flist),3);

hbar=1;%1.05457173e-34; %m^2*kg/s

%% First define the matrices for Fplus, Fminus, Fz
for f=1:numel(Flist);
F=Flist(f);    
Fp=zeros(2*F+1);
Fm=zeros(2*F+1);
Fz=zeros(2*F+1);

vp=2*F+2:2*F+2:(2*F+1).^2;
mfp=F-1:-1:-F;
Fp(vp)=hbar*sqrt(F*(F+1)-mfp.*(mfp+1));

vm=2:2*F+2:(2*F+1).^2;
mfm=F:-1:-F+1;
Fm(vm)=hbar*sqrt(F*(F+1)-mfm.*(mfm-1));

mf=F:-1:-F;
vz=1:2*F+2:(2*F+1).^2;
Fz(vz)=hbar*mf;

%% Vector components

Fmat_vector{f,1}=(1/(F*2))*(Fp+Fm)/2; %u=x
Fmat_vector{f,2}=(1/(F*2))*(Fp-Fm)/(2*1i); %u=y
Fmat_vector{f,3}=(1/(F*2))*(Fz); %u=z
end

function Fmat_tensor=F_tensor(Flist);
%% This function calculates the FuFv components for certain f that will be an
%% input.
%% Inputs:  F is the hyperfine quantum number
%% Outputs: Fmat_tensor is a 1x9 cell where each element is a matrix corresponding to
%%          the tensor component that multiplies the electric field components in cartesian basis
%% Each matrix writen in the |F,mF> basis, where the eigenvectors are arranged such that |F,F>=(1 0 0 ...)'

% The cell Fmat_tensor contains the combinations (uv)=(xx,xy,xz,yz,yy,yz,zx,zy,zz) 

Fmat_tensor=cell(numel(Flist),9);

hbar=1;%1.05457173e-34; %m^2*kg/s

%% First define the matrices for Fplus, Fminus, Fz, F^2
for f=1:numel(Flist);
F=Flist(f);    
Fp=zeros(2*F+1);
Fm=zeros(2*F+1);
Fz=zeros(2*F+1);
F2=(hbar^2)*F*(F+1)*eye(2*F+1);

vp=2*F+2:2*F+2:(2*F+1).^2;
mfp=F-1:-1:-F;
Fp(vp)=hbar*sqrt(F*(F+1)-mfp.*(mfp+1));

vm=2:2*F+2:(2*F+1).^2;
mfm=F:-1:-F+1;
Fm(vm)=hbar*sqrt(F*(F+1)-mfm.*(mfm-1));

mf=F:-1:-F;
vz=1:2*F+2:(2*F+1).^2;
Fz(vz)=hbar*mf;

%% Tensor components

Fmat_tensor{f,1}=(3/(F*(2*F-1)))*((1/2)*(F2-Fz^2+(1/2)*(Fp^2+Fm^2))-(1/3)*F2); %first component u=x,v=x
Fmat_tensor{f,2}=(3/(F*(2*F-1)))*(1/2)*(1/(1i*2))*(Fp^2-Fm^2); %u=x,v=y
Fmat_tensor{f,3}=(3/(F*(2*F-1)))*(1/2)*(1/2)*(Fp*Fz+Fz*Fp+Fm*Fz+Fz*Fm); %u=x,v=z

Fmat_tensor{f,4}=Fmat_tensor{f,1}; %u=y,v=x
Fmat_tensor{f,5}=(3/(F*(2*F-1)))*((1/2)*(F2-Fz^2-(1/2)*(Fp^2+Fm^2))-(1/3)*F2); %u=y,v=y
Fmat_tensor{f,6}=(3/(F*(2*F-1)))*(1/2)*(1/(1i*2))*(Fp*Fz-Fz*Fp-Fm*Fz+Fz*Fm); %u=y,v=z

Fmat_tensor{f,7}=Fmat_tensor{f,3};%u=z,v=x
Fmat_tensor{f,8}=Fmat_tensor{f,6};%u=z,v=y
Fmat_tensor{f,9}=(3/(F*(2*F-1)))*(Fz^2-(1/3)*F2); %u=z,v=z

end
