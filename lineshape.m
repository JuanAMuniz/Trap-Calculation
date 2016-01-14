function [lsh,n] = lineshape(delta,r,z,Ug,Ue,T)

% Constants and parameters
hbar = 1.0545718*1e-34;
mCs = 2.20695*1e-25;
kB = 1.38065*1e-23;
kHz = 1e3;
MHz = 1e6;
Gamma = 2*pi*5.2*MHz;
uK = 1e-6;

%Calculate partition functions
dr = r(2)-r(1);
dz = z(2)-z(1);

%Z = 2*pi*dr*dz*sum(exp(-Ug./(kB*T))*r);

for i=1:numel(delta);
    det = delta(i) ;
    Ndet(i) = sum((1/(sum(exp(-Ug./(kB*T))*r))).*exp(-Ug./(kB*T)).*(1./(1+4.*((hbar*det+Ug-Ue)./(hbar*Gamma)).^2))*r);
end

n = exp(-Ug./(kB*T));
N = sum((1/sum(exp(Ug./(kB*T))*r))*(1./(-1+exp(Ug/(kB*T))))*r)
lsh = Ndet;

end
