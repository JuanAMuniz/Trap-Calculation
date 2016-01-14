clc
close all
clear all
format long g

% Constants
hbar = 1.0545718*1e-34;
mCs = 2.20695*1e-25;
kB = 1.38065*1e-23;
kHz = 1e3;
MHz = 1e6;
Gamma = 2*pi*5.2*MHz;
uK = 1e-6;

% Trap
U0 = kB*300*uK;
W0 = 120*1e-6;
lambda = 854.7*1e-9;
Wr = (1/W0)*sqrt(4*U0/mCs);
Wz = (2*pi/lambda)*sqrt(2*U0/mCs);

%Wr = 2*pi*1*kHz;
%Wz = 2*pi*350*kHz;
Tr = 50*uK;
Tz = 50*uK;

% Integral
delta = 2*pi*[-50:1:50]*MHz;

Zr = (2*pi*kB*Tr/(mCs*Wr^2));
Zz = (2*pi*kB*Tz/(mCs*Wz^2))^(1/2);
Lr = sqrt(hbar/(mCs*Wr))*sqrt((1+exp(-hbar*Wr/(kB*Tr)))/(1-exp(-hbar*Wr/(kB*Tr))))
Lz = sqrt(hbar/(mCs*Wz))*sqrt((1+exp(-hbar*Wz/(kB*Tz)))/(1-exp(-hbar*Wz/(kB*Tz))))
fr = Wr/(2*pi)
fz = Wz/(2*pi)

T = [10:10:200]*uK;

for k = 1:numel(T);
    Tr = T(k); Tz = T(k);
    for i=1:numel(delta);
        det = delta(i);
        potx = @(x) mCs*Wr^2.*(x.^2)./2;
        poty = @(y) mCs*Wz^2.*(y.^2)./2;
        dN = @(x,y) x.*exp(-potx(x)./(kB*Tr)).*exp(-poty(y)./(kB*Tz)).*(1./(1+4.*((hbar*det-U0+potx(x)+poty(y))./(hbar*Gamma)).^2));
        dN2 = @(x,y) (x.*exp(-(mCs*Wr^2.*(x.^2))./(2*kB*Tr))).*(exp(-(mCs*Wz^2.*(y.^2))./(2*kB*Tz))).*(1./(1+4.*((hbar*det-U0+mCs*0.5.*((Wr.*x).^2+(Wz.*y).^2))./(hbar*Gamma)).^2));
        %Ndet(i) = 2*pi*integral2(dN,0,2000e-6,-40e-6,40e-6,'Method','iterated','AbsTol',0,'RelTol',1e-8)/(Zr*Zz);
        Zr = (2*pi*kB*Tr/(mCs*Wr^2));
        Zz = (2*pi*kB*Tz/(mCs*Wz^2))^(1/2);
        Ndet(i,k) = 2*pi*integral2(dN,0,3000e-6,-40e-6,40e-6,'Method','iterated','AbsTol',1e2,'RelTol',1e-12)/(Zr*Zz);
        Ndet2(i,k) = 2*pi*integral2(dN2,0,3000e-6,-40e-6,40e-6,'Method','iterated','AbsTol',1e2,'RelTol',1e-12)/(Zr*Zz);
    end
    k
end

figure
subplot(2,2,1)
plot(delta/(2*pi*MHz),Ndet)
Nnorm = bsxfun(@rdivide,Ndet,max(Ndet))
subplot(2,2,2)
plot(delta/(2*pi*MHz),Nnorm)
subplot(2,2,3)
plot(delta/(2*pi*MHz),Ndet2)
Nnorm2 = bsxfun(@rdivide,Ndet2,max(Ndet2))
subplot(2,2,4)
plot(delta/(2*pi*MHz),Nnorm2)


figure
pcolor(Nnorm2)
shading interp

% 
% hold on
% Tr = 100*uK;
% Tz = 100*uK;
% 
% Zr = (2*pi*kB*Tr/(mCs*Wr^2));
% Zz = (2*pi*kB*Tz/(mCs*Wz^2))^(1/2);
% 
% for i=1:numel(delta)
%     det = delta(i);
%     potx = @(x) mCs*Wr^2.*(x.^2)./2;
%     poty = @(y) mCs*Wz^2.*(y.^2)./2;
%     dN = @(x,y) x.*exp(-potx(x)./(kB*Tr)).*exp(-poty(y)./(kB*Tz)).*(1./(1+4.*((hbar*det-U0+potx(x)+poty(y))./(hbar*Gamma)).^2));
%     %dN = @(x,y) (x.*exp(-(mCs*Wr^2.*x.^2)./(2*kB*Tr))).*(exp(-(mCs*Wz^2.*y.^2)./(2*kB*Tz))).*(1./(1+4.*((hbar*det-U0+mCs*0.5.*((Wr.*x).^2+(Wz.*y).^2))./(hbar*Gamma)).^2));
%     Ndet(i) = 2*pi*integral2(dN,0,100e-6,-20e-6,20e-6,'Method','iterated','AbsTol',0,'RelTol',1e-8)/(Zr*Zz);
% end
% intN = sum(Ndet)*(delta(2)-delta(1))
% 
% plot(delta/(2*pi*MHz),Ndet/max(Ndet),'r')
% xlabel('\Delta (MHz)')
% ylabel('Ntot')
% 
% 




