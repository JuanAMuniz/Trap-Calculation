%% Free space lattice for AC stark shifts calculation

clc
close all
clear all
format long g

% Define parameters

w0_red = 42e-6;
det_red = -160e9 ;
Pred = 120e-3; %in W
theta = (pi/180)*0;

fCs = 351.72196e12;
c = 299792458;
eta = 376.73;
e0 = 8.854187817e-12;
kb = 1.3806488e-23;
f_red = fCs+det_red; lambda_red = c./f_red; zR_red = pi*w0_red^2/lambda_red;

r.x = [0:100e-9:w0_red*3];
r.y = r.x;
r.z = [-lambda_red/4:lambda_red/2/40:lambda_red/4];
[X,Z] = meshgrid(r.x,r.z);

% Lattice electric field
E0_red = sqrt(4*eta/(pi*w0_red^2)); %amplitude
E_red = E0_red*1./(sqrt(1+(Z./zR_red).^2)).*exp(-(X.^2)./(w0_red^2*(1+(Z./zR_red).^2))); % gaussian field
Etot_redx = E_red.*exp(1i*(2*pi/lambda_red).*Z)+cos(theta)*E_red.*exp(-1i*(2*pi/lambda_red).*Z); %X component
Etot_redy = sin(theta)*E_red.*exp(-1i*(2*pi/lambda_red).*Z); %Y cmponent

% Polarizability
plotpolarizability = 0;
if plotpolarizability == 1;
    for i = 1:1;
    f_pol = [c./0.9e-6:1e9:c/0.84e-6];
    polout = CesiumDLinePol(f_pol,'lightshift','au');

    figure
    subplot(3,3,1)
    plot(c./f_pol/1e-6,polout(1:3:size(polout,1),1))
    hold on
    plot(c./f_pol/1e-6,polout(1:3:size(polout,1),2))
    ylim([-1e5 1e5])
    xlabel('\lambda (um)');ylabel('\alpha_{scalar}(\lambda)');legend('6S_{1/2} F=3','6S_{1/2} F=4');title('\alpha_{scalar}(\lambda) for 6S_{1/2} F=3 and 6S_{1/2} F=4')

    subplot(3,3,2)
    plot(c./f_pol/1e-6,polout(2:3:size(polout,1),1))
    hold on
    plot(c./f_pol/1e-6,polout(2:3:size(polout,1),2))
    ylim([-1e5 1e5])
    xlabel('\lambda (um)');ylabel('\alpha_{vector}(\lambda)');legend('6S_{1/2} F=3','6S_{1/2} F=4');title('\alpha_{vector}(\lambda) for 6S_{1/2} F=3 and6S_{1/2} F=4')

    subplot(3,3,3)
    plot(c./f_pol/1e-6,polout(3:3:size(polout,1),1))
    hold on
    plot(c./f_pol/1e-6,polout(3:3:size(polout,1),2))
    ylim([-0.5e2 0.5e2])
    xlabel('\lambda (um)');ylabel('\alpha_{tensor}(\lambda)');legend('6S_{1/2} F=3','6S_{1/2} F=4');title('\alpha_{tensor}(\lambda) for 6S_{1/2} F=3 and 6S_{1/2} F=4')

    subplot(3,3,4)
    plot(c./f_pol/1e-6,polout(1:3:size(polout,1),2))
    hold on
    plot(c./f_pol/1e-6,polout(1:3:size(polout,1),4))
    ylim([-1e5 1e5])
    xlabel('\lambda (um)');ylabel('\alpha_{scalar}(\lambda)');legend('6S_{1/2} F=4','6P_{1/2} F=4');title('\alpha_{scalar}(\lambda) for 6S_{1/2} F=4 and 6P_{1/2} F=4')

    subplot(3,3,5)
    plot(c./f_pol/1e-6,polout(2:3:size(polout,1),2))
    hold on
    plot(c./f_pol/1e-6,polout(2:3:size(polout,1),4))
    ylim([-1e5 1e5])
    xlabel('\lambda (um)');ylabel('\alpha_{vector}(\lambda)');legend('6S_{1/2} F=4','6P_{1/2} F=4');title('\alpha_{vector}(\lambda) for 6S_{1/2} F=4 and 6P_{1/2} F=4')

    subplot(3,3,6)
    plot(c./f_pol/1e-6,polout(3:3:size(polout,1),2))
    hold on
    plot(c./f_pol/1e-6,polout(3:3:size(polout,1),4))
    ylim([-0.5e2 0.5e2])
    xlabel('\lambda (um)');ylabel('\alpha_{tensor}(\lambda)');legend('6S_{1/2} F=4','6P_{3/2} F=5');title('\alpha_{tensor}(\lambda) for 6S_{1/2} F=4 and 6P_{1/2} F=4')

    subplot(3,3,7)
    plot(c./f_pol/1e-6,polout(1:3:size(polout,1),2))
    hold on
    plot(c./f_pol/1e-6,polout(1:3:size(polout,1),8))
    ylim([-1e5 1e5])
    xlabel('\lambda (um)');ylabel('\alpha_{scalar}(\lambda)');legend('6S_{1/2} F=4','6P_{3/2} F=5');title('\alpha_{scalar}(\lambda) for 6S_{1/2} F=4 and 6P_{3/2} F=5' )

    subplot(3,3,8)
    plot(c./f_pol/1e-6,polout(2:3:size(polout,1),2))
    hold on
    plot(c./f_pol/1e-6,polout(2:3:size(polout,1),8))
    ylim([-1e5 1e5])
    xlabel('\lambda (um)');ylabel('\alpha_{vector}(\lambda)');legend('6S_{1/2} F=4','6P_{3/2} F=5');title('\alpha_{vector}(\lambda) for 6S_{1/2} F=4 and 6P_{3/2} F=5' )

    subplot(3,3,9)
    plot(c./f_pol/1e-6,polout(3:3:size(polout,1),2))
    hold on
    plot(c./f_pol/1e-6,polout(3:3:size(polout,1),8))
    ylim([-0.5e2 0.5e2])
    xlabel('\lambda (um)');ylabel('\alpha_{tensor}(\lambda)');legend('6S_{1/2} F=4','6P_{1/2} F=4');title('\alpha_{tensor}(\lambda) for 6S_{1/2} F=4 and 6P_{3/2} F=5' )

    figure
    subplot(2,1,1)
    plot(c./f_pol/1e-6,polout(1:3:size(polout,1),1))
    hold on
    plot(c./f_pol/1e-6,polout(1:3:size(polout,1),2))
    ylim([-1e5 1e5])
    xlabel('\lambda (um)');ylabel('\alpha_{scalar}(\lambda)');legend('6S_{1/2} F=3','6S_{1/2} F=4');

    subplot(2,1,2)
    plot(c./f_pol/1e-6,(polout(1:3:size(polout,1),2)-polout(1:3:size(polout,1),1))./polout(1:3:size(polout,1),2))
    xlabel('\lambda (um)');ylabel('(\alpha_{scalar}^{6S_{1/2}F=4}(\lambda)-\alpha_{scalar}^{6S_{1/2}F=3}(\lambda))/\alpha_{scalar}^{6S_{1/2}F=4}(\lambda)');
    ylim([-0.5 0.5])
    end
end

Etot2 = sqrt(abs(Etot_redx).^2+abs(Etot_redy).^2);

figure
contourf(r.x/1e-6,r.z/1e-6,Etot2,50)

% Add the right polarization, polarizability, power for real case

E.x(:,:,1) = Etot_redx;
E.y(:,:,1) = Etot_redy;
E.z(:,:,1) = 0*Etot_redx;

rr.x = r.z;
rr.y = r.x;

Flist = [3 4 3 4 2 3 4 5];
P = [Pred];
f = [f_red];

% Magnetic field components and hyperfine Lande g-factor
B.B = [0,0,0];
B.gf = [-1/4,1/4,-1/12,1/12,-2/3,0,4/15,2/5];

% Before executing go back to original directory
H = Hstark2D_B(E,rr,f,P,Flist,B);

% Allocate space
U_starkF3_yz = zeros(numel(rr.x),numel(rr.y),7);
U_starkF4_yz = zeros(numel(rr.x),numel(rr.y),9);
U_starkF4p_yz = zeros(numel(rr.x),numel(rr.y),9);
U_starkF5p_yz = zeros(numel(rr.x),numel(rr.y),11);

tic
for i = 1:numel(rr.x);
    for j = 1:numel(rr.y);
        U_starkF3_yz(i,j,:) = H{i,j,1};  %is an array for F=3, S 1/2, it has the different "mf" from the diagonalization
        U_starkF4_yz(i,j,:) = H{i,j,2};  %is an array for F=4, S 1/2
        U_starkF4p_yz(i,j,:) = H{i,j,4}; %is an array for F'=4, P 1/2
        U_starkF5p_yz(i,j,:) = H{i,j,8}; %is an array for F'=5, P 3/2,
    end
end
toc

figure
subplot(3,3,1);pcolor(rr.x/1e-6,rr.y/1e-6,U_starkF3_yz(:,:,4)');shading interp;axis tight;colorbar
subplot(3,3,4);plot(rr.x/1e-6,squeeze(U_starkF3_yz(:,round(size(rr.y,2)/2),:)));shading interp;axis tight
subplot(3,3,7);plot(rr.y/1e-6,squeeze(U_starkF3_yz(round(size(rr.x,2)/2),:,:)));shading interp;axis tight
subplot(3,3,2);pcolor(rr.x/1e-6,rr.y/1e-6,U_starkF4_yz(:,:,5)');shading interp;axis tight;colorbar
subplot(3,3,5);plot(rr.x/1e-6,squeeze(U_starkF4_yz(:,round(size(rr.y,2)/2),:)));shading interp;axis tight
subplot(3,3,8);plot(rr.y/1e-6,squeeze(U_starkF4_yz(round(size(rr.x,2)/2),:,:)));shading interp;axis tight
subplot(3,3,3);pcolor(rr.x/1e-6,rr.y/1e-6,U_starkF5p_yz(:,:,5)');shading interp;axis tight;colorbar
subplot(3,3,6);plot(rr.x/1e-6,squeeze(U_starkF5p_yz(:,round(size(rr.y,2)/2),:)));shading interp;axis tight
subplot(3,3,9);plot(rr.y/1e-6,squeeze(U_starkF5p_yz(round(size(rr.x,2)/2),:,:)));shading interp;axis tight

% Probe 4-5'
figure
subplot(3,2,1);plot(rr.x/1e-6,squeeze(U_starkF4_yz(:,1,:))); xlabel('z (\mu m)');ylabel('U_{6S_{1/2},F=4} (mK)');title('U_{6S_{1/2},F=4}(0,0,z) (mK)')
subplot(3,2,2);plot(rr.y/1e-6,squeeze(U_starkF4_yz(1,:,:))); xlabel('x (\mu m)');ylabel('U_{6S_{1/2},F=4} (mK)');title('U_{6S_{1/2},F=4}(x,0,0) (mK)')
subplot(3,2,3);plot(rr.x/1e-6,squeeze(U_starkF5p_yz(:,1,:))); xlabel('z (\mu m)');ylabel('U_{6P_{3/2},F=5} (mK)');title('U_{6P_{3/2},F=5}(0,0,z) (mK)')
subplot(3,2,4);plot(rr.y/1e-6,squeeze(U_starkF5p_yz(1,:,:))); xlabel('x (\mu m)');ylabel('U_{6P_{3/2},F=5} (mK)');title('U_{6P_{3/2},F=5}(x,0,0) (mK)')
subplot(3,2,5);plot(rr.x/1e-6,squeeze(U_starkF4p_yz(:,1,:))); xlabel('z (\mu m)');ylabel('U_{6P_{1/2},F=4} (mK)');title('U_{6P_{1/2},F=4}(0,0,z) (mK)')
subplot(3,2,6);plot(rr.y/1e-6,squeeze(U_starkF4p_yz(1,:,:))); xlabel('x (\mu m)');ylabel('U_{6P_{1/2},F=4} (mK)');title('U_{6P_{1/2},F=4}(x,0,0) (mK)')

% Probe 4-5'
figure
subplot(3,2,1);plot(rr.x/1e-6,squeeze(U_starkF4_yz(:,1,:))/(kb*1e-3)); xlabel('z (\mu m)');ylabel('U_{6S_{1/2},F=4} (mK)');title('U_{6S_{1/2},F=4}(0,0,z) (mK)')
subplot(3,2,2);plot(rr.y/1e-6,squeeze(U_starkF4_yz(1,:,:))/(kb*1e-3)); xlabel('x (\mu m)');ylabel('U_{6S_{1/2},F=4} (mK)');title('U_{6S_{1/2},F=4}(x,0,0) (mK)')
subplot(3,2,3);plot(rr.x/1e-6,squeeze(U_starkF5p_yz(:,1,:))/(kb*1e-3)); xlabel('z (\mu m)');ylabel('U_{6P_{3/2},F=5} (mK)');title('U_{6P_{3/2},F=5}(0,0,z) (mK)')
subplot(3,2,4);plot(rr.y/1e-6,squeeze(U_starkF5p_yz(1,:,:))/(kb*1e-3)); xlabel('x (\mu m)');ylabel('U_{6P_{3/2},F=5} (mK)');title('U_{6P_{3/2},F=5}(x,0,0) (mK)')
subplot(3,2,5);plot(rr.x/1e-6,squeeze(U_starkF4p_yz(:,1,:))/(kb*1e-3)); xlabel('z (\mu m)');ylabel('U_{6P_{1/2},F=4} (mK)');title('U_{6P_{1/2},F=4}(0,0,z) (mK)')
subplot(3,2,6);plot(rr.y/1e-6,squeeze(U_starkF4p_yz(1,:,:))/(kb*1e-3)); xlabel('x (\mu m)');ylabel('U_{6P_{1/2},F=4} (mK)');title('U_{6P_{1/2},F=4}(x,0,0) (mK)')



Ues = mean(U_starkF5p_yz(:,:,:),3);
Ugs = mean(U_starkF3_yz(:,:,:),3);




%% Integrate

% Constants
hbar = 1.0545718*1e-34;
mCs = 2.20695*1e-25;
kB = 1.38065*1e-23;
kHz = 1e3;
MHz = 1e6;
Gamma = 2*pi*5.2*MHz;
uK = 1e-6;

delta = 2*pi*[-50:0.5:50]*MHz;
rcord = (rr.y)';
zcord= rr.x;

figure
subplot(2,1,1)
pcolor(rcord/1e-6,zcord/1e-9,Ues/(kb*1e-3));shading flat;colorbar
xlabel('r (um)');ylabel('z (nm)');title('Excited state potential (mK)')
subplot(2,1,2)
pcolor(rcord/1e-6,zcord/1e-9,Ugs/(kb*1e-3));shading flat;colorbar
xlabel('r (um)');ylabel('z (nm)');title('Ground state potential (mK)')

%% Expected shifts

zind = round(size(Ugs,1)/2);
U0e = squeeze(U_starkF5p_yz(zind,1,:));
U0g = Ugs(zind,1)*ones(size(U0e));
DeltaU = (U0e - U0g)/(hbar*2*pi)/(1e6);
DeltaU = unique(round(DeltaU*100)/100);
figure
plot(DeltaU,'--o')
ylabel('Shift (MHz)');title('Expected shift for current trap ')

% 
% figure('units','normalized','outerposition',[0 0 1 1])
% filename = 'trapatchip_lowT2_BE.gif';
% subplot(4,2,[1:2]);pcolor(rcord/1e-6,zcord/1e-9,Ugs/(kb*1e-3));shading flat;colorbar
% xlabel('r (um)');ylabel('z (nm)');title('Ground state potential (mK)')
% 
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d'
% startZ = [1 0 30e-9 0];
% startR = [1 0 5e-6 0];
% zind = round(size(Ugs,1)/2);
% T=[5:5:250]*uK;
% %T=[5:20:-Ugs(zind,1)/kB/uK]*uK;
% 
% for i=1:numel(T);
%     i
%     [lshT(i,:),nden] = lineshape(delta,rcord,zcord,Ugs,Ues,T(i));
%     subplot(4,2,[3:4]);pcolor(rcord/1e-6,zcord/1e-9,nden);shading flat; colorbar;xlabel('r (um)');ylabel('z (um)');title(['Density at T = ',num2str(T(i)/uK),'uK'])
%     subplot(4,2,5);plotyy(zcord/1e-9,Ugs(:,1)/uK/kB,zcord/1e-9,nden(:,1)/max(nden(:,1)));nZfit = fit(zcord',nden(:,1)/max(nden(:,1)),gaussEqn,'Start', startZ); Zwidth(i) = nZfit.c;xlabel('z (nm)');legend('Trap potential along z (uK)','Density profile along z');title('Trap and density along z')
%     subplot(4,2,7);plot(nZfit,zcord,nden(:,1)/max(nden(:,1)));xlabel('z ');title('Gaussian fit')
%     subplot(4,2,6);plotyy(rcord/1e-6,Ugs(zind,:)/uK/kB,rcord/1e-6,nden(zind,:)/max(nden(zind,:)));nRfit = fit(rcord,nden(zind,:)'/max(nden(zind,:)),gaussEqn,'Start', startR); Rwidth(i) = nRfit.c;xlabel('r (um)');legend('Trap potential along r (uK)','Density profile along r');title('Trap and density along r')
%     subplot(4,2,8);plot(nRfit,rcord,nden(zind,:)/max(nden(zind,:)));xlabel('r ');title('Gaussian fit')
%     pause(1)
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1;
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
% end
% 
% lshT(i+1,:) = lineshape(delta,rcord,zcord,0*Ugs,0*Ues,T(end));
% 
% figure
% plot(delta/(2*pi*MHz),lshT)
% xlabel('Detuning (MHz)'),ylabel('Lineshape');title(['Lineshape for different temperatures'] )
% figure
% pcolor(delta/(2*pi*MHz),T/uK,lshT(1:end-1,:));shading interp;colorbar
% xlabel('Detuning (MHz)'),ylabel('Temperature (uK)');title('Lineshapes for different temperatures')
% 
% % Harmonic fitting
% 
% zlim = 10;
% pz = polyfit(zcord(zind-zlim:zind+zlim)',Ugs(zind-zlim:zind+zlim,1),2);
% rlim = 150;
% pr = polyfit(rcord(1:rlim)',Ugs(zind,1:rlim),2);
% 
% figure
% subplot(2,1,1)
% plot(zcord,Ugs(:,1))
% hold on
% plot(zcord(zind-zlim:zind+zlim),pz(1)*zcord(zind-zlim:zind+zlim).^2+pz(2)*zcord(zind-zlim:zind+zlim)+pz(3),'ro')
% fz = sqrt(pz(1)*2/mCs)/2/pi
% xlabel('z (m)');ylabel('Ugs(z,0) (J)');title(['Harmonic fitting near trap bottom: fz = ',num2str(fz/1e3),'kHz'])
% subplot(2,1,2)
% plot(rcord,Ugs(zind,:))
% hold on
% plot(rcord(1:rlim),pr(1)*rcord(1:rlim).^2+pr(2)*rcord(1:rlim)+pr(3),'ro')
% fr = sqrt(pr(1)*2/mCs)/2/pi
% xlabel('r (m)');ylabel('Ugs(z0,r) (J)');title(['Harmonic fitting near trap bottom: fr = ',num2str(fr/1e3),'kHz'])
% 
% Wr = 2*pi*fr;
% Wz = 2*pi*fz;
% U0 = Ugs(zind,1);
% 
% for k = 1:numel(T);
%     for i=1:numel(delta);
%         det = delta(i);
%         dN = @(x,y) (x.*exp(-(mCs*Wr^2.*(x.^2))./(2*kB*T(k)))).*(exp(-(mCs*Wz^2.*(y.^2))./(2*kB*T(k)))).*(1./(1+4.*((hbar*det+U0+mCs*0.5.*((Wr.*x).^2+(Wz.*y).^2))./(hbar*Gamma)).^2));
%         Zr = (2*pi*kB*T(k)/(mCs*Wr^2));
%         Zz = (2*pi*kB*T(k)/(mCs*Wz^2))^(1/2);
%         Ndet(i,k) = 2*pi*integral2(dN,0,2000e-6,-40e-6,40e-6,'Method','iterated','AbsTol',1e2,'RelTol',1e-12)/(Zr*Zz);
%     end
%     Lr(k) = sqrt(hbar/(mCs*Wr))*sqrt((1+exp(-hbar*Wr/(kB*T(k))))/(1-exp(-hbar*Wr/(kB*T(k)))));
%     Lz(k) = sqrt(hbar/(mCs*Wz))*sqrt((1+exp(-hbar*Wz/(kB*T(k))))/(1-exp(-hbar*Wz/(kB*T(k)))));
%     k
% end
% 
% figure
% pcolor(delta/(2*pi*MHz),T/uK,Ndet');shading interp;colorbar
% 
% Ndetnorm = bsxfun(@rdivide,Ndet,max(Ndet)); % Normalized (scaled) matrix by column
% NlshTnorm = bsxfun(@rdivide,lshT(1:end-1,:)',max(lshT(1:end-1,:)'));
% 
% figure
% subplot(2,1,1)
% pcolor(delta/(2*pi*MHz),T/uK,Ndetnorm');shading interp;colorbar
% xlabel('Detuning (MHz)'),ylabel('Temperature (uK)');title('Lineshapes for different temperatures - HO')
% subplot(2,1,2)
% pcolor(delta/(2*pi*MHz),T/uK,NlshTnorm');shading interp;colorbar
% xlabel('Detuning (MHz)'),ylabel('Temperature (uK)');title('Lineshapes for different temperatures - OL')
% 
% figure
% plot(T/uK,Lr/1e-6)
% hold on
% plot(T/uK,Lz/1e-9,'r')
% plot(T/uK,Rwidth/1e-6,'g')
% plot(T/uK,Zwidth/1e-9,'k')
% xlabel(' T(uK)')
% ylabel('Thermal width')
% legend('HO Lr (um)','HO Lz (nm)','OL Lr (um)','OL Lz (nm)')
% title('Thermal length for HO and OL model')
% 
% %% Expected shifts
% U0e = squeeze(U_starkF5p_yz(zind,1,:));
% U0g = Ugs(zind,1)*ones(size(U0e));
% DeltaU = (U0e - U0g)/(hbar*2*pi)/(1e6);
% DeltaU = unique(round(DeltaU*100)/100);
% figure
% plot(DeltaU,'--o')
% ylabel('Shift (MHz)');title('Expected shift for current trap ')

%% Scattering rate

Gamma_sc = (Gamma/det_red/2/pi)*(mean(U0g)/hbar);
Trec = 198e-9;
t_trap = -(mean(U0g)/kb/2)/((1/3)*Trec*Gamma_sc)

