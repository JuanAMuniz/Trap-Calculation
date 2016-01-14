clc
close all
clc


open('enfigtmp (9).fig')
h = gcf; %current figure handle
axesObjs2 = get(h, 'Children');  %axes handles
dataObjs2 = get(axesObjs2, 'Children'); %handles to low-level graphics objects in axes

objTypes2 = get(dataObjs2, 'Type');  %type of low-level graphics object
xdata1 = get(dataObjs2, 'XData');  %data from low-level grahics objects
ydata1 = get(dataObjs2, 'YData');

%objTypes2 = get(dataObjs2{2}, 'Type');  %type of low-level graphics object
%xdata1 = get(dataObjs2{2}, 'XData');  %data from low-level grahics objects
%ydata1 = get(dataObjs2{2}, 'YData');
%xdata1 = xdata1{3}; ydata1 = ydata1{3};
[x,y]=ginput(6)
xdiff1 = x(2:end)-x(1);
mf = [5:-1:0]
figure
plot(xdiff1,'--or')
title('Peak detuning for moving lattice (45ms) at 853.5nm')
xlabel('Peak index (mf)')
ylabel('Detuning from FS resonance (MHz)')

open('853p1lattice.fig')
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
xth = get(dataObjs, 'YData');
hold on
plot(xdiff1+24,'--or')

diffx = diff(x);
mfmean = (mf(1:end-1)+mf(2:end))/2;

figure
subplot(2,1,1)
plot(xdata1,ydata1,'--o')
xlabel('Detuning (MHz)')
ylabel('Flou signal')
title('Flou signal for moving lattice (45ms) at 853.1nm')
subplot(2,1,2)
plot(mfmean,diffx,'--o')
pexp = polyfit(mfmean,diffx',1);
hold on
plot(mfmean,diff(xth),'--or')
pth = polyfit(mfmean,diff(xth),1);
plot(mfmean,pexp(1)*mfmean+pexp(2),'-.k')
xlabel('Mean mf')
ylabel('Peak difference (MHz)')
title('Peak difference for moving lattice (45ms) at 853.1nm')
legend(['Experiment p=',num2str(pexp(1))],['Calculation p=',num2str(pth(1))])







