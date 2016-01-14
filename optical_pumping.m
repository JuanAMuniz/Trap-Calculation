function optical_pumping(initialState,intensityVector,evolutionTime)
% out = optical_pumping(initialState,intensityVector,evolutionTime) returns
% the ground state populations as a function of time. The initial state is
% expressed as a 16x1 vector, and the pump light intensities (mW/cm^2) are
% given in a 1x12 vector: [33',34',43',44'], ordered sig-,pi,sig+. The 
% times are given in seconds.
% Ex.
%
% optical_pumping(ones([16,1])/16,[0,0.01,0.01,0],0:1e-6:1000e-6)
%
% runs the program starting in a uniform superposition, with 34' and 43'
% light, for 1 ms. Currently the program only runs for the D2 3' and 4'
% excited states.
%% define constants
a = 1.87052e8; % prefactor for pump matrix
inty = intensityVector;
%% clean up input
initialState = initialState/sum(initialState);
%% create matrices
A = defineTransitionMatrix;
D = (A.^2)./repmat(sum(A.^2),[16,1]); % normalized decay matrix
P = a*abs(A.^2).*(...
    inty(1)*padarray(padarray(padarray(eye(6),[0,1],'pre'),[1,0],...
    'post'),[9,9],'post')... % F=3 => F'=3 sigma-
    + inty(2)*padarray(eye(7),[9,9],'post')... % F=3 => F'=3 pi
    + inty(3)*padarray(padarray(padarray(eye(6),[0,1],'post'),[1,0],...
    'pre'),[9,9],'post')... % F=3 => F'=3 sigma+
    + inty(10)*padarray(circshift(eye(9),[1,0]),[7,7],...
    'pre')...% F=4 => F'=4 sigma-
    + inty(11)*padarray(eye(9),[7,7],'pre')...% F=4 => F'=4 pi
    + inty(12)*padarray(circshift(eye(9),[0,1]),[7,7],...
    'pre')...% F=4 => F'=4 sigma+
    + inty(4)*padarray(padarray(padarray(eye(7),[2,0],'post'),[7,0],...
    'pre'),[0,9],'post')... % F=3 => F'=4 sigma-
    + inty(5)*padarray(padarray(padarray(eye(7),[1,0]),[7,0],'pre'),...
    [0,9],'post')... % F=3 => F'=4 pi
    + inty(6)*padarray(padarray(padarray(eye(7),[2,0],'pre'),[7,0],...
    'pre'),[0,9],'post')... % F=3 => F'=4 sigma+
    + inty(7)*padarray(padarray(padarray(eye(7),[0,2],'pre'),[9,0],...
    'post'),[0,7],'pre')... % F=4 => F'=3 sigma-
    + inty(8)*padarray(padarray(padarray(eye(7),[0,1]),[9,0],'post'),...
    [0,7],'pre')... % F=4 => F'=3 pi
    + inty(9)*padarray(padarray(padarray(eye(7),[0,2],'post'),[9,0],...
    'post'),[0,7],'pre')... % F=4 => F'=3 sigma+
    ); % pump matrix
C = (D*P).*(~eye(16)); % mask off diagonal elements
C = C - repmat(sum(C),[16,1]).*eye(16);
%% generate plot
t = evolutionTime;
0;
W = zeros(16,length(t));
x = initialState;
for j = 1:length(t)
    W(:,j) = expm(t(j)*C)*x;
end
set(0,'defaultaxesfontname','Helvetica');
set(0,'defaulttextfontname','Helvetica');
for j = 1:16
    if j < 8
        subplot(2,9,j+1);
    else
        subplot(2,9,j+2);
    end
    area(1e3*t,W(j,:),'linewidth',0.1);
    set(gca,'YLim',[0,1]);
    [f,m] = ixfun(j);
    title(['$$|',num2str(f)',',',num2str(m),'\rangle$$'],...
        'interpreter','latex')
    xlabel('time (ms)');
end
end

function A = defineTransitionMatrix()
A = zeros(16);
% F = 3 to F' = 3
A(1,1) = -sqrt(9/32); % pi |3,-3> to |3,-3>
A(1,2) = -sqrt(3/32); % sigma- |3,-2> to |3,-3>
A(2,1) = sqrt(3/32);  % sigma+ |3,-3> to |3,-2>
A(2,2) = -sqrt(1/8);  % pi |3,-2> to |3,-2>
A(2,3) = -sqrt(5/32); % sigma- |3,-1> to |3,-2>
A(3,2) = sqrt(5/32);  % sigma+ |3,-2> to |3,-1>
A(3,3) = -sqrt(1/32); % pi |3,-1> to |3,-1>
A(3,4) = -sqrt(3/16); % sigma- |3,0> to |3,-1>
A(4,3) = sqrt(3/16);  % sigma+ |3,-1> to |3,0>
A(4,4) = 0;           % pi |3,0> to |3,0>
A(4,5) = -sqrt(3/16); % sigma- |3,1> to |3,0>
A(5,4) = sqrt(3/16);  % sigma+ |3,0> to |3,1>
A(5,5) = sqrt(1/32);  % pi |3,1> to |3,1>
A(5,6) = -sqrt(5/32); % sigma- |3,2> to |3,1>
A(6,5) = sqrt(5/32);  % sigma+ |3,1> to |3,2>
A(6,6) = sqrt(1/8);   % pi |3,2> to |3,2>
A(6,7) = -sqrt(3/32); % sigma- |3,3> to |3,2>
A(7,6) = sqrt(3/32);  % sigma+ |3,2> to |3,3>
A(7,7) = sqrt(9/32);  % pi |3,3> to |3,3>
% F = 3 to F' = 4
A(8,1) = sqrt(5/24);     % sigma- |3,-3> to |4,-4>
A(9,1) = -sqrt(5/96);    % pi |3,-3> to |4,-3>
A(9,2) = sqrt(5/32);     % sigma- |3,-2> to |4,-3>
A(10,1) = sqrt(5/672);   % sigma+ |3,-3> to |4,-2>
A(10,2) = -sqrt(5/56);   % pi |3,-2> to |4,-2>
A(10,3) = sqrt(25/224);  % sigma- |3,-1> to |4,-2>
A(11,2) = sqrt(5/224);   % sigma+ |3,-2> to |4,-1>
A(11,3) = -sqrt(25/224); % pi |3,-1> to |4,-1>
A(11,4) = sqrt(25/336);  % sigma- |3,0> to |4,-1>
A(12,3) = sqrt(5/112);   % sigma+ |3,-1> to |4,0>
A(12,4) = -sqrt(5/42);   % pi |3,0> to |4,0>
A(12,5) = sqrt(5/112);   % sigma- |3,1> to |4,0>
A(13,4) = sqrt(25/336);  % sigma+ |3,0> to |4,1>
A(13,5) = -sqrt(25/224); % pi |3,1> to |4,1>
A(13,6) = sqrt(5/224);   % sigma- |3,2> to |4,1>
A(14,5) = sqrt(25/224);  % sigma+ |3,1> to |4,2>
A(14,6) = -sqrt(5/56);   % pi |3,2> to |4,2>
A(14,7) = sqrt(5/672);   % sigma- |3,3> to |4,2>
A(15,6) = sqrt(5/32);    % sigma+ |3,2> to |4,3>
A(15,7) = -sqrt(5/96);   % pi |3,3> to |4,3>
A(16,7) = sqrt(5/24);    % sigma+ |3,3> to |4,4>
% F = 4 to F' = 3
A(1,8) = sqrt(7/72);  % sigma+ |4,-4> to |3,-3>
A(2,9) = sqrt(7/96);  % sigma+ |4,-3> to |3,-2>
A(3,10) = sqrt(5/96); % sigma+ |4,-2> to |3,-1>
A(4,11) = sqrt(5/144);% sigma+ |4,-1> to |3,0>
A(5,12) = sqrt(1/48); % sigma+ |4,0> to |3,1>
A(6,13) = sqrt(1/96); % sigma+ |4,1> to |3,2>
A(7,14) = sqrt(1/288);% sigma+ |4,2> to |3,3>
A(1,9) = sqrt(7/288); % pi |4,-3> to |3,-3>
A(2,10) = sqrt(1/24); % pi |4,-2> to |3,-2>
A(3,11) = sqrt(5/96); % pi |4,-1> to |3,-1>
A(4,12) = sqrt(1/18); % pi |4,0> to |3,0>
A(5,13) = sqrt(5/96); % pi |4,1> to |3,1>
A(6,14) = sqrt(1/24); % pi |4,2> to |3,2>
A(7,15) = sqrt(7/288);% pi |4,3> to |3,3>
A(1,10) = sqrt(1/288);% sigma- |4,-2> to |3,-3>
A(2,11) = sqrt(1/96);% sigma- |4,-1> to |3,-2>
A(3,12) = sqrt(1/48);% sigma- |4,0> to |3,-1>
A(4,13) = sqrt(5/144);% sigma- |4,1> to |3,0>
A(5,14) = sqrt(5/96);% sigma- |4,2> to |3,1>
A(6,15) = sqrt(7/96);% sigma- |4,3> to |3,2>
A(7,16) = sqrt(7/72);% sigma- |4,4> to |3,3>
% F = 4 to F' = 4
A(8,8) = -sqrt(7/30);     % pi |4,-4> to |4,-4>
A(8,9) = -sqrt(7/120);    % sigma- |4,-3> to |,4,-4>
A(9,8) = sqrt(7/120);     % sigma+ |4,-4> to |,4,-3>
A(9,9) = -sqrt(21/160);   % pi |4,-3> to |,4,-3>
A(9,10) = -sqrt(49/480);  % sigma- |4,-2> to |,4,-3>
A(10,9) = sqrt(49/480);   % sigma+ |4,-3> to |4,-2>
A(10,10) = -sqrt(7/120);  % pi |4,-2> to |4,-2>
A(10,11) = -sqrt(21/160); % sigma- |4,-1> to |4,-2>
A(11,10) = sqrt(21/160);  % sigma+ |4,-2> to |4,-1>
A(11,11) = sqrt(7/480);   % pi |4,-1> to |4,-1>
A(11,12) = -sqrt(7/48);   % sigma- |4,0> to |4,-1>
A(12,11) = sqrt(7/48);    % sigma+ |4,-1> to |4,0>
A(12,12) = 0;             % pi |4,0> to |4,0>
A(12,13) = -sqrt(7/48);   % sigma- |4,1> to |4,0>
A(13,12) = sqrt(7/48);    % sigma+ |4,0> to |4,1>
A(13,13) = sqrt(7/480);   % pi |4,1> to |4,1>
A(13,14) = -sqrt(21/160); % sigma- |4,2> to |4,1>
A(14,13) = sqrt(21/160);  % sigma+ |4,1> to |4,2>
A(14,14) = sqrt(7/120);   % pi |4,2> to |4,2>
A(14,15) = -sqrt(49/480); % sigma- |4,3> to |4,2>
A(15,14) = sqrt(49/480);  % sigma+ |4,2> to |4,3>
A(15,15) = sqrt(21/160);  % pi |4,3> to |4,3>
A(15,16) = -sqrt(7/120);  % sigma- |4,4> to |4,3>
A(16,15) = sqrt(7/120);   % sigma+ |4,3> to |4,4>
A(16,16) = sqrt(7/30);    % pi |4,4> to |4,4>
end

function [f,m] = ixfun(ind)
if ind < 8
    f = 3;
    zeeman = -3:3;
    m = zeeman(ind);
else
    f = 4;
    ind = ind - 7;
    zeeman = -4:4;
    m = zeeman(ind);
end
end
