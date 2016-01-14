function out = CesiumDLinePol(frq,varargin)
% CesiumDLinePol(frq) returns the scalar, vector (axial) and tensor
% polarizabilities for the hyperfine states of the 6S1/2, 6P1/2, and 6P3/2
% states of Cesium-133 for a list of frequencies (Hz). The output is
% formatted into 25 rows:
%
%  1: frequency (Hz)
%  2: 6S1/2 F=3 scalar polarizability
%  3: 6S1/2 F=3 vector polarizability
%  4: 6S1/2 F=3 tensor polarizability
%  5: 6S1/2 F=4 scalar polarizability
%  6: 6S1/2 F=4 vector polarizability
% ...
% 25: 6P3/2 F=5 tensor polarizability
%
% CesiumDLinePol takes two optional arguments. Specifying 'lightshift'
% returns the polarizabilities in a 3x8 matrix. Specifying 'au' returns the
% polarizabilities in atomic units. Ex.
%
% >> CesiumDLinePol(351.72e12,'lightshift','au')
%
% ans =
%
% 1.0e+07 *
%
%    0.4047    2.2331   -0.0001   -0.0001    0.2435    0.9470   -0.3453   -0.1261
%   -0.3038    2.2335   -0.0001    0.0001   -0.4870         0   -0.5524   -0.3781
%         0         0         0         0   -0.0696    0.7892   -0.1381    0.1260
%
% Cesium transition data are taken from Ref. [1] and [2], hyperfine energy
% corrections are calculated using Ref. [3], polarizabilities are calulated
% using Eq. (3) of Ref. [4], and hyperfine energy states are incorporated
% using equation (2.18) from Ref. [5], but modified for the Wigner-Eckart
% normalization from Ref. [4]. More detailed notes are in the Evernote
% notebook, entitled "Cesium D Line Polarizability Calculation."
%
% REFERENCES
% [1] F. Le Kien, P. Schneeweiss, and A. Rauschenbeutel, Eur. Phys. J. D
% 67, 92 (2013).
% [2] B. Arora, M. Safronova, and C. Clark, Phys. Rev. A 76, 052509 (2007).
% [3] D. A. Steck, Cesium D Line Data, available online at
% http://steck.us/alkalidata (revision 2.1.4, 23 December 2010).
% [4] P. Rosenbusch, S. Ghezali, V. Dzuba, V. Flambaum, K. Beloy, and
% A. Derevianko, Phys. Rev. A 79, 013404 (2009).
% [5] M. J. Martin, Quantum Metrology and Many-Body Physics: Pushing the
% Frontier of the Optical Lattice Clock, University of Colorado, 2013.

%% Read varargin
numvarargs = length(varargin);
lsmode = 0; aumode = 0;
if numvarargs > 2
    error('CesiumDLinePol:TooManyInputs',...
        'Takes at most two optional inputs.');
end

if sum(strcmp(varargin,'lightshift'))
    lsmode = 1;
end

if sum(strcmp(varargin,'au'))
    aumode = 1;
end

recognizedArguments = 0;
for ii = 1:length(varargin)
    recognizedArguments = recognizedArguments ...
        + sum(strcmp(varargin{ii},{'au','lightshift'}));
end

if recognizedArguments < length(varargin)
    error('CesiumDLinePol:UnrecognizedArguments',...
        'One or more of the optional arguments are not recognized.')
end

%% Specify constants
% (SI unless otherwise specified)
const.eps0 = 625000/(pi*22468879468420441); %vacuum permittivity
% const.eps0 = 8.854e-12; % permittivity of free space
const.c = 299792458; % speed of light
const.hbar = 1.05457148e-34; % reduced Planck constant
const.I = 7/2; % Cesium-133 nuclear spin
const.a0 = 5.29177211e-11; % Bohr radius
const.qe = 1.6021766e-19; % elementary charge
const.auPol = (4*pi*const.eps0*const.a0^3); % atomic unit polarizability

%% Import transition data
% from Ref. [1] and [2]; table is formatted
% n' L' J' n L J \nu (Hz) d (a.u.)
load('CsDipoleTransitions.mat');
transitions = CsDipoleTransitions;
% partition into three manifolds: 6S1/2, 6P1/2, and 6P3/2
S1hTransitions = transitions(transitions(:,4)==6 ...
    & transitions(:,5) == 0 & transitions (:,6) ==1/2,:);
P1hTransitions = transitions(transitions(:,4)==6 ...
    & transitions(:,5) == 1 & transitions (:,6) ==1/2,:);
P3hTransitions = transitions(transitions(:,4)==6 ...
    & transitions(:,5) == 1 & transitions (:,6) ==3/2,:);

%% Calculate polarizability coefficients
% Outputs formatted
% n' L' J' n L J F \nu d a0coeff a1coeff a2coeff
S1hPol = pol(S1hTransitions,const);
P1hPol = pol(P1hTransitions,const);
P3hPol = pol(P3hTransitions,const);

%% Calculate propagators and polarizabilities
% 6S1/2
[g0,g1] = propagator(S1hPol,2*pi*frq,const);
S1hAlpha0 = repmat(S1hPol(:,10),[1,length(frq)]).*g0;
S1hAlpha1 = repmat(S1hPol(:,11),[1,length(frq)]).*g1;
S1hAlpha2 = repmat(S1hPol(:,12),[1,length(frq)]).*g0;

% 6P1/2
[g0,g1] = propagator(P1hPol,2*pi*frq,const);
P1hAlpha0 = repmat(P1hPol(:,10),[1,length(frq)]).*g0;
P1hAlpha1 = repmat(P1hPol(:,11),[1,length(frq)]).*g1;
P1hAlpha2 = repmat(P1hPol(:,12),[1,length(frq)]).*g0;

% 6P3/2
[g0,g1] = propagator(P3hPol,2*pi*frq,const);
P3hAlpha0 = repmat(P3hPol(:,10),[1,length(frq)]).*g0;
P3hAlpha1 = repmat(P3hPol(:,11),[1,length(frq)]).*g1;
P3hAlpha2 = repmat(P3hPol(:,12),[1,length(frq)]).*g0;

%% Sum over transitions
[unq,~,ic] = unique(S1hPol(:,5:8),'rows','stable');
[row,~] = size(unq);
[unq2,~,ic2] = unique(P1hPol(:,5:8),'rows','stable');
[row2,~] = size(unq2);
[unq3,~,ic3] = unique(P3hPol(:,5:8),'rows','stable');
[row3,~] = size(unq3);
rowptr = 1;
out = zeros(3*row + 3*row2 + 3*row3 + 1,length(frq));
out(1,:) = frq;
for ii = 1:row
    out(rowptr + 1,:) = sum(S1hAlpha0(ic==ii,:));
    out(rowptr + 2,:) = sum(S1hAlpha1(ic==ii,:));
    out(rowptr + 3,:) = sum(S1hAlpha2(ic==ii,:));
    rowptr = rowptr + 3;
end

for ii = 1:row2
    out(rowptr + 1,:) = sum(P1hAlpha0(ic2==ii,:));
    out(rowptr + 2,:) = sum(P1hAlpha1(ic2==ii,:));
    out(rowptr + 3,:) = sum(P1hAlpha2(ic2==ii,:));
    rowptr = rowptr + 3;
end

for ii = 1:row3
    out(rowptr + 1,:) = sum(P3hAlpha0(ic3==ii,:));
    out(rowptr + 2,:) = sum(P3hAlpha1(ic3==ii,:));
    out(rowptr + 3,:) = sum(P3hAlpha2(ic3==ii,:));
    rowptr = rowptr + 3;
end

if aumode
    out(2:end,:) = out(2:end,:)/const.auPol;
end

if lsmode
    [~,cols] = size(out); out2 = zeros([3*cols,8]);
    for ii = 1:cols
        out2((1:3)+3*(ii-1),:) = reshape(out(2:end,ii),[3,8]);
    end
    out = out2;
end
end

function out = pol(transitions,const)
% pol(transitions,const) returns polarizability coefficients for the
% input transitions.
%% Add hyperfine states and splittings
tnsns = hyperfine(transitions,const);
tnsns = hfs(tnsns,const);

%% calculate scalar, vector and tensor prefactors (Ref. [4] & [5])
[unq2,~,ic] = unique(tnsns(:,[3:4,7:8]),'rows','stable');
a0u = zeros(length(unq2),1); a1u = a0u; a2u = a0u;
for ii = 1:length(unq2)
    Jp = unq2(ii,1); Fp = unq2(ii,2); J = unq2(ii,3); F = unq2(ii,4);
    K = 0; % scalar coefficient
    a0u(ii) = sqrt(2*K+1)*(-1)^(K+2*F)...
        *Wigner6jcoeff(1,1,K,F,F,Fp)...
        *(2*F+1)*(2*Fp+1)...
        *(-1)^(Fp+F+Jp+J+2*const.I+Jp-J)...
        *Wigner6jcoeff(J,Jp,1,Fp,F,const.I)^2;
    a0u(ii) = a0u(ii)*1/sqrt(3)/sqrt(2*F+1);
    K = 1; % vector coefficient
    a1u(ii) = sqrt(2*K+1)*(-1)^(K+2*F)...
        *Wigner6jcoeff(1,1,K,F,F,Fp)...
        *(2*F+1)*(2*Fp+1)...
        *(-1)^(Fp+F+Jp+J+2*const.I+Jp-J)...
        *Wigner6jcoeff(J,Jp,1,Fp,F,const.I)^2;
    a1u(ii) = a1u(ii)*(-1)/sqrt(2)/sqrt(2*F+1)*2*F/sqrt(F*(F+1));
    K = 2; % tensor coefficient
    a2u(ii) = sqrt(2*K+1)*(-1)^(K+2*F)...
        *Wigner6jcoeff(1,1,K,F,F,Fp)...
        *(2*F+1)*(2*Fp+1)...
        *(-1)^(Fp+F+Jp+J+2*const.I+Jp-J)...
        *Wigner6jcoeff(J,Jp,1,Fp,F,const.I)^2;
    a2u(ii) = a2u(ii)*(-2)/sqrt(6)*2*F*(2*F-1)...
        *sqrt(factorial(2*F-2)/factorial(2*F+3));
end
a0 = a0u(ic); a1 = a1u(ic); a2 = a2u(ic);
djj = tnsns(:,10)*(const.qe*const.a0);
% Our convention is opposite, requiring negation
tnsns(:,10) = -a0.*djj.^2;
tnsns(:,11) = -a1.*djj.^2;
tnsns(:,12) = -a2.*djj.^2;
out = tnsns;
end

function out = hyperfine(in,const)
% hyperfine(in,const) adds hyperfine transitions to transition matrix.
%% select initial state lj
[Ci,~,ic] = unique(in(:,2:3),'rows','stable');
[row, ~] = size(Ci);
% determine number of F states per manifold, preallocate space
numFi = abs(abs(Ci(:,2)+const.I)-abs(Ci(:,2)-const.I)) + 1;
Fi = zeros(sum(numFi),1); id = Fi';
% loop over initial states
rwptr = 1;
for ii = 1:row
    J = Ci(ii,2);
    ext = [abs(J-const.I),abs(J+const.I)]; % extrema
    rng = min(ext):1:max(ext); % range
    Fi(rwptr:rwptr+length(rng)-1) = rng; % fill in Fs
    id(rwptr:rwptr+length(rng)-1) = ii; % for indexing
    rwptr = rwptr + numFi(ii); % point row to next block
end
%% select final state lj
[Cf,~,fc] = unique(in(:,5:6),'rows','stable');
[row, ~] = size(Cf);
% determine number of F states per manifold, preallocate space
numFf = abs(abs(Cf(:,2)+const.I)-abs(Cf(:,2)-const.I)) + 1;
Ff = zeros(sum(numFf),1); fd = Ff';
% loop over initial states
rwptr = 1;
for ii = 1:row
    J = Cf(ii,2);
    ext = [abs(J-const.I),abs(J+const.I)]; % extrema
    rng = min(ext):1:max(ext); % range
    Ff(rwptr:rwptr+length(rng)-1) = rng; % fill in Fs
    fd(rwptr:rwptr+length(rng)-1) = ii; % for indexing
    rwptr = rwptr + numFf(ii); % point row to next block
end
%% populate output matrix
[~, col] = size(in);
% determine number of rows we're going to need
numRows = numFi(ic).*numFf(fc)';
% preallocate output matrix (adding two columns)
out = zeros(sum(numRows),col+2);
% kludgey way to generalize to other atoms?
numFi = vertcat(1,numFi)-1;
numFf = vertcat(1,numFf)-1;
rwptr = 1;
for ii=1:length(numRows)
    % fill in input values
    out(rwptr:rwptr+numRows(ii)-1,[1:3,5:7,9:10]) = ...
        repmat(in(ii,:),[numRows(ii),1]);
    % fill in F, F'
    fi = Fi(ic(ii)+sum(numFi(1:ic(ii))):ic(ii)+sum(numFi(1:ic(ii)))+numFi(ic(ii)+1));
    ff = Ff(fc(ii)+sum(numFf(1:fc(ii))):fc(ii)+sum(numFf(1:fc(ii)))+numFf(fc(ii)+1));
    out(rwptr:rwptr+numRows(ii)-1,4) = repmat(fi,[numel(ff),1]);
    tmp = repmat(ff,[1,numel(fi)]);
    out(rwptr:rwptr+numRows(ii)-1,8) = reshape(tmp',[numel(tmp),1]);
    % step row pointer
    rwptr = rwptr + numRows(ii);
end
end

function out = hfs(transitions,const)
% out = hfs(transitions,constants) calculates the hyperfine energy shift
% for the specified transition in Joules. This formula is given in Dan
% Steck's cesiumnumbers.pdf in Eq. (16) (Ref. [3]), and the values are
% taken from table 5.

%% local constants
Ag = 2.2981579425e9; % magnetic dipole constant for 6S1/2, in Hz
A1 = 291.920175e6; % magnetic dipole constant for 6P1/2, in Hz
A2 = 50.2882723e6; % magnetic dipole constant for 6P3/2, in Hz
B = -0.493417e6; % electric quadrupole constant for 6P3/2, in Hz
C = 0.56070e3; % magnetic octupole constant for 6P3/2, in Hz

%% subselect D1 and D2 transitions (n=6)
idx = transitions(:,1) == 6 & transitions(:,5) == 6;
tnsns = transitions(idx,:);

%% calculate K, K'
Fp = tnsns(:,4);
Jp = tnsns(:,3);
F = tnsns(:,8);
J = tnsns(:,7);
II = const.I;
K = F.*(F+1)-II*(II+1)-J.*(J+1);
Kp = Fp.*(Fp+1)-II*(II+1)-Jp.*(Jp+1);

%% calculate ground state shift
groundStateShift = Ag.*K/2;

%% calculate excited state shift
excitedStateShift = (A1*(Jp<=1/2)+A2*(Jp>1/2)).*Kp/2;
tmp = (Jp>1/2).*((3/2)*(Kp.*(Kp+1)-2*II*(II+1)*Jp.*(Jp+1))...
    ./(4*II*(2*II-1)*Jp.*(2*Jp-1))*B +...
    (5*Kp.^2.*(Kp/4+1) + Kp.*(II*(II+1)+Jp.*(Jp+1)...
    +3-3*II*(II+1)*Jp.*(Jp+1))...
    -5*II*(II+1)*Jp.*(Jp+1))...
    ./(II*(II-1)*(2*II-1)*Jp.*(Jp-1).*(2*Jp-1))*C);
tmp(isnan(tmp)) = 0;
excitedStateShift = excitedStateShift + tmp;

%% modify appropriate transitions
transitions(idx,9) = transitions(idx,9) ...
    + excitedStateShift - groundStateShift;
out = transitions;
end

function [g0, g1] = propagator(states,omega,const)
% [g0 g1] = propagator(states,omega,const) takes the polarizability
% tables and uses them to calculate the propagator for the states
% tabulated.

[Unq,~,ic] = unique(states(:,9),'rows','stable');
g0u = zeros(length(Unq),1); g1u = g0u;
g0 = zeros(length(ic),length(omega)); g1 = g0;
for jj = 1:length(omega);
    for ii = 1:length(Unq)
        omgij = 2*pi*Unq(ii);
        g0u(ii) = const.hbar^(-1)*((omgij-omega(jj))^(-1)...
            +(omgij+omega(jj))^(-1));
        g1u(ii) = const.hbar^(-1)*((omgij-omega(jj))^(-1)...
            -(omgij+omega(jj))^(-1));
    end
    g0(:,jj) = g0u(ic); g1(:,jj) = g1u(ic);
end
end
