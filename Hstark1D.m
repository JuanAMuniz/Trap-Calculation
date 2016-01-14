function [eH,H]=Hstark1D(E,r,fl,Pl,Flist)

%%  This program calculates the eigenenergies for each hyperfine state of interest at each position. This progrm assumes the electric filed is defined along one line.
%% Inputs:  E: electric field, in a 3D structure, in cartesian components. Can contain different profiles, according to how many modes are incorporated. The first component refers to the position index and the second to the diferent wavelengths.
%%          r: structure containing the different positions
%%          fl: array containing the different frequencies at which we want to calculate the trap potential
%%          Pl: power correspondig to each frequency
%%          Flist: lists the hyperfine states of interest
%% Output: H is a 2D cell array, where each cell contains the eigenenergies for each position (row) and each hyperfine state of interest in order (column).
%% Outline: * in case more than one electric field mode is used for the calculation, each quantity (E,wl,Pl) needs to be order in the same way.
%%          * we calculate the scalar-vector-tensor dynamical polarixabilities at the beggining and the hyperfine dependent terms of each Hamiltonian at the beggining
%%          * The matrices are written in the |F,mF> basis with |F,F>=(1,0,...) and going down

%% Initialize cell

H=cell(numel(r),numel(Flist));
eH=cell(numel(r),numel(Flist));

%% Calculate polarizabilities

alpha=CesiumDLinePol(fl,'lightshift');%Generates polariZabilities using Andrew's code

%By defult we calculate it for the ground states F=3 and 4, the D1
%excited states F'=3 and 4, and the D2 excited states F'=2,3,4and 5

alpha_scalar=alpha(1:3:3*numel(fl),:);    % scalar polarizability
alpha_vector=alpha(2:3:3*numel(fl),:);    % vector polarizability
alpha_tensor=alpha(3:3:3*numel(fl),:);   % tensor polarizability
%Each alpha contains as many rows as frequencies we are interested and each column corresponds to an hyperfne state as the ones mentioned before (in that order!) 

%% Calculate angular momentum quantities

[Fvec,Ften]=Fmat(Flist);

%Initialization
for k=1:numel(Flist);
    Hinit{k}=zeros(2*Flist(k)+1);
end

%% For loop over positions
tic
for ix=1:numel(r);
            % Auxiliary variables 
            Hmat=cell(1,numel(Flist));
            
            for w=1:numel(fl);  %loop over wavelengths
                %Reorganize electric field components in positive and
                %negative frequency components
                Ep=sqrt(Pl(w))*(1/2)*[E.x(ix,w),E.y(ix,w),E.z(ix,w)];
                Em=sqrt(Pl(w))*(1/2)*[conj(E.x(ix,w)),conj(E.y(ix,w)),conj(E.z(ix,w))];

                E2scalar=sum(Em.*Ep);
                Evector=cross(Em,Ep);
                Etensor=[Em(1)*Ep(1),Em(1)*Ep(2),Em(1)*Ep(3),Em(2)*Ep(1),Em(2)*Ep(2),Em(2)*Ep(3),Em(3)*Ep(1),Em(3)*Ep(2),Em(3)*Ep(3)];

                %Create identity matrices times E2vector and E2tensor
                %components. Ordered in F and cartesian basis style
                
                %[E2vector,E2tensor]=EIdentities(Evector,Etensor,Flist,Fmatrices_vector,Fmatrices_tensor);
                [FEv,FEt]=multiplyFandE(Flist,Fvec,Ften,Evector,Etensor);
                
                for k=1:numel(Flist);   %start loop over possible hyperfine levels
                    if isempty(Hmat{k})==1  %Start matrix as zeros
                        Hmat{k}=Hinit{k};%zeros(2*Flist(k)+1);
                    end
                    %comupte the hamiltonian matrix for each k in the |F,mF>
                    %basis defined by the total angular momentum operator
                    Hmat{k}=Hmat{k}-1i*alpha_vector(w,k)*FEv{k};%-alpha_scalar(w,k)*E2scalar*eye(2*Flist(k)+1)-1i*alpha_vector(w,k)*FEv{k}-alpha_tensor(w,k)*FEt{k};
                    if w==numel(fl);
                        [A,B]=eig(Hmat{k});
                        H{ix,k}=diag(real(B));   %final diagonalization
                        eH{ix,k}=A;
                    end
                    
                end
            end
            clear Hmat
            
        %Loop over positions     
end
toc

% This function multiplies the tensorial anglular momentum components (vector and tensor) times
% the electeic field vector and tensor components
function [FEv,FEt]=multiplyFandE(Flist,Fv,Ft,Ev,Et);

for i=1:numel(Flist);
    FEv{i}=0;
    FEt{i}=0;
    for k=1:3;
        FEv{i}=FEv{i}+Ev(k)*Fv{i,k};
    end
    for k=1:9;
        FEt{i}=FEt{i}+Et(k)*Ft{i,k};
    end
end
