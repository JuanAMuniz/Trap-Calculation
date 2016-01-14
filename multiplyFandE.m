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
