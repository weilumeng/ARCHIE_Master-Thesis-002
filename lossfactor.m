%Net injections
pn=Ag*results.x-pd;                                         %net injections/withdrawals
pn=pn./baseMVA;                                             %p.u. values for the lineflow

%Line flow
pijk=(PTDF*Ag*results.x-PTDF*pd)./baseMVA;                  %lineflow in p.u.
pijk_quad=pijk.^2;                                          %line transfer quadriert (lineflow)

%Resistance of the lines
R=1./G;                                                     %Resistance of the lines
for i=1:length(branchdata(:,1))                             %R=0, wenn G=0 und ansonsten R=Inf annimmt 
    if R(i)==Inf
        R(i)=0;
    end
end

%Aggregated loss of the system
Ploss=R.*(pijk_quad);                                       %Ploss p.u. for each branch
Ploss=Ploss*100;                                            %Ploss MW for each branch
Ptotalloss=ones(1,length(Ploss))*(Ploss);                   %Total Ploss in MW



%Loss Factor (not working)
LF=zeros(nb,1);
for i=1:nb
    for m=1:nl
        if full(Cft(m,i))~=0
                LF(i)=LF(i)+2*R(m)*PTDF(m,i)*pijk(m);
        end
    end
end

%Working loss factor
LFwhat=zeros(nb,1);
for i=1:nb
LFwhat(i)=sum(2*pijk.*PTDF(:,i).*R);
end

%Delivery Factor
DF=1-LFwhat;


%Fictinous Nodal Demand
E=zeros(nb,1);
for i=1:nb
    for m=1:nl
        if full(Cft(m,i))~=0
        E(i)=E(i)+0.5*Ploss(m);
        end
    end
end

%Loss Distribution Factor
LDF=zeros(nb,1);
for i=1:nb
    LDF(i)=E(i)/sum(E);
end

%Estimated LF, DF, E and Ploss
Ploss_est=0;
LF_est=0;
DF_est=1;
E_est=0;
