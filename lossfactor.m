%% Resistance of the lines
R=1./G;                                                     %Resistance of the lines
for i=1:length(branchdata(:,1))                             %R=0, wenn G=0 und ansonsten R=Inf annimmt 
    if R(i)==Inf
        R(i)=0;
    end
end

%% Aggregated loss of the system


%Line flow in p.u. quadriert
lineflow_base=lineflow./baseMVA;
lineflow_quad=lineflow_base.*lineflow_base;
pn_base=pn./baseMVA;


%Ploss for each branch in MW
Ploss=lineflow_base.*R;
Ptotalloss_base=(lineflow_quad'*R);
Ptotalloss=Ptotalloss_base.*baseMVA;


%% Loss and Delivery Factor

        
LF=zeros(nb,1);
for i=1:nb
    for k=1:nl
        for j=1:nb
            LF(i)=LF(i)+(2*R(k)*PTDF(k,i)*((PTDF(k,j)*pn_base(j))));
        end
    end
end



%Delivery Factor
DF=1-LF;



%% Fictinous Nodal Demand
E=zeros(nb,1);
for i=1:nb
    for k=1:nl
        if full(Cft(k,i))~=0
        E(i)=E(i)+0.5*baseMVA*lineflow_quad(k)*R(k);
        end
    end
end



%% Abbruchkriterien


%Sollte gegen Null gehen
DifferencePloss=abs(Ptotalloss-Ploss_est);
 

%Estimated LF, DF, E and Ploss
Ploss_est=Ptotalloss;
LF_est=LF;
DF_est=DF;
E_est=E;


%Das sollte gleich losses sein
scheduled_gen=sum(results.x);
scheduled_losses=sum(results.x-pd);
scheduled_RESULT=ones(1,nb)*(results.x-pd)-Ploss_est;