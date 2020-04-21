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
Ploss=lineflow_quad.*R*baseMVA;
Ptotalloss_base=(lineflow_quad'*R);
Ptotalloss=Ptotalloss_base*baseMVA;


%% Loss and Delivery Factor


LFx=(2.*lineflow_base.*PTDF)'*R;


%Delivery Factor
DF=1-LFx;


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
Ploss_est=w*Ploss_est+(1-w)*Ptotalloss;
LF_est=LFx;
DF_est=w*DF_est+(1-w)*DF;
E_est_old=w*E_est_old+(1-w)*E_est;
E_est=w*E_est+(1-w)*E;


%Das sollte gleich losses sein
scheduled_gen=sum(results.x);
scheduled_losses=sum(Ag*results.x-pd);
scheduled_RESULT=ones(1,nb)*(Ag*results.x-pd)-Ploss_est;