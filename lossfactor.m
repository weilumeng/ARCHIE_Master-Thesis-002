%% Resistance of the lines
R=1./G;                                                     
for i=1:length(branchdata(:,1))                              
    if R(i)==Inf
        R(i)=0;
    end
end

%% Aggregated loss of the system

% Line flow in p.u. 
lineflow_base=lineflow./baseMVA;
lineflow_quad=lineflow_base.*lineflow_base;
pn_base=pn./baseMVA;


% Net System loss and Ploss for each branch in MW
Ploss=lineflow_quad.*R*baseMVA;
Ptotalloss_base=(lineflow_quad'*R);
Ptotalloss=Ptotalloss_base*baseMVA;


%% Loss and Delivery Factor

% Loss Factor
LFx=(2.*lineflow_base.*PTDF)'*R;
%Delivery Factor
DF=1-LFx;

%% Fictitious Nodal Demand
E=zeros(nb,1);
for i=1:nb
    for k=1:nl
        if full(Cft(k,i))~=0
        E(i)=E(i)+0.5*baseMVA*lineflow_quad(k)*R(k);
        end
    end
end


%% Termination criterion

% Should converge towards 0
DifferencePloss=abs(Ptotalloss-Ploss_est);
 
%% Updating the Iteration parameters
% Estimated LF, DF, E and Ploss
Ploss_est=w*Ploss_est+(1-w)*Ptotalloss;
LF_est=LFx;
DF_est=w*DF_est+(1-w)*DF;
E_est_old=w*E_est_old+(1-w)*E_est;
E_est=w*E_est+(1-w)*E;

%% More termination criterions

% Injected power after convergence
scheduled_gen=sum(results.x);
% Losses after convergence
scheduled_losses=sum(Ag*results.x-pd);
% Should convergence towards 0
scheduled_RESULT=ones(1,nb)*(Ag*results.x-pd)-Ploss_est;