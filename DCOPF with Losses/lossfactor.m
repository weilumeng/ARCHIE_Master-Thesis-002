%% Resistance of the lines
R=status.*branch(:,3);
R=R.*tap;
%% Aggregated loss of the system

% Line flow in p.u. 
lineflow_base=lineflow./baseMVA;
lineflow_quad=lineflow_base.*lineflow_base;


% Net System loss and Ploss for each branch in MW
Ploss=lineflow_quad.*R;
Ptotalloss=sum(Ploss)*baseMVA;


%% Loss and Delivery Factor

% Loss Factor
LF=(2.*lineflow_base.*PTDF)'*R;
%Delivery Factor
DF=1-LF;

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

DF_est=w*DF_est+(1-w)*DF;

E_est=E;

%% More termination criterions

% Injected power after convergence
scheduled_gen=sum(results.x);
% Losses after convergence
scheduled_losses=sum(Ag*results.x-pd);
% Should convergence towards 0
scheduled_RESULT=ones(1,nb)*(Ag*results.x-pd)-Ploss_est;