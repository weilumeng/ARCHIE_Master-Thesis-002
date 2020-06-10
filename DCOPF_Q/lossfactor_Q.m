%% NOT SURE : Active and reactive power losses lead to the most problems so I'm not quite sure
%% if everything is being handled correctly or if the equations in the paper are wrong 

%% Resistance and Reactance of the lines with tap ratios and status check for inactive lines
R=status.*branchdata(:,3);
X_res=status.*branchdata(:,4);
R=R.*tap;
X_res=X_res.*tap;

%% Aggregated loss of the system
%based on Equation (11) of "A Linear LMP Model for Active and Reactive Power with Power Loss"

% Line flow in p.u. 
P_lineflow_base=lineflow(:,1)./baseMVA;

% Line flow in p.u. squared for the loss formulation
P_lineflow_quad=P_lineflow_base.*P_lineflow_base;


% Active Power loss for each line in p.u.
Ploss=P_lineflow_quad.*R;  
% Active power system loss in MW
Ptotalloss=sum(Ploss)*baseMVA;

% Reactive Power loss for each line in p.u.
Qloss=(P_lineflow_quad.*X_res);
% Reactive Power system loss in MVar
Qtotalloss=sum(Qloss)*baseMVA;

%% Loss and Delivery Factor

% Active Power Loss Factor
LF_P=(2.*P_lineflow_base.*GSF_PP)'*R;

% Reactive Power Loss Factor
LF_Q=(2.*P_lineflow_base.*GSF_PQ)'*X_res;

%Delivery Factor for Active and Reactive Power
DF_P=1-LF_P;
DF_Q=1-LF_Q;

%% Fictitious Nodal Demand (Active power)
%% Cft is used to check if branches are connected to a bus or not to allocate the losses of the branches to the right buses
E_P=zeros(nb,1);
for i=1:nb
    for k=1:nl
        if full(Cft(k,i))~=0
        E_P(i)=E_P(i)+0.5*baseMVA*P_lineflow_quad(k)*R(k);
        end
    end
end
%% Fictitious Nodal Demand (Reactive power)
E_Q=zeros(nb,1);
for i=1:nb
    for k=1:nl
        if full(Cft(k,i))~=0
        E_Q(i)=E_Q(i)+0.5*baseMVA*P_lineflow_quad(k)*X_res(k);
        end
    end
end

%% Convergence criterion

% Should converge towards 0
DifferencePloss=abs(Ptotalloss-Ploss_est);
DifferenceQloss=abs(Qtotalloss-Qloss_est);

%% Updating the Iteration parameters
% Estimated LF, DF, E and Ploss updated with the damping method

%System losses
Ploss_est=w*Ploss_est+(1-w)*Ptotalloss;
Qloss_est=w*Qloss_est+(1-w)*Qtotalloss;

%Delivery Factors
DF_P_est=w*DF_P_est+(1-w)*DF_P;
DF_Q_est=w*DF_Q_est+(1-w)*DF_Q;

%Fictinous nodal demand (active power)
E_P_est_old=w*E_P_est_old+(1-w)*E_P_est;
E_P_est=w*E_P_est+(1-w)*E_P;

%Fictinous nodal demand (reactive power)
E_Q_est_old=w*E_Q_est_old+(1-w)*E_Q_est;
E_Q_est=w*E_Q_est+(1-w)*E_Q;

%% More termination criterions

% Injected power after convergence
scheduled_gen=sum(results.x(1:ng));
scheduled_gen_Q=sum(results.x(ng+1:end));
% Losses after convergence
scheduled_losses=sum(Ag*results.x(1:ng)-pd);
scheduled_losses_Q=sum(Ag*results.x(ng+1:end)-qd);
% Should convergence towards 0
scheduled_RESULT=ones(1,nb)*(Ag*results.x(1:ng)-pd)-Ploss_est;
scheduled_RESULTQ=ones(1,nb)*(Ag*results.x(ng+1:end)-qd)-Qloss_est;