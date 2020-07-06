%% Case Data
clc;clear;
tic
casedata='IEEE118bus_constrained.m';
run(fullfile(casedata))

%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('ptdf_matrix_Q.m'))

%% Generator Bus Connection Matrix

Ag=zeros(nb,ng);
for i=1:nb
    for j=1:ng
        if gendata(j,1)==i
            Ag(i,j)=1;
        end
    end
end

%% Setting up Matrices for the Quadratic Programming Solver
% Cost function to minimize
p_cost=gencostdata(:,6)';

% Left side of A*x=b ; Equality Matrix as in A*x=b
Aeq=ones(1,ng);

% Shunt injection and demand of the whole system 
charging_injection=(BD)*ones(nb,1).*baseMVA;

% Right side of A*x=b
p_demand=ones(1,nb)*busdata(:,3);
q_demand=ones(1,nb)*busdata(:,4)-sum(charging_injection);

%% Setting up Matrices for GSF Formulation
%Lineflow limit constraints based on the GSF Formulation e.g. [A1 A2]*x<[b1] 
% with x as the optimization vector consisting of [Pg;Qg] Pg=Active Power Injection Qg=Reactive Power
% Injection

A1=GSF_PP*Ag;
A2=GSF_PQ*Ag;
pd=busdata(:,3);
qd=busdata(:,4)-charging_injection;

b1=prat(1:nl)+GSF_PP*(pd)+GSF_PQ*(qd);
b2=-prat(1:nl)+GSF_PP*(pd)+GSF_PQ*(qd);

%% Voltage formulation
%Voltage constraints based on the Nodal Injection e.g. [AP_Voltage AQ_Voltage]*x<[B_voltage_max] 
% with x as the optimization vector consisting of [Pg;Qg] Pg=Active Power Injection Qg=Reactive Power
% Injection

vbase=1;            %Sets the base voltage of all buses to 1 p.u.
                    %only affects case118 and is needed for the linear 
                    %assumptions to hold


% Active power voltage coupling
AP_Voltage=X(nb+1:end,1:nb)*Ag./baseMVA;

% Reactive power voltage coupling
AQ_Voltage=X(nb+1:end,nb+1:end)*Ag./baseMVA;

% Right side of the voltage constraint
B_voltage_max=X(nb+1:end,nb+1:end)*qd./baseMVA+X(nb+1:end,1:nb)*pd./baseMVA+(vmax-vbase);

B_voltage_min=X(nb+1:end,nb+1:end)*qd./baseMVA+X(nb+1:end,1:nb)*pd./baseMVA+(vmin-vbase);

%% Generation limit for each generator

p_lb=gendata(:,10);
p_ub=gendata(:,9);
q_lb=gendata(:,5);
q_ub=gendata(:,4);

%% Quadratic cost functions for each generator
Qmatrix=zeros(2*size(A1,2),2*size(A1,2));
for i=1:length(gencostdata(:,5))
    Qmatrix(i,i)=gencostdata(i,5);
end

for i=1:length(gencostdata(:,5))
    Qmatrix(ng+i,ng+i)=gencostdata(i,5);
end

%% Quadratic Programming with Gurobi Solver
params.numericfocus=3;


names={num2str(zeros(1,2*ng))};
for i=1:ng
names(i) = {num2str("P"+num2str(i))};
end
for k=ng+1:2*ng
    names(k) = {num2str("Q"+num2str(k-ng))};
end
model.varnames = names;
model.obj = [p_cost zeros(1,length(p_cost))];
model.Q = sparse(Qmatrix);
model.A = [sparse(A1) sparse(A2);sparse(A1) sparse(A2);sparse(AP_Voltage) sparse(AQ_Voltage);sparse(AP_Voltage) sparse(AQ_Voltage);sparse(Aeq) sparse(zeros(1,size(Aeq,2))); sparse(zeros(1,size(Aeq,2))) sparse(Aeq)];
model.sense = [repmat('<',size(A1,1),1); repmat('>',size(A1,1),1); repmat('<',size(AP_Voltage,1),1); repmat('>',size(AQ_Voltage,1),1);repmat('=',size(Aeq,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = full([b1(:); b2(:); B_voltage_max(:); B_voltage_min(:); p_demand(:); q_demand(:)]);
model.lb = [p_lb; q_lb];
model.ub = [p_ub; q_ub];

gurobi_write(model, 'DCOPF_Q.lp');

results = gurobi(model,params);

%% Lagrange multipliers

lambda.ineqlin = -results.pi(1:2*size(A1,1));
lambda.ineqlin_voltage= -results.pi(2*size(A1,1)+1:end-2);
lambda.eqlin = results.pi(end-1:end);

%% Generation Cost, Congestion Cost and LMP for each bus
generationcost=lambda.eqlin;
p_congestioncost=(lambda.ineqlin'*[-GSF_PP ; -GSF_PP])';
q_congestioncost=(lambda.ineqlin'*[-GSF_PQ ; -GSF_PQ])';
p_lmp=generationcost(1)+p_congestioncost;
q_lmp=generationcost(2)+q_congestioncost;

%% Net injections
pn=Ag*results.x(1:ng)-pd;
qn=Ag*results.x(ng+1:end)-qd;


Voltageatbus=vbase+(X(nb+1:end,1:nb)*pn+X(nb+1:end,nb+1:end)*qn)./baseMVA;


Q_Voltagecostmin=(-lambda.ineqlin_voltage(1:nb)'*X(nb+1:end,nb+1:end))';
Q_Voltagecostmax=(-lambda.ineqlin_voltage(nb+1:end)'*X(nb+1:end,nb+1:end))';
Q_Voltagecosttotal=(Q_Voltagecostmin+Q_Voltagecostmax)./baseMVA;           


P_Voltagecostmin=(-lambda.ineqlin_voltage(1:nb)'*X(nb+1:end,1:nb))';
P_Voltagecostmax=(-lambda.ineqlin_voltage(nb+1:end)'*X(nb+1:end,1:nb))';
P_Voltagecosttotal=(P_Voltagecostmin+P_Voltagecostmax)./baseMVA;

%% Lineflow based on GSF_PP and net injections
P_lineflow=GSF_PP*pn;
PQ_lineflow=GSF_PQ*qn;
Q_lineflow=GSF_QQ*qn;
QP_lineflow=GSF_QP*pn;

lineflow=[P_lineflow+PQ_lineflow, Q_lineflow+QP_lineflow];


q_lmp=q_lmp+Q_Voltagecosttotal;
p_lmp=p_lmp+P_Voltagecosttotal;
toc

%% ACOPF MATPOWER
% AC Optimal Power flow
mpopt = mpoption('model','ACOPF');
mpc.bus(:,VM)=1;

resultAC = runopf(mpc,mpopt);
BusVolAC = resultAC.bus(:,VM);

BranchFlowAC = ((resultAC.branch(:,PF)-resultAC.branch(:,PT))/2);
BranchFlowAC_Q= ((resultAC.branch(:,QF)-resultAC.branch(:,QT))/2);

ALMPAC=(resultAC.bus(:,LAM_P));
RLMPAC=(resultAC.bus(:,LAM_Q));

AEA=sum(abs((p_lmp-ALMPAC)./ALMPAC))/nb;

mpopt = mpoption('model','DCOPF');
resultDC = rundcopf(mpc,mpopt);
ALMPDC=(resultDC.bus(:,LAM_P));
AEADC=sum(abs((ALMPDC-ALMPAC)./ALMPAC))/nb;

%% PLOTTING
figure
plot(1:nb,BusVolAC,'Marker','*','LineWidth',0.6,'Color','#A2142F','MarkerFaceColor','#A2142F');hold on;
plot(1:nb,Voltageatbus,'Marker','o','LineWidth',0.6, 'Color', '#EDB120','MarkerFaceColor','#EDB120'); hold on;
plot(1:nb,vmax,'Marker','v','LineWidth',1.0,'Color','#0072BD', 'MarkerSize', 4, 'MarkerFaceColor', '#0072BD'); hold on;
plot(1:nb,vmin,'Marker','^','LineWidth',1.0,'Color','#0072BD', 'MarkerSize', 4, 'MarkerFaceColor', '#0072BD')
legend('ACOPF','My DCOPF', 'Vmax', 'Vmin');
xlabel('Node');
ylabel('Voltage Magnitude');
title('Voltage Magnitude Comparison');

figure
plot(1,AEA,'Marker','o','LineWidth',0.6,'Color','#A2142F','MarkerFaceColor','#A2142F');hold on;
plot(1,AEADC,'Marker','*','LineWidth',0.6,'Color','#A2142F','MarkerFaceColor','#A2142F');hold on;
legend('Average Error : My DCOPF', 'Matpower DCOPF');
xlabel('Loadlevel');
ylabel('Error Percentage');
title('Average Error Analysis');

figure
plot(1:nl,BranchFlowAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nl,lineflow(:,1),'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','My DCOPF');
xlabel('Branch');
ylabel('MW');
title('Active Power Flow Comparison');

figure
plot(1:nl,BranchFlowAC_Q,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nl,lineflow(:,2),'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','My DCOPF');
xlabel('Branch');
ylabel('MVar');
title('Reactive Power Flow Comparison');

figure
plot(1:nb,ALMPAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nb,p_lmp,'Marker','o','LineWidth',0.5, 'Color', '#EDB120'); hold on;
plot(1:nb,ALMPDC, 'Marker','.','LineWidth',0.5, 'Color', '#0072BD')
legend('ACOPF','My DCOPF', 'DCOPF');
xlabel('Bus');
ylabel('$/MWh');
title('ALMP Comparison');

figure
plot(1:nb,RLMPAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nb,q_lmp,'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','My DCOPF');
xlabel('Bus');
ylabel('$/MVar');
title('RLMP Comparison');