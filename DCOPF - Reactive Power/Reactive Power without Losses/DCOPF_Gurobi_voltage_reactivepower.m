%% Case Data
clc;clear;
tic
casedata='IEEE33bus.m';
run(fullfile(casedata))


%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('ptdf_matrix_Q.m'))

%% Setting up Matrices for the Quadratic Programming Solver

p_cost=zeros(1,nb);
for i=1:nb
    for j=1:length(gendata(:,6))
        if i==gendata(j,1)
            p_cost(i)=gencostdata(j,6); 
        end
    end
end


Aeq=zeros(1,nb);
for i=1:nb
    for j=1:length(gendata(:,1))
        if i==gendata(j,1)
            Aeq(i)=1;
        end
    end
end


% Demand P & Q

charging_injection=zeros(1,nb)*diag(BD)*baseMVA;

p_demand=ones(1,nb)*busdata(:,3);
q_demand=ones(1,nb)*busdata(:,4)-charging_injection;
%% Setting up Matrices for GSF Formulation
A1=GSF_PP;
A2=GSF_PQ;
pd=busdata(:,3);
qd=busdata(:,4);%-diag(BD)*baseMVA;


b1=prat(1:nl)+GSF_PP*(pd)+GSF_PQ*(qd);
b2=prat(1:nl)-GSF_PP*(pd)-GSF_PQ*(qd);

A_P=[A1; -A1];
A_Q=[A2; -A2];
b=[b1;b2];

%% Voltage formulation
% Active power voltage coupling

A_Voltage_P=(X(nb+1:end,1:nb))./baseMVA;

AP_Voltage=[A_Voltage_P; -A_Voltage_P];

% Reactive power voltage coupling
A_Voltage_Q=(X(nb+1:end,nb+1:end))./baseMVA;

AQ_Voltage=[A_Voltage_Q; -A_Voltage_Q];
%Voltage constraints
B_Voltage_max=(vmax-vbase)+(X(nb+1:end,1:nb)*pd+X(nb+1:end,nb+1:end)*qd)./baseMVA;

B_Voltage_min=(vmin-vbase)+(X(nb+1:end,1:nb)*pd+X(nb+1:end,nb+1:end)*qd)./baseMVA;

B_Voltage=[B_Voltage_max; -B_Voltage_min];


%% Generation limit for each generator
p_lb=zeros(1,nb);
p_ub=zeros(1,nb);
q_lb=zeros(1,nb);
q_ub=zeros(1,nb);
for i=1:nb
    for j=1:length(gendata(:,6))
        if i==gendata(j,1)
            p_lb(i)=gendata(j,10);
            p_ub(i)=gendata(j,9);
            q_lb(i)=gendata(j,5);
            q_ub(i)=gendata(j,4);
        end
    end
end

%% Quadratic cost functions for each generator
Qmatrix=zeros(2*size(A_P,2),2*size(2*A_P,2));
for i=1:nb
    for k=1:1*length(gencostdata(:,5))
        if gencostdata(k,1)~=0
            if gendata(k,1)==i
                Qmatrix(i,i)=gencostdata(k,5);
            end
        end
    end
end
for i=1:nb
    for k=1:1*length(gencostdata(:,5))
        if gencostdata(k,1)~=0
            if gendata(k,1)==i
                Qmatrix(nb+i,nb+i)=gencostdata(k,5);
            end
        end
    end
end


%% Quadratic Programming with Gurobi Solver
params.numericfocus=3;
params.FeasibilityTol=1e-5;
params.presolve=0;
params.OptimalityTol=1e-5;
params.Dualreductions=0;
names={num2str(zeros(1,2*nb))};
for i=1:nb
names(i) = {num2str("P"+num2str(i))};
end
for k=nb+1:2*nb
    names(k) = {num2str("Q"+num2str(k-nb))};
end
model.varnames = names;
model.obj = [p_cost zeros(1,length(p_cost))];
model.Q = sparse(Qmatrix);
model.A = [sparse(A_P) sparse(A_Q); sparse(AP_Voltage) sparse(AQ_Voltage); sparse(Aeq) sparse(zeros(1,size(Aeq,2))); sparse(zeros(1,size(Aeq,2))) sparse(Aeq)];
model.sense = [repmat('<',size(A_P,1),1); repmat('<',size(AP_Voltage,1),1); repmat('=',size(Aeq,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = full([b(:); B_Voltage(:); p_demand(:); q_demand(:)]);
model.lb = [p_lb, q_lb];
model.ub = [p_ub, q_ub];

gurobi_write(model, 'DCOPF_QP_Q1.lp');

results = gurobi(model,params);

%% Lagrange multipliers

lambda.lower = max(0,results.rc);
lambda.upper = -min(0,results.rc);
lambda.ineqlin = -results.pi(1:size(A_P,1));
lambda.ineqlin_voltage= -results.pi(size(A_P,1)+1:end-2);
lambda.eqlin = -results.pi(end-1:end);


%% Generation Cost, Congestion Cost and LMP for each bus
generationcost=-lambda.eqlin;
p_congestioncost=(lambda.ineqlin'*[-GSF_PP ; GSF_PP])';
q_congestioncost=(lambda.ineqlin'*[-GSF_PQ ; GSF_PQ])';
p_lmp=generationcost(1)+p_congestioncost;
q_lmp=generationcost(2)+q_congestioncost;


%% Net injections
pn=results.x(1:nb)-pd;
qn=results.x(nb+1:end)-qd;

%% Voltage at each bus
Voltageatbus=vbase+((X(nb+1:end,1:nb)*pn)+(X(nb+1:end,nb+1:end)*(qn)))./baseMVA;


%% Lambda multiplier for voltage constraints
P_Voltagecostmax=(-lambda.ineqlin_voltage(1:nb)'*X(nb+1:end,1:nb))';
P_Voltagecostmin=(lambda.ineqlin_voltage(nb+1:end)'*X(nb+1:end,1:nb))';
P_Voltagecost_total=(P_Voltagecostmax+P_Voltagecostmin)./baseMVA;

Q_Voltagecostmax=(-lambda.ineqlin_voltage(1:nb)'*X(nb+1:end,nb+1:end))';
Q_Voltagecostmin=(lambda.ineqlin_voltage(nb+1:end)'*X(nb+1:end,nb+1:end))';
Q_Voltagecost_total=(Q_Voltagecostmax+Q_Voltagecostmin)./baseMVA;

%% Lineflow based on GSF_PP and net injections
P_lineflow=GSF_PP*pn;
PQ_lineflow=GSF_PQ*qn;
Q_lineflow=GSF_QQ*qn;
QP_lineflow=GSF_QP*pn;


lineflow=[P_lineflow+PQ_lineflow, Q_lineflow+QP_lineflow];

ldispatch=[results.x(1:nb), results.x(nb+1:end)];


p_lmp=p_lmp+P_Voltagecost_total;
q_lmp=q_lmp+Q_Voltagecost_total;

toc

%% ACOPF MATPOWER
% AC Power flow
define_constants;
mpc = ext2int(loadcase('case33bw'));
%mpc.bus(:,VMAX)=1.025;
%mpc.bus(:,VMIN)=0.975;
mpopt = mpoption('model','ACOPF');
resultAC = runopf(mpc,mpopt);
BusVolAC = resultAC.bus(:,VM);
BusAglAC = resultAC.bus(:,VA);
BranchFlowAC = (resultAC.branch(:,PF)-resultAC.branch(:,PT))/2;
BranchFlowAC_Q= (resultAC.branch(:,QF)-resultAC.branch(:,QT))/2;
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
plot(1:nb,p_lmp,'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','My DCOPF');
xlabel('Bus');
ylabel('$/h');
title('ALMP Comparison');

figure
plot(1:nb,RLMPAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nb,q_lmp,'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','My DCOPF');
xlabel('Bus');
ylabel('$/h');
title('RLMP Comparison');