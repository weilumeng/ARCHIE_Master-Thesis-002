%% Case Data
clc;clear;
tic
casedata='IEEE33bus.m';
run(fullfile(casedata))



%% Loadcurve

%loadprofiledata=readtable('202002110000.xlsx', 'PreserveVariableNames', true, 'Range', 'C8:C103', 'ReadVariableNames', true);
%loadcurve=table2array(loadprofiledata(:,1));
%loadcurve=str2double(loadcurve);
%loadcurve_base=loadcurve./mean(loadcurve);

%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('gsf_matrix.m'))


%% Iteration Setup
%step=0.25;
iter=1;
maxiter=100;
%minute=0;

%% Generator Bus Connection Matrix (sets the generators at the corresponding buses)

Ag=zeros(nb,ng);
for i=1:nb
    for j=1:ng
        if gendata(j,1)==i
            Ag(i,j)=1;
        end
    end
end

%% Loadcurve timesteps

%while step<24
%run(fullfile(casedata))
%busdata(:,3)=busdata(:,3)*loadcurve_base(step*4);
%iter=1;

%% Setup for estimated losses
Ploss_est=0;                %estimated active power losses
DF_P_est=ones(nb,1);        %Active Power Delivery Factor

Qloss_est=0;                %estimated reactive power losses
DF_Q_est=ones(nb,1);        %Reactive Power Delivery Factor

E_P_est=zeros(nb,1);        %Active Power Fictinous Nodal Demand
                            %Active Power Fictinous Nodal Demand of iteration-1   (used for iteration damping)

E_Q_est=zeros(nb,1);        %Reactive Power Fictinous Nodal Demand
                            %Reactive Power Fictinous Nodal Demand of iteration-1 (used for iteration damping)

gendispatch=zeros(2*ng,maxiter);                    %Generator dispatch for each iteration
mismatchdispatch=ones(ng,maxiter);                  %Mismatch of generator dispatch after each iteration as convergence criterion


%% Setting up Matrices for the Quadratic Programming Solver
% Active power linear cost function to minimize
p_cost=gencostdata(1:ng,6)';

%% Iteration for the loss formulation starts here

while iter~=maxiter

% Assigning the Active power delivery factors to the generation units
Aeq_P=ones(1,ng);
for i=1:nb
   for j=1:ng
        if gendata(j,1)==i
            Aeq_P(j)=DF_P_est(i);
        end
    end
end

% Assigning the Reactive power delivery factors to the generation units
Aeq_Q=ones(1,ng);
for i=1:nb
   for j=1:ng
        if gendata(j,1)==i
            Aeq_Q(j)=DF_Q_est(i);
        end
    end
end

% Active/Reactive Power Demand with Delivery Factor and Shunt injection
charging_injection=(BD)*ones(nb,1).*baseMVA;                %NOT SURE IF RIGHT : 
                                                            %Shunt injection solely based on Susceptance Shunt elements for
                                                            %comparisson check Equation (3) of "A Linear LMP Model for Active and Reactive Power with Power Loss"

%Right side of the Active and reactive Power balance equation with Delivery Factor and Loss
%offset (check OPF formulation of "A Linear LMP Model for Active and
%Reactive Power with Power Loss" Equation (15))

p_demand=ones(1,nb)*(busdata(:,3).*DF_P_est)-Ploss_est;
q_demand=ones(1,nb)*(busdata(:,4).*DF_Q_est)-Qloss_est-sum(charging_injection);

%% Setting up Matrices for GSF Formulation
A1=GSF_PP*Ag;                                       %Left side of line transfer limit Equation (active power)
A2=GSF_PQ*Ag;                                       %Left side of line transfer limit Equation (reactive power)
pd=busdata(:,3)+E_P_est;                            %Nodal active power demand

qd=busdata(:,4)+E_Q_est-charging_injection;         %NOT SURE IF RIGHT:
                                                    %Nodal reactive power
                                                    %demand - charging
                                                    %injection of shunt
                                                    %elements at each bus
                                                    
b1=prat(1:nl)+GSF_PP*(pd)+GSF_PQ*(qd);              %right side of line transfer limit Equation
b2=prat(1:nl)-GSF_PP*(pd)-GSF_PQ*(qd);              %right side of line transfer limit Equation

A_P=[A1; -A1];                                      %Ordering of the transfer limit equation for gurobi formulation
A_Q=[A2; -A2];                                      %for +Plim and -Plim (line transfer limit) 
b=[b1;b2];

%% Voltage formulation
%% NOT SURE if im handling the voltage constraints correct
%% Based on the paper "A Linear LMP Model for Active and Reactive Power with Power Loss"
%% What im doing is A*x=b with A being the left side of the Voltage Equation and "x" being computed by Gurobi
%% then the left side "b" is for Vmax e.g.
%% b = VoltageMax(taken from the casedata) - Voltage basepoint of each bus (taken from the casedata) ... 
%% + Partial X matrix * (Active Power Demand+Active Power Nodal Losses) + Partial X matrix * (Reactive Power Demand + Reactive Power Nodal Losses)
%% I'm dividing the left and right side by the baseMVA to have the Power Injection/Withdrawal in p.u.
%% by doing that the Voltage differences are being calculated on p.u.

%% In the paper the Nodal losses are being addedd as a positive injection at each bus, which makes absolutely no sense and has to be a error on their side
%% other papers that do the same kind of nodal losses always have the losses as a negative injection at each bus

% Active power voltage coupling
A_Voltage_P=X(nb+1:end,1:nb)*Ag./baseMVA;       %Left side of Voltage Equation (active power)

% Reactive power voltage coupling
A_Voltage_Q=X(nb+1:end,nb+1:end)*Ag./baseMVA;   %Left side of Voltage Equation (reactive power)

%Voltage constraints (Right side of Voltage Equation)
B_Voltage_max=(vmax-vbase)+(X(nb+1:end,1:nb)*(pd)+X(nb+1:end,nb+1:end)*(qd))./baseMVA;  
B_Voltage_min=(vmin-vbase)+(X(nb+1:end,1:nb)*(pd)+X(nb+1:end,nb+1:end)*(qd))./baseMVA;

B_Voltage=[B_Voltage_max; -B_Voltage_min];           %Ordering for Gurobi for <Vmax and >Vmin

%% Generation limit for each generator

p_lb=gendata(:,10);             %Active power lower generation constraint
p_ub=gendata(:,9);              %Active power upper generation constraint
q_lb=gendata(:,5);              %Reactive power lower generation constraint
q_ub=gendata(:,4);              %Reactive power upper generation constraint

%% Quadratic cost functions for each generator
Qmatrix=zeros(2*size(A_P,2),2*size(2*A_P,2));               %Active power quadratic costs
for i=1:length(gencostdata(:,5))
    Qmatrix(i,i)=0*gencostdata(i,5);
end
for i=1:length(gencostdata(:,5))
    Qmatrix(ng+i,ng+i)=gencostdata(i,5);                    %Reactive power quadratic costs
end

%% Quadratic Programming with Gurobi Solver
params.numericfocus=3;                                      %Gurobi parameter -> slows down convergence but is more precise
%params.FeasibilityTol=1e-9;
%params.OptimalityTol=1e-2;
%params.Crossover=0;
%params.presolve=0;
names={num2str(zeros(1,2*length(gendata(:,1))))};           %Gurobi names of the optimization vector "x"
for i=1:length(gendata(:,1))
names(i) = {num2str("P"+num2str(i))};
end
for k=length(gendata(:,1))+1:2*length(gendata(:,1))
    names(k) = {num2str("Q"+num2str(k-length(gendata(:,1))))};
end
model.varnames = names;
model.obj = [p_cost zeros(1,length(p_cost))];               %Gurobi linear objective costs
model.Q = sparse(Qmatrix);                                  %Gurobi quadratic objective costs
%% Left side of A*x=b
model.A = [sparse(A_P) sparse(A_Q); sparse(A_Voltage_P) sparse(A_Voltage_Q); sparse(-A_Voltage_P) sparse(-A_Voltage_Q); sparse(Aeq_P) sparse(zeros(1,size(Aeq_P,2))); sparse(zeros(1,size(Aeq_Q,2))) sparse(Aeq_Q)];
%% Equality sense as in <,> or =
model.sense = [repmat('<',size(A_P,1),1); repmat('<',size(A_Voltage_P,1),1); repmat('<',size(A_Voltage_Q,1),1); repmat('=',size(Aeq_P,1),1); repmat('=',size(Aeq_Q,1),1)];
%% Right side of A*x=b
model.rhs = full([b(:); B_Voltage(:); p_demand(:); (q_demand(:))]);
%% Generation constraints of the gurobi model
model.lb = [p_lb, q_lb];
model.ub = [p_ub, q_ub];

%% Gurobi model gets written in a file called 'DCOPF_QP_Q.lp' to check if constraints are being handled in the right way
gurobi_write(model, 'DCOPF_QP_Q.lp');
%% Gurobi solving the given model
results = gurobi(model,params);

%% Saving the generation dispatch results at each iteration
for i=iter:iter
gendispatch(:,i)=results.x;
end
%% Calculating the squared mismatch of the generation dispatches at iter -> iter+1 as a convergence criterion
if iter~=1
mismatchdispatch=(gendispatch(1:ng,iter-1)-gendispatch(1:ng,iter)).^2+(gendispatch(ng+1:end,iter-1)-gendispatch(ng+1:end,iter)).^2;
end

%% Bus injections of Active and Reactive Power 
%% = Generation-Demand-Nodal losses
%% For reactive power the shunt injection are being subtracted of the nodal demand
pn=Ag*results.x(1:ng)-pd;
qn=Ag*results.x(ng+1:end)-qd;

%% Damping parameter (necessary for >100 busses) [im using w=0.87 for case 118 (higher values lead to slower convergence but lower values might diverge]
%% the idea is to take a part of the iteration-1 to slow down a fast convergence that may lead to the model diverging
w=0.0;
%if iter~=1
%pn_old=Ag*gendispatch(1:ng,iter-1)-pd;
%pn=w*pn_old+(1-w)*pn;

%qn_old=Ag*gendispatch(ng+1:end,iter-1)-qd;
%qn=w*qn_old+(1-w)*qn;
%end

%% Lineflow based on GSF Matrix and bus injections

P_lineflow=GSF_PP*pn;
PQ_lineflow=GSF_PQ*qn;
Q_lineflow=GSF_QQ*qn;
QP_lineflow=GSF_QP*pn;

%Damped update of lineflows used for the loss calculation
%if iter~=1
%    P_lineflow_old=GSF_PP*pn_old;
%    P_lineflow=w*P_lineflow_old+(1-w)*P_lineflow;
    
%    PQ_lineflow_old=GSF_PQ*qn_old;
%    PQ_lineflow=w*PQ_lineflow_old+(1-w)*PQ_lineflow;
%end

% Matrix of active/reactive power generation dispatches
generation_dispatch=[Ag*results.x(1:ng), Ag*results.x(ng+1:end)];

% lineflow of active/reactive power
lineflow=[P_lineflow+PQ_lineflow, Q_lineflow+QP_lineflow];
%% Loss formulation
run(fullfile('lossfactor_Q.m'));

%% Criterion to stop further iterations (convergence criterion)
if any(abs(mismatchdispatch) > 1e-9)
    iter=iter+1;
else
    iter=maxiter;
end
end
toc

%% Voltage at each bus
Voltageatbus=vbase+((X(nb+1:end,1:nb)*pn)+(X(nb+1:end,nb+1:end)*qn))./baseMVA;

%% Lagrange multipliers (based on the constraints of the gurobi model for calculation of the LMPs)

lambda.lower = max(0,results.rc);
lambda.upper = -min(0,results.rc);
lambda.ineqlin = -results.pi(1:size(A_P,1));
lambda.ineqlin_voltage= -results.pi(size(A_P,1)+1:end-2);
lambda.eqlin = results.pi(end-1:end);


%% Generation Cost, Congestion Cost and LMP for each bus
generationcost=lambda.eqlin;

P_congestioncost=(lambda.ineqlin'*[-GSF_PP ; GSF_PP])';
Q_congestioncost=(lambda.ineqlin'*[-GSF_PQ ; GSF_PQ])';

P_losscost=lambda.eqlin(1)*(DF_P_est-1);
Q_losscost=lambda.eqlin(2)*(DF_Q_est-1);

%% LMP Calculation
%% Lambda multiplier for voltage constraints for LMP calculation
P_Voltagecostmax=(-lambda.ineqlin_voltage(1:nb)'*X(nb+1:end,1:nb))';
P_Voltagecostmin=(-lambda.ineqlin_voltage(nb+1:end)'*X(nb+1:end,1:nb))';
P_Voltagecost_total=(P_Voltagecostmax+P_Voltagecostmin);

Q_Voltagecostmax=(-lambda.ineqlin_voltage(1:nb)'*X(nb+1:end,nb+1:end))';
Q_Voltagecostmin=(-lambda.ineqlin_voltage(nb+1:end)'*X(nb+1:end,nb+1:end))';
Q_Voltagecost_total=(Q_Voltagecostmax+Q_Voltagecostmin);


p_lmp=generationcost(1)+P_congestioncost+P_losscost+P_Voltagecost_total;
q_lmp=generationcost(2)+Q_congestioncost+Q_losscost+Q_Voltagecost_total;

%% ACOPF MATPOWER
% AC Power flow

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