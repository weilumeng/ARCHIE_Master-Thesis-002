%% Case Data
clear;
tic
casedata='IEEE30bus.m';
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
q_cost=0.0*gencostdata(:,6)'; %% for now 0% of active power costs


% Assigning the delivery factors to the generation units
Aeq=ones(1,ng);
for i=1:nb
   for j=1:ng
        if gendata(j,1)==i
            Aeq(j)=1;
        end
    end
end

% Demand with Delivery Factor
p_demand=ones(1,nb)*busdata(:,3);
q_demand=ones(1,nb)*busdata(:,4)-ones(1,nb)*Bs;

%% Setting up Matrices for GSF Formulation
A1=GSF_PP*Ag;
A2=GSF_PQ*Ag;
pd=busdata(:,3);
qd=busdata(:,4)-Bs;

b1=prat(1:length(branchdata(:,1)))+GSF_PP*(pd)+GSF_PQ*(qd);
b2=prat(1:length(branchdata(:,1)))-GSF_PP*(pd)-GSF_PQ*(qd);

A_P=[A1; -A1];
A_Q=[A2; -A2];
b=[b1;b2];

%% Voltage formulation
% Active power voltage coupling
AP_Voltage=(Xmat(nb+1:end,1:nb)./baseMVA)*Ag;

% Reactive power voltage coupling
AQ_Voltage=(Xmat(nb+1:end,nb+1:end)./baseMVA)*Ag;


B_voltage_max=vmax+Xmat(nb+1:end,1:nb)*pd./baseMVA+Xmat(nb+1:end,nb+1:end)*qd./baseMVA;

B_voltage_min=vmin+Xmat(nb+1:end,1:nb)*pd+Xmat(nb+1:end,nb+1:end)*qd;

%% Generation limit for each generator
p_lb=gendata(:,10);
p_ub=gendata(:,9);
q_lb=gendata(:,5);
q_ub=gendata(:,4);

%q_ub=[-5.4; 1.8; 34.5; 32.5; 7; 37.5];


%% Quadratic cost functions for each generator
Qmatrix=zeros(size(A_P,2),size(A_P,2));
for i=1:length(gencostdata(:,5))
    Qmatrix(i,i)=gencostdata(i,5);
end
for k=length(gencostdata(:,5))+1:2*length(gencostdata(:,5))
    Qmatrix(k,k)=q_cost(k-length(gencostdata(:,5)));
end

%% Quadratic Programming with Gurobi Solver
names={num2str(zeros(1,2*length(gendata(:,1))))};
for i=1:length(gendata(:,1))
names(i) = {num2str("P"+num2str(i))};
end
for k=length(gendata(:,1))+1:2*length(gendata(:,1))
    names(k) = {num2str("Q"+num2str(k-length(gendata(:,1))))};
end
model.varnames = names;
model.obj = [p_cost zeros(1,length(p_cost))];
model.Q = sparse(Qmatrix);
model.A = [sparse(A_P) sparse(A_Q);sparse(AP_Voltage) sparse(AQ_Voltage);sparse(AP_Voltage) sparse(AQ_Voltage);sparse(Aeq) sparse(zeros(1,size(Aeq,2))); sparse(zeros(1,size(Aeq,2))) sparse(Aeq)];
model.sense = [repmat('<',size(A_P,1),1); repmat('<',size(AP_Voltage,1),1); repmat('<',size(AQ_Voltage,1),1);repmat('=',size(Aeq,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = full([b(:);B_voltage_min(:); B_voltage_max(:); p_demand(:);q_demand(:)]);
if ~isempty(p_lb)
    model.lb = [p_lb; q_lb];
else
    model.lb = -inf(size(model.A,2),1); % default lb for MATLAB is -inf
end
if ~isempty(p_ub)
    model.ub = [p_ub; q_ub];
else
    model.ub= inf(size(model.A,2),1);
end

gurobi_write(model, 'DCOPF_QP_Q.lp');

results = gurobi(model);

%% Lagrange multipliers

lambda.lower = max(0,results.rc);
lambda.upper = -min(0,results.rc);
lambda.ineqlin = -results.pi(1:size(A_P,1));
lambda.ineqlin_voltage= -results.pi(size(A_P,1)+1:end-2);
lambda.eqlin = results.pi(end-2:end);


%% Generation Cost, Congestion Cost and LMP for each bus
generationcost=lambda.eqlin;
p_congestioncost=(lambda.ineqlin'*[-GSF_PP ; GSF_PP])';
q_congestioncost=(lambda.ineqlin'*[-GSF_PQ ; GSF_PQ])';
p_lmp=generationcost(1)+p_congestioncost;
q_lmp=generationcost(2)+q_congestioncost;


%% Net injections
pn=Ag*results.x(1:ng)-pd;
qn=Ag*results.x(ng+1:end)-qd;


Voltageatbus=Xmat(nb+1:end,1:nb)*pn+Xmat(nb+1:end,nb+1:end)*qn;

V1=zeros(1,1);
for k=1:nb
    V1=V1+Xmat(nb+1,k)*pn(k)+Xmat(nb+1,nb+k)*qn(k);
end
V13=zeros(1,1);
for k=1:nb
    V13=V13+Xmat(nb+13,k)*pn(k)+Xmat(nb+13,nb+k)*qn(k);
end

Voltagecostmin=(lambda.ineqlin_voltage(1:nb)'*Xmat(nb+1:end,nb+1:end))';

Voltagecostmax=(lambda.ineqlin_voltage(nb+1:end)'*Xmat(nb+1:end,nb+1:end))';


%% Lineflow based on GSF_PP and net injections
P_lineflow=GSF_PP*pn;
PQ_lineflow=GSF_PQ*qn;
Q_lineflow=GSF_QQ*qn;
QP_lineflow=GSF_QP*pn;

lineflow=[P_lineflow+PQ_lineflow, Q_lineflow+QP_lineflow];
toc