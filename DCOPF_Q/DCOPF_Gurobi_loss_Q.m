%% Case Data
clear;
tic
casedata='IEEE30bus.m';
run(fullfile(casedata))

%% Loadcurve

%loadprofiledata=readtable('202002110000.xlsx', 'PreserveVariableNames', true, 'Range', 'C8:C103', 'ReadVariableNames', true);
%loadcurve=table2array(loadprofiledata(:,1));
%loadcurve=str2double(loadcurve);
%loadcurve_base=loadcurve./mean(loadcurve);

%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('ptdf_matrix_Q.m'))

%% Iteration Setup

%step=0.25;
iter=1;
maxiter=2;
%minute=0;

%% Generator Bus Connection Matrix

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
   

%% Estimated Losses
Ploss_est=0;
Qloss_est=0;
LF_est=zeros(nb,1);
DF_est=ones(nb,1);
E_est=zeros(nb,1);
E_est_old=zeros(nb,1);
gendispatch=zeros(ng,maxiter);
mismatchdispatch=ones(ng,maxiter);

%% Setting up Matrices for the Quadratic Programming Solver
% Cost function to minimize
p_cost=gencostdata(:,6)';
q_cost=0.0*gencostdata(:,6)'; %% for now 10% of active power costs

%% Iteration for the loss formulation starts here

while iter~=maxiter

% Assigning the delivery factors to the generation units
Aeq=ones(1,ng);
for i=1:nb
   for j=1:ng
        if gendata(j,1)==i
            Aeq(j)=DF_est(i);
        end
    end
end

% Demand with Delivery Factor
p_demand=ones(1,nb)*(busdata(:,3).*DF_est);
q_demand=ones(1,nb)*(busdata(:,4).*DF_est);
% Demand + Losses
p_demand=p_demand-Ploss_est;
q_demand=q_demand-Qloss_est;

%% Setting up Matrices for GSF Formulation
A1=GSF_PP*Ag;
A2=GSF_PQ*Ag;
pd=busdata(:,3);
qd=busdata(:,4);

b1=prat(1:length(branchdata(:,1)))+GSF_PP*(pd+E_est)+GSF_PQ*(qd);
b2=prat(1:length(branchdata(:,1)))-GSF_PP*(pd+E_est)-GSF_PQ*(qd);

A_P=[A1; -A1];
A_Q=[A2; -A2];
b=[b1;b2];


%% Generation limit for each generator
p_lb=gendata(:,10);
p_ub=gendata(:,9);
q_lb=gendata(:,5);
q_ub=gendata(:,4);


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
model.A = [sparse(A_P) sparse(A_Q);sparse(Aeq) sparse(zeros(1,size(Aeq,2))); sparse(zeros(1,size(Aeq,2))) sparse(Aeq)]; % A must be sparse
model.sense = [repmat('<',size(A_P,1),1); repmat('=',size(Aeq,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = full([b(:);p_demand(:);q_demand(:)]); % rhs must be dense
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
lambda.eqlin = results.pi((size(A_P,1)+1):end);


%% Generation Cost, Congestion Cost and LMP for each bus
generationcost=lambda.eqlin;
p_congestioncost=(lambda.ineqlin'*[-GSF_PP ; GSF_PP])';
q_congestioncost=(lambda.ineqlin'*[-GSF_PQ ; GSF_PQ])';
p_lmp=generationcost(1)+p_congestioncost;
q_lmp=generationcost(2)+q_congestioncost;

%% Saving the generation dispatch results at each iteration
for i=iter:iter
gendispatch(:,i)=results.x(1:ng);
end
if iter~=1
mismatchdispatch=(gendispatch(:,iter-1)-gendispatch(:,iter)).^2;
end

%% Net injections
pn=Ag*results.x(1:ng)-pd-E_est;
qn=Ag*results.x(ng+1:end)-qd;


%V1=Xmat(6,1:5)*pn./baseMVA+Xmat(6,6:end)*qn./baseMVA;
%V2=Xmat(7,1:5)*pn./baseMVA+Xmat(7,6:end)*qn./baseMVA;
%V3=Xmat(36,1:33)*pn;

%% Damping parameter (necessary for >100 busses)
w=0.0;
if iter>1
pn_old=Ag*gendispatch(:,iter-1)-pd-E_est_old;
pn=w*pn_old+(1-w)*pn;
end

%% Lineflow based on GSF_PP and net injections
P_lineflow=GSF_PP*pn+GSF_PQ*qn;
Q_lineflow=GSF_QQ*qn+GSF_QP*pn;

if iter>1
lineflow_old=GSF_PP*pn_old;
P_lineflow=w*lineflow_old+(1-w)*P_lineflow;
end

%% Loss formulation
run(fullfile('lossfactor_Q.m'));

%% Checking at which iteration the loss formulation converges
checkiter=iter;

%% Adding loss costs to LMPs
if iter~=1
losscost=lambda.eqlin*(DF_est-1);
p_lmp=p_lmp+losscost;
end

%% Criterion to stop further iterations
if any(abs(mismatchdispatch) > 1e-8)
    iter=iter+1;
else
    iter=maxiter;
end

end


%% Printing out the solutions

%if iter==maxiter || iter==2
if iter>maxiter
%Printing out the Solution
for v=1:length(names)
fprintf('\n %s = %2.2f MW\n', names{v}, results.x(v));
end

fprintf('\n Generation cost = %g $/MWh\n\n', generationcost);

busnumber=busdata(:,1);

for v=1:length(p_congestioncost)
    if p_congestioncost(v)~=0
    fprintf(' Congestion Cost @ Bus %g = %4.2f $/MWh\n', busnumber(v) , p_congestioncost(v));
    end
end

fprintf('\n\n   LMP -- Bus Number ||    Generation Cost ||  Congestion Cost \n \n')
for v=1:length(p_lmp)
        if p_lmp(v)~=0
         fprintf(' The LMP @ Bus %2.4g is %4.2f $/MWh generation cost and %4.2f $/MWh congestion cost\n', busnumber(v) , generationcost, p_congestioncost(v));
        end
end


fprintf('\n              Line flow Table \n')
fprintf('   From Bus  ||    To Bus   ||   Line Flow \n\n')
for v=1:length(P_lineflow)
        fprintf('%8.0f     ||%8.0f     ||%13.6f MW  \n', fb(v) , tb(v) , P_lineflow(v));
end

fprintf('\n The objective value is %4.2f $/h\n', results.objval);
end

%% Saving the results in .mat files
%if ~exist(fullfile('loadprofile', casedata), 'dir')
%       mkdir('loadprofile', casedata)
%end
%f1=num2str(minute*60);
%f5=num2str((minute+0.25)*60);
%f2='load at timestep - hour';
%f3=num2str(floor(step-0.25));
%f4='- minute -';
%filename=strcat(f2,f3,f4,f1,f4,f5);
%f=fullfile('loadprofile', casedata , filename);
%save(f);

%% Step increase for 24 hour loadflow analysis
%minute=minute+0.25;
%if minute > 0.75
%    minute=0;
%end
%step=step+0.25;
%end
toc