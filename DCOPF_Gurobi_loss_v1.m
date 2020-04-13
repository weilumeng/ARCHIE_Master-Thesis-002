%% Loading Data
clear;
IEEE118bus;

%% Generation Shift Factor/Power Transfer Distribution Matrix
ptdf_matrix 

%% Setting up the Data
maxiter=9;

% Generator Bus Connection Matrix
Ag=zeros(nb,nb);
for i=1:nb
    for k=1:ng
            if gendata(k,1)==i
                Ag(i,k)=1;
            end
    end
end

%Setting up Matrices for the Quadratic Programming Solver
f=zeros(1,nb);
for i=1:nb
    for j=1:length(gendata(:,6))
        if i==gendata(j,1)
            f(i)=gencostdata(j,6);  %Cost function to minimize
        end
    end
end

%Estimated Losses
Ploss_est=0;
LF_est=zeros(nb,1);
DF_est=ones(nb,1);
E_est=zeros(nb,1);

%% Iteration starts here

for iter=1:maxiter

% Generation Side of Pg=Pd Equation with Delivery Factor
Aeq=zeros(1,nb);
for i=1:nb
    for j=1:length(gendata(:,1))
        if i==gendata(j,1)
            Aeq(i)=1;
        end
    end
end

Aeq=DF_est'.*Aeq;

% Demand Side of Pg=Pd Equation with Delivery Factor and Active Power
% Losses
pd=busdata(:,3);
demand=ones(1,nb)*(DF_est.*pd);             
demand=demand-Ploss_est;                       

%Setting up Matrices for PTDF Formulation
b1=prat(1:length(branchdata(:,1)))+PTDF*(pd);
b2=prat(1:length(branchdata(:,1)))-PTDF*(pd);
A=[PTDF; -PTDF];
b=[b1; b2];


% Generation limit for each generator
lb=zeros(1,nb);
ub=zeros(1,nb);
for i=1:nb
    for j=1:length(gendata(:,6))
        if i==gendata(j,1)
            lb(i)=gendata(j,10);
            ub(i)=gendata(j,9); %Cost function to minimize
        end
    end
end


% Quadratic cost functions for each generator
Qmatrix=zeros(size(A,2),size(A,2));
for i=1:nb
        for k=1:length(gencostdata(:,5))
            if gencostdata(k,1)~=0
                if gendata(k,1)==i
                Qmatrix(i,i)=gencostdata(k,5);
                end
            end
        end
end

% Quadratic Programming with Gurobi Solver
names=num2cell(zeros(1,nb));
for i=1:nb
names(i) = {num2str("P"+num2str(i))};
end
model.varnames = names;
model.obj = f; 
model.Q = sparse(Qmatrix);
model.A = [sparse(A); sparse(Aeq)]; % A must be sparse
model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = full([b(:); demand(:)]); % rhs must be dense
if ~isempty(lb)
    model.lb = lb;
else
    model.lb = -inf(size(model.A,2),1); % default lb for MATLAB is -inf
end
if ~isempty(ub)
    model.ub = ub;
end

gurobi_write(model, 'DCOPF_QP.lp');

results = gurobi(model);
lambda.lower = max(0,results.rc);
lambda.upper = -min(0,results.rc);
lambda.ineqlin = -results.pi(1:size(A,1));
lambda.eqlin = results.pi((size(A,1)+1):end);



%Generation Cost, Congestion Cost and Powerflow
generationcost=lambda.eqlin;
congestioncost=(lambda.ineqlin'*[-PTDF; PTDF])';
lmp=generationcost+congestioncost;

% Net injection/withdrawl
pn=results.x-pd;

lineflow=PTDF*pn;

%Loss Factor Script
lossfactor;

losscost=lambda.eqlin*(DF-1);

end


%Printing out the Solution
if iter==maxiter
for v=1:length(names)
fprintf('\n %s = %2.2f MW\n', names{v}, results.x(v));
end

for v=1:length(generationcost)
fprintf('\n Generation cost = %g $/MWh\n\n', generationcost);
end

busnumber=busdata(:,1);

for v=1:length(congestioncost)
    if congestioncost(v)~=0
    fprintf(' Congestion Cost @ Bus %g = %4.2f $/MWh\n', busnumber(v) , congestioncost(v));
    end
end

fprintf('\n\n   LMP -- Bus Number ||    Generation Cost ||  Congestion Cost \n \n')
for v=1:length(lmp)
        if lmp(v)~=0
         fprintf(' The LMP @ Bus %2.4g is %4.2f $/MWh generation cost and %4.2f $/MWh congestion cost\n', busnumber(v) , generationcost, congestioncost(v));
        end
end


fprintf('\n              Line flow Table \n')
fprintf('   From Bus  ||    To Bus   ||   Line Flow \n\n')
for v=1:length(lineflow)
        fprintf('%8.0f     ||%8.0f     ||%13.6f MW  \n', fb(v) , tb(v) , lineflow(v));
end

fprintf('\n The objective value is %4.2f $/h\n', results.objval);
end