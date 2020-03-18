clear;
IEEE118bus;
ptdf_matrix %Generation Shift Factor/Power Transfer Distribution Matrix


%Generator Bus Connection Matrix
Ag=zeros(nb,ng);

for i=1:nb
    for j=1:ng
        if gendata(j,1)==i
            Ag(i,j)=1;
        end
    end
end

%Setting up Matrices for the Quadratic Programming Solver
f=gencostdata(:,6)';    %Cost function to minimize
Aeq=-ones(1,length(gendata(:,1)));  %Generation Side of Pg=Pd Equation
beq=-sum(busdata(:,3));             %Demand Side of Pg=Pd Equation

%Setting up Matrices for PTDF Formulation
A1=PTDF*Ag;    
pd=busdata(:,3);      
b1=prat(1:length(branchdata(:,1)))+PTDF*pd;
b2=prat(1:length(branchdata(:,1)))-PTDF*pd;
A=[A1; -A1];
b=[b1;b2];

%Generation Limit for each generator
lb=gendata(:,10);
ub=gendata(:,9);

%Quadratic cost functions for each generator
Qmatrix=zeros(size(A,2),size(A,2));
for i=1:length(gencostdata(:,5))
    if gencostdata(i,5)~=0
        Qmatrix(i,i)=gencostdata(i,5);
    end
end

%Quadratic Programming with Gurobi Solver

for i=1:length(gendata(:,1))
names(i) = {num2str("P"+num2str(i))};
end
model.varnames = names;
model.obj = f; 
model.Q = [sparse(Qmatrix)];
model.A = [sparse(A); sparse(Aeq)]; % A must be sparse
model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = full([b(:); beq(:)]); % rhs must be dense
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
lambda.eqlin = -results.pi((size(A,1)+1):end);

genergy=lambda.eqlin;
conjcost=(lambda.ineqlin'*[-PTDF ; PTDF])';
lmp=genergy+conjcost;
lineflow=PTDF*(Ag*results.x-pd);

%Printing out the Solution

for v=1:length(names)
fprintf('\n %s = %2.2f\n', names{v}, results.x(v));
end

for v=1:length(genergy)
fprintf('\n Generation cost = %g $/MWh\n\n', genergy);
end

busnumber=busdata(:,1);

for v=1:length(conjcost)
    if conjcost(v)~=0
    fprintf(' Congestion Cost @ Bus %g = %4.2f $/MWh\n', busnumber(v) , conjcost(v));
    end
end

fprintf(' LMP -- Bus Number || Generation Cost || Congestion Cost \n \n')
for v=1:length(lmp)
        if lmp(v)~=0
         fprintf(' The LMP @ Bus %2.2g is %4.2f $/MWh generation cost and %4.2f $/MWh congestion cost\n', busnumber(v) , genergy, conjcost(v));
    end
end


fprintf('\n               Line flow Table \n')
fprintf('   From Bus  ||    To Bus   ||   Line Flow \n\n')
for v=1:length(lineflow)
        fprintf('%8.0f     ||%8.0f     ||%10.5f MW  \n', fb(v) , tb(v) , lineflow(v));
end

fprintf('\n The objective value is %4.2f $/h\n', results.objval);
