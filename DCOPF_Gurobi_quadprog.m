clear;
IEEE33bus_modified;
ptdf_matrix %Generation Shift Factor/Power Transfer Distribution Matrix
        
Ag=zeros(nb,ng);

for i=1:nb
    for j=1:ng
        if gendata(j,1)==i
            Ag(i,j)=1;
        end
    end
end


f=gencostdata(:,6)';
Aeq=-ones(1,length(gendata(:,1)));
beq=-sum(busdata(:,3));

A1=PTDF*Ag;
       
pd=busdata(:,3);
        
b1=prat(1:length(branchdata(:,1)))+PTDF*pd;
b2=prat(1:length(branchdata(:,1)))-PTDF*pd;
A=[A1; -A1];
b=[b1;b2];

lb=gendata(:,10);
ub=gendata(:,9);

Qmatrix=zeros(5,5);

%Quadratic Programming with Gurobi

names = {'P', 'P1', 'P2', 'P3', 'P4'};
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


for v=1:length(names)
fprintf('\n %s = %g\n', names{v}, results.x(v));
end

for v=1:length(genergy)
fprintf('\n Generation cost is %g\n\n', genergy);
end

busnames=busdata(:,1);

for v=1:length(conjcost)
    if conjcost(v)~=0
    fprintf('Congestion Cost @ Bus %g is %4.2f $/h\n', busnames(v) , conjcost(v));
    end
end

fprintf('LMP -- Bus Number || Generation Cost || Congestion Cost \n \n')
for v=1:length(lmp)
        if lmp(v)~=0
         fprintf('The LMP @ Bus %g is %4.2f $/MWh generation cost and %4.2f $/MWh congestion cost\n', busnames(v) , genergy, conjcost(v));
    end
end

fprintf('\nThe objective value is %4.2f $/MWh\n', results.objval);

fprintf('Line flow Table \n')
fprintf('   From Bus  ||    To Bus   ||   Line Flow \n\n')
for v=1:length(lineflow)
        fprintf('%8.0f     ||%8.0f     ||%10.5f MW  \n', fb(v) , tb(v) , lineflow(v));
end
