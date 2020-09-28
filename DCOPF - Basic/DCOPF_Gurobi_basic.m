%% Case Data
clc;clear;
mpc = (loadcase('case118_constrained'));

%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('ptdf_matrix_basic.m'))


%% Generator Bus Connection Matrix

Ag=zeros(nb,ng);
for i=1:nb
    for j=1:ng
        if gen(j,1)==i
            Ag(i,j)=1;
        end
    end
end

%% Loadcurve timesteps

%while step<24
%run(fullfile(casedata))
%busdata(:,3)=busdata(:,3)*loadcurve_base(step*4);
   


%% Setting up Matrices for the Quadratic Programming Solver
% Cost function to minimize
f=gencost(:,6)';


% Assigning the delivery factors to the generation units
Aeq=ones(1,ng);


% Demand Side of Pg=Pd Equation
beq=ones(1,nb)*bus(:,3);                                  


%% Setting up Matrices for PTDF Formulation
A1=PTDF*Ag;    
pd=bus(:,3);

b1=prat(1:nl)+PTDF*(pd);
b2=prat(1:nl)-PTDF*(pd);

A=[A1; -A1];
b=[b1;b2];


%% Generation limit for each generator
lb=gen(:,10);
ub=gen(:,9);


%% Quadratic cost functions for each generator
Qmatrix=zeros(size(A,2),size(A,2));
for i=1:length(gencost(:,5))
    if gencost(i,5)~=0
        Qmatrix(i,i)=gencost(i,5);
    end
end

%% Quadratic Programming with Gurobi Solver
names={num2str(zeros(1,length(gen(:,1))))};
for i=1:ng
names(i) = {num2str("P"+num2str(i))};
end
model.varnames = names;
model.obj = f; 
model.Q = sparse(Qmatrix);
model.A = [sparse(A); sparse(Aeq)]; % A must be sparse
model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = full([b(:); beq(:)]); % rhs must be dense

gurobi_write(model, 'DCOPF_basic.lp');

results = gurobi(model);

%% Lagrange multipliers

lambda.lower = max(0,results.rc);
lambda.upper = -min(0,results.rc);
lambda.ineqlin = -results.pi(1:size(A,1));
lambda.eqlin = results.pi((size(A,1)+1):end);


%% Generation Cost, Congestion Cost and LMP for each bus
generationcost=lambda.eqlin;
congestioncost=(lambda.ineqlin'*[-PTDF ; PTDF])';
lmp=generationcost+congestioncost;


%% Net injections
pn=Ag*results.x-pd;

%% Lineflow based on PTDF and net injections
lineflow=PTDF*pn;



%% Printing out the solutions

%Printing out the Solution
for v=1:length(names)
fprintf('\n %s = %2.2f MW\n', names{v}, results.x(v));
end

for v=1:length(generationcost)
fprintf('\n Generation cost = %g $/MWh\n\n', generationcost);
end

busnumber=bus(:,1);

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


%% ACOPF MATPOWER
define_constants;
mpopt = mpoption('model','ACOPF');
resultAC = runopf(mpc,mpopt);
BusVolAC = resultAC.bus(:,VM);
BusAglAC = resultAC.bus(:,VA);
BranchFlowAC = (resultAC.branch(:,PF)-resultAC.branch(:,PT))/2;
BranchFlowAC_Q= (resultAC.branch(:,QF)-resultAC.branch(:,QT))/2;
ALMPAC=(resultAC.bus(:,LAM_P));
RLMPAC=(resultAC.bus(:,LAM_Q));

AEA=sum(abs((lmp-ALMPAC)./ALMPAC))/nb;

mpopt = mpoption('model','DCOPF');
resultDC = rundcopf(mpc,mpopt);
ALMPDC=(resultDC.bus(:,LAM_P));
AEADC=sum(abs((ALMPDC-ALMPAC)./ALMPAC))/nb;

%% PLOTTING
figure
plot(1:nl,BranchFlowAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nl,lineflow,'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','Generischer DCOPF');
xlabel('Leitung');
ylabel('Wirkleistungsfluss (MW)');
title('Wirkleistungsfluss Vergleich');


figure
plot(1:nb,ALMPAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nb,lmp,'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','Generischer DCOPF');
xlabel('Knoten');
ylabel('Wirkleistungskosten ($/MWh)');
title('ALMP Vergleich');
