%% Case Data
clc;clear;
mpc = (loadcase('case118'));


%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('ptdf_matrix_loss.m'))

%% Iteration Setup

iter=1;
maxiter=100;

%% Generator Bus Connection Matrix

Ag=zeros(nb,ng);
for i=1:nb
    for j=1:ng
        if gen(j,1)==i
            Ag(i,j)=1;
        end
    end
end


%% Estimated Losses
Ploss_est=0;

DF_est=ones(nb,1);
E_est=zeros(nb,1);

gendispatch=zeros(ng,maxiter);
mismatchdispatch=ones(ng,maxiter);

%% Setting up Matrices for the Quadratic Programming Solver
% Cost function to minimize
f=gencost(:,6)';

%% Iteration for the loss formulation starts here

while iter~=maxiter

% Assigning the delivery factors to the generation units
Aeq=ones(1,ng);
for i=1:nb
   for j=1:ng
        if gen(j,1)==i
            Aeq(j)=DF_est(i);
        end
    end
end

% Demand Side of Pg=Pd Equation with Delivery Factor
beq=ones(1,nb)*(bus(:,3).*DF_est);
%Demand Side of Pg=Pd Equation + Losses
beq=beq-Ploss_est;                                      


%% Setting up Matrices for PTDF Formulation
A1=PTDF*Ag;    
pd=bus(:,3);

b1=prat(1:length(branch(:,1)))+PTDF*(pd+E_est);
b2=prat(1:length(branch(:,1)))-PTDF*(pd+E_est);

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
for i=1:length(gen(:,1))
names(i) = {num2str("P"+num2str(i))};
end
model.varnames = names;
model.obj = f; 
model.Q = sparse(Qmatrix);
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
else
    model.ub= inf(size(model.A,2),1);
end

gurobi_write(model, 'DCOPF_loss.lp');

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

%% Saving the generation dispatch results at each iteration
for i=iter:iter
gendispatch(:,i)=results.x;
end
if iter~=1
mismatchdispatch=(gendispatch(:,iter-1)-gendispatch(:,iter)).^2;
end

%% Net injections
pn=Ag*results.x-pd-E_est;

%% Damping parameter (necessary for >100 busses)
if nb>100
w=0.75;
else
    w=0;
end

%% Lineflow based on PTDF and net injections
lineflow=PTDF*pn;


%% Loss formulation
run(fullfile('lossfactor.m'));

%% Checking at which iteration the loss formulation converges
checkiter=iter;

%% Adding loss costs to LMPs
if iter~=1
losscost=lambda.eqlin*(DF_est-1);
lmp=lmp+losscost;
end

%% Criterion to stop further iterations
if any(abs(mismatchdispatch) > 1e-8)
    iter=iter+1;
else
    iter=maxiter;
end

end


%% Printing out the solutions

if iter==maxiter || iter==2
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
end

%% MATPOWER Reference Values (ACOPF AND DCOPF)
define_constants;
mpopt = mpoption('model','ACOPF');
resultAC = runopf(mpc,mpopt);
BusVolAC = resultAC.bus(:,VM);
BusAglAC = resultAC.bus(:,VA);
BranchFlowAC = (resultAC.branch(:,PF)-resultAC.branch(:,PT))/2;
BranchFlowAC_Q= (resultAC.branch(:,QF)-resultAC.branch(:,QT))/2;
ALMPAC=(resultAC.bus(:,LAM_P));
RLMPAC=(resultAC.bus(:,LAM_Q));

Fehler=(abs((lmp-ALMPAC)./ALMPAC))*100;
AEA=sum(abs((lmp-ALMPAC)./ALMPAC))/nb;

mpopt = mpoption('model','DCOPF');
resultDC = rundcopf(mpc,mpopt);
BranchFlowDC = (resultDC.branch(:,PF)-resultDC.branch(:,PT))/2;
ALMPDC=(resultDC.bus(:,LAM_P));
AEADC=sum(abs((ALMPDC-ALMPAC)./ALMPAC))/nb;
FehlerDC=(abs((ALMPDC-ALMPAC)./ALMPAC))*100;



%% PLOTTING
fh1=figure('Name','Leistungsflüsse','Menu','none','ToolBar','none', 'Position', [   488   342   560   420]);
plot(1:nl,BranchFlowAC,'Marker','*','LineWidth',1,'Color','#A2142F', 'MarkerSize', 4,'LineStyle', '-');hold on;
plot(1:nl,lineflow,'Marker','o','LineWidth',1, 'Color', '#EDB120','MarkerSize', 4,'LineStyle', '-'); hold on;
plot(1:nl,BranchFlowDC,'Marker','.','LineWidth',0.8, 'Color', '#0072BD','MarkerSize', 6,'LineStyle', '-')
xlim([1 nl+1]);
%ylim([-0.5 4]);
legend('MATPOWER ACOPF','DCOPF-L', 'MATPOWER DCOPF');
xlabel('Leitung','FontSize',10,'FontName','Cambria');
ylabel('Wirkleistungsfluss (MW)','FontSize',10,'FontName','Cambria');
%title('Vergleich der Wirkleistungsflüsse','FontSize',10,'FontName','Cambria');
set(gca,'FontSize',10, 'FontName','Cambria')
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('Leistungsflüsse Verlust' ,'-dsvg')


energycost=(1:nb)';
energycost(1:nb)=generationcost;

fh2=figure('Name','LMP Zerlegung','Menu','none','ToolBar','none');
set(gca,'FontSize',10, 'FontName','Cambria')
t=tiledlayout(3,1);
nexttile
plot(1:nb,energycost,'Marker','.','LineWidth',0.8,'Color','#1b9e77','MarkerSize', 9,'LineStyle', '-');hold on;
legend('Energiekosten','FontSize',10, 'FontName','Cambria');
xlim([1 nb+1]);
%ylim([19 21]);
nexttile([2 1])
plot(1:nb,congestioncost,'Marker','.','LineWidth',0.8, 'Color', '#d95f02','MarkerSize', 9,'LineStyle', '-'); hold on;
plot(1:nb,losscost, 'Marker','.','LineWidth',0.8, 'Color', '#7570b3','MarkerSize', 9,'LineStyle', '-')
%ylim([0 2.5]);
xlim([1 nb+1]);
legend('Engpasskosten' , 'Verlustkosten','FontSize',10, 'FontName','Cambria');
xlabel(t,'Knoten','FontSize',10,'FontName','Cambria');
ylabel(t,'Wirkleistungskosten ($/MWh)','FontSize',10,'FontName','Cambria');
%title(t,'LMP Zerlegung des vorgeschlagenen DCOPF mit Verlusten','FontSize',10,'FontName','Cambria', 'Fontweight', 'bold');
t.TileSpacing = 'none';
t.Padding = 'none';
print('LMP Zerlegung' ,'-dsvg')


colorbars='#EDB120';
colorbarsDC='#0072BD';

fh3= figure('Name','ALMP Dif', 'Menu','none','ToolBar','none');
set(gca,'FontSize',10, 'FontName','Cambria')
t=tiledlayout(3,1);
nexttile([2 1]);
plot(1:nb,ALMPAC,'Marker','*','LineWidth',0.8,'Color','#A2142F','MarkerSize', 5,'LineStyle', '-');hold on;
plot(1:nb,lmp,'Marker','o','LineWidth',0.8, 'Color', '#EDB120','MarkerSize', 5,'LineStyle', '-'); hold on;
plot(1:nb,ALMPDC, 'Marker','.','LineWidth',0.8, 'Color', '#0072BD','MarkerSize', 6,'LineStyle', '-'); hold on;
ylabel('Wirkleistungskosten ($/MWh)','FontSize',10,'FontName','Cambria');
legend('MATPOWER ACOPF','DCOPF-L', 'MATPOWER DCOPF');
%ylim([18.75 23]);
xlim([1 nb+1]);
nexttile
bars=bar(1:nb,[Fehler,FehlerDC]); hold on;
bars(1).FaceColor=colorbars;
bars(2).FaceColor=colorbarsDC;
bars(1).EdgeColor=colorbars;
bars(2).EdgeColor=colorbarsDC;
%ylim([0 25]);
xlim([1 nb+1]);
ylabel('rel. Fehler (%)','FontSize',10,'FontName','Cambria');
%ytickformat('percentage')
legend('rel. Fehler - DCOPF-L','rel. Fehler - MATPOWER DCOPF');
xlabel(t,'Knoten','FontSize',10,'FontName','Cambria');
%title('Vergleich der ALMPs','FontSize',10,'FontName','Cambria');
t.TileSpacing = 'none';
t.Padding = 'none';
print('ALMP Verlust' ,'-dsvg')