%% Case Data
clear;
tic
casedata='IEEE33bus_modified.m';
run(fullfile(casedata))

%% Loadcurve

%loadprofiledata=readtable('202002110000.xlsx', 'PreserveVariableNames', true, 'Range', 'C8:C103', 'ReadVariableNames', true);
%loadcurve=table2array(loadprofiledata(:,1));
%loadcurve=str2double(loadcurve);
%loadcurve_base=loadcurve./mean(loadcurve);

%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('ptdf_matrix_basic.m'))

%% Loadcurve Setup

%step=0.25;
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
   


%% Setting up Matrices for the Quadratic Programming Solver
% Cost function to minimize
f=gencostdata(:,6)';


% Assigning the delivery factors to the generation units
Aeq=ones(1,ng);


% Demand Side of Pg=Pd Equation
beq=ones(1,nb)*busdata(:,3);                                  


%% Setting up Matrices for PTDF Formulation
A1=PTDF*Ag;    
pd=busdata(:,3);

b1=prat(1:nl)+PTDF*(pd);
b2=prat(1:nl)-PTDF*(pd);

A=[A1; -A1];
b=[b1;b2];


%% Generation limit for each generator
lb=gendata(:,10);
ub=gendata(:,9);


%% Quadratic cost functions for each generator
Qmatrix=zeros(size(A,2),size(A,2));
for i=1:length(gencostdata(:,5))
    if gencostdata(i,5)~=0
        Qmatrix(i,i)=gencostdata(i,5);
    end
end

%% Quadratic Programming with Gurobi Solver
names={num2str(zeros(1,length(gendata(:,1))))};
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

%% ACOPF MATPOWER
% AC Power flow
define_constants;
mpc = ext2int(loadcase('case33bwm'));
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
%figure
%plot(1:nb,BusVolAC,'Marker','*','LineWidth',0.6,'Color','#A2142F','MarkerFaceColor','#A2142F');hold on;
%plot(1:nb,Voltageatbus,'Marker','o','LineWidth',0.6, 'Color', '#EDB120','MarkerFaceColor','#EDB120'); hold on;
%plot(1:nb,vmax,'Marker','v','LineWidth',1.0,'Color','#0072BD', 'MarkerSize', 4, 'MarkerFaceColor', '#0072BD'); hold on;
%plot(1:nb,vmin,'Marker','^','LineWidth',1.0,'Color','#0072BD', 'MarkerSize', 4, 'MarkerFaceColor', '#0072BD')
%legend('ACOPF','My DCOPF', 'Vmax', 'Vmin');
%xlabel('Node');
%ylabel('Voltage Magnitude');
%title('Voltage Magnitude Comparison');

figure
plot(1,AEA,'Marker','o','LineWidth',0.6,'Color','#A2142F','MarkerFaceColor','#A2142F');hold on;
plot(1,AEADC,'Marker','*','LineWidth',0.6,'Color','#A2142F','MarkerFaceColor','#A2142F');hold on;
legend('Average Error : My DCOPF', 'Matpower DCOPF');
xlabel('Loadlevel');
ylabel('Error Percentage');
title('Average Error Analysis');

figure
plot(1:nl,BranchFlowAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nl,lineflow,'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','My DCOPF');
xlabel('Branch');
ylabel('MW');
title('Active Power Flow Comparison');

%figure
%plot(1:nl,BranchFlowAC_Q,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
%plot(1:nl,lineflow(:,2),'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
%legend('ACOPF','My DCOPF');
%xlabel('Branch');
%ylabel('MVar');
%title('Reactive Power Flow Comparison');

figure
plot(1:nb,ALMPAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
plot(1:nb,lmp,'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
legend('ACOPF','My DCOPF');
xlabel('Bus');
ylabel('$/MWh');
title('ALMP Comparison');

%figure
%plot(1:nb,RLMPAC,'Marker','*','LineWidth',0.5,'Color','#A2142F');hold on;
%plot(1:nb,q_lmp,'Marker','o','LineWidth',0.5, 'Color', '#EDB120');
%legend('ACOPF','My DCOPF');
%xlabel('Bus');
%ylabel('$/MVARh');
%title('RLMP Comparison');