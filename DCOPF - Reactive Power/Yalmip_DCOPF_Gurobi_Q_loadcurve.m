%% This algorithm is used for the 24 hour analysis

%% Case Data
clc;clear;
define_constants;
casedata='case69m';
mpc = (loadcase(casedata));

% for cases without branch limits (such as case33bw)
mpc.branch(:,RATE_A)=9999;


%Voltage constraints
mpc.bus(:,VMAX)=1.03;
mpc.bus(:,VMIN)=0.97;

%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('Yalmip_ptdf_matrix_Q.m'))


%% Loadcurve setup based on 10. February 2020

loadprofiledata=readtable('202002100000.xlsx', 'PreserveVariableNames', true, 'Range', 'C8:C104', 'ReadVariableNames', true);
loadcurve=table2array(loadprofiledata(:,1));
loadcurve=str2double(loadcurve);
loadcurve_base=loadcurve./mean(loadcurve);

%% Iteration Setup for the 24 hour loadcurve
step=0;             %needed for the loadcurve
minute=0;           %needed for the loadcurve
iter=0;

gendispatch=zeros(2*ng,100);                    %Generator dispatch for each iteration
mismatchdispatch=zeros(ng,100);

while step<24
mpc = (loadcase(casedata));

% for cases without branch limits (such as case33bw)
mpc.branch(:,RATE_A)=9999;

%Voltage constraints
mpc.bus(:,VMAX)=1.03;
mpc.bus(:,VMIN)=0.97;

mpc.bus(:,PD)=mpc.bus(:,PD)*loadcurve_base(1+step*4);
mpc.bus(:,QD)=mpc.bus(:,QD)*loadcurve_base(1+step*4);

[baseMVA, bus, gen, branch, gencost] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
iter=0;
maxiter=15;
E_P_est_old=0.01;
while (sum(abs(E_P_est-E_P_est_old))+sum(abs(E_Q_est-E_Q_est_old))>0.00000001 && iter <= maxiter  )
%% Model construction with yalmip toolbox
%variables
P=sdpvar(ng,1);
Q=sdpvar(ng,1);
Pl=sdpvar(nl,1);
Ql=sdpvar(nl,1);
V=sdpvar(nb,1);
%constraints
constraint=[];
%% Power balance
constraint=[constraint,(sum(DF_P_est(gen(:,GEN_BUS)).*P)==sum(DF_P_est.*bus(:,PD))./baseMVA+sum(E_P_est)):'Pbala'];
%constraint=[constraint,(sum(DF_Q_est(gen(:,GEN_BUS)).*Q)==sum(DF_Q_est.*bus(:,QD))./baseMVA-sum(sum(Bbus))+sum(E_Q_est)):'Qbala'];
% The constraint of reactive power balance should be removed, since it is
% actually reflected by the voltage constraint below.

%% Power flow limits
constraint=[constraint,(GSF_PP(:,gen(:,GEN_BUS))*P-GSF_PP(:,:)*(bus(:,PD)./baseMVA+E_P_est) + GSF_PQ(:,gen(:,GEN_BUS))*Q-GSF_PQ(:,:)*(bus(:,QD)./baseMVA+E_Q_est)==Pl):'Pl'];
constraint=[constraint,(GSF_QP(:,gen(:,GEN_BUS))*P-GSF_QP(:,:)*(bus(:,PD)./baseMVA+E_P_est) + GSF_QQ(:,gen(:,GEN_BUS))*Q-GSF_QQ(:,:)*(bus(:,QD)./baseMVA+E_Q_est)==Ql):'Ql'];
constraint=[constraint,((-branch(:,RATE_A)/baseMVA)<=Pl(:)):'Plmin'];
constraint=[constraint,(Pl(:)<=(branch(:,RATE_A)/baseMVA)):'Plmax'];
%% Node voltage
constraint=[constraint,(X(nb+1:2*nb,gen(:,GEN_BUS))*P-X(nb+1:2*nb,1:nb)*((bus(:,PD)./baseMVA)+E_P_est) + X(nb+1:2*nb,nb+gen(:,GEN_BUS))*Q-X(nb+1:2*nb,nb+1:2*nb)*((bus(:,QD)./baseMVA)+E_Q_est)==V):'V'];
constraint=[constraint,(bus(1:nb,VMIN)<=V(1:nb)):'Vmin'];
constraint=[constraint,(V(1:nb)<=bus(1:nb,VMAX)):'Vmax'];

ang=X(1:nb,gen(:,GEN_BUS))*P-X(1:nb,1:nb)*((bus(:,PD)./baseMVA)-E_P_est) + X(1:nb,nb+gen(:,GEN_BUS))*Q-X(1:nb,nb+1:2*nb)*(bus(:,QD)./baseMVA+E_Q_est);
%% Generation output limits
constraint=[constraint,((gen(:,PMIN)./baseMVA)<=P<=(gen(:,PMAX)./baseMVA)):'Plim'];
constraint=[constraint,((gen(:,QMIN)./baseMVA)<=Q<=(gen(:,QMAX)./baseMVA)):'Qlim'];

%% Objective function
if(length(mpc.gencost(:,1))==2*ng) % both active and reactive power prices
    if(mpc.gencost(1,MODEL)==POLYNOMIAL)
        if(mpc.gencost(1,NCOST)==2)
            objective=baseMVA*(P'*mpc.gencost(1:ng,COST)+sum(mpc.gencost(1:ng,COST+1))+ ...
                Q'*mpc.gencost(ng+1:2*ng,COST)+sum(mpc.gencost(ng+1:2*ng,COST+1)));
        elseif(mpc.gencost(1,NCOST)==3)
            objective=((P'*baseMVA).^2)*(mpc.gencost(1:ng,COST))+((Q'*baseMVA).^2)*(mpc.gencost(ng+1:2*ng,COST))+ ...
                baseMVA*(P'*mpc.gencost(1:ng,COST+1)+sum(mpc.gencost(1:ng,COST+2))+ ...
            Q'*mpc.gencost(ng+1:2*ng,COST+1)+sum(mpc.gencost(ng+1:2*ng,COST+2)));
        else
            error('only linear or second-order polynomial accepted')
        end
    elseif(mpc.gencost(1,MODEL)==PW_LINEAR)
        error('piecewise linearization not considered')
    end
elseif(length(mpc.gencost(:,1))==ng) % only active power prices
    if(mpc.gencost(1,MODEL)==POLYNOMIAL)
        if(mpc.gencost(1,NCOST)==2)
            objective=baseMVA*(P'*mpc.gencost(1:ng,COST)+sum(mpc.gencost(1:ng,COST+1)));
        elseif(mpc.gencost(1,NCOST)==3)
            objective=((P'*baseMVA).^2)*(mpc.gencost(1:ng,COST))+ baseMVA*(P'*mpc.gencost(1:ng,COST+1)+sum(mpc.gencost(1:ng,COST+2)));
        else
            error('only linear or second-order polynomial accepted')
        end
    elseif(mpc.gencost(1,MODEL)==PW_LINEAR)
        error('piecewise linearization not considered')
    end
else
    error('price data error')
end
%% Solve by gurobi
options=sdpsettings('solver','gurobi');
solution = optimize(constraint,objective,options);

%% Return values
busAng=value(ang)-sum(value(ang))/nb;busVol=value(V);
genP=value(P);genQ=value(Q);branchP=(value(Pl)*baseMVA);branchQ=(value(Ql)*baseMVA);
objectiveValue=value(objective);
          
gendispatch(:,iter+1)=[genP;genQ];
if iter~=0
mismatchdispatch=(gendispatch(1:ng,iter)-gendispatch(1:ng,iter+1)).^2+(gendispatch(ng+1:end,iter)-gendispatch(ng+1:end,iter+1)).^2;
end


run(fullfile('Yalmip_lossfactor_Q.m'));
iter=iter+1;

end
%% LMP calculation
    EnegyP=-dual(constraint('Pbala'))./baseMVA;
  % EnegyQ=-dual(constraint('Qbala'))./baseMVA;
    CongP=-dual(constraint('Pl'))./baseMVA;
    CongQ=-dual(constraint('Ql'))./baseMVA;
    LimV=-dual(constraint('V'))./baseMVA;
    %LMP decomposition
    PLMPE=zeros(nb,1);
    PLMPC=zeros(nb,1);
    PLMPV=zeros(nb,1);
    QLMPE=zeros(nb,1);
    QLMPC=zeros(nb,1);
    QLMPV=zeros(nb,1);
    for i=1:nb
        PLMPE(i)=EnegyP*DF_P_est(i);
        PLMPC(i)=CongP'*GSF_PP(:,i)+CongQ'*GSF_QP(:,i);
        PLMPV(i)=LimV'*X(nb+1:2*nb,i);
       %  QLMPE(i)=EnegyQ*DF_Q_est(i);
        %QLMPE(i)=0;
        QLMPC(i)=CongP'*GSF_PQ(:,i)+CongQ'*GSF_QQ(:,i);
        QLMPV(i)=LimV'*X(nb+1:2*nb,nb+i);
    end
    LMP_P=PLMPE+PLMPC+PLMPV;
    LMP_Q=QLMPE+QLMPC+QLMPV;
    


%% Clearing Yalmip decision variables and constraints - they cannot be saved in binary format    
constraint=[];objective=[];ang=[];V=[];P=[];Q=[];Pl=[];Ql=[]; 


%% MATPOWER Reference Values (ACOPF AND DCOPF)
% AC Optimal Power flow
mpopt = mpoption('model','ACOPF');

resultAC = runopf(mpc,mpopt);
BusVolAC = resultAC.bus(:,VM);

BranchFlowAC = ((resultAC.branch(:,PF)-resultAC.branch(:,PT))/2);
BranchFlowAC_Q= ((resultAC.branch(:,QF)-resultAC.branch(:,QT))/2);

ALMPAC=(resultAC.bus(:,LAM_P));
RLMPAC=(resultAC.bus(:,LAM_Q));

% Error Plots DCOPF-QV and MATPOWER ACOPF

Fehler_ALMP=(abs((LMP_P-ALMPAC)./ALMPAC))*100;
AEA_ALMP=sum(abs((LMP_P-ALMPAC)./ALMPAC))/nb;

Fehler_RLMP=(abs((LMP_Q-RLMPAC)./RLMPAC))*100;
AEA_RLMP=sum(abs((LMP_Q-RLMPAC)./RLMPAC))/nb;

% DC Optimal Power flow
mpopt = mpoption('model','DCOPF');
resultDC = rundcopf(mpc,mpopt);
BranchFlowDC = (resultDC.branch(:,PF)-resultDC.branch(:,PT))/2;
ALMPDC=(resultDC.bus(:,LAM_P));
AEADC=sum(abs((ALMPDC-ALMPAC)./ALMPAC))/nb;
FehlerDC=(abs((ALMPDC-ALMPAC)./ALMPAC))*100;
BusVolDC = resultDC.bus(:,VM);

%% Saving the results in .mat files
day='_10_february_2020';
if ~exist(fullfile('loadprofile', casedata), 'dir')
       mkdir('loadprofile', casedata)
end
f7='very_tight_voltage';
f1='_loadcurve after_';
f2=num2str(floor(step));
f3='_hours from minute_';
f4=num2str(minute*60);
f5='_up to minute_';
f6=num2str((minute+0.25)*60);
filename=strcat(casedata,f7,f1,f2,f3,f4,f5,f6,day);
f=fullfile('loadprofile', casedata , filename);
save(f);

%% Step increase for 24 hour loadflow analysis
minute=minute+0.25;
if minute > 0.75
    minute=0;
end
step=step+0.25;
end

%% Plots are handled by the file loading_data_from_loadcurve_files