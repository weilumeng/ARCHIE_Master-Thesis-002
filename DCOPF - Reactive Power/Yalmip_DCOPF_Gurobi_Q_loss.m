%% Case Data
clc;clear;
define_constants;
mpc = (loadcase('case118m'));

% for cases without branch limits (such as case33bw)
%mpc.branch(:,RATE_A)=9999;

%Scaling the load level (if needed)
mpc.bus(:,[PD QD])=mpc.bus(:,[PD QD])*1.0;

%Voltage constraints
mpc.bus(:,VMAX)=1.03;
mpc.bus(:,VMIN)=0.97;
  

%% Generation Shift / Power Transfer Distribution Matrix

run(fullfile('Yalmip_ptdf_matrix_Q.m'))

%% Iteration Setup
gendispatch=zeros(2*ng,100);                   
mismatchdispatch=ones(ng,100); 
maxiter=15;
iter=0;
E_P_est_old=0.01;

%% Iteration starts here
while (sum(abs(E_P_est-E_P_est_old))+sum(abs(E_Q_est-E_Q_est_old))>0.00000001 && iter<=maxiter )
%% Model construction with yalmip toolbox
%decission variables
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
% The constraint of reactive power balance is inactive, because its
% reflected by the voltage constraint below.

%% Power flow limits
constraint=[constraint,(GSF_PP(:,gen(:,GEN_BUS))*P-GSF_PP(:,:)*(bus(:,PD)./baseMVA+E_P_est) + GSF_PQ(:,gen(:,GEN_BUS))*Q-GSF_PQ(:,:)*(bus(:,QD)./baseMVA+E_Q_est)==Pl):'Pl'];
constraint=[constraint,(GSF_QP(:,gen(:,GEN_BUS))*P-GSF_QP(:,:)*(bus(:,PD)./baseMVA+E_P_est) + GSF_QQ(:,gen(:,GEN_BUS))*Q-GSF_QQ(:,:)*(bus(:,QD)./baseMVA+E_Q_est)==Ql):'Ql'];
constraint=[constraint,((-branch(:,RATE_A)/baseMVA)<=Pl(:)):'Plmin'];
constraint=[constraint,(Pl(:)<=(branch(:,RATE_A)/baseMVA)):'Plmax'];

%% Node voltage
constraint=[constraint,(X(nb+1:2*nb,gen(:,GEN_BUS))*P-X(nb+1:2*nb,1:nb)*((bus(:,PD)./baseMVA)+E_P_est) + X(nb+1:2*nb,nb+gen(:,GEN_BUS))*Q-X(nb+1:2*nb,nb+1:2*nb)*((bus(:,QD)./baseMVA)+E_Q_est)==V):'V'];
constraint=[constraint,(bus(1:nb,VMIN)<=V(1:nb)):'Vmin'];
constraint=[constraint,(V(1:nb)<=bus(1:nb,VMAX)):'Vmax'];

%% Node voltage angle
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
                baseMVA*(P'*mpc.gencost(1:ng,COST+1)+sum(mpc.gencost(1:ng,COST+2)));
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
options=sdpsettings('solver','gurobi');%,'showprogress',1,'savesolverinput',1,'savesolveroutput',1);%'cplex.Diagnostics','on','cplex.display','on');
solution = optimize(constraint,objective,options);

%% Return values
busAng=value(ang)-sum(value(ang))/nb;busVol=value(V);
genP=value(P);genQ=value(Q);branchP=(value(Pl)*baseMVA);branchQ=(value(Ql)*baseMVA);


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


%% MATPOWER Reference Values (ACOPF AND DCOPF)
% AC Optimal Power flow
mpopt = mpoption('model','ACOPF');

resultAC = runopf(mpc,mpopt);
BusVolAC = resultAC.bus(:,VM);

BranchFlowAC = ((resultAC.branch(:,PF)-resultAC.branch(:,PT))/2);
BranchFlowAC_Q= ((resultAC.branch(:,QF)-resultAC.branch(:,QT))/2);

% Error Plots DCOPF-QV and MATPOWER ACOPF
ALMPAC=(resultAC.bus(:,LAM_P));
RLMPAC=(resultAC.bus(:,LAM_Q));

Fehler_ALMP=(abs((LMP_P-ALMPAC)./ALMPAC))*100;
AEA_ALMP=sum(abs((LMP_P-ALMPAC)./ALMPAC))/nb;

Fehler_RLMP=(abs((LMP_Q-RLMPAC)./RLMPAC))*100;
AEA_RLMP=sum(abs((LMP_Q-RLMPAC)./RLMPAC))/nb;

% DC Optimal Power flow

mpopt = mpoption('model','DCOPF');
resultDC = rundcopf(mpc,mpopt);
BranchFlowDC = (resultDC.branch(:,PF)-resultDC.branch(:,PT))/2;

% Error Plots MATPOWER DCOPF and MATPOWER ACOPF
ALMPDC=(resultDC.bus(:,LAM_P));
AEADC=sum(abs((ALMPDC-ALMPAC)./ALMPAC))/nb;
FehlerDC=(abs((ALMPDC-ALMPAC)./ALMPAC))*100;
BusVolDC = resultDC.bus(:,VM);


%% PLOTTING

fh1=figure('Name','Spannung','Menu','none','ToolBar','none');
set(gca,'FontSize',10, 'FontName','Cambria')
plot(1:nb,BusVolAC,'Marker','*','LineWidth',1,'Color','#A2142F','MarkerFaceColor','none','MarkerSize', 5,'LineStyle', '-');hold on;
plot(1:nb,busVol,'Marker','o','LineWidth',1, 'Color', '#EDB120','MarkerFaceColor','none','MarkerSize', 5,'LineStyle', '-'); hold on;
plot(1:nb,BusVolDC,'Marker','.','LineWidth',0.9, 'Color', '#0072BD','MarkerSize', 6,'LineStyle', '-')
plot(1:nb,mpc.bus(:,VMAX),'Marker','v','LineWidth',1.0,'Color','#5e3c99', 'MarkerSize', 4, 'MarkerFaceColor', '#5e3c99'); hold on;
plot(1:nb,mpc.bus(:,VMIN),'Marker','^','LineWidth',1.0,'Color','#5e3c99', 'MarkerSize', 4, 'MarkerFaceColor', '#5e3c99')
legend('MATPOWER ACOPF','DCOPF-QV', 'MATPOWER DCOPF', 'Spannungsmaximum', 'Spannungsminimum');
ylim([0.96 1.06]);
xlim([1 nb+1]);
xlabel('Knoten','FontSize',10,'FontName','Cambria');
ylabel('Spannung (p.u.)','FontSize',10,'FontName','Cambria');
%title('Vergleich der Spannungsamplitude','FontSize',10,'FontName','Cambria');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width) ax_height];
print('Spannung' ,'-dsvg')



fh2=figure('Name','Wirkleistung','Menu','none','ToolBar','none');
set(gca,'FontSize',10, 'FontName','Cambria')
plot(1:nl,BranchFlowAC,'Marker','*','LineWidth',1,'Color','#A2142F','MarkerSize', 5,'LineStyle', '-');hold on;
plot(1:nl,branchP,'Marker','o','LineWidth',1, 'Color', '#EDB120','MarkerSize', 5,'LineStyle', '-');hold on;
plot(1:nl,BranchFlowDC,'Marker','.','LineWidth',0.2, 'Color', '#0072BD','MarkerSize', 5,'LineStyle', '-')
xlim([1 nl+1])
%ylim([-1.5 2.5])
legend('MATPOWER ACOPF','DCOPF-QV', 'MATPOWER DCOPF');
xlabel('Leitung','FontSize',10,'FontName','Cambria');
ylabel('Wirkleistungsfluss (MW)','FontSize',10,'FontName','Cambria');
title('Vergleich der Wirkleistungsflüsse','FontSize',10,'FontName','Cambria');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('Wirkleistung' ,'-dsvg')



fh3=figure('Name','Blindleistung','Menu','none','ToolBar','none');
set(gca,'FontSize',10, 'FontName','Cambria')
plot(1:nl,BranchFlowAC_Q,'Marker','*','LineWidth',1,'Color','#A2142F','MarkerSize', 5,'LineStyle', '-');hold on;
plot(1:nl,branchQ,'Marker','o','LineWidth',1, 'Color', '#EDB120','MarkerSize', 5,'LineStyle', '-');
legend('MATPOWER ACOPF','DCOPF-QV');
xlabel('Leitung','FontSize',10, 'FontName','Cambria');
ylabel('Blindleistungsfluss (MVar)','FontSize',10, 'FontName','Cambria');
title('Vergleich der Blindleistungsflüsse','FontSize',10, 'FontName','Cambria');
%ylim([-0.3 1])
xlim([1 nl+1]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('Blindleistung' ,'-dsvg')

% Colors for the bar chart
colorbars='#EDB120';
colorbarsDC='#0072BD';
%

fh4=figure('Name','ALMP_Dif','Menu','none','ToolBar','none');
set(gca,'FontSize',10, 'FontName','Cambria')
t=tiledlayout(3,1);
nexttile([2,1])
plot(1:nb,ALMPAC,'Marker','*','LineWidth',0.8,'Color','#A2142F','MarkerSize', 4,'LineStyle', '-');hold on;
plot(1:nb,LMP_P,'Marker','o','LineWidth',0.8, 'Color', '#EDB120','MarkerSize', 4,'LineStyle', '-'); hold on;
plot(1:nb,ALMPDC, 'Marker','.','LineWidth',0.8, 'Color', '#0072BD','MarkerSize', 4,'LineStyle', '-')
ylabel('Wirkleistungskosten ($/MWh)','FontSize',10, 'FontName','Cambria','Color','k');
legend('MATPOWER ACOPF','DCOPF-QV', 'MATPOWER DCOPF');
%ylim([16 32])
xlim([1 nb+1])
nexttile
barP=bar(1:nb,[Fehler_ALMP,FehlerDC]); hold on;
barP(1).FaceColor=colorbars;
barP(2).FaceColor=colorbarsDC;
barP(1).EdgeColor=colorbars;
barP(2).EdgeColor=colorbarsDC;
ylabel('rel. Fehler (%)','FontSize',10,'FontName','Cambria','Color','k');
%ylim([0 50])
yticks(0:10:50)
xlim([1 nb+1])
legend('rel. Fehler - DCOPF-QV','rel. Fehler - MATPOWER DCOPF');
xlabel(t,'Knoten','FontSize',10,'FontName','Cambria');
%title(t,'Vergleich der ALMPs','FontSize',10, 'FontName','Cambria','Fontweight', 'bold');
t.TileSpacing = 'none';
t.Padding = 'none';
print('ALMP_Dif' ,'-dsvg')


fh5=figure('Name','RLMP','Menu','none','ToolBar','none');
set(gca,'FontSize',10, 'FontName','Cambria')
t=tiledlayout(3,1);
nexttile([2 1])
plot(1:nb,RLMPAC,'Marker','*','LineWidth',0.8,'Color','#A2142F','MarkerSize', 4,'LineStyle', '-');hold on;
plot(1:nb,LMP_Q,'Marker','o','LineWidth',0.8, 'Color', '#EDB120','MarkerSize', 4,'LineStyle', '-');
ylabel('Blindleistungskosten ($/MVar)','FontSize',10, 'FontName','Cambria','Color','k');
%ylim([0 15]);
xlim([1 nb+1]);
legend('ACOPF','DCOPF-QV');
nexttile
barQ=bar(1:nb,Fehler_RLMP,0.4);
barQ.FaceColor=colorbars;
barQ.EdgeColor=colorbars;
ylabel('rel. Fehler (%)','FontSize',10,'FontName','Cambria','Color','k');
%ylim([0 50]);
%yticks(0:10:50);
xlim([1 nb+1]);
legend('rel. Fehler - DCOPF-QV');
xlabel(t,'Knoten','FontSize',10, 'FontName','Cambria');
%title(t,'Vergleich der RLMPs','FontSize',10, 'FontName','Cambria','Fontweight', 'bold');
t.TileSpacing = 'none';
t.Padding = 'none';
print('RLMP' ,'-dsvg')