clear;clc;
%% Iteration setup to go through the calculated files at each timestep
step=0;
minute=0;
Loadcurvestruct=[];LMP_P=[];LMP_Q=[];ALMPAC=[];RLMPAC=[];ALMPDC=[];BusVolAC=[];BusVolDC=[];busVol=[];Fehler_ALMP=[];Fehler_RLMP=[];AEA_ALMP=[];AEA_RLMP=[];AEADC=[];FehlerDC=[];
while step<24
%% select which casefiles to load (names have to match) 
%% results of the 24 hour loadcurve are stored in the loadprofile folder and are accessed here to plot the graphs
casedata='case69m';
day='_10_february_2020';
f7='very_tight_voltage';             
f1='_loadcurve after_';
f2=num2str(floor(step));
f3='_hours from minute_';
f4=num2str(minute*60);
f5='_up to minute_';
f6=num2str((minute+0.25)*60);
filename=strcat(casedata,f7,f1,f2,f3,f4,f5,f6,day);
f=fullfile('loadprofile', casedata , filename);
%% Datamatrix in which all results are stored
Loadcurvestruct=[Loadcurvestruct,load(f)];
minute=minute+0.25;
if minute > 0.75
    minute=0;
end
step=step+0.25;
end
timeline=0.25:0.25:24;
%% Accessing the data from the datamatrix at each timepoint
for i=1:96
LMP_P=[LMP_P,Loadcurvestruct(i).LMP_P];
LMP_Q=[LMP_Q,Loadcurvestruct(i).LMP_Q];
ALMPAC=[ALMPAC,Loadcurvestruct(i).ALMPAC];
RLMPAC=[RLMPAC,Loadcurvestruct(i).RLMPAC];
ALMPDC=[ALMPDC,Loadcurvestruct(i).ALMPDC];
BusVolAC=[BusVolAC,Loadcurvestruct(i).BusVolAC];
BusVolDC=[BusVolDC,Loadcurvestruct(i).BusVolDC];
busVol=[busVol,Loadcurvestruct(i).busVol];
Fehler_ALMP=[Fehler_ALMP,Loadcurvestruct(i).Fehler_ALMP];
Fehler_RLMP=[Fehler_RLMP,Loadcurvestruct(i).Fehler_RLMP];
FehlerDC=[FehlerDC,Loadcurvestruct(i).FehlerDC];
AEA_ALMP=[AEA_ALMP,Loadcurvestruct(i).AEA_ALMP];
AEA_RLMP=[AEA_RLMP,Loadcurvestruct(i).AEA_RLMP];
AEADC=[AEADC,Loadcurvestruct(i).AEADC];
end

%% Defining and renaming variables for the 24 hour plot
loadcurve_base=Loadcurvestruct(1).loadcurve_base;
[maxload,maxloadindex]=max(loadcurve_base);
[minload,minloadindex]=min(loadcurve_base);
mpc=Loadcurvestruct(1).mpc;
nb=size(mpc.bus(:,1),1);
busnumber=1:1:nb;
timeline=timeline';
AEA_ALMP_24_hour=mean(AEA_ALMP*100);
AEADC_24_hour=mean(AEADC*100);
AEA_RLMP_24_hour=mean(AEA_RLMP*100);

%% Clearing the datamatrix to increase the performance for further plotting
clearvars Loadcurvestruct;
%% 24 hour Loadcurve Plotting in 3D
% ALMP of the DCOPF-QV

fh1=figure('Name','ALMP_24');
colormap parula
surf(timeline,busnumber,LMP_P, 'EdgeColor', 'flat', 'FaceColor', 'none')
set(gca,'FontSize',10, 'FontName','Cambria')
xlabel('Zeit (h)','FontSize',10, 'FontName','Cambria');
ylabel('Knoten','FontSize',10, 'FontName','Cambria')
zlabel('Wirkleistungskosten ($/MWh)','FontSize',10, 'FontName','Cambria');
view([-60 20]);
xticks(0:4:24)
ylim([1 nb])
%zlim([24 50])
%view([45 30]);
%title('ALMP 24 Stunden Lastkurve des DCOPF-QV');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('ALMP_24' ,'-dsvg')

% RLMP of the DCOPF-QV

fh2=figure('Name','RLMP_24');
surf(timeline,busnumber,LMP_Q, 'EdgeColor', 'flat', 'FaceColor', 'none')
set(gca,'FontSize',10, 'FontName','Cambria')
xlabel('Zeit (h)','FontSize',10, 'FontName','Cambria');
ylabel('Knoten','FontSize',10, 'FontName','Cambria')
zlabel('Blindleistungskosten ($/MVarh)','FontSize',10, 'FontName','Cambria');
view([-60 20]);
%view([45 30]);
xticks(0:4:24)
ylim([1 nb])
%zlim([2.5 30])
%title('RLMP 24 Stunden Lastkurve des DCOPF-QV');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('RLMP_24' ,'-dsvg')


% ALMP of the ACOPF

fh3=figure('Name','ACOPF ALMP_24','Menu','none','ToolBar','none');
surf(timeline,busnumber,ALMPAC, 'EdgeColor', 'flat', 'FaceColor', 'none')
set(gca,'FontSize',10, 'FontName','Cambria')
xlabel('Zeit (h)','FontSize',10, 'FontName','Cambria');
ylabel('Knoten','FontSize',10, 'FontName','Cambria')
zlabel('Wirkleistungskosten ($/MWh)','FontSize',10, 'FontName','Cambria');
view([-60 20]);
xticks(0:4:24)
ylim([1 nb])
%zlim([14 32])
%view([45 30]);
%title('ALMP 24 Stunden Lastkurve des MATPOWER ACOPF');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('ACOPF ALMP_24' ,'-dsvg')

% RLMP of the ACOPF

fh4=figure('Name','ACOPF RLMP_24');
surf(timeline,busnumber,RLMPAC, 'EdgeColor', 'flat', 'FaceColor', 'none')
set(gca,'FontSize',10, 'FontName','Cambria')
xlabel('Zeit (h)','FontSize',10, 'FontName','Cambria');
ylabel('Knoten','FontSize',10, 'FontName','Cambria')
zlabel('Blindleistungskosten ($/MVarh)','FontSize',10, 'FontName','Cambria');
view([-60 20]);
%view([45 30]);
xticks(0:4:24)
ylim([1 nb])
%zlim([4.5 18])
%title('RLMP 24 Stunden Lastkurve des MATPOWER ACOPF');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('ACOPF RLMP_24' ,'-dsvg')

% ALMP Error Plot of the DCOPF-QV

fh5=figure('Name', 'Error ALMP 24');
colormap(summer)
surf(timeline,busnumber,Fehler_ALMP, 'EdgeColor', 'flat', 'FaceColor', 'none');
set(gca,'FontSize',10, 'FontName','Cambria')
xlabel('Zeit (h)','FontSize',10, 'FontName','Cambria');
ylabel('Knoten','FontSize',10, 'FontName','Cambria')
zlabel('rel. Fehler (%)','FontSize',10, 'FontName','Cambria');
%title('Relative Abweichung in [%] der DCOPF-QV ALMP 24 Stunden Lastkurve');
view([-60 20]);
xticks(0:4:24)
ylim([1 nb])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('Error ALMP 24' ,'-dsvg')

% RLMP Error Plot of the DCOPF-QV

fh6=figure('Name', 'Error RLMP 24');
colormap((summer))
surf(timeline,busnumber,Fehler_RLMP, 'EdgeColor', 'flat', 'FaceColor', 'none');
set(gca,'FontSize',10, 'FontName','Cambria')
xlabel('Zeit (h)','FontSize',10, 'FontName','Cambria');
ylabel('Knoten','FontSize',10, 'FontName','Cambria')
zlabel('rel. Fehler (%)','FontSize',10, 'FontName','Cambria');
%title('Relative Abweichung in [%] der DCOPF-QV RLMP 24 Stunden Lastkurve');
view([-60 20]);
xticks(0:4:24)
ylim([1 nb])
%zlim([0 25]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom (ax_width-0.01) ax_height];
print('Error RLMP 24' ,'-dsvg')
