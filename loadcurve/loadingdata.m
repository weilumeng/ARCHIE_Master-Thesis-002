clear;clc;
step=0;
minute=0;
Loadcurvestruct=[];LMP_P=[];LMP_Q=[];ALMPAC=[];RLMPAC=[];ALMPDC=[];
while step<23.75
casedata='case69m';
f7='tight_voltage';             %%ACTIVE ONLY WHEN NEEDED
f1='_loadcurve after_';
f2=num2str(floor(step));
f3='_hours from minute_';
f4=num2str(minute*60);
f5='_up to minute_';
f6=num2str((minute+0.25)*60);
filename=strcat(casedata,f7,f1,f2,f3,f4,f5,f6);
f=fullfile('loadprofile', casedata , filename);
Loadcurvestruct=[Loadcurvestruct,load(f)];
minute=minute+0.25;
if minute > 0.75
    minute=0;
end
step=step+0.25;
end
timeline=0.25:0.25:23.75;
for i=1:95
LMP_P=[LMP_P,Loadcurvestruct(i).LMP_P];
LMP_Q=[LMP_Q,Loadcurvestruct(i).LMP_Q];
ALMPAC=[ALMPAC,Loadcurvestruct(i).ALMPAC];
RLMPAC=[RLMPAC,Loadcurvestruct(i).RLMPAC];
ALMPDC=[ALMPDC,Loadcurvestruct(i).ALMPDC];
end

loadcurve_base=Loadcurvestruct(1).loadcurve_base;
[maxload,maxloadindex]=max(loadcurve_base);
[minload,minloadindex]=min(loadcurve_base);
mpc=Loadcurvestruct(1).mpc;
nb=size(mpc.bus(:,1),1);
busnumber=1:1:nb;
timeline=timeline';

clearvars Loadcurvestruct;
%% Plot 3D
figure
surf(timeline,busnumber,LMP_P, 'EdgeColor', 'flat', 'FaceColor', 'none')
xlabel('Timesteps');
ylabel('Bus')
zlabel('$/MWh');
view([-60 15]);
%view([45 30]);
colormap parula;
cb=colorbar;
pos=get(cb,'Position');
set(cb,'Position',pos+[0.1,0,0.0,0.0]);
title('ALMP Loadcurve');

figure
surf(timeline,busnumber,LMP_Q, 'EdgeColor', 'flat', 'FaceColor', 'none')
xlabel('Timesteps');
ylabel('Bus')
zlabel('$/MVarh');
view([-60 15]);
%view([45 30]);
colormap parula;
cb=colorbar;
pos=get(cb,'Position');
set(cb,'Position',pos+[0.1,0,0.0,0.0]);
title('RLMP Loadcurve');

figure
surf(timeline,busnumber,ALMPAC, 'EdgeColor', 'flat', 'FaceColor', 'none')
xlabel('Timesteps');
ylabel('Bus')
zlabel('$/MWh');
view([-60 15]);
%view([45 30]);
colormap parula;
cb=colorbar;
pos=get(cb,'Position');
set(cb,'Position',pos+[0.1,0,0.0,0.0]);
title('ALMP ACOPF Loadcurve');

figure
surf(timeline,busnumber,RLMPAC, 'EdgeColor', 'flat', 'FaceColor', 'none')
xlabel('Timesteps');
ylabel('Bus')
zlabel('$/MVarh');
view([-60 15]);
%view([45 30]);
colormap parula;
cb=colorbar;
pos=get(cb,'Position');
set(cb,'Position',pos+[0.1,0,0.0,0.0]);
title('RLMP ACOPF Loadcurve');

%% Plot 2d
set_linewidth = 1.5;
set_markersize = 5;
color_AC_min=[0.6350 0.0780 0.1840];
color_AC_max = [0.4660 0.6740 0.1880];
color_LMPL_max = [0.4940 0.1840 0.5560];
color_LMPL_min = [0.9290 0.6940 0.1250];

figure('Name','ALMP')
plot(1:nb,ALMPAC(:,minloadindex),'Marker','*','LineWidth',set_linewidth,'Color',color_AC_min, 'MarkerSize',set_markersize,'MarkerEdgeColor',color_AC_min,'MarkerFaceColor',color_AC_min);hold on;
plot(1:nb,LMP_P(:,minloadindex),'Marker','o','LineWidth',set_linewidth, 'Color', color_LMPL_min, 'MarkerSize',set_markersize,'MarkerEdgeColor',color_LMPL_min,'MarkerFaceColor',color_LMPL_min); hold on;
plot(1:nb,ALMPDC(:,minloadindex), 'Marker','.','LineWidth',set_linewidth, 'Color', '#0072BD', 'MarkerSize',set_markersize);hold on;
plot(1:nb,ALMPAC(:,maxloadindex),'Marker','s','LineWidth',set_linewidth, 'MarkerSize',set_markersize,'Color',color_AC_max,'MarkerEdgeColor',color_AC_max,'MarkerFaceColor',color_AC_max);hold on;
plot(1:nb,LMP_P(:,maxloadindex),'Marker','+','LineWidth',set_linewidth, 'MarkerSize',set_markersize,'Color',color_LMPL_max,'MarkerEdgeColor',color_LMPL_max,'MarkerFaceColor',color_LMPL_max); hold on;
plot(1:nb,ALMPDC(:,maxloadindex),'LineWidth',set_linewidth, 'MarkerSize',set_markersize,'Color','#0072BD');
legend('ACOPF min','My DCOPF min', 'DCOPF min','ACOPF max','My DCOPF max', 'DCOPF max','Location','best','Orientation','vertical');
axis([1,nb,min(ALMPAC(:,maxloadindex))-5,max(ALMPAC(:,maxloadindex))+5]);
xlabel('Bus','FontSize',12,'FontName','Cambria');
ylabel('$/MWh','FontSize',12,'FontName','Cambria');
title('ALMP Comparison','FontSize',12,'FontName','Cambria');

figure('Name','RLMP')
plot(1:nb,RLMPAC(:,minloadindex),'Marker','*','LineWidth',set_linewidth,'Color',color_AC_min, 'MarkerSize',set_markersize,'MarkerEdgeColor',color_AC_min,'MarkerFaceColor',color_AC_min);hold on;
plot(1:nb,LMP_Q(:,minloadindex),'Marker','o','LineWidth',set_linewidth, 'Color', color_LMPL_min, 'MarkerSize',set_markersize,'MarkerEdgeColor',color_LMPL_min,'MarkerFaceColor',color_LMPL_min); hold on;
plot(1:nb,RLMPAC(:,maxloadindex),'Marker','s','LineWidth',set_linewidth,'Color',color_AC_max, 'MarkerSize',set_markersize,'MarkerEdgeColor',color_AC_max,'MarkerFaceColor',color_AC_max);hold on;
plot(1:nb,LMP_Q(:,maxloadindex),'Marker','+','LineWidth',set_linewidth,'Color',color_LMPL_max, 'MarkerSize',set_markersize,'MarkerEdgeColor',color_LMPL_max,'MarkerFaceColor',color_LMPL_max); hold on;
axis([1,nb,min(RLMPAC(:,minloadindex))-5,max(RLMPAC(:,maxloadindex))+5]);
legend('ACOPF min','My DCOPF min','ACOPF max','My DCOPF max','Location','best','Orientation','vertical');
xlabel('Bus','FontSize',12,'FontName','Cambria');
ylabel('$/MVarh','FontSize',12,'FontName','Cambria');
title('RLMP Comparison','FontSize',12,'FontName','Cambria');