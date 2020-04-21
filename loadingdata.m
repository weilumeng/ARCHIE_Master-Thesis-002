clear;
step=1;
while step<24
f1=num2str(step);
f2='timestep';
filename=strcat(f2,f1);
f=fullfile('loadprofile',filename);
load(f)
figure
ylim1=(min(lmp)-0.01*min(lmp));
ylim2=(max(lmp)+0.01*min(lmp));
plot(busnumber,lmp, 'Marker', '*', 'Color', 'r');
title('LMP Plot')
xlabel('Bus')
ylabel('$/MW')
xticks((linspace(1,max(busnumber),max(busnumber))));
ylim([ylim1 ylim2]);
xlim([min(busnumber) max(busnumber)])
ytickformat('usd')
step=step+1;
end