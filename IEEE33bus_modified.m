%% system MVA base
baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
busdata = [  %% (Pd and Qd are specified in kW & kVAr here, converted to MW & MVAr below)
	1	3	0	0	0	0	1	1	0	1	1	1	1;
	2	1	180	0	0	0	1	1	0	10	1	1.1	0.9;
	3	1	320	0	0	0	1	1	0	10	1	1.1	0.9;
	
	
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
gendata = [
	1	0	0	10	-10	1	100	1	600	0	0	0	0	0	0	0	0	0	0	0	0;
    3	0	0   10	-10	1	100	1	300	0	0	0	0	0	0	0	0	0	0	0	0;
    
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
branchdata = [  %% (r and x specified in ohms here, converted to p.u. below)
	1	2	0	0.25	0	200	0	0	0	0	1	-360	360;
	2	3	0	0.25	0	200	0	0	0	0	1	-360	360;
    1	3	0	0.25	0	200	0	0	0	0	1	-360	360;
		
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencostdata = [
	2	0	0	3	0	10	0;
    2	0	0	3	0	20	0;
    

];


%% convert branch impedances from Ohms to p.u.

%Vbase = busdata(1, 10) * 1e3;      %% in Volts
%Sbase = baseMVA * 1e6;              %% in VA
%branchdata(:, [3 4]) = branchdata(:, [3 4]) / (Vbase^2 / Sbase);

%% convert loads from kW to MW
%busdata(:, [3, 4]) = busdata(:, [3, 4]) / 1e3;