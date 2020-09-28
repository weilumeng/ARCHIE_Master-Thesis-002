define_constants;                    %needed for the mpc function of matpower       

[baseMVA, bus, gen, branch, gencost] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);

%% Acquiring Data

nl=size(branch(:,1),1);         %Number of lines/branches
nb=size(bus(:,1),1);            %Number of buses
ng=size(gen(:,1),1);            %Number of generators
status = branch(:,11);          %Checks the status of active/inactive lines
fb=branch(:,1);                 %From Bus ...
tb=branch(:,2);                 %To Bus ...
G= status ./ branch(:,3);       %Conductance of the line
B= status ./ branch(:,4);       %Susceptance of the line
prat=branch(:,6);               %Rated line flow limit
slack = find(bus(:, 2) == 3);   %Finding the reference bus
                                    %Demand of each bus


% Tap Ratios
tap=ones(nl,1);
i=find(branch(:,9));
tap(i)=branch(i,9);
B= B ./ tap;
G= G ./ tap;




%PTDF Matrix Construction

iline = [(1:nl)'; (1:nl)'];
Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);
Bf = sparse(iline, [fb; tb], [B; -B], nl, nb);
Bbus = Cft' * Bf;
noslack = find((1:nb)' ~= slack);
PTDF = zeros(nl, nb);
PTDF(:, noslack) = full(Bf(:, noslack) / Bbus(noslack, noslack));
%PTDF(:, :) = full(Bf(:, :) / Bbusnew(:, :));

%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Inf

for i=1:length(prat)
    if prat(i)==0
        prat(i)=Inf;
    end
end
