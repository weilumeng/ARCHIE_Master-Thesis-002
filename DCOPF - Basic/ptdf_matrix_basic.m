%% Acquiring Data

nl=size(branchdata(:,1),1);         %Number of lines/branches
nb=size(busdata(:,1),1);            %Number of buses
ng=size(gendata(:,1),1);            %Number of generators
status = branchdata(:,11);          %Checks the status of active/inactive lines
fb=branchdata(:,1);                 %From Bus ...
tb=branchdata(:,2);                 %To Bus ...
G= status ./ branchdata(:,3);       %Conductance of the line
B= status ./ branchdata(:,4);       %Susceptance of the line
prat=branchdata(:,6);               %Rated line flow limit
slack = find(busdata(:, 2) == 3);   %Finding the reference bus
                                    %Demand of each bus


% Tap Ratios
tap=ones(nl,1);
i=find(branchdata(:,9));
tap(i)=branchdata(i,9);
B= B ./ tap;
G= G ./ tap;




%PTDF Matrix Construction

iline = [(1:nl)'; (1:nl)'];
Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);
Bf = sparse(iline, [fb; tb], [B; -B], nl, nb);
Bbusnew = Cft' * Bf;
noslack = find((1:nb)' ~= slack);
PTDF = zeros(nl, nb);
PTDF(:, noslack) = full(Bf(:, noslack) / Bbusnew(noslack, noslack));
%PTDF(:, :) = full(Bf(:, :) / Bbusnew(:, :));

%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Inf

for i=1:length(prat)
    if prat(i)==0
        prat(i)=Inf;
    end
end
