%% Acquiring Data

nl=size(branchdata(:,1),1);         %Number of lines/branches
nb=size(busdata(:,1),1);            %Number of buses
ng=size(gendata(:,1),1);            %Number of generators
status = branchdata(:,11);          %Checks the status of active/inactive lines
fb=branchdata(:,1);                 %From Bus ...
tb=branchdata(:,2);                 %To Bus ...
iline = [(1:nl)'; (1:nl)'];

G= status ./ branchdata(:,3);       %Conductance of the line
B= status ./ branchdata(:,4);       %Susceptance of the line

prat=branchdata(:,6);               %Rated line flow limit

slack = find(busdata(:, 2) == 3);   %Finding the reference bus
noslack = find((1:nb)' ~= slack);

% Tap Ratios
tap=ones(nl,1);
i=find(branchdata(:,9));
tap(i)=branchdata(i,9);
B= B ./ tap;
G= G ./ tap;


%PTDF Matrix Construction

Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);
Bf = sparse(iline, [fb; tb], [B; -B], nl, nb);
Bbus = Cft' * Bf;

PTDF = zeros(nl, nb);
PTDF(:, noslack) = full(Bf(:, noslack) / Bbus(noslack, noslack));
%PTDF(:, :) = full(Bf(:, :) / Bbus(:, :));


%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Inf

for i=1:nl
    if prat(i)==0
        prat(i)=Inf;
    end
end
