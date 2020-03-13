%Acquiring Data

nl=size(branchdata(:,1),1);         %Number of lines/branches
nb=size(busdata(:,1),1);            %Number of buses
ng=size(gendata(:,1),1);            %Number of generators
status = branchdata(:, 11);         %Checks the status of active/inactive lines
fb=branchdata(:,1);                 %From Bus ...
tb=branchdata(:,2);                 %To Bus ...
r= status ./ branchdata(:,3);       %Resistance of the line
X= status ./ branchdata(:,4);       %Reactance of the line
prat=branchdata(:,6);               %Rated line flow limit
slack = find(busdata(:, 2) == 3);   %Finding the reference bus
    % Slack bus Distribution NOT WORKING PROPERLY FOR MORE THAN 1 SLACK


%Checking the Data for number of slack buses; used for distribution of
%weights

if length(slack)==1
    ns=slack;
else 
    ns=1;
end


%PTDF Matrix Construction

iline = [(1:nl)'; (1:nl)'];
Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);
Bf = sparse(iline, [fb; tb], [X; -X], nl, nb);
Bbusnew = Cft' * Bf;
noref   = (2:nb)';      %% use bus 1 for voltage angle reference
noslack = find((1:nb)' ~= ns);
PTDF = zeros(nl, nb);
PTDF(:, noslack) = full(Bf(:, noref) / Bbusnew(noslack, noref));


%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Pmax of the slack bus generator

for i=1:length(prat)
    if prat(i)==0
        prat(i)=gendata(ns,9);
    end
end


%Slack distribution - may be needed for different cases // currently not
%working properly
if length(slack) ~= 1
    if size(slack, 2) == 1  %% slack is a vector of weights
        slack=slack/sum(slack);
PTDF = PTDF *(eye(nb,nb)-slack*ones(1,nb));
    else
        PTDF = PTDF * slack;
    end
end