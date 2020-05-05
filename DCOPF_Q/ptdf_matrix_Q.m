%% Acquiring Data

nl=size(branchdata(:,1),1);         %Number of lines/branches
nb=size(busdata(:,1),1);            %Number of buses
ng=size(gendata(:,1),1);            %Number of generators
status = branchdata(:,11);          %Checks for active/inactive lines
fb=branchdata(:,1);                 %From Bus ...
tb=branchdata(:,2);                 %To Bus ...
G= status ./ branchdata(:,3);       %Conductance of the line
B= status ./ branchdata(:,4);       %Susceptance of the line
Bc = status .* branchdata(:, 5);    %total line charging susceptance
Gs = busdata(:,5);
Bs = busdata(:,6);
prat=branchdata(:,6);               %Rated line flow limit
vmin=busdata(:,13);                 %max voltage in p.u.
vmax=busdata(:,12);                 %min voltage in p.u.
slack = find(busdata(:, 2) == 3);   %Finding the reference bus


%% Tap Ratios
tap=ones(nl,1);
i=find(branchdata(:,9));
tap(i)=branchdata(i,9);
B= B ./ tap;
G= G ./ tap;
                                                    
for i=1:length(branchdata(:,1))                              
    if G(i)==Inf
        G(i)=0;
    end
end
%% PTDF Matrix Construction

iline = [(1:nl)'; (1:nl)'];
Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);
Bf = sparse(iline, [fb; tb], [B; -B], nl, nb);
Bbus = Cft' * Bf;
noref   = (2:nb)';      %% use bus 1 for voltage angle reference
noslack = find((1:nb)' ~= slack);
noslack2= find((1:(2*nb))' ~= slack);

%% recommended version
PTDF = zeros(nl, nb);
PTDF(:, noslack) = full(Bf(:, noref) / Bbus(noslack, noref));


%% Constructing the Bus Matrices w/ shunts
Bf_shunt = sparse(iline, [fb; tb], [B; -B], nl, nb);
Bbus_shunt = Cft' * Bf_shunt;

%Adding shunts (not sure if thats correct)
for i=1:nb
    for j=1:nb
        if i==branchdata(j,1)
            if branchdata(j,5)~=0
                Bbus_shunt(i,i)=Bbus_shunt(i,i)+(1/(branchdata(j,5)/2));
           end
        end
    end
end

Gf = sparse(iline, [fb; tb], [G; -G], nl, nb);
Gbus = Cft' * Gf;

Gf_shunt = sparse(iline, [fb; tb], [G; -G], nl, nb);
Gbus_shunt = Cft' * Gf_shunt;

% Adding shunts (not sure if thats correct)
for i=1:nb
    for j=1:nb
        if i==branchdata(j,1)
            if branchdata(j,5)~=0
                Gbus_shunt(i,i)=Gbus_shunt(i,i)+(1/(branchdata(j,5)/2));
            end
        end
    end
end

%% Constructing the Inverse of C
C=[Bbus -Gbus_shunt;Gbus Bbus_shunt];


% inverted C --> Xmat filled with zeros at slack bus row and column

Xmat=zeros(2*nb,2*nb);
Xmat(noslack2,noslack2)=-(C(noslack2,noslack2)\eye(2*nb-1));


%% Constructing new partial GSF matrices 

GSF_PP =  full(Gf * Xmat(nb+1:end,1:nb))    -full(Bf * Xmat(1:nb,1:nb));
GSF_PQ =  full(Gf * Xmat(nb+1:end,nb+1:end))-full(Bf * Xmat(1:nb,nb+1:end));
GSF_QP = -full(Gf * Xmat(1:nb,1:nb))        -full(Bf * Xmat(nb+1:end,1:nb));
GSF_QQ = -full(Gf * Xmat(1:nb,nb+1:end))    -full(Bf * Xmat(nb+1:end,nb+1:end));


%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Inf

for i=1:length(prat)
    if prat(i)==0
        prat(i)=Inf;
    end
end



G= status ./ branchdata(:,3);
G= G ./ tap;