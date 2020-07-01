% Load case
define_constants;                          %currently not needed

mpc = ext2int(loadcase('case33bw'));        %used for makeYbus function my cases have a different name so here 
                                           %you have to select the corresponding matpower casename (will fix this later on)

%% Acquiring Data
nl=size(branchdata(:,1),1);         %Number of lines/branches
nb=size(busdata(:,1),1);            %Number of buses
ng=size(gendata(:,1),1);            %Number of generators
status = branchdata(:,11);          %Checks for active/inactive lines
fb=branchdata(:,1);                 %From Bus ...
tb=branchdata(:,2);                 %To Bus ...
prat=branchdata(:,6);               %Rated line flow limit
vmin=busdata(:,13);                 %min voltage in p.u.
vmax=busdata(:,12);                 %max voltage in p.u.
slack = find(busdata(:, 2) == 3);   %finding the reference/slack bus
noslack = find((1:nb)' ~= slack);   %bus numbering except the slack bus
noslack2= find((1:(2*nb))' ~= slack); %bus numbering from 1 to 2*number of buses except slack bus (used for C matrix with transmission networks)
noslack3= [noslack;nb+noslack];     %bus numbering from 1 to 2*number of buses except for 2 slack bus rows/columns (used for C matrix with distribution networks)
iline = [(1:nl)'; (1:nl)'];         %used for incidence matrix with sparse function
vbase=busdata(:,8);                 %basevalue for the voltages in p.u.



%% Tap ratios
tap = ones(nl, 1);                                      %% default tap ratio = 1
i=find(branchdata(:,9));                                %% tap ratio based on makeYbus function of matpower
tap(i)=branchdata(i,9);                                 
tap = tap .* exp(1j*pi/180 * branchdata(:, 10));

%% Conductance and Susceptance of the series elements of the lines
G= status.*branchdata(:,3)./(branchdata(:,4).^2+branchdata(:,3).^2);        %conductance and susceptance taken from the paper
B= -status.*branchdata(:,4)./(branchdata(:,4).^2+branchdata(:,3).^2);       %"A State-Independent Linear Power Flow Model With Accurate Estimation of Voltage Magnitude"
G= G./tap;
B= B./tap;

%% Ybus created with matpower function
Ybus = makeYbus(mpc);
Gbus = real(Ybus);
Bbus = imag(Ybus);


%% Creating the Submatrices for the C Matrix based on "A Linear LMP Model for Active and Reactive Power with Power Loss"
GP = -Gbus;                        %C = -[BP -GP;GQ BQ]              
BD = diag(sum(Bbus));             %Susceptance shunt elements
GD = diag(sum(Gbus));             %Conductance shunt elements  
BP = Bbus - BD;
BQ = Bbus;
GQ = Gbus - GD;                   %GQ approxiately equals -Gbus


%% Incidence Matrix either 1 if branches connect from bus a to b or -1 if branches connect from bus b to a 
%% or 0 if no connection between buses

Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);

%% Incidence Matrix applied to series conductance and susceptance
Bf = sparse(iline, [fb; tb], [B; -B], nl, nb);
Gf = sparse(iline, [fb; tb], [G; -G], nl, nb);


%% Constructing the Inverse of C
C=[BP GP;GQ BQ];

if nb==118                                                              %removing 1 or 2 rows/columns currently based on number of buses
Cinverse=full(C(noslack2,noslack2)\eye(2*nb-1));                        %removes 1 row and column for transmission networks
else
Cinverse=full(C(noslack3,noslack3)\eye(size(C(noslack3,noslack3))));    %removes 2 rows and columns for distribution networks
end
% inverted C --> Xmat filled with zeros at slack bus row and column
X=zeros(2*nb,2*nb);
if nb==118
    X(noslack2,noslack2)=Cinverse;                                      %1 removed row/column filled with zeros
else
X(noslack3,noslack3)=Cinverse;                                          %2 removed rows/columns filled with zeros
end

C(:,slack)=[];
C(slack,:)=[];
Cinverse=pinv(full(C));
X=zeros(2*nb);
X(noslack2,noslack2)=-Cinverse;

%X=zeros(2*nb,2*nb);
%Cinverse=pinv(full(C(noslack2,noslack2)));
%X(noslack2,noslack2)=Cinverse;
%% Constructing partial GSF matrices based on "A Linear LMP Model for Active and Reactive Power with Power Loss"
GSF_PP =  ( Gf *  X(nb+1:end,1:nb)      -Bf * X(1:nb,1:nb));
GSF_PQ =  ( Gf *  X(nb+1:end,nb+1:end)  -Bf * X(1:nb,nb+1:end));
GSF_QP =  (-Gf *  X(1:nb,1:nb)          -Bf * X(nb+1:end,1:nb));
GSF_QQ =  (-Gf *  X(1:nb,nb+1:end)      -Bf * X(nb+1:end,nb+1:end));


%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Inf
for i=1:length(prat)
    if prat(i)==0
        prat(i)=Inf;
    end
end