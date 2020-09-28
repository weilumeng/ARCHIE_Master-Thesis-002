%% Load case
define_constants;                    %needed for the mpc function of matpower       

[baseMVA, bus, gen, branch, gencost] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);

%% Acquiring Data

nl=size(branch(:,1),1);         %Number of lines/branches
nb=size(bus(:,1),1);            %Number of buses
ng=size(gen(:,1),1);            %Number of generators
status = branch(:,11);          %Checks for active/inactive lines
fb=branch(:,1);                 %From Bus ...
tb=branch(:,2);                 %To Bus ...
iline = [(1:nl)'; (1:nl)'];


prat=branch(:,6);               %Rated line flow limit

slack = find(bus(:, 2) == 3);   %finding the reference/slack bus
noslack = find((1:nb)' ~= slack);
noslack2= find((1:(2*nb))' ~= slack);
noslack3 =[noslack;nb+noslack];


%% Tap Ratios
tap=ones(nl,1);
i=find(branch(:,9));
tap(i)=branch(i,9);
tap = tap .* exp(1j*pi/180 * branch(:, 10));

G= status.*branch(:,3)./(branch(:,3).^2+branch(:,4).^2);
B= -status.*branch(:,4)./(branch(:,3).^2+branch(:,4).^2);
G=G./tap;
B=B./tap;


%% Ybus created with matpower function
Ybus = makeYbus(mpc);
Gbus = real(Ybus);
Bbus = imag(Ybus);

%% Constructing the Bus Matrices

GP = Gbus;
BD = diag(sum(Bbus'));             %shunt elements
GD = diag(sum(Gbus'));
BP = -Bbus + BD;
BQ = -Bbus;
GQ = -Gbus + GD;                   %GQ approxiately equals -Gbus
%% GSF Matrix Construction

Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);

Bf = sparse(iline, [fb; tb], [B; -B], nl, nb);
Gf = sparse(iline, [fb; tb], [G; -G], nl, nb);


%% Constructing the Inverse of C
C=[BP,GP;GQ,BQ];

%Cinverse=pinv(full(C(noslack2,noslack2)));
Cinverse=(C(noslack2,noslack2)\speye(size(C(noslack2,noslack2))));
%Cinverse=lsqminnorm(C(noslack2,noslack2),speye(2*nb-1));              %removes 1 row and column for transmission networks
%Cinverse=(C(noslack3,noslack3)\speye(size(C(noslack3,noslack3))));    %removes 2 rows and columns for distribution networks

%% Pseudo inversion with 1 row and column removed

%Cinverse=pinv(full(C));                              

% inverted C --> X filled with zeros at slack bus row and column
X=zeros(2*nb);

                                    
%X(noslack3,noslack3)=Cinverse;     %1 removed row/column filled with zeros
X(noslack2,noslack2)=Cinverse;      %2 removed rows/columns filled with zeros


%% Constructing new partial GSF matrices 
GSF_PP =  full( Gf *  X(nb+1:end,1:nb)      -Bf * X(1:nb,1:nb));
GSF_PQ =  full( Gf *  X(nb+1:end,nb+1:end)  -Bf * X(1:nb,nb+1:end));
GSF_QP =  full(-Gf *  X(1:nb,1:nb)          -Bf * X(nb+1:end,1:nb));
GSF_QQ =  full(-Gf *  X(1:nb,nb+1:end)      -Bf * X(nb+1:end,nb+1:end));


E_P_est=0*ones(nb,1);
E_Q_est=0*ones(nb,1);
E_P_est_old=0*ones(nb,1);
E_Q_est_old=0*ones(nb,1);
DF_P_est=ones(nb,1);
DF_Q_est=ones(nb,1);

%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Inf
for i=1:nl
    if prat(i)==0
        prat(i)=Inf;
    end
end