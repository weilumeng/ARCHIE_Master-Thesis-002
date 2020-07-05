%% Load case
define_constants;                    %needed for the mpc function of matpower       

mpc = ext2int(loadcase('case118m')); %casedata for matpower

%% Acquiring Data

nl=size(branchdata(:,1),1);         %Number of lines/branches
nb=size(busdata(:,1),1);            %Number of buses
ng=size(gendata(:,1),1);            %Number of generators
status = branchdata(:,11);          %Checks for active/inactive lines
fb=branchdata(:,1);                 %From Bus ...
tb=branchdata(:,2);                 %To Bus ...
iline = [(1:nl)'; (1:nl)'];


prat=branchdata(:,6);               %Rated line flow limit

slack = find(busdata(:, 2) == 3);   %finding the reference/slack bus
noslack = find((1:nb)' ~= slack);
noslack2= find((1:(2*nb))' ~= slack);
noslack3 =[noslack;nb+noslack];

vbase=busdata(:,8);
vmin=busdata(:,13);                 %max voltage in p.u.
vmax=busdata(:,12);                 %min voltage in p.u.

%% Tap Ratios
tap=ones(nl,1);
i=find(branchdata(:,9));
tap(i)=branchdata(i,9);
tap = tap .* exp(1j*pi/180 * branchdata(:, 10));

G= status.*branchdata(:,3)./(branchdata(:,4).^2+branchdata(:,3).^2);
B= -status.*branchdata(:,4)./(branchdata(:,3).^2+branchdata(:,4).^2);
B= B ./ tap;
G= G ./ tap;


%% Ybus created with matpower function
Ybus = makeYbus(mpc);
Gbus = real(Ybus);
Bbus = imag(Ybus);

%% Constructing the Bus Matrices

GP = Gbus;
BD = diag(sum(Bbus));             %shunt elements
GD = diag(sum(Gbus));
BP = -Bbus + BD;
BQ = -Bbus;
GQ = -Gbus + GD;                       %GQ approxiately equals -Gbus


%% GSF Matrix Construction

Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);

Bf = sparse(iline, [fb; tb], [B; -B], nl, nb);
Gf = sparse(iline, [fb; tb], [G; -G], nl, nb);



%% Constructing the Inverse of C
C=[BP GP;GQ BQ];

Cinverse=(C(noslack2,noslack2)\speye(2*nb-1));                         %removes 1 row and column for transmission networks
%Cinverse=(C(noslack3,noslack3)\speye(size(C(noslack3,noslack3))));    %removes 2 rows and columns for distribution networks


%% Pseudo inversion with 1 row and column removed
%%while it inverts the matrix the Pseudo Inverse of C is unprecise and shouldnt be used
%%which can easily be checked by Cinverse*C which should be close to the
%%identity matrix and fails to do so

%Cinverse=pinv(full(C(noslack2,noslack2)));                              


% inverted C --> X filled with zeros at slack bus row and column
X=zeros(2*nb);

X(noslack2,noslack2)=Cinverse;                                      %1 removed row/column filled with zeros
%X(noslack3,noslack3)=Cinverse;                                     %2 removed rows/columns filled with zeros

%% Constructing new partial GSF matrices 
GSF_PP =  full( Gf *  X(nb+1:end,1:nb)      -Bf * X(1:nb,1:nb));
GSF_PQ =  full( Gf *  X(nb+1:end,nb+1:end)  -Bf * X(1:nb,nb+1:end));
GSF_QP =  full(-Gf *  X(1:nb,1:nb)          -Bf * X(nb+1:end,1:nb));
GSF_QQ =  full(-Gf *  X(1:nb,nb+1:end)      -Bf * X(nb+1:end,nb+1:end));


%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Inf
for i=1:nl
    if prat(i)==0
        prat(i)=Inf;
    end
end