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
Gs = busdata(:,5);                  %Conductance Shunt
Bs = busdata(:,6);                  %Susceptance Shunt
prat=branchdata(:,6);               %Rated line flow limit
vmin=busdata(:,13);                 %max voltage in p.u.
vmax=busdata(:,12);                 %min voltage in p.u.
slack = find(busdata(:, 2) == 3);   %Finding the reference/slack bus
noslack = find((1:nb)' ~= slack);
noslack2= find((1:(2*nb))' ~= slack);
iline = [(1:nl)'; (1:nl)'];
vbase=busdata(:,8);



%% Tap Ratios
tap=ones(nl,1);
i=find(branchdata(:,9));
tap(i)=branchdata(i,9);
tap = tap .* exp(1j*pi/180 * branchdata(:, 10));

G=status.*branchdata(:,3)./(branchdata(:,4).^2+branchdata(:,3).^2);
B=-status.*branchdata(:,4)./(branchdata(:,3).^2+branchdata(:,4).^2);
B= B ./ tap;
G= G ./ tap;

Y_res= status ./ (branchdata(:,3) + 1j*branchdata(:,4));

%%YBUS
Ys = status ./ (branchdata(:, 3) + 1j * branchdata(:, 4));
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;


Ysh = (busdata(:, 5) + 1j * busdata(:, 6)) / baseMVA; %% vector of shunt admittances

Yf = sparse(iline, [fb; tb], [Yff; Yft], nl, nb);
Yt = sparse(iline, [fb; tb], [Ytf; Ytt], nl, nb);

%% build Ybus
Ybus = sparse([fb;fb;tb;tb], [fb;tb;fb;tb], [Yff;Yft;Ytf;Ytt], nb, nb) + ... %% branch admittances
            sparse(1:nb, 1:nb, Ysh, nb, nb);        %% shunt admittance

%% GSF Matrix Construction

Cft = sparse(iline, [fb;tb], [ones(nl, 1); -ones(nl, 1)], nl, nb);

Bf = sparse(iline, [fb; tb], [B; -B], nl, nb);
Gf = sparse(iline, [fb; tb], [G; -G], nl, nb);
Y_series = sparse(iline, [fb; tb], [Y_res; -Y_res], nl, nb);


Bf_new=imag(Y_series);
Gf_new=real(Y_series);


%% Constructing the Bus Matrices
Gbus = real(Ybus);
Bbus = imag(Ybus);

GP = Gbus;
BD = diag(sum(Bbus));             %shunt elements
BP = -Bbus + BD;
BQ = -Bbus;
GQ = -Gbus;                       %GQ approxiately equals -Gbus

%% Constructing the Inverse of C
C=[BP GP;GQ BQ];
Cinverse=C(noslack2,noslack2)\eye(2*nb-1);

% inverted C --> Xmat filled with zeros at slack bus row and column
X=zeros(2*nb,2*nb);
X(noslack2,noslack2)=Cinverse;


%% Constructing new partial GSF matrices 
GSF_PP =  full( Gf *  X(nb+1:end,1:nb)      -Bf * X(1:nb,1:nb));
GSF_PQ =  full( Gf *  X(nb+1:end,nb+1:end)  -Bf * X(1:nb,nb+1:end));
GSF_QP =  full(-Gf *  X(1:nb,1:nb)          -Bf * X(nb+1:end,1:nb));
GSF_QQ =  full(-Gf *  X(1:nb,nb+1:end)      -Bf * X(nb+1:end,nb+1:end));


%Checking for unrestricted line flows - if line flow at line i is  
%unrestricted (=0) set it to Inf
for i=1:length(prat)
    if prat(i)==0
        prat(i)=Inf;
    end
end