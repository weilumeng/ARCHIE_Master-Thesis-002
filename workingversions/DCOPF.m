clear;
IEEE33bus;
ptdf_matrix %calculate Generation Shift Factor Matrix
        
       
nd=length(busdata(:,3));
ng=length(gendata(:,1));
        
Ag=zeros(nb,ng);

for i=1:nb
    for j=1:ng
        if gendata(j,1)==i
            Ag(i,j)=1;
        end
    end
end


f=gencostdata(:,6)';
Aeq=-ones(1,length(gendata(:,1)));
Beq=-sum(busdata(:,3));

A1=H*Ag;
       
pd=busdata(:,3);
        
B1=prat(1:length(branchdata(:,1)))+H*pd;
B2=prat(1:length(branchdata(:,1)))-H*pd;
A=[A1; -A1];
B=[B1;B2];

lb=gendata(:,10);
ub=gendata(:,9);
        

[x,fval,exitflag,output,lamda] = linprog(f,A,B,Aeq,Beq,lb,ub)
genergy=lamda.eqlin
conjcost=(lamda.ineqlin'*[-H ; H])'
lmp=genergy+conjcost
lineflow=H*(Ag*x-pd)

        