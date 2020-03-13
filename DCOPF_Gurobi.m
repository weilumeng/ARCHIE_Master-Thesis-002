clear;
IEEE33bus_modified;
ptdf_matrix %Generation Shift Factor/Power Transfer Distribution Matrix
        
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
beq=-sum(busdata(:,3));

A1=PTDF*Ag;
       
pd=busdata(:,3);
        
b1=prat(1:length(branchdata(:,1)))+PTDF*pd;
b2=prat(1:length(branchdata(:,1)))-PTDF*pd;
A=[A1; -A1];
b=[b1;b2];

lb=gendata(:,10);
ub=gendata(:,9);

[x,fval,exitflag,output,lambda] = gurobi_linprog(f,A,b,Aeq,beq,lb,ub)


genergy=lambda.eqlin
conjcost=(lambda.ineqlin'*[-PTDF ; PTDF])'
lmp=genergy+conjcost
lineflow=PTDF*(Ag*x-pd)

