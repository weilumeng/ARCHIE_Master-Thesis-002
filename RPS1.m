        %Program to calualte locational marginal price for the problem
        %given in Sahidapur book page 395
        %This program calculates the LMP using DCOPf model
        %Program developed by Dr.Arunachalam Sundaram
        clear
        basemva = 100;  accuracy = 0.001; maxiter = 50;



        gexdata=[1 10 600 1 0 600
                 2 20 200 3 0 300 ];
        dexdata=[1 25 320 3 0 320
                 2 35 180 2 0 180];
     
        %        No  code Mag.    Degree  MW    Mvar  MW  Mvar Qmin Qmax     Mvar
        busdata=[1   1    1.0    0.0     0.0   0.0    0     0.0    0   0       0
                 2   0    1.0    0.0     180   0.0    0.0   0.0    0   0       0
                 3   2    1.0    0.0     320   0.0    0.0   0.0    0   0       0
                ];
        linedata =[...
                 1            2      0       0.25         0     1          200
                 2            3      0       0.25         0     1          200
                 1            3      0       0.25         0     1          200];

        slack=1;
        mkgsfmatrix %Program to calculate Generation Shift Factor Matrix
        nbus=length(busdata(:,1));
        nd=length(dexdata(:,1));
        ng=length(gexdata(:,1));
        prat=linedata(:,7);
        Ag=zeros(nbus,ng);
        Ad=zeros(nbus,nd);
        for i=1:nbus
            for j=1:ng
                if gexdata(j,4)==i
                    Ag(i,j)=1;
                end
            end
        end

        for i=1:nbus
            for j=1:nd
                if dexdata(j,4)==i
                    Ad(i,j)=1;
                end
            end
        end

        f=gexdata(:,2)' 
        Aeq=[-1 -1] ;
        Beq=-500;

        A1=sf*Ag ;


        pd=busdata(:,5);
        B1=prat(1:3)+sf*pd
        B2=prat(1:3)-sf*pd
        A=[A1; -A1];
        B=[B1;B2];

        lb=gexdata(:,5) ;
        ub=gexdata(:,6);



        [x,fval,exitflag,output,lamda] = linprog(f,A,B,Aeq,Beq,lb,ub)
        genergy=lamda.eqlin
        conjcost=(lamda.ineqlin'*[-sf ; sf])'
        lmp=genergy+conjcost
        lineflow=sf*(Ag*x-pd)

        