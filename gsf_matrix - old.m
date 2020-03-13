%Acquiring Data
NL=length(branchdata(:,1));
nbus=length(busdata(:,1));
sn=branchdata(:,1);
rn=branchdata(:,2);
r=branchdata(:,3);
X=branchdata(:,4);
prat=branchdata(:,6);
ns=slack;

%Checking for unrestricted line flows
%if line flow at line i is unrestricted (=0) set it to Pmax of the slack bus
%generator
for i=1:length(prat)
    if prat(i)==0
        prat(i)=gendata(ns,9);
    end
end

%Creating the Admitance Matrix Y / Susceptance Matrix
y1=zeros(nbus,nbus);
for i=1:NL
    y1(sn(i),rn(i))=-1/X(i);
    y1(rn(i),sn(i))=-1/X(i);
    y1(sn(i),sn(i))=y1(sn(i),sn(i))+1/X(i);
    y1(rn(i),rn(i))=y1(rn(i),rn(i))+1/X(i);
end
y=y1;

%removing row and column corresponding to slack bus
y(ns,:)=[];
y(:,ns)=[];

%inverting y to form Bbus Matrix
Bbus= y \ eye(size(y));

%Adding zeros corresponding to slack bus row and column
Xbus=[zeros(size(Bbus,1),1) Bbus];
Xbus_new=[zeros(1,size(Xbus,2)); Xbus];

%Creating Bline Matrix
Bline=zeros(NL,nbus);
for i=1:NL
    Bline(i,sn(i))=1/X(i);
    Bline(i,rn(i))=-1/X(i);
end

%Calculating the GSF/PTDF Matrix
ptdf=Bline*Xbus_new;

