% Clearing Workspace and Loading Data
clear;
IEEE33bus;
% Assessing the Data
NL=length(branchdata(:,1));
nbus=length(busdata(:,1));
sn=branchdata(:,1);
rn=branchdata(:,2);
r=branchdata(:,3);
X=branchdata(:,4);
prat=branchdata(:,6);
slack=1;
ns=slack;

y1=zeros(nbus,nbus);

%********************formation of  negative Ybus[Y]------------------
for i=1:NL
    y1(sn(i),rn(i))=-1/X(i);
    y1(rn(i),sn(i))=-1/X(i);
    y1(sn(i),sn(i))=y1(sn(i),sn(i))+1/X(i);
    y1(rn(i),rn(i))=y1(rn(i),rn(i))+1/X(i);
end
y1;
y=y1;
%------------to form Bus Suspectence matrix [B']---------------
%------shifting of Nth row to corresponding Slack row
for j=1:nbus
    y(ns,j)=y(nbus,j);
end
%-------shifting of Nth Column to corresponding Slack column-------------
for i=1:nbus
    y(i,ns)=y(i,nbus);
end
%----removing last row and column we get [B']--------------
n=nbus-1;

for i=1:n
    for j=1:n
        Bp(i,j)=y(i,j);
    end
end
Bp;
rec=inv(Bp);

xx=zeros(NL,NL);
node=zeros(NL,nbus);
for i=1:NL
    xx(i,i)=1/X(i);
    sn(i,1)=branchdata(i,1);
    rn(i,1)=branchdata(i,2);
    node(i,sn(i,1))=1;
    node(i,rn(i,1))=-1;
end
node1=node;
node(:,ns)=node(:,nbus);
node(:,nbus)=[];

sf=xx*node*rec;
sf=[sf zeros(NL,1)];
sf(:,nbus)=sf(:,ns);
sf(:,ns)=zeros;
sf;