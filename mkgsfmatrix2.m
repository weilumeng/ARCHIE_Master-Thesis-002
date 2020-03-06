% program to calculate generation shift factor matrix developed by
% Dr.Arunachalam Sundaram

NL=length(linedata(:,1));
nbus=length(busdata(:,1));
sn=linedata(:,1);
rn=linedata(:,2);
r=linedata(:,3);
X=linedata(:,4);
prat=linedata(:,7);
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

%----removing last row and column we get [B']--------------
n=nbus-1;

for i=1:n
    for j=1:n
        Bp(i,j)=y(i,j);
    end
end
Bp;
Bbus=inv(Bp);


Xbus=[zeros(size(Bbus,1),1) Bbus];
Xbus_new=[zeros(1,size(Xbus,2)); Xbus];


Bline=zeros(NL,nbus);
for i=1:NL
    Bline(i,sn(i))=1/X(i);
    Bline(i,rn(i))=-1/X(i);
end
  
ptdf=Bline*Xbus_new;

%xx=zeros(NL,NL);
%node=zeros(NL,nbus);
%for i=1:NL
%    xx(i,i)=1/X(i);
 %   sn(i,1)=linedata(i,1);
  %  rn(i,1)=linedata(i,2);
   % node(i,sn(i,1))=1;
    %node(i,rn(i,1))=-1;
%end
%node1=node;
%node(:,ns)=node(:,nbus);
%node(:,nbus)=[];

%sf=xx*node*Bbus;
%sf=[sf zeros(NL,1)];
%sf(:,nbus)=sf(:,ns);
%sf(:,ns)=zeros;
%sf