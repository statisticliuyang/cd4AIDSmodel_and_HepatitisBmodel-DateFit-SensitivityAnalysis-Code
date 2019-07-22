function ss = AIDSCMs(k,data)
y00=data.ydata(1,:)/9194;
y0=[y00,k(1:4),1-sum(k(1:4))-sum(y00)];
   
l=length(data.ydata(:,1));
L=1:l;

[~,y]=ode45(@AIDSCM,L,y0,[],k(5:15));
y=y*9194;

ssp=zeros(4,1);
for j = 1:4
  ssp(j) = sum((y(2:l,j)-data.ydata(2:l,j)).^2 );
end
ss = sum(ssp);
end