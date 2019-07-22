function ymod = AIDSMCpreD(data,k)
t=data.ydata(:,1);
y0=[data.ydata(1,2:5)'./9194 ; k(1:4) ; (1-sum(k(1:4))-sum(data.ydata(1,2:5)/9194))];
[~,y]=ode45(@AIDSCM,t,y0,[],k(5:15));
ymod=y(:,5:8)*9194;
end