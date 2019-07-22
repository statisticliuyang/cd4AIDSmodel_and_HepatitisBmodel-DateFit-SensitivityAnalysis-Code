function ss = YGs(k,data)
y0=[k(5),k(6)];
l=length(data.ydata(:,1));
L=1:l+1;
cs =[798.4097,0.0321];
po=cs(1)./(1+(cs(1)/557-1)*exp(-cs(2)*(L-1)));
[~,y]=ode45(@YGmodel,L,y0,[],k(1:4));
y=y(:,2).*po';
y=y-[0,y(1:l)']';
ss=sum((y(2:l+1)-data.ydata(:,2)).^2);
end
