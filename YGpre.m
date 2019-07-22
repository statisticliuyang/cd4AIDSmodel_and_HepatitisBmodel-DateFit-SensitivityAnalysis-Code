function ymod = YGpre(data,theta)
y0=[theta(5),theta(6)];
l=length(data.ydata(:,1));
L=1:l+1;
cs =[798.4097,0.0321];
po=cs(1)./(1+(cs(1)/557-1)*exp(-cs(2)*(L-1)));
[~,y]=ode45(@YGmodel,2:l+1,y0,[],theta(1:4));
ymod=y.*po(2:l+1)';
end