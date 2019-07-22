%logistics回归（人口模型）
x=[557 563 568 573 578 583 588 593 598 603]';
t=[1:10]';
t0=t(1);
x0=x(1);
fun=@(cs,td)cs(1)./(1+(cs(1)/x0-1)*exp(-cs(2)*(td-t0)));%S型曲线方程 这里为增长速率越来越小
%最后为零 故函数一直增长直到固定值（不会下降）
cs=lsqcurvefit(fun,rand(2,1),t(2:end),x(2:end),zeros(2,1))%确定参数
year=[1:10];%%横坐标范围
xhat=fun(cs,year);
plot(year,xhat);


clear model data parama options
%data.ydata(:,1) = (1:1:13)';
%data.ydata(:,2) = [9821	11117 11958	14656 17986	22909 14909	13118 13456	13249 8509 8536	9557]';
data.ydata(:,1) = (1:1:8)';
data.ydata(:,2) = [22909 14909	13118 13456	13249 8509 8536	9557]'/10000;
l=length(data.ydata(:,1));
L=1:l+1;
cs =[798.4097,0.0321];
po=cs(1)./(1+(cs(1)/557-1)*exp(-cs(2)*(L-1)));



k00 = [0.1,0.5,0,0,0,0]';
[k0,ss0] = fminsearch(@YGs,k00,[],data)
mse = ss0/(length(data.ydata)-5)

% figure(1); clf
% y0=[k0(5) ,k0(6)];
% [t,y] = ode45(@YGmodel,L,y0,[],k0(1:4));
% y=y.*po';
% plot(t,y(:,1),':.',t,y(:,2),'--','linewidth',2.0)
% hold on
% fi=k0(6)*557+cumsum(data.ydata(:,2));
% plot(L(2:l+1),fi,'*')
% hold off
% title('Data and fitted model');



k0 = [0.36,0.36,0.3,0.02,0.05,0.05]';
params = {
    {'q',  k0(1), 0, 0.72,  k0(1), Inf,   1,      0}
    {'\lambda_i',  k0(2), 0, 1,  k0(2), Inf,   1,      0}
    {'\lambda_e',  k0(3), 0, 0.6,  k0(3), Inf,   1,      0}
    {'\mu',  k0(4), 0, 0.04,  k0(3), Inf,   1,      0}
    {'e_1(0)', k0(5) , 0.02,0.4,  k0(4),   0.1,   1,      1}
    {'e_2(0)', k0(6) , 0.02,0.4,  k0(4),   0.1,   1,      1}
    };

model.ssfun = @YGs;
model.sigma2 = mse ;

options.nsimu = 1000000;
options.updatesigma = 1;

[results,chain,s2chain] = mcmcrun(model,data,params,options);


figure(2); clf
mcmcplot(chain,[],results,'chainpanel')

figure(3); clf
mcmcplot(sqrt(s2chain),[],[],'dens',2)
title('error std')
chainstats(chain,results)

figure(4); clf
y0=results.mean(5:6);
[t,y] = ode45(@YGmodel,L,y0,[],results.mean(1:4));
y=y.*po';
plot(t,y(:,1),':.',t,y(:,2),'--','linewidth',2.0)
hold on
fi=results.mean(6)*557+cumsum(data.ydata(:,2));
plot(L(2:l+1),fi,'*')
hold off
title('Data and fitted model');

figure(5); clf
res = YGpre(data,results.mean);
data.ydata(:,3) = fi;
data.ydata(:,2) = res(:,1);
out = mcmcpred(results,chain,[],data,@YGpre,1024);
h = mcmcpredplot(out,data,1);


datd = load('Dyg.dat');
date = load('Eyg.dat'); 

figure(5); clf
subplot(1,3,1)
plot(1:8,date(:,1),'o-','linewidth',1.5,'color','black')
box off
xlim([1,8]);
xlabel('Number of Variables');ylabel(' Importance');
title('Variable Importance(Random Forest)');

subplot(1,3,2)
plot(1:8,date(:,2),'o-','linewidth',1.5,'color','black')
box off
xlim([1,8]);
xlabel('Number of Variables');ylabel(' Importance');
title('Variable Importance(Random Forest)');

sr=partialcorr(datd,'Type','Spearman');
subplot(1,3,3)
bar(sr(2:9,1))
box off
colormap(white)
xlim([1,8]);
xlabel('Number of Variables');ylabel(' Importance');
title('Variable Importance(PRCCs)');
 









