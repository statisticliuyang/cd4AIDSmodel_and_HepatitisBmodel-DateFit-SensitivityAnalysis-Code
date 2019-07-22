datah = load('dat.dat');
data.ydata(:,1) = datah(:,1)/10000;
data.ydata(:,2) = datah(:,2)/10000;
data.ydata(:,3) = datah(:,3)/10000;
data.ydata(:,4) = datah(:,4)/10000;
pp=9194;
y00=data.ydata(1,:)/pp;

k0=[0.00004,0.00008,0.00009,0.00008,0.95,0.725,0.725,0.1,0.5,0.5,0.004,0.025,0.7,0.2,0.2];
params = {
    {'D_{10}',  k0(1), 0, 0.00008,  NaN, Inf,   1,      0}
    {'D_{20}',  k0(2), 0, 0.00016,  NaN, Inf,   1,      0}
    {'D_{30}',  k0(3), 0, 0.00018,  NaN, Inf,   1,      0}
    {'D_{40}',  k0(4), 0, 0.00016,  NaN, Inf,   1,      0}
    {'q_0',  k0(5), 0, 1,  NaN, Inf,  1,      0}
    {'\lambda_{i1}',  k0(7), 0.45, 1,  NaN, Inf,   1,      0}
    {'\lambda_{e1}',  k0(6), 0.45, 1,  NaN, Inf,   1,      0}
    {'f_1',  k0(8), 0, 0.2,  NaN, Inf,   1,      0}
    {'f_2',  k0(9), 0.2, 0.6,  NaN, Inf,   1,      0}
    {'qua',  k0(10), 0.4, 0.6,  NaN, Inf,   1,      0}
    {'qub',  k0(11), 0, 0.008,  NaN, Inf,   1,      0}
    {'quc',  k0(12), 0, 0.05,  NaN,Inf,   1,      0}
    {'qda',  k0(13), 0.4, 1,  NaN, Inf,   1,      0}
    {'qdb',  k0(14), 0.1, 0.3,  NaN, Inf,   1,      0}
    {'qdc',  k0(15), 0.1, 0.3,  NaN, Inf,   1,      0}
    };

model.ssfun = @AIDSCMs;
model.sigma2 = 0.01 ;

options.nsimu = 2000000;
options.updatesigma = 1;

[results,chain,s2chain] = mcmcrun(model,data,params,options);
chainstats(chain,results)

figure(1); clf
mcmcplot(chain,[],results,'chainpanel')



figure(2); clf
y0=[y00,results.mean(1:4) ,1-sum(y00)-sum(results.mean(1:4))];
[t,y] = ode45(@AIDSCM,1:length(data.ydata(:,1)),y0,[],results.mean(5:15));
y=y*pp;
sum(y(10,1:8))
sum(y(10,1:4))
sum(y(10,5:8))
for i=1:4
    subplot(2,3,i)
    plot(t,y(:,i),'k--',t,data.ydata(:,i),'ko','linewidth',1.5);
    box off
    ylim([0,2.5]);xlim([1,10.2]);
    xlabel('time(year)');ylabel('infected persons(ten thousand)');
    legend(['fitting curve of U_',num2str(i)],['real data of U_',num2str(i)],'Location','NorthWest');
    title(['fitting result of U_',num2str(i)])
end
subplot(2,3,5)
plot(t,y(:,5),'k:.',t,y(:,6),'k--',t,y(:,7),'k:',t,y(:,8),'k-','linewidth',1.5);
box off
    xlabel('time(year)');ylabel(' infected persons(ten thousand)');
    xlim([1,10.2]);
    legend('D_1','D_2','D_3','D_4','Location','NorthWest');
    

for i=1:4
    tt(:,i)=y(:,i)+y(:,i+4);
    dr(:,i)=y(:,i)./tt(:,i);
end
   subplot(2,3,6)
   plot(t,dr(:,1),'k:.',t,dr(:,2),'k--',t,dr(:,3),'k:',t,dr(:,4),'k-','linewidth',1.5)
    xlabel('time(year)');ylabel('diagnostic rate');
    legend('CD4\geq500','500>CD4\geq350','350>CD4\geq200','CD4<200','Location','NorthWest');
    xlim([1,10.2]);
   box off

   
figure(3); clf
l=length(data.ydata(:,1));
data2.ydata(:,1) = 1:l;
data2.ydata(:,2) = data.ydata(1:l,1);
data2.ydata(:,3) = data.ydata(1:l,2);
data2.ydata(:,4) = data.ydata(1:l,3);
data2.ydata(:,5) = data.ydata(1:l,4);
out=mcmcpred(results,chain,[],data2,@AIDSMCpre,500);
mcmcpredplot(out,data2,1);
figure(4); clf
out=mcmcpred(results,chain,[],data2,@AIDSMCpreD,500);
mcmcpredplot(out);
   
   
   
datd = load('D.dat');
date = load('E.dat'); 
figure(5); clf
subplot(1,3,1)
plot(1:32,date(:,1),'o-','linewidth',1.5,'color','black')
box off
xlim([1,33]);
xlabel('Number of Variables');ylabel(' Importance');
title('Variable Importance(Random Forest)');
subplot(1,3,2)
plot(1:32,date(:,2),'o-','linewidth',1.5,'color','black')
box off

sr=partialcorr(datd,'Type','Spearman');
xlim([1,33]);
xlabel('Number of Variables');ylabel(' Importance');
title('Variable Importance(Random Forest)');
subplot(1,3,3)
bar(sr(2:33,1))
box off
colormap(white)
xlim([1,33]);
xlabel('Number of Variables');ylabel('Importance');
title('Variable Importance(PRCCs)');

figure(6); clf
subplot(1,2,1)
y0=[y00,results.mean(1:4) ,1-sum(y00)-sum(results.mean(1:4))];
[t,y] = ode45(@AIDSCM,1:length(data.ydata(:,1)),y0,[],results.mean(5:15));
y1=y*9194;
th=[results.mean(5),0.4,results.mean(7:15)];
[t,y] = ode45(@AIDSCM,1:length(data.ydata(:,1)),y0,[],th);
y2=y*9194;
th=[results.mean(5),0.3,results.mean(7:15)];
[t,y] = ode45(@AIDSCM,1:length(data.ydata(:,1)),y0,[],th);
y3=y*9194;
plot(t,y1(:,1),'k-',t,y2(:,1),'k--',t,y3(:,1),'k:.',t,y1(:,2),'k-',t,y2(:,2),'k--',t,y3(:,2),'k:.','linewidth',2.0)
xlim([1,10.2]);
legend('\lambda_b=0.101','\lambda_b=0.05','\lambda_b=0.15','Location','NorthWest');
title('fitting result of U1 and U2');
xlabel('time(year)');ylabel(' infected persons(ten thousand)');
subplot(1,2,2)
plot(t,y1(:,3),'k-',t,y2(:,3),'k--',t,y3(:,3),'k:.',t,y1(:,4),'k-',t,y2(:,4),'k--',t,y3(:,4),'k:.','linewidth',2.0)
xlim([1,10.2]);
legend('\lambda_b=0.101','\lambda_b=0.05','\lambda_b=0.15','Location','NorthWest');
title('fitting result of U3 and U4');
xlabel('time(year)');ylabel(' infected persons(ten thousand)');









