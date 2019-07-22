function ydot = YGmodel(t,y,k)
q=k(1);
lambda_i=k(2);
lambda_e=k(3);
beta_2=0.1761/10000;
beta_1=1.5*beta_2;
mu=k(4);
ydot = [
     q*(1-y(1)-y(2))*(y(1)*lambda_i+y(2)*lambda_e)-y(1)*(beta_1+mu);
     (1-q)*(1-y(1)-y(2))*(y(1)*lambda_i+y(2)*lambda_e)-y(2)*beta_2;
];
end

