clc
clear

% Monte Carlo
E= 0.02;
V=12;
R_inf= 2*(1-E); 
R_sup= 2*(1+E);
n_s=1000000;
minimo = inf;
maximo = -inf;

for k=1:n_s
nrand = rand;
if nrand < minimo
   minimo = nrand;
end
if nrand > maximo
   maximo = nrand;
end

R(:,k)= R_inf +(R_sup - R_inf)*nrand;
I(:,k)= V/R(k);
end


I_inf=min(I);
I_max=max(I);
I_MC=[I_inf I_max];
figure(1)
boxplot(I);

figure(2)
hist(I)


% % AI (Matematica Intervalar)
% R=infsup(1.96,2.04);
% I_AI=V/R;
% 
% 
% 
% disp('----------------------------------')
% disp('        Corrente Intervalar ')
% disp('----------------------------------')
% fprintf('%s'   , '| MC')
% fprintf('%s'    ,'|  ');
% fprintf('[%4.4f; %4.4f]' ,I_MC);
% fprintf('%s \n' ,' |  ');
% fprintf('%s'   , '| AI')
% fprintf('%s'    ,'|  ');
% fprintf('[%4.4f; %4.4f]' ,[inf(I_AI) sup(I_AI)]);
% fprintf('%s \n' ,' |  ');
% disp('---------------------------------')

