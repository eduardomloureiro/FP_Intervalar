clc
clear
close all

% Monte Carlo
num_simulacoes = 1000000;
E = 0.05;
f_inf = (1-E); 
f_sup = (1+E);
R = 2*ones(num_simulacoes,1);
V = 12;

nrand = zeros(1,num_simulacoes);

for k = 1:num_simulacoes
nrand(k) = rand;

R(k)= R(k)*(f_inf +(f_sup - f_inf)*rand);
I(:,k) = V/R(k);
end


I_inf = min(I);
I_max = max(I);
I_MC = [I_inf I_max];
figure
boxplot(I);

figure
hist(I)

figure
hist(R)

figure
hist(nrand)



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

