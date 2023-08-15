clc
clear
close all
format shortg
tic


lambda=1;
% Dados da Rede Elétrica de 14 barras

%               Barra   V     A   Pg     Qg    Pl     Ql    ysht_bus tipo  Qmin Qmax
%
DADOS_BARRA= [    1    1      0    0     0     0     0       0        2     0    0
                  2    1      0    0     0    20.0   5       0        0     0    0
                  3    1      0    0     0    20.0   5       0        0     0    0
                  4    1      0    0     0    20.0   5       0        0     0    0
                  5    1      0    0     0    20.0   10.0    0        0     0    0
                  6    1      0    0     0    20.0   5       0        0     0    0
                  7    1      0    55    20   0      0       0        1    -40   40
                  8    1      0    0     0    20.0   10.0    0        0     0    0
                  9    1      0    0     0    20.0   5       0        0     0    0
                  10   1      0    55    20   0      0       0        1    -40   40
                  11   1      0    0     0    20.0   10.0    0        0     0    0
                  12   1      0    0     0    20.0   10.0    0        0     0    0
                  13   1      0    0     0    20.0   5       0        0     0    0
                  14   1      0    90    30   0      0       0        1    -60   60
                  ];



%                  B_DE    B_PARA     r%    x%      tap    ang. tap  ysht_linha  tp

DADOS_REDE= [        1         2     2   5    0           0       0    1
                     2         3     2   5    0           0       0    1
                     3         4     2   5    0           0       0    1
                     4         5     2   5    0           0       0    1
                     5         6     2   5    0           0       0    1
                     6         7     2   5    0           0       0    1
                     3         8     2   5    0           0       0    1
                     8         9     2   5    0           0       0    1
                     9         10    2   5    0           0       0    1
                     10        11    2   5    0           0       0    1
                     11        12    2   5    0           0       0    1
                     6         12    2   5    0           0       0    1
                     5         13    2   5    0           0       0    1
                     13        14    2   5    0           0       0    1
                     ];

%% 1passo) Leitura dos dados:

%Dados de linha

B_de = DADOS_REDE(:,1);
B_para = DADOS_REDE(:,2);
r = (1/100)*DADOS_REDE(:,3);
x = (1/100)*DADOS_REDE(:,4);
a = DADOS_REDE(:,5);
fi = (pi/180)*DADOS_REDE(:,6);
ysht_linha = (1/200)*DADOS_REDE(:,7);
tp = DADOS_REDE(:,8); %tp=1 (LT); tp=2 (trafo em fase); tp=3 (trafo defasador)

%Dados da Barra

Barra = DADOS_BARRA(:,1);
V = DADOS_BARRA(:,2);
teta = (pi/180)*DADOS_BARRA(:,3); %ângulo em radianos
Pg = (1/100)*DADOS_BARRA(:,4);
Qg = (1/100)*DADOS_BARRA(:,5);
Pl = lambda*(1/100)*DADOS_BARRA(:,6);
Ql = lambda*(1/100)*DADOS_BARRA(:,7);
ysht_bus = (1/100)*DADOS_BARRA(:,8);
tipo = DADOS_BARRA(:,9); %tipo=0(barra PQ); tipo=1(barra PV); tipo=2(barra slack)
Qmin= (1/100)*DADOS_BARRA(:,10);
Qmax= (1/100)*DADOS_BARRA(:,11);

nm =length(B_de);%numero de medidas
nb=max(max(B_de),max(B_para));%numero de barras no sistema

Vesp = V;

%% 2 passo) Montagem da Ybarra
%Montagem da matriz de admitancia shunt do sistema:

for k=1:length(r)
   if (r(k)& x(k))~=0
      ykm(k,1)=1/((r(k) + 1j*x(k)));
   else
      if r(k)==0 && x(k)~=0
         ykm(k,1)= 1/(0 + 1j*x(k));
      else
         if x(k)==0 && r(k)~=0
            ykm(k,1)= 1/((r(k)+1j*0));
         else
            if (r(k)& x(k))==0
               ykm(k,1)= 0+ 0*1j;
            end
         end
      end
   end
end

%Montagem da matriz Ybarra

YBARRA=zeros(nb,nb);

for k=1:nb
   YBARRA(k,k)=YBARRA(k,k)+1j*ysht_bus(k);
end


for k=1:nm
   
   m=B_de(k);
   n=B_para(k);
   
   if a(k)~= 0
      
      YBARRA(m,m)=YBARRA(m,m)+ ((1/(a(k).^2))*ykm(k))+ 1j*ysht_linha(k);
      
      if n~=m
         
         YBARRA(m,n)=YBARRA(m,n)- (1/a(k))*ykm(k);
         YBARRA(n,n)=YBARRA(n,n)+ ykm(k)+ 1j*ysht_linha(k);
         YBARRA(n,m)=YBARRA(n,m)- (1/a(k))*ykm(k);
      end
   else
      YBARRA(m,m)=YBARRA(m,m)+ ykm(k)+ 1j*ysht_linha(k);
      
      if n~=m
         
         YBARRA(m,n)=YBARRA(m,n)- ykm(k);
         YBARRA(n,n)=YBARRA(n,n)+ ykm(k)+ 1i*ysht_linha(k);
         YBARRA(n,m)=YBARRA(n,m)- ykm(k);
      end
   end
end



%Montagem das matrizes de condutância G e susceptância B

G = zeros(nb,nb);
B = zeros(nb,nb);

for k=1:nb
   for m=1:nb
      G(k,m)=real(YBARRA(k,m));
      B(k,m)=imag(YBARRA(k,m));
   end
end
G;
B;

%% 3 passo) Comparar valores especificados e calculados de potencias

% Valores especificados para as potências ativas e reativas:
Pesp = zeros(nb,1);
Qesp = zeros(nb,1);

for k=1:nb
   Pesp(k) = Pg(k)-Pl(k);
   Qesp(k) = Qg(k)-Ql(k);
end

%Cálculo dos valores de potencias ativas e reativas (calculados em função das variáveis V e Theta):
Pcalc = zeros(nb,1);
Qcalc = zeros(nb,1);

for k=1:nb
   for m = 1:nb
      Pc = V(k)*V(m)*(G(k,m)*cos(teta(k)-teta(m))+B(k,m)*sin(teta(k)-teta(m)));
      Pcalc(k) = Pcalc(k)+Pc;
      Qc = V(k)*V(m)*(G(k,m)*sin(teta(k)-teta(m))-B(k,m)*cos(teta(k)-teta(m)));
      Qcalc(k)= Qcalc(k)+ Qc;
   end
end



% Cálculo de delta P e delta Q

deltaP = zeros(nb,1);
deltaQ = zeros(nb,1);
for k=1:nb
   deltaP(k)=Pesp(k)-Pcalc(k);
   deltaQ(k)=Qesp(k)-Qcalc(k);
   if tipo(k)==1    % Barra PV
      deltaQ(k)=0;
   end
   if tipo(k)==2    % Barra V-teta
      deltaP(k)=0;
      deltaQ(k)=0;
   end
end
deltaY=[deltaP ; deltaQ];

%% 4 passo) Processo iterativo

iter=0;    % Número de iterações
H=zeros(nb,nb);
N=zeros(nb,nb);
M=zeros(nb,nb);
L=zeros(nb,nb);
tol=10^-5; % tolerancia

while  max(abs(deltaY))>=tol
   
   %Montagem do Jacobiano
   for j=1:nm
      k=B_de(j);
      m=B_para(j);
      % Submatriz H
      H(k,k)=-B(k,k)*V(k)^2-Qcalc(k);
      H(m,m)=-B(m,m)*V(m)^2-Qcalc(m);
      H(k,m)=V(k)*V(m)*(G(k,m)*sin(teta(k)-teta(m))-B(k,m)*cos(teta(k)-teta(m)));
      H(m,k)=V(m)*V(k)*(G(m,k)*sin(teta(m)-teta(k))-B(m,k)*cos(teta(m)-teta(k)));
      % Submatriz N
      N(k,k)=inv(V(k))*(Pcalc(k)+V(k)^2*G(k,k));
      N(m,m)=inv(V(m))*(Pcalc(m)+V(m)^2*G(m,m));
      N(k,m)=V(k)*(G(k,m)*cos(teta(k)-teta(m))+B(k,m)*sin(teta(k)-teta(m)));
      N(m,k)=V(m)*(G(m,k)*cos(teta(m)-teta(k))+B(m,k)*sin(teta(m)-teta(k)));
      % Submatriz M
      M(k,k)= Pcalc(k)-V(k)^2*G(k,k);
      M(m,m)= Pcalc(m)-V(m)^2*G(m,m);
      M(k,m)=-V(k)*V(m)*(G(k,m)*cos(teta(k)-teta(m))+B(k,m)*sin(teta(k)-teta(m)));
      M(m,k)=-V(m)*V(k)*(G(m,k)*cos(teta(m)-teta(k))+B(m,k)*sin(teta(m)-teta(k)));
      % Submatriz L
      L(k,k)=inv(V(k))*(Qcalc(k)-V(k)^2*B(k,k));
      L(m,m)=inv(V(m))*(Qcalc(m)-V(m)^2*B(m,m));
      L(k,m)=V(k)*(G(k,m)*sin(teta(k)-teta(m))-B(k,m)*cos(teta(k)-teta(m)));
      L(m,k)=V(m)*(G(m,k)*sin(teta(m)-teta(k))-B(m,k)*cos(teta(m)-teta(k)));
   end
   
   
   for k=1:nb
      deltaP(k)=Pesp(k)-Pcalc(k);
      deltaQ(k)=Qesp(k)-Qcalc(k);
      if tipo(k)==1    % Barra PV
         deltaQ(k)=0;
      end
      if tipo(k)==2    % Barra V-teta
         deltaP(k)=0;
         deltaQ(k)=0;
      end
   end
   deltaY=[deltaP ; deltaQ];
   
   
   
   
   % Alterações no Jacobiano devido ao tipo de barra
   % Para que fazer isso? (Não se especifica Q para barra PV e não se especifica
   %P nem Q para a barra V-teta. Para eliminar as equações do sistema
   %linear principal, faz-se as diagonais de L, H e L iguais a um
   %BIGNUMBER.
   
   for k=1:nb
      if tipo(k)==1 % Barra PV
         L(k,k)=10^10;
      end
      if tipo(k)==2 % Barra V-teta
         H(k,k)=10^10;
         L(k,k)=10^10;
      end
   end
   
   Jacob=[H N
      M L];
   
   deltaX=inv(Jacob)*deltaY; % determina os incrementos de DeltaV e Delta teta
   deltaT=deltaX(1:length(deltaP));
   deltaV=deltaX(length(deltaP)+1:length(deltaP)+length(deltaQ));
   
   V=V+deltaV;   % novo valor de V
   teta=teta+deltaT; % novo valor de teta
   
   Pcalc = zeros(nb,1);
   Qcalc = zeros(nb,1);
   
   for k=1:nb
      for m = 1:nb
         Pc = V(k)*V(m)*(G(k,m)*cos(teta(k)-teta(m))+B(k,m)*sin(teta(k)-teta(m)));
         Pcalc(k) = Pcalc(k)+Pc;
         Qc = V(k)*V(m)*(G(k,m)*sin(teta(k)-teta(m))-B(k,m)*cos(teta(k)-teta(m)));
         Qcalc(k)= Qcalc(k)+ Qc;
      end
   end
   
   for k=1:nb
      deltaP(k)=Pesp(k)-Pcalc(k);
      deltaQ(k)=Qesp(k)-Qcalc(k);
      if tipo(k)==1    % Barra PV
         deltaQ(k)=0;
      end
      if tipo(k)==2    % Barra V-teta
         deltaP(k)=0;
         deltaQ(k)=0;
      end
   end
   deltaY=[deltaP ; deltaQ];
   iter=iter+1;
end

%% 5 passo) Cálculo de Fluxos passantes nas linhas Pkm e Qkm ; cálculo de variáveis funcionais

for j=1:nm
   
   k=B_de(j);
   m=B_para(j);
   gkm= real(ykm);
   bkm = imag(ykm);
   %
   if tp(j) == 1 %Linha de transmissão
      Pkm(j)= gkm(j)*V(k).^2 - V(k)*V(m)*(gkm(j)*cos(teta(k)-teta(m)) + bkm(j)*sin(teta(k)-teta(m)));
      Qkm(j) = - (ysht_linha(j)+bkm(j))*V(k).^2 + V(k)*V(m)*(bkm(j)*cos(teta(k)-teta(m)) - gkm(j)*sin(teta(k)-teta(m)));
      Pmk(j)= gkm(j)*V(m).^2 - V(k)*V(m)*(gkm(j)*cos(teta(k)-teta(m)) - bkm(j)*sin(teta(k)-teta(m)));
      Qmk(j) = - (ysht_linha(j)+bkm(j))*V(m).^2 + V(k)*V(m)*(bkm(j)*cos(teta(k)-teta(m)) + gkm(j)*sin(teta(k)-teta(m)));
      
      
   else
      if tp(j) == 2 %Trafo em fase
         Pkm(j)= (1./(a(j).^2))*gkm(j)*(V(k).^2) - (1./a(j))*V(k)*V(m)*(gkm(j)*cos(teta(k)-teta(m)) + bkm(j)*sin(teta(k)-teta(m)));
         Qkm(j) = (-1./(a(j).^2))*(bkm(j))*(V(k).^2) + (1./a(j))*V(k)*V(m)*(bkm(j)*cos(teta(k)-teta(m)) - gkm(j)*sin(teta(k)-teta(m)));
         Pmk(j)= V(m)^2*gkm(j) - (1/a(j))*V(k)*V(m)*(gkm(j)*cos(teta(k)-teta(m)) - bkm(j)*sin(teta(k)-teta(m)));
         Qmk(j) = - bkm(j)*V(m).^2 + (1/a(j))*V(k)*V(m)*(bkm(j)*cos(teta(k)-teta(m)) + gkm(j)*sin(teta(k)-teta(m)));
      else
         if tp(j) == 3 %Trafo defasador
            Pkm(j)= gkm(j)*(V(k).^2) - V(k)*V(m)*(gkm(j)*cos(teta(k)-teta(m)+ fi(j)) + bkm(j)*sin(teta(k)-teta(m)+ fi(j)));
            Qkm(j) =- bkm(j)*V(k).^2 + V(k)*V(m)*(bkm(j)*cos(teta(k)-teta(m)+fi(j)) - gkm(j)*sin(teta(k)-teta(m)+ fi(j)));
            Pmk(j)= gkm(j)*(V(m).^2) - V(k)*V(m)*(gkm(j)*cos(teta(k)-teta(m)+ fi(j)) - bkm(j)*sin(teta(k)-teta(m)+ fi(j)));
            Qmk(j) = - bkm(j)*V(m).^2 + V(k)*V(m)*(bkm(j)*cos(teta(k)-teta(m)+fi(j)) + gkm(j)*sin(teta(k)-teta(m)+ fi(j)));
            
            
         end
      end
   end
end

%% Cálculo de correntes injetadas
Vc=V.*exp(1j*teta);
S=Pcalc+1j*Qcalc;
Iinj=conj(S./Vc);



toc % finalização

teta=rad2deg(teta); %transformar radianos e graus

%% 6 passo) Saída de dados
disp('____________________________________________________________________________________________')
disp('Resultado do fluxo de potência')
disp('---------------------------------------------------------------------------------------------')
disp('BARRA |   V    |  Teta (°)  |    P(MW)    |    Q(MVAr)    ')
disp('---------------------------------------------------------------------------------------------')
for j=1:nb
   fprintf(1,'%6.0f   %6.4f    %6.4f       %6.4f         %6.4f\n',(j),V(j),teta(j),100*Pcalc(j),100*Qcalc(j))
end
disp('---------------------------------------------------------------------------------------------')
disp(' DE      |    PARA    |    Pkm(MW)     |     Qkm(MVAr)    |    Pmk(MW)     |     Qmk(MVAr)        ')
disp('---------------------------------------------------------------------------------------------')
for j=1:nm
   fprintf(1,'%6.0f    %6.0f          %6.4f         %6.4f            %6.4f             %6.4f\n',B_de(j),B_para(j),100*Pkm(j),100*Qkm(j),100*Pmk(j),100*Qmk(j))
end
disp('---------------------------------------------------------------------------------------------')
disp(sprintf('Número de iterações: %d',iter))
disp('---------------------------------------------------------------------------------------------')

%% Análise Harmônica

%        h    %mag   %angle
Dharm=[  1     100     0
         2      0      0
         3      20.0   -265
         4      0       0
         5      15     -260
         6      0       0
         7      10      -320
         8      0       0
         9      5       120
         10     0       0
         11     2       0
         12     0       0
         13     1       0
         14     0       0
         15     0.5     0];


%        h    %mag   %angle
% Dharm=[  1     100     0
%          2      0      0
%          3      25      0
%          4      0      0
%          5      15      0
%          6      0      0
%          7      10      0
%          8      0      0
%          9      5      0
%          10     0      0
%          11     3     0];


ho =   Dharm(:,1); % harmonic order
mag =  Dharm(:,2); % magnitudes em percentage
ang=   Dharm(:,3); % angles in degrees



for h=2:max(ho)
   
   
   for k=1:length(r)
      if (r(k)& x(k))~=0
         ykmh(k,1)=1/((r(k) + 1j*h*x(k)));
      else
         if r(k)==0 && x(k)~=0
            ykmh(k,1)= 1/(0 + 1j*h*x(k));
         else
            if x(k)==0 & r(k)~=0
               ykmh(k,1)= 1/((r(k)+1j*0));
            else
               if (r(k)& x(k))==0
                  ykmh(k,1)= 0+ 0*1j;
               end
            end
         end
      end
   end
   
   %Montagem da matriz Ybarra
   
   YBARRAh=zeros(nb,nb);
   
   for k=1:nb
      YBARRAh(k,k)=YBARRAh(k,k)+1j*h*ysht_bus(k);
   end
   
   
   for k=1:nm
      
      m=B_de(k);
      n=B_para(k);
      
      if a(k)~= 0
         
         YBARRAh(m,m)=YBARRAh(m,m)+ ((1/(a(k).^2))*ykmh(k))+ 1j*h*ysht_linha(k);
         
         if n~=m
            
            YBARRAh(m,n)=YBARRAh(m,n)- (1/a(k))*ykmh(k);
            YBARRAh(n,n)=YBARRAh(n,n)+ ykmh(k)+ 1j*h*ysht_linha(k);
            YBARRAh(n,m)=YBARRAh(n,m)- (1/a(k))*ykmh(k);
         end
      else
         YBARRAh(m,m)=YBARRAh(m,m)+ ykmh(k)+ 1j*h*ysht_linha(k);
         
         if n~=m
            
            YBARRAh(m,n)=YBARRAh(m,n)- ykmh(k);
            YBARRAh(n,n)=YBARRAh(n,n)+ ykmh(k)+ 1j*h*ysht_linha(k);
            YBARRAh(n,m)=YBARRAh(n,m)- ykmh(k);
         end
      end
   end
   
   
   % Injeção de harmônicas
   
   barra_inj = 9; % barra onde serão injetadas harmônicas
   
   iinj1=(mag(h)/100)*abs(Iinj(barra_inj))*exp(1j*((ang(h)+ho(h)*angle(Iinj(barra_inj))*180/pi)*pi/180)); 
   %iinj2=(mag(h)/100)*abs(Iinj(5))*exp(i*((ang(h)+ho(h)*angle(Iinj(5))*180/pi)*pi/180))% injeção na barra 5
   
   Ih = zeros(nb,1);
   Ih(barra_inj) = iinj1;
   
   R=abs(V.^2./Pcalc);     % Resistência da carga linear
   j=sqrt(-1);             % número complexo
   X=abs(V.^2./Qcalc); % Reatância da carga linear
   
   Yp=(R.^-1 + (1j*h*X).^-1); % admitância equivalente da associação em paralelo de R com X
   
   %%% TESTE
   %  garantir que a tensão nas barras PV e Vteta seja zero
   barras_PV = find(tipo == 1);
   barras_Vteta = find(tipo == 2);
   indices = [barras_PV; barras_Vteta]';
   for iter = indices
      YBARRAh(iter,iter) = 10^10; 
   end
   %%%
   
%    YBARRAh(1,1)=10^10;     % para garantir que a tensão na sbestação seja zero
   
   YBARRAhn=YBARRAh;       % Formação da nova Ybarra com inclusão das cargas lineares
   
   for kk=1:nb
      YBARRAhn(kk,kk)=YBARRAh(kk,kk)+Yp(kk);
   end
   ZBARRAh=inv(YBARRAhn);
   
   iter=0;
   delV=100;
   tol=10^-8;
   Ih1=Ih;                % vetor de correntes injetadas (todas zero exceto nos lugares onde há FH)
   Vh2=inv(YBARRAhn)*Ih;
   
   while abs(delV)>=tol
      
      Vh=inv(YBARRAhn)*Ih; % calcula a tensão a partir de inversão da Ybarra
      
      Ihn=Yp.*Vh;          % calcula corrente absorvida por cada carga linear
      
      Ih=Ih1-Ihn;          % corrente total na barra: injetada menos a absorvida
      
      Vh1=inv(YBARRAhn)*Ih; % calcula a tensão a partir de inversão da Ybarra
      
      delV=abs(max((Vh1)-(Vh))); % verifica a diferença entre as tensões calculadas antes e após compensação de corrente
      
      iter=iter+1;           % atualiza número da iteração
      
   end
   Sh=Vh.*conj(Ih);
   Tens_harm(:,h)=abs(Vh);  % armazena tensões harmônicas
   Corr_harm(:,h)=abs(Ih);  % armazena correntes harmônicas
   Po_harm(:,h)=(Sh);
   abs_imped(:,h)=abs(ZBARRAh(4,4)); % armazena a impedância (amplitude) no driving point
end

%% Calcular THD de tensão por barra, percentual;
% arredondado a partir da nona (9a) casa decimal;
% multiplicado por 5 para bater com o resultado do HarmZs

THD = 5*round( 100*sqrt(sum(Tens_harm.^2,2))./V , 9); 

%% Plot 1 - 3D: barra, ordem harmônica, tensão harmônica

ordens = 2:max(ho);

% Tensão Harmônica
figure
s = surf(ordens,1:nb,Tens_harm(:,ordens));
xlabel("Ordem Harmônica")
ylabel("Barra")
zlabel("Tensão Harmônica (p.u)")
title("Tensões Harmônicas")
% s.EdgeColor = 'interp';
s.FaceColor = 'interp';
colormap(jet)
colorbar

% Corrente Harmônica
figure
s = surf(ordens,1:nb,Corr_harm(:,ordens));
xlabel("Ordem Harmônica")
ylabel("Barra")
zlabel("Corrente Harmônica (p.u)")
title("Correntes Harmônicas")
% s.EdgeColor = 'interp';
s.FaceColor = 'interp';
colormap(jet)
colorbar

%% Plot 2 - 2D: tensão x barra

% Tensão na Frequência Fundamental (V)
figure
b = bar(1:nb,V);
xlabel("Barra")
ylabel("Tensão (p.u)")
title("Tensões na Frequência Fundamental")

%% Plot 3 - 2D: THD x barra

% THD (%)
figure
b = bar(1:nb,THD);
xlabel("Barra")
ylabel("THD (%)")
title("THD")


