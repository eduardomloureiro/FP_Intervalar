clc  
format short
tic
clear all

% Carregamento
lambda=1;
% Dados da Rede Elétrica 

%                 Barra   V     A   Pg     Qg    Pl     Ql    ysht_bus tipo  Qmin Qmax
               
DADOS_BARRA= [      1    1      0    0     0     0     0       0        2     0    0
                    2    1      0    0     0    20.0   5.0     0        0     0    0
                    3    1      0    0     0    20.0   5.0     0        0     0    0
                    4    1      0    0     0    20.0   5.0     0        0     0    0
                    5    1      0    0     0    20.0   10.0    0        0     0    0
                    6    1      0    0     0    20.0   5.0     0        0     0    0
                    7    1      0    55    20      0    0      0        1     -40  -40
                    8    1      0    0     0    20.0   10.0    0        0     0    0
                    9    1      0    0     0    20.0   5.0     0        0     0    0
                   10    1      0    55    20     0      0     0        1     -40  -40
                   11    1      0    0     0    20.0   10.0    0        0     0    0
                   12    1      0    0     0    20.0   10.0    0        0     0    0
                   13    1      0    0     0    20.0   5.0     0        0     0    0
                   14    1      0    90    30    0     0    0        1     -60  -60
                    ];
                
                

%                  B_DE    B_PARA  r%   x%  tap    ang. tap  ysht_linha  tp
            
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
 
 B_de=DADOS_REDE(:,1);
 B_para=DADOS_REDE(:,2);
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
           ykm(k,1)=1/((r(k) + i*x(k)));
       else
           if r(k)==0 & x(k)~=0
               ykm(k,1)= 1/(0 + i*x(k));
           else
               if x(k)==0 & r(k)~=0
                   ykm(k,1)= 1/((r(k)+i*0));
               else
                   if (r(k)& x(k))==0
                       ykm(k,1)= 0+ 0*i;
                   end
               end
           end
       end
   end
 
 %Montagem da matriz Ybarra
 
 YBARRA=zeros(nb,nb);

 for k=1:nb
    YBARRA(k,k)=YBARRA(k,k)+i*ysht_bus(k);
end

   
for k=1:nm
    
    m=B_de(k);
    n=B_para(k);
    
    if a(k)~= 0
   
    YBARRA(m,m)=YBARRA(m,m)+ ((1/(a(k).^2))*ykm(k))+ i*ysht_linha(k);
    
        if n~=m
      
        YBARRA(m,n)=YBARRA(m,n)- (1/a(k))*ykm(k);
        YBARRA(n,n)=YBARRA(n,n)+ ykm(k)+ i*ysht_linha(k);
        YBARRA(n,m)=YBARRA(n,m)- (1/a(k))*ykm(k);
        end
    else
        YBARRA(m,m)=YBARRA(m,m)+ ykm(k)+ i*ysht_linha(k);
    
        if n~=m
      
        YBARRA(m,n)=YBARRA(m,n)- ykm(k);
        YBARRA(n,n)=YBARRA(n,n)+ ykm(k)+ i*ysht_linha(k);
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
Vc=V.*exp(i*teta);
S=Pcalc+i*Qcalc;
Iinj=conj(S./Vc);
    
teta=rad2deg(teta); %transformar radianos e graus
