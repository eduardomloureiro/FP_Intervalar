clc, clear, close all
format shortg
tic


Fluxo_Potencia_Determ;

% Solução Determinística
solucao_V = V;
solucao_Ang = teta;

% Jacobiana Determinística
Jacob;

% Incerteza
incerteza = 0.05;

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
                    
% Incerteza associada as cargas

P_incert = -(DADOS_BARRA(:,6)*(1-incerteza));
Q_incert = -(DADOS_BARRA(:,7)*(1-incerteza));

% Derivada Primeira

derivada_prim = inv(Jacob)*([P_incert;Q_incert]);

% % Atualização das cargas na matriz DADOS_BARRA
% 
% DADOS_BARRA(:,6) = DADOS_BARRA(:,6)*(1-incerteza);
% DADOS_BARRA(:,7) = DADOS_BARRA(:,7)*(1-incerteza);
% 
% delta = abs(incerteza)*ones(length(DADOS_BARRA(:,6))*2,1);

% 
        incerteza = -incerteza;         % Incerteza superior no algoritmo
        
        V_int = solucao_V + derivada_prim(1:nb); % Vetor de solução de primeira ordem  
        Ang_int = solucao_Ang + derivada_prim(nb+1:end);
                
        incerteza = incerteza;         % Incerteza inferior no algoritmo
        Vsol= solucao_V + derivada_prim(1:nb); % Vetor de solução de primeira ordem  
        Int(:,2)= Vsol;     

%Int.V = [min(Int(:,2)')' max(Int(:,2)')'];
Int_V.V = [min(([Int(:,2)])')' max(([Int_V.V])')'];
Int.Ang = [min([Int(:,2)]')' max([Int(:,2)]')'];

toc; % finalização

