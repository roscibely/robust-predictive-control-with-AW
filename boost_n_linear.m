% FUNÇÃO REPRESENTATIVA DO CONVERSOR
function [A, B, C]=boost_n_linear(Pot,Dcycle,Vo,L,Co,Rco)
% Modelo do conversor boost não linear contínuo no tempo
% Usar funções de resolução por equação diferencial - RK4 por exemplo
% Formulação do no modelo no espaço de estados
% [A, B, C]=boost_n_linear(Pot,Dcycle,Vo,L,Co,Rco)
% 
% x_dot = A.x+B.Vg;
% y = C.x;
% 
% Usando resolução por RK4
% 
%     % Modelo RK4
%     k1=A*x(:,i)+B*Vg(i);
%     k2=A*(x(:,i)+0.5*Ts*k1)+B*Vg(i);
%     k3=A*(x(:,i)+0.5*Ts*k2)+Bp*Vg(i);
%     k4=A*(x(:,i)+Ts*k3)+B*Vg(i);
%     x(:,i+1)=x(:,i)+(Ts/6)*(k1+2*(k2+k3)+k4);
%     % Equação da saida
%     y(i)=C*x(:,i); 
% 
%     % sendo i o passo de iteração do loop e Ts o tempo de amostragem

    
    
Ro=Vo*Vo/Pot;
Dcycle1=1-Dcycle;


RcoRo=Rco*Ro/(Rco+Ro);
R1=Dcycle1*Dcycle1*Ro+Dcycle*Dcycle1*RcoRo;


%Comportamento do sistema quando;
%A) a chave estiver ligada
A1=[0 0;0 -1/(Co*(Rco+Ro))];
B1=[1/L 0]'; 
C1=[0 Ro/(Ro+Rco)];

%B) a chave estiver desligada
A2=[-(RcoRo)/L  -Ro/(L*(Ro+Rco)); Ro/(Co*(Ro+Rco))  -1/(Co*(Ro+Rco))];
B2=B1; C2=[(RcoRo) Ro/((Rco+Ro))];  

%matrizes linearizadas
A=A1*Dcycle+A2*Dcycle1;
B=B1*Dcycle+B2*Dcycle1;
C=C1*Dcycle+C2*Dcycle1;