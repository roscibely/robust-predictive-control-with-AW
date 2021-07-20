function [Ah,Bh,A,B,C,D]=boost_discreto(Pot,Vg,Vo,L,Co,Rco,Ts,g,h)
% Doutorado em Engenharia Elétrica - UFC
% Função matemática do conversor boost discreto
% Aluno:Marcus Vinicius Silvério Costa
Ro=Vo*Vo/Pot;
Dcycle=1-(Vg/Vo);
Dcycle_=1-Dcycle;

RcoRo=Rco*Ro/(Rco+Ro);
R_=Dcycle_*Dcycle_*Ro+Dcycle*Dcycle_*RcoRo;

At=[-Dcycle_*RcoRo/L -Dcycle_*Ro/(L*(Rco+Ro)); 
    Dcycle_*Ro/(Co*(Rco+Ro)) -1/(Co*(Rco+Ro)) ];

Bt1=(Ro/L)*(Dcycle_*Ro+Rco)/(Rco+Ro);
Bt2=-(Ro)/(Co*(Rco+Ro));
Bt=(Vg/R_)*[Bt1 Bt2]';

Ct=[Dcycle_*RcoRo (Ro)/(Rco+Ro)];

Dt=-Vg*RcoRo/R_;

[m,n]=size(At);
[m,nb]=size(Bt);
s = expm([[At Bt]*Ts; zeros(nb,n+nb)]);
A = s(1:n,1:n);
B = s(1:n,n+1:n+nb);
% A=eye(2)+Ts*At;
% B=Ts*Bt;
% [A,B]=c2d(At,Bt,Ts);
C=Ct; D=Dt;


Ah=[A zeros(2,1);-h*C g];
Bh=[B;-h*D];