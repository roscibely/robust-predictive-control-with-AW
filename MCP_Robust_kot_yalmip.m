function [gammas,Fs]=MCP_Robust_kot_yalmip(A,B,Qc,Rc,umax,xk)
[m,n]=size(B{1});
L=size(A,2);
% inicialização das LMIS
gamma=sdpvar(1);
Q=sdpvar(m,m,'symmetric');
Y=sdpvar(n,m);
x_k=sdpvar(m,1);
LMIs=[Q>=0];
for i=1:L
M11=Q;
M21=(A{i}*Q+B{i}*Y); M22=Q;
M31=(sqrt(Qc)*Q); M32=zeros(m); M33=gamma*eye(m);
M41=(sqrt(Rc)*Y); M42=zeros(n,m);M43=zeros(n,m); M44=gamma*eye(n);
M{i}=[M11 M21' M31' M41';M21 M22 M32' M42';M31 M32 M33 M43'; M41 M42 M43 M44];
LMIs=[LMIs,M{i}>=0];
end
LMIs=[LMIs, [1 x_k';x_k Q]>0];
LMIs=[LMIs, [umax*umax*eye(n) Y;Y' Q]>0];

% ops=sdpsettings('solver','sedumi','sedumi.eps',1e-5);
% ops=sdpsettings('solver','lmilab');

ops=sdpsettings;

% % Solvers
ops.solver='sedumi';    
% ops.solver='sdpt3';
% ops.solver='lmilab';

ops.tol=1e-5;
ops.verbose=1;  % Essa opção oculta os cálculos realizados pelo solver

model = optimizer(LMIs, gamma,ops,x_k,{Y,Q,gamma});
out=model{xk};
Ys=out{1};
Qs=out{2};
gammas=out{3};
Fs=Ys*inv(Qs);

