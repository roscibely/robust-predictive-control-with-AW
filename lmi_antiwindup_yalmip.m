function [Fs,Qs,Ls,gamma]=lmi_antiwindup_yalmip(A,B,C,D)
[m,n]=size(B{1});
mu=sdpvar(1);
Q=sdpvar(m,m,'symmetric');
U=sdpvar(n,n,'symmetric');
L=sdpvar(n,m);
Xa=sdpvar(n,n, 'full');
politopo = length(A);
% Usando o YALMIP
LMIs=[Q>=0,U>=0, mu>=0];


% Processo de Restrição
for i=1:politopo
M11=-(Xa+Xa'-Q);
M12=-L'; M21=M12';
M13=zeros(m,n); M31=M13';
M14=(C{i}*Xa+D{i}*L)'; M41=M14';
M15=(A{i}*Xa+B{i}*L)'; M51=M15';

M22=-2*U; 
M23=eye(n); M32=M23';
M24=U*D{i}'; M42=M24';
M25=U*B{i}'; M52=M25';

M33=-mu*eye(n);
M34=zeros(n,n); M43=M34';
M35=zeros(n,m); M53=M35';

M44=-eye(n);
M45=zeros(n,m); M54=M45';

M55=-(Xa+Xa'-Q);

Maw=[M11 M12 M13 M14 M15;
     M21 M22 M23 M24 M25;
     M31 M32 M33 M34 M35;
     M41 M42 M43 M44 M45;
     M51 M52 M53 M54 M55];
end
LMIs=[LMIs, Maw<0];
ops=sdpsettings;

% % Solvers
ops.solver='sedumi';    
% ops.solver='sdpt3';
% ops.solver='lmilab';

ops.tol=1e-5;
ops.verbose=1;  % Essa opção oculta os cálculos realizados pelo solver

optimize(LMIs, mu,ops);
gamma=sqrt(value(mu));
Qs = value(Q);
Xs=value(Xa);
Ls = value(L);
Fs = Ls*inv(Xs);

checkset(LMIs)