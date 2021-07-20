%  Marcus Vinicius - Versao final 3
close all, clear all, clc
Pot=[1000 380]; Vo=48; Vg=[36 26]; L=36e-6; Co=4400e-6; Rco=26.7e-3;
Ts=1e-3; Ta=Ts;
g=1;h=1;
Q1=diag([0.1 0.1 0.1],0); R1=0.1;
umax=.5;
% Formulação da resposta impulso para potência nominal
t_=0:Ts:20e-3; %gera uma curva de 20ms
[Ahn,Bhn,An,Bn,Cn,Dn]=boost_discreto(Pot(1),Vg(1),Vo,L,Co,Rco,Ts,g,h);
u=0;
xset=[(Pot(1)/Vg(1)) Vo]';
for k=1:length(t_)
xset(:,k+1)=An*xset(:,k)+Bn*u;
end
% Armazena-se 20 pontos
x1set=xset(2,1:end-1);
for i=1:length(Pot)
    for j=1:length(Vg)
        [Ah{2*(i-1)+j},Bh{2*(i-1)+j},A{2*(i-1)+j},B{2*(i-1)+j},C{2*(i-1)+j},D{2*(i-1)+j}]=boost_discreto(Pot(i),Vg(j),Vo,L,Co,Rco,Ts,g,h);
    end
end
for k=1:length(t_)
xhset(:,k)=[(Pot(1)/Vg(1)) x1set(k) 0]';
[gamma(k),F(:,:,k)]=MCP_Robust_kot_yalmip(Ah,Bh,Q1,R1,umax,xhset(:,k));
end
% Ganho do controlador AW relaxado
[Fs,Qs,Ls,gamma]=lmi_antiwindup_yalmip(A,B,C,D);
fprintf('Ganho AW Relaxado: %.10f %.10f \n', Fs(1), F(2));
% Find Normalized Coprime factors of system [a,b,c,d].
a_=An+Bn*Fs;
b_=Bn;
c_=Cn+Dn*Fs;
d_=Dn;
Mz_Ia=a_;Mz_Ib=b_; Mz_Ic=Fs; Mz_Id=0;
Nza=a_;Nzb=b_; Nzc=c_; Nzd=d_;
% Condições iniciais
t=0:Ts:300e-3;
R=[Vo*ones(1,ceil(length(t)*0.5)) Vo*ones(1,ceil(length(t)*0.5))];
Ppot=[Pot(2) Pot(2)*ones(1,ceil((length(t)-1)*0.25)) Pot(2)*ones(1,ceil((length(t)-1)*0.25)) Pot(1)*ones(1,ceil((length(t)-1)*0.25)) Pot(1)*ones(1,ceil((length(t)-1)*0.25))];
Vvg=[Vg(1) Vg(1)*ones(1,ceil((length(t)-1)*0.25)) Vg(2)*ones(1,ceil((length(t)-1)*0.25)) Vg(2)*ones(1,ceil((length(t)-1)*0.25)) Vg(1)*ones(1,ceil((length(t)-1)*0.25))];
%% Modelo Por MPC-AW-RELAX
N=20;
Kmpc=-F(:,:,N);
fprintf('MPC: %.10f \n', Kmpc);
x=[Pot(1)/Vg(2); Vg(2)]; y=0; u=0; v=0; Vo=48; xh=[0 0]'; x1=xh;x2=xh;
eh=0; ud=0;yd=0; u_d=0;y_d=0;
for i=1:length(t)    
    % Sinal de controle LMI
    u(i)=(-Kmpc(1:end-1)*x(:,i)-Kmpc(end)*v)-ud;
    u_ant=u(i);
    if u(i)<0, u(i)=0;end
    if u(i)>umax, u(i)=umax;end
    u_dep=u(i);
    % Aplicação do Antiwind up
    u_til(i)=u_ant-u_dep; 
    % M(z)-I
    x1(:,i+1)=(An+Bn*Fs)*x1(:,i)+Bn*u_til(i);
    ud=Fs*x1(:,i);
    % Nz(z) = G2M(z)
    yd=(Cn+Dn*Fs)*x1(:,i)+Dn*u_til(i);    
    % Modelo da planta
    [A_, B_, C_]=boost_n_linear(Ppot(i),u(i),Vo,L,Co,Rco);    
    % Saturação da corrente no indutor
    if x(1,i)<0, x(1,i)=0;end
    if x(1,i)>105, x(1,i)=105;end     
    % Saturação da tensão no Capacitor
    if x(2,i)<0, x(2,i)=0;end
    if x(2,i)>63, x(2,i)=63;end    
    % Modelo RK4
    k1=A_*x(:,i)+B_*Vvg(i);
    k2=A_*(x(:,i)+0.5*Ts*k1)+B_*Vvg(i);
    k3=A_*(x(:,i)+0.5*Ts*k2)+B_*Vvg(i);
    k4=A_*(x(:,i)+Ts*k3)+B_*Vvg(i);
    x(:,i+1)=x(:,i)+(Ts/6)*(k1+2*(k2+k3)+k4);
    % Equação da saida
    y(i)=C_*x(:,i);         
    if y(i)<0, y(i)=0;end
    y_in=y(i)+yd;    
    % Integrador
    v=g*v+h*(R(i)-y_in);    
    Vo=y(i);    
    y_d(i)=yd;
    u_d(i)=ud;
end
ympc_aw=y; umpc_aw=u; xmpc_aw=x;
%% Modelo MPC
x=[Pot(1)/Vg(2); Vg(2)]; y=0; u=0; v=0; Vo=48; xh=[0 0]'; x1=xh;x2=xh;
eh=0; ud=0;yd=0; u_d=0;y_d=0;
for i=1:length(t)
    % Sinal de controle LMI
    u(i)=(-Kmpc(1:end-1)*x(:,i)-Kmpc(end)*v);
    u_ant=u(i);
    if u(i)<0, u(i)=0;end
    if u(i)>umax, u(i)=umax;end
    % Modelo da planta
    [A_, B_, C_]=boost_n_linear(Ppot(i),u(i),Vo,L,Co,Rco);
    % Saturação da corrente no indutor
    if x(1,i)<0, x(1,i)=0;end
    if x(1,i)>105, x(1,i)=105;end
    % Saturação da tensão no Capacitor
    if x(2,i)<0, x(2,i)=0;end
    if x(2,i)>63, x(2,i)=63;end
    % Modelo RK4
    k1=A_*x(:,i)+B_*Vvg(i);
    k2=A_*(x(:,i)+0.5*Ts*k1)+B_*Vvg(i);
    k3=A_*(x(:,i)+0.5*Ts*k2)+B_*Vvg(i);
    k4=A_*(x(:,i)+Ts*k3)+B_*Vvg(i);
    x(:,i+1)=x(:,i)+(Ts/6)*(k1+2*(k2+k3)+k4);
    % Equação da saida
    y(i)=C_*x(:,i);      
    if y(i)<0, y(i)=0;end
    y_in=y(i)+yd;
    % Integrador
    v=g*v+h*(R(i)-y_in); Vo=y(i);
end
ympc=y; umpc=u; xmpc=x;
%%
[Fs,Qs1,Ls,gamma]=lmi_antiwindup_yalmip_n_relax(A,B,C,D);
fprintf('Não relaxado: %.10f \n', Fs);
% Find Normalized Coprime factors of system [a,b,c,d].
a_=An+Bn*Fs;
b_=Bn;
c_=Cn+Dn*Fs;
d_=Dn;
Mz_Ia1=a_;Mz_Ib1=b_; Mz_Ic1=Fs; Mz_Id1=0;
Nza1=a_;Nzb1=b_; Nzc1=c_; Nzd1=d_;
% Condições iniciais
t=0:Ts:300e-3;
R=[Vo*ones(1,ceil(length(t)*0.5)) Vo*ones(1,ceil(length(t)*0.5))];
Ppot=[Pot(2) Pot(2)*ones(1,ceil((length(t)-1)*0.25)) Pot(2)*ones(1,ceil((length(t)-1)*0.25)) Pot(1)*ones(1,ceil((length(t)-1)*0.25)) Pot(1)*ones(1,ceil((length(t)-1)*0.25))];
Vvg=[Vg(1) Vg(1)*ones(1,ceil((length(t)-1)*0.25)) Vg(2)*ones(1,ceil((length(t)-1)*0.25)) Vg(2)*ones(1,ceil((length(t)-1)*0.25)) Vg(1)*ones(1,ceil((length(t)-1)*0.25))];
%% Modelo Por MPC-AW
N=20;
Kmpc=-F(:,:,N);
x1=[Pot(1)/Vg(2); Vg(2)]; y1=0; u1=0; v1=0; Vo=48; xh1=[0 0]'; x11=xh1; x2=xh1;
eh=0; ud1=0;yd1=0; u_d1=0;y_d1=0;
for i=1:length(t)  
    % Sinal de controle LMI
    u1(i)=(-Kmpc(1:end-1)*x1(:,i)-Kmpc(end)*v1)-ud1;
    u_ant=u1(i);
    if u1(i)<0, u1(i)=0;end
    if u1(i)>umax, u1(i)=umax;end
    u_dep=u1(i);    
    % Aplicação do Antiwind up
    u_til1(i)=u_ant-u_dep; 
    % M(z)-I
    x11(:,i+1)=(An+Bn*Fs)*x11(:,i)+Bn*u_til1(i);
    ud1=Fs*x11(:,i);
    % Nz(z) = G2M(z)
    yd1=(Cn+Dn*Fs)*x11(:,i)+Dn*u_til1(i);   
    % Modelo da planta
    [A_, B_, C_]=boost_n_linear(Ppot(i),u1(i),Vo,L,Co,Rco);
    % Saturação da corrente no indutor
    if x1(1,i)<0, x1(1,i)=0;end
    if x1(1,i)>105, x1(1,i)=105;end 
    % Saturação da tensão no Capacitor
    if x1(2,i)<0, x1(2,i)=0;end
    if x1(2,i)>63, x1(2,i)=63;end  
    % Modelo RK4
    k1=A_*x1(:,i)+B_*Vvg(i);
    k2=A_*(x1(:,i)+0.5*Ts*k1)+B_*Vvg(i);
    k3=A_*(x1(:,i)+0.5*Ts*k2)+B_*Vvg(i);
    k4=A_*(x1(:,i)+Ts*k3)+B_*Vvg(i);
    x1(:,i+1)=x1(:,i)+(Ts/6)*(k1+2*(k2+k3)+k4);
    % Equação da saida
    y1(i)=C_*x1(:,i);
    if y1(i)<0, y1(i)=0;end
    y_in=y1(i)+yd1;
    % Integrador
    v1=g*v1+h*(R(i)-y_in);
    Vo=y1(i);
    y_d1(i)=yd1;
    u_d1(i)=ud1;
end
ympcaw=y1; umpcaw=u1; xmpcaw=x1;
%%
figure(3)
plot(t,ympc_aw,'ro-',t,ympc,'g*-',t,ympcaw,'c.-',t,R(1:end-1),'k-.','linewidth',1.5), legend('MPC_{AW-RELAX}','MPC','MPC_{AW}','Ref'), grid
axis([t(1) t(end) 25 65])
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('Output Voltage (V)','fontsize',25,'fontname','Times New Roman','fontangle','normal')

figure(4)
plot(t,umpc_aw,'ro-',t,umpc,'g*-',t,umpcaw,'c.-','linewidth',1.5), legend('MPC_{AW-RELAX}','MPC','MPC_{AW}'), grid
axis([t(1) t(end) .15 0.55])
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('Duty Cycle','fontsize',25,'fontname','Times New Roman','fontangle','normal')

figure(5)
plot(t,Ppot./ympc_aw,'ro-',t,Ppot./ympc,'g*-',t,Ppot./ympcaw,'c.-','linewidth',1.5), legend('MPC_{AW-RELAX}','MPC','MPC_{AW}'), grid
axis([t(1) t(end) 5 25])
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('Current Io (A)','fontsize',25,'fontname','Times New Roman','fontangle','normal')

figure(6)
plot(t,xmpc_aw(1,1:end-1),'ro-',t,xmpc(1,1:end-1),'g*-',t,xmpcaw(1,1:end-1),'c.-','linewidth',1.5), legend('MPC_{AW-RELAX}','MPC','MPC_{AW}'), grid
axis([t(1) t(end) 0 105])
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('Current I_L (A)','fontsize',25,'fontname','Times New Roman','fontangle','normal')
%%


