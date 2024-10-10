 function Yd=tenseg_dyn_x_xdot_CTS_membrane(t,Y,data_in)
% Input:
%   t: current time
%   Y0: current X, Xd values:   Y0=[X;Xd];
% Output:
%   Yd=[Xd,Xdd]

global l_tc l_l t_tc t_l E_tc n n_d l0_tc l0_l A_p  sigma_l
C_t=data_in.C_t;
C_l=data_in.C_l;
Ia=data_in.Ia;
Ib=data_in.Ib;
S_tc=data_in.S_tc;
E_tc=data_in.E_tc;
A_tc=data_in.A_tc;
A_p=data_in.A_p;
t_p=data_in.t_p;  %板的厚度
B_epsilon=data_in.B_epsilon;
C_pl_bar=data_in.C_pl_bar;
D=data_in.D;  %转换矩阵

l0_tc_t=data_in.l0_tc_t;            
l0_l_t=data_in.l0_l_t; 

n0=data_in.N(:);
w_t=data_in.w_t;           % external force
dnb_t=data_in.dnb_t;       % forced node displacement
dnb_d_t=data_in.dnb_d_t;    %velocity of forced moved node
dnb_dd_t=data_in.dnb_dd_t;    %velocity of forced moved node

M=data_in.M;
Cd=data_in.Cd;   %阻尼矩阵
dt=data_in.dt;

nf=numel(Y)/2;
na=Y(1:nf,:);              %free node cooridnate
na_d=Y(nf+1:end,:);         %free node velocity
%%
% tspan=0:dt:tf;
ind=floor(t/dt)+1;
% dnb=dnb_t(:,ind);
% nb_d=dnb_d_t(:,ind); %this is the velocity of fixed node
% nb_dd=dnb_dd_t(:,ind); %this is the acceleration of fixed node

% Get current pinned nodes displacement 
if ischar(dnb_t)
    run(dnb_t) % if dnb is a string, run that script
elseif size(dnb_t,2)==1
    dnb = dnb_t; % dnb can be constant
else
    dnb = dnb_t(:,ind); % or dnb can be time-varying
end

% Get current pinned nodes velocity 
if ischar(dnb_d_t)
    run(dnb_d_t) % if dnb is a string, run that script
elseif size(dnb_d_t,2)==1
    nb_d = dnb_d_t; % dnb can be constant
else
    nb_d = dnb_d_t(:,ind); % or dnb can be time-varying
end

% Get current pinned nodes acceleration 
if ischar(dnb_dd_t)
    run(dnb_dd_t) % if dnb is a string, run that script
elseif size(dnb_dd_t,2)==1
    nb_dd = dnb_dd_t; % dnb can be constant
else
    nb_dd = dnb_dd_t(:,ind); % or dnb can be time-varying
end

% Get current external forces
if ischar(w_t)
    run(w_t) % if w_t is a string, run that script
elseif size(w_t,2)==1
    w = w_t; % w_t can be constant
else
    w = w_t(:,ind); % or W can be time-varying
end

% Get current rest length
if ischar(l0_tc_t)
    run(l0_tc_t) % if w_t is a string, run that script
elseif size(l0_tc_t,2)==1
    l0_tc = l0_tc_t; % w_t can be constant
else
    l0_tc = l0_tc_t(:,ind); % or W can be time-varying
end

% Get current rest angle
if ischar(l0_l_t)
    run(l0_l_t) % if w_t is a string, run that script
elseif size(l0_l_t,2)==1
    l0_l = l0_l_t; % w_t can be constant
else
    l0_l = l0_l_t(:,ind); % or W can be time-varying
end
%%   calculate stiffness matrix
% w=w_t(:,ind);
nb=Ib'*n0+dnb;
n=Ia*na+Ib*nb;
n_d=Ia*na_d+Ib*nb_d;

l_t=sqrt(sum((reshape(n,3,[])*C_t').^2))';  % elements' length of truss
l_tc=S_tc*l_t;
l_l=sqrt(sum((reshape(n,3,[])*C_l').^2))';  % elements' length of panel lines
    
H_t=reshape(n,3,[])*C_t';
H_l=reshape(n,3,[])*C_l';
Delta_l_l=l_l-l0_l;
Delta_l_tc=l_tc-l0_tc;
Cell_H_t=mat2cell(H_t,3,ones(1,size(H_t,2)));          % transfer matrix H into a cell: Cell_H
Cell_H_l=mat2cell(H_l,3,ones(1,size(H_l,2)));          % transfer matrix H into a cell: Cell_H
A_2t=kron(C_t',eye(3))*blkdiag(Cell_H_t{:})*diag(l_t.^-1);
A_2tc=A_2t*S_tc';
A_2l=kron(C_l',eye(3))*blkdiag(Cell_H_l{:})*diag(l_l.^-1);
t_tc=diag(E_tc)*diag(A_tc)*diag(l0_tc.^-1)*Delta_l_tc;
t_l=(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*Delta_l_l;
sigma_l=D*inv(B_epsilon)*C_pl_bar*Delta_l_l;        
%% calculate accerlation
a=Ia'*(w-M*Ib*nb_dd-Cd*n_d-A_2tc*t_tc-A_2l*t_l);
b=inv(Ia'*M*Ia);
na_dd=(Ia'*M*Ia)\(Ia'*(w-M*Ib*nb_dd-Cd*n_d-A_2tc*t_tc-A_2l*t_l));      %dynamic equation
Yd=[na_d;na_dd];
