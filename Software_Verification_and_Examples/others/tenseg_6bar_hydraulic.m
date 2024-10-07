clc;clear;close all;
%% material
[consti_data,~,~,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties:'linear_elastic' 'multielastic' 'plastic'
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)
Eb=1;
Es=1;
%% cross section design cofficient
thick=6e-3;              % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid11
%% set parameter
c_b=0.5;           % coefficient of safty of bars 0.5
c_s=0.3;           % coefficient of safty of strings 0.3

substep=1;             % 荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=0;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
limit=1e-6;
%% N C of the structure
a=1;
b=a*sqrt((5+sqrt(5))/10);
c=a*sqrt((5-sqrt(5))/10);
r=0.5*b+c;
idx=(0:9).';
z=exp(idx/10*2*pi*1i);
p=[b*real(z),b*imag(z),0.5*b*(-1).^idx;0,0,r;0,0,-r];
N=p';
% Manually specify connectivity indices.bar（指定连接关系）
C_b_in = [1 12;2 5;3 9;4 8;6 11;7 10];   % This is indicating the bar connection
% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);

% Manually specify connectivity indices.string
C_s_in = [1 2;2 3;4 5;5 6;6 7;7 8;9 10;1 10;2 4;4 6;6 8;8 10;1 3;3 5;7 9;1 9;2 12;4 12;8 12;10 12;3 11;5 11;7 11;9 11];  % This is indicating the string connection
% Convert the above matrices into full connectivity matrices.
C_s = tenseg_ind2C(C_s_in,N);

C=[C_b;C_s];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% Boundary constraints
pinned_X=[]; pinned_Y=[]; pinned_Z=[];
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% Group information
%generate group index
% gr={1:6;7:30}; % group of 1 vertical bar(out) ...2 vertical bar(in) 3 top string 4 diagonal string 5 circlur string 6 top string 7 diagonal string 8 circlur string bottom 9 circluar string top
gr={};
Gp=tenseg_str_gp(gr,C);    %generate group matrix

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

%SVD of equilibrium matrix
[U1,U2,V1,V2,S]=tenseg_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
% index_gp=1;                   % number of groups with designed force
% fd=-10000;                        % force in bar is given as -1000
% [q_gp,t_gp,q,t,z]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
% t(1:6)=-100;   % to specify the compression members
t=zeros(ne,1);
%% cross sectional design
index_b=find(t<0);          % index of bar in compression
index_b=1:6;
index_s=setdiff(1:ne,index_b);	% index of strings
A=0.01*ones(ne,1);
E=1e9*ones(ne,1);
l0=l;
l0(index_b)=l(index_b)+0.1;
% l0=E.*A.*l./(t+E.*A);
rho=ones(ne,1);
mass=rho.*A.*l0;
%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[]; % exert force on top node of cable dome, exertnal node
ind_dnb=[]; dnb0=[];
ind_dl0=[]; dl0=[];   %extent rest length of bar
[w_t,dnb_t,l0_t,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0,dl0,l0,b,gravity,[0;0;9.8],C,mass);

% case 'force' 'displacement' 'rest_length'

%% equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;
data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_t;% forced movement of pinned nodes
data.substep=substep;    % substep
data.limit=limit;
% for func stress_strain cal
data.Eb=Eb;
data.Es=Es;
data.f_int=t;
data.l_int=l;

% nonlinear analysis
% tic
data_out=static_solver(data);        %solve equilibrium using mNewton method

% output final configuration
t_t=data_out.t_out;          %member force in every step
n_t=data_out.n_out;          %nodal coordinate in every step
N_out=data_out.N_out{1,1};
tenseg_plot(N_out,C_b,C_s);
N=N_out;


%%
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
%%
q_index=data_out.q_out(1,1)/data_out.q_out(7,1);
l_index=data_out.l_out(1,1)/data_out.l_out(7,1);

%% Hydraulic analysis
%% Hagen Poiseuille's law for laminar flow
rho_f=900;miu_f=0.03;g=9.8;      % 流体材料参数
r_s=0.004;r_n=0.0045;ts=0.001;r_b=0.001;
A_b=pi*r_b^2;A_s=pi*r_s^2;A_n=pi*r_n^2;
n_in=1;n_out=6;
v_in=2;p_out=0;
v=[zeros(1,(n_in-1)),v_in,zeros(1,(n_out-n_in-1)),-v_in,zeros(1,nn-n_out)]';

Re=rho_f*v_in*2*r_s/miu_f;

H_s=N*C_s';                     % strings's direction matrix
l_s=sqrt(diag(H_s'*H_s));       % strings' length
L=C_s'*diag(1./l_s)*C_s;        % 系数矩阵
k=pi*r_s^4/8/miu_f;
Q=A_s*v;

% 得把系数矩阵拆分 
% P=[zeros(1,n_out-1),p_out,zeros(1,nn-n_out)]';                           % 如果将节点压强初始化，除出口定压节点外，其余置零
% k*L*P=Q;
% k*[L(:,1:n_out-1),L(:,n_out),L(:,nn-n_out)]*[P1,p_out,P2]'=Q;
% [L(:,1:n_out-1),L(:,nn-n_out)]*[P1,P2]'=Q/k-L(:,n_out)*p_out
L0=[L(:,1:n_out-1),L(:,(n_out+1):nn)];
P0=[L(:,1:n_out-1),L(:,(n_out+1):nn)]\(Q/k-L(:,n_out)*p_out);
P1=P0(1:n_out-1,1);P2=P0(n_out:nn-1,1);
P=[P1',p_out,P2']';          % 每个节点的压强
delta_P=-C_s*P;              % 每根管段中的压强差
q=k*(delta_P./l_s);          % 每根管段中的流量
V_s=q/A_s;                   % 每根管段中的流速
Re_s=rho_f*abs(V_s*2)*r_s/miu_f;  % 每根管段中的雷诺数
%% simplified FSI calculation   transform the action of the fluid into external forces
% V_s 每根管段中的流速  P 每个节点的压强
%与V_s和q有关的系数矩阵的计算
ns=size(C_s,1);
kv=C_s'*diag(V_s);
kq=pi*r_s^2*abs(kv);
kqv=rho_f*kv.*kq;
kp0=zeros(nn,ns);
%与p有关的系数矩阵的计算
for i=1:nn
    As=C_s';
    col=find(As(i,:));                % 与第i个节点相连的管段编号
    for j=1:size(col,2) 
        col1=find(C_s(col(:,j),:));     % 第j个管段连接的节点编号
        col2=find(col1~=i);                 % 管段j上，找和第i节点相连的节点编号
        kp0(i,col(:,j))=P(col1(:,col2),1);
    end  
end
kp=A_s*kp0.*C_s';
kF=kqv+kp;

H_s=N*C_s';                     % strings's direction matrix
l_s=sqrt(diag(H_s'*H_s));       % strings' length
n=H_s./(ones(3,1)*l_s');        % 索的单位方向向量

F_n=kF*n';     % n为方向向量，包含进出口导管的方向向量
w_f=[];       % fluid pressure acts as external force

%%
%% Boundary constraints
pinned_X=[]; pinned_Y=[]; pinned_Z=[];
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% Group information
%generate group index
% gr={1:6;7:30}; % group of 1 vertical bar(out) ...2 vertical bar(in) 3 top string 4 diagonal string 5 circlur string 6 top string 7 diagonal string 8 circlur string bottom 9 circluar string top
gr={};
Gp=tenseg_str_gp(gr,C);    %generate group matrix

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N_out,C,Gp,Ia);

%SVD of equilibrium matrix
[U1,U2,V1,V2,S]=tenseg_svd(A_1ag);
diag(S)
A_1a*t_t