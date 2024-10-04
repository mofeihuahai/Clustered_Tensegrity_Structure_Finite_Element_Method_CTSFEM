%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%计算加载后的结构节点的位移%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% [1] structure design(define N, C, boundary constraints, clustering,
% calculate equilibrium matrix,
% group matrix,prestress mode, minimal mass design)
% [2] calculate tangent stiffness matrix, material
% stiffness, geometry stiffness,
% [3] dynamic simulation

%EXAMPLE
clc; clear all; close all;
% Global variable、
pully =1;                      %pully influence 1 or 0
if pully 
  [consti_datap,Ep,sigmap,rho_p]=material_lib_p('Steel_Q345');  
end
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','nylon_rope');

% %写function
% [consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib_multi('material-1','material-2','material-3',''...);
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=1; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% static analysis set
substep=1;                                     %荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 
%% N C of the structure
% Manually specify node positions of double layer prism.
R1=0.28;R2=0.24;R3=0.24;R4=0.195;          %radius
hl=0.05;     %滑轮长
p=12;          %complexity for cable dome
% generate node in one unit
beta1=2*pi/p;beta2=pi/p;
T1=[cos(beta1) -sin(beta1) 0
    sin(beta1) cos(beta1) 0
    0 0 1];
T2=[cos(beta2) -sin(beta2) 0
    sin(beta2) cos(beta2) 0
    0 0 1];
N0=[R1*T2*[0;0;0],R2*T2*[0;0;0],R3*T1*[1;0;0],R4*T1*[1;0;0]];      %initial N
N=[];
Nc1=zeros(3,p);
Nc2=zeros(3,p);
for i=1:p
    alpha=beta2+(i-1)*beta1;
    if i<=3
            r=2^(1/2)*sin(pi/4)/sin(3*pi/4-alpha)/2;
            Nc1(1,i)=r*cos(alpha);
            Nc1(2,i)=r*sin(alpha);
            Nc2(1,i)=(r-hl)*cos(alpha);
            Nc2(2,i)=(r-hl)*sin(alpha);
    end
        if i>3 && i<=6
            r=2^(1/2)*sin(pi/4)/sin(5*pi/4-alpha)/2;
            Nc1(1,i)=-1*r*sin(alpha-pi/2);
            Nc1(2,i)=r*cos(alpha-pi/2);
            Nc2(1,i)=-1*(r-hl)*sin(alpha-pi/2);
            Nc2(2,i)=(r-hl)*cos(alpha-pi/2);
        end
         if i>6 && i<=9
            r=2^(1/2)*sin(pi/4)/sin(7*pi/4-alpha)/2;
            Nc1(1,i)=-1*r*cos(alpha-pi);
            Nc1(2,i)=-1*r*sin(alpha-pi);
            Nc2(1,i)=-1*(r-hl)*cos(alpha-pi);
            Nc2(2,i)=-1*(r-hl)*sin(alpha-pi);
         end
        if i>9 && i<=12
            r=1.9208^(1/2)*sin(pi/4)/sin(9*pi/4-alpha)/2;
            Nc1(1,i)=r*sin(alpha-3*pi/2);
            Nc1(2,i)=-1*r*cos(alpha-3*pi/2);
            Nc2(1,i)=(r-hl)*sin(alpha-3*pi/2);
            Nc2(2,i)=-1*(r-hl)*cos(alpha-3*pi/2);
        end
end
for i=1:p    %rotate nodes
 N=[N,T1^(i-1)*N0];
end
for i=1:p
    N(:,1+4*(i-1))=Nc1(:,i);
    N(:,2+4*(i-1))=Nc2(:,i);
end

C_b_in=[[1:4:4*p-3]',[2:4:4*p-2]';[3:4:4*p-1]',[4:4:4*p]';];
C_s_in=[[2:4:4*p-2]',[3:4:4*p-1]';[3:4:4*p-1]',[6:4:4*p-2,2]';[4:4:4*p]',[8:4:4*p,4]'];

C_b = tenseg_ind2C(C_b_in,N);%%将上面输入c_b_in转化为关系矩阵
C_s = tenseg_ind2C(C_s_in,N);
C=[C_s;C_b];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
% Plot the structure to make sure it looks right
tenseg_plot(N,[],C);%绘图
title('Cable dome');
%% Boundary constraints
pinned_X=([1:4:4*p-3])'; pinned_Y=([1:4:4*p-3])'; pinned_Z=([1:4:4*p-3])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%% Group/Clustered information 
%generate group index
% gr=[];
gr={[1:2*p]',[2*p+1:3*p]'};
Gp=tenseg_str_gp2(gr,C);    %generate group matrix
% S=eye(ne);                  % no clustering matrix
S=Gp';                      % clustering matrix is group matrix
%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

l_c=S*l;                % length vector CTS  %同组构件的总长
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);%A_lag*V2=0   U2'*A_lag=0  V2即为自应力模态矩阵

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;   %w0 72*1  w0a  36*1

%prestress design
index_gp=[2];                   % number of groups with designed force %加载初始力的索引，由此力计算其他部分的力
fd=50;                        % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
t_c=pinv(S')*t;    %力（分组后） 2*1
q_c=pinv(S')*q;    %力密度（分组后） 2*1
%% cross sectional design 
index_b=find(t_c<0);              % index of bar in compression  查找杆件
index_s=setdiff(1:size(S,1),index_b);	% index of strings       查找索件 
[A_b,A_s,A_c,A,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
%redefine young's modulus
%E_c_p=2e11*ones(24,1);%+
%E_c=[E_c_s;E_c_p];%+
%E=S'*E_c;     %Young's modulus CTS
I=eye(length(E_c));
index_p=(3:26);
E_c=I*[Es*ones(2,1);Ep*ones(numel(index_p),1)];  
E=S'*E_c;
A_c_s=[0.3e-6;0.3e-6];%+
A_c_p=1e-6*ones(24,1);
A_c=[A_c_s;A_c_p];
A=S'*A_c;     % Cross sectional area CTS
l0=(t+E.*A).\E.*A.*l;%+
l0_c=S*l0;
rho=[1700;1700;7870*ones(24,1)];%+
mass=S'*rho.*A.*l0;
%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
% ind_w=[3:3:6*p];w=-2*ones(2*p,1);
% f1=-4;   %外环节点节点力
% f2=-6;   %内环节点节点力
% ind_w=[3*[2:4:46] 3*[3:4:47] 3*[4:4:48]];w=[f1*ones(p,1);f2*ones(p,1);f2*ones(p,1)];
ind_w=[];w=[];
ind_dnb=[3*[1:4:4*p-3]']; 
%dnb0=0.2*sin(linspace(0,1,p)'*4*pi+1/2*pi);
dnb0=zeros(length(ind_dnb),1);
% m=2/3;                     %%形状系数，控制高差
% for i=1:p
%     dnb0(i)=(N(1,4*i-3)^2-N(2,4*i-3)^2)/m;
% end
dnb0=[0 0.16 0.32 0.32 0.16 0 0 0.16 0.32 0.32 0.16 0];    %抬升高度
% dnb0=1.5*[0.0784 0.0392 -0.0392 -0.0784 -0.0392 0.0392 0.0784 0.0392 -0.0392 -0.0784 -0.0392 0.0392];
ind_dl0_c=[]; dl0_c=[];
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;0;9.8],C,mass);
%% Step1: equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_ct;
data.substep=substep;    % substep

% nonlinear analysis
% data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
data_out1=static_solver_CTS(data);
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
N_out1=data_out1.N_out;
N_out1=N_out1{1};
tenseg_plot(N_out1,C_b,C_s);
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N_out1,C,Gp,Ia);
l0=(t_t+E.*A).\E.*A.*l;
%% 加外部荷载前的索力调节
l_c=S*l;                % length vector CTS  %同组构件的总长
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_2ag);%A_lag*V2=0   U2'*A_lag=0  V2即为自应力模态矩阵
fd=75;                                     %constant force
w_ta=Ia'*w_t;
t1=pinv(A_2ag)*w_ta;
I=eye(length(l_c));
index_h=2;                                %index of string with constant force
e_d=I(:,index_h);
z=(e_d'*V2)\(fd-e_d'*t1);        %self-stress coefficient
t2=z.*V2;
t=t1+t2;
l0_c=(t+E_c.*A_c).\E_c.*A_c.*l_c; 
%% 索力调节后的二次平衡态计算
data.l0_t=l0_c;
data.N=N;
data_out2=static_solver_CTS(data);
t_t=data_out2.t_out;          %member force in every step
n_t=data_out2.n_out;          %nodal coordinate in every step
N_out1=data_out2.N_out{1};
%% 加上滑轮自重后的平衡态计算
% ind_w=[3:3:6*p];w=-2*ones(2*p,1);
f1=-0.5;   %外环节点节点力
f2=-0.8;   %内环节点节点力
f3=-0.8;
ind_w=[3*[2:4:46] 3*[3:4:47] 3*[4:4:48]];w=[f1*ones(p,1);f2*ones(p,1);f3*ones(p,1)];
ind_dnb=[]; dnb0=[];
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;0;9.8],C,mass);
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_ct;% forced movement of pinned nodes
data.N=N_out1;
data_out2=static_solver_CTS(data);
t_t2=data_out2.t_out;          %member force in every step
n_t2=data_out2.n_out;          %nodal coordinate in every step
N_out2=data_out2.N_out;
tenseg_plot(N_out2{1},C_b,C_s);
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N_out2{1},C,Gp,Ia);
l02=(t_t2+E.*A).\E.*A.*l;
l0_c2=S*l02;
%% 单个节点加载后的平衡态计算
substep=6;                %分级加载数量   
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;0;9.8],C,mass);
ind_wp=3*[12:4:32];wp=[-4.9];    %加载节点索引，加载节点力的大小
w_p=zeros(3*nn,1);
w_p(ind_wp)=wp;
w_p=w_p*linspace(0,1,substep);
w_t=repmat(w_t,1,substep);
w_p=w_t+w_p;
data.w_t=w_p;              % 节点自身重力加上加载力
data.N=N_out2{1};
data.dnb_t=dnb_t;
data.substep=substep;
data_out3=static_solver_CTS(data);
t_t3=data_out3.t_out;          
n_t3=data_out3.n_out;          
N_out=data_out3.N_out;
N_out3=N_out{end};
tenseg_plot(N_out3,C_b,C_s);
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N_out3,C,Gp,Ia);
l03=(t_t3(:,end)+E.*A).\E.*A.*l;
l0_c3=S*l03;
%% 输出分级加载时的节点位移图
ind_j=[24 48];                      %所需计算的节点索引
d_j=plot_nodal_displacement(substep,{'dz'},ind_j,N_out2{1},N_out);
%d_j表示节点的位移，各列的四个元素分别表示dx、dy、dz以及节点的位移
%% 保存数据
save ('fd=75半.mat');