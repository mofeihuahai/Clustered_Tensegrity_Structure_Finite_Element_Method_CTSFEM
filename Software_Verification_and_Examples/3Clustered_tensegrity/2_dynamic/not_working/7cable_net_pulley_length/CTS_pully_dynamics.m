%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%A Clustered Cable Net(deployable)%%%%%%%%
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
% Global variable
pully =1;                      %pully influence 1 or 0
if pully 
  [consti_datap,Ep,sigmap,rho_p]=material_lib_p('Steel_Q345');  
end
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','nylon_rope');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

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
gravity=1;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 

%dynamic analysis set
dt=0.00001;               % time step in dynamic simulation
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=3;                   % final time of dynamic simulation
out_dt=0.2;            % output data interval(approximately, not exatly)

amplitude=0;            % amplitude of external force of ground motion 
period=0.5;             %period of seismic
%% N C of the structure
% Manually specify node positions of double layer prism.
R1=0.28;R2=0.24;R3=0.134;R4=0.084;          %radius
p=12;          %complexity for cable dome
% generate node in one unit
beta1=2*pi/p;beta2=pi/p';
T1=[cos(beta1) -sin(beta1) 0
    sin(beta1) cos(beta1) 0
    0 0 1];
T2=[cos(beta2) -sin(beta2) 0
    sin(beta2) cos(beta2) 0
    0 0 1];
N0=[R1*[1;0;0],R2*[1;0;0],R3*T2*[1;0;0],R4*T2*[1;0;0]];      %initial N
N=[];
for i=1:p    %rotate nodes
 N=[N,T1^(i-1)*N0];
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
A_1ac=A_1a*S';          %equilibrium matrix CTS
A_2ac=A_2a*S';          %equilibrium matrix CTS
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
I=eye(length(E_c));
index_p=(3:26);
E_c=I*[Es*ones(2,1);Ep*ones(numel(index_p),1)];  
E=S'*E_c;
A_c_s=[0.5e-6;0.5e-6];%+
A_c_p=1e-6*ones(24,1);
A_c=[A_c_s;A_c_p];
A=S'*A_c;     % Cross sectional area CTS
l0=(t+E.*A).\E.*A.*l;%+
l0_c=S*l0;
rho=[1700;1700;7870*ones(24,1)];%+
mass=S'*rho.*A.*l0;
%% tangent stiffness matrix
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
% plot the mode shape of tangent stiffness matrix
num_plt=1:4;
plot_mode(K_mode,k,N,Ia,C_b,C_s,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.8,saveimg);
%% input file of ANSYS
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping
%% mode analysis
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode  Kt_aa*V_mode=Ia'*M*Ia  *  V_mode  *  D1
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
plot_mode(V_mode,omega,N,Ia,C_b,C_s,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.8,saveimg);
%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
%ind_w=[3:3:6*p];w=-2*ones(2*p,1);
ind_w=[];w=[];
ind_dnb=[3*[1:4:4*p-3]']; 
%dnb0=0.2*sin(linspace(0,1,p)'*4*pi+1/2*pi);
dnb0=1.5*[0.0784 0.0392 -0.0392 -0.0784 -0.0392 0.0392 0.0784 0.0392 -0.0392 -0.0784 -0.0392 0.0392];
ind_dl0_c=[]; dl0_c=[];
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;0;9.8],C,mass);
%% Step1: statics: equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_ct;% forced movement of pinned nodes
data.substep=substep;    % substep

% nonlinear analysis
% data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
data_out1=static_solver_CTS_z(data);
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
N_out=data_out1.N_out;
N=N_out{end};   %+
H=N*C';          %++
l=sqrt(diag(H'*H));  %++
l0=(t_t+E.*A).\E.*A.*l;%++
l0_c=S*l0;           %++
%% Step 2: dynamics:change rest length of strings
% time step
if auto_dt
dt=pi/(1/8*max(omega)); 	% time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and 
ind_w=[];w=[];
ind_dl0_c=[1]; dl0_c=[-2];
% ind_dl0_c=[2]; dl0_c=[-120];
[w_t,l0_ct]=tenseg_load_prestress_CTS(tspan,ind_w,w,'ramp',ind_dl0_c,dl0_c,l0_c,gravity,[0;0;0],C,mass);

% boundary node motion info
[~,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'vib_force',gravity,[0;0;9.8],C,mass,[1,2],amplitude,period);

% give initial speed of free coordinates
n0a_d=zeros(numel(a),1);                    %initial speed in X direction
    %% dynamics calculation

% input data
data.N=N_out{end}; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia; data.Ib=Ib;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;           % external force
data.dnb_t=dnb_t; data.dnb_d_t=dnb_d_t;  data.dnb_dd_t=dnb_dd_t; % forced movement of pinned nodes
data.l0_t=l0_ct;         % forced movement of pinned nodes
data.n0a_d=n0a_d;        %initial speed of free coordinates
data.M=M;data.D=D;
data.rho=S'*rho;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;

%% dynamic analysis
% solve dynamic equation
data_out=dynamic_solver_CTS(data);        %solve ODE of dynamic equation
% time history of structure
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate 
l_t=data_out.l_t;   %time history of members' length 
nd_t=data_out.nd_t;   %time history of nodal coordinate


%% plot member force 
tenseg_plot_result(out_tspan,t_t([1:3,2*p+1:2*p+3],:),{'1','2','3','4','5','6'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t([3*4-2,3*4],:),{'4X','4Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);

%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,90])

%% save output data
if savedata==1
    save (['CTS_pully_dynamics',material{1},'.mat']);
end
%% make video of the dynamic
name=['CTS_pully_dynamics',material{1},'_slack_',num2str(material{2})];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t,C_b,C_s,[],min(numel(out_tspan),50),name,savevideo,material{2})

%output data to tecplot
tenseg_tecplot(C,n_t,t_t,interp1([min(radius),max(radius)],[0.2,0.8],radius));
