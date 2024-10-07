%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A CTS T bar example: nonlinear shape control %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% [1] structure design(calculate equilibrium matrix,
% group matrix,prestress mode, minimal mass design)
% [2] modal analysis(calculate tangent stiffness matrix, material
% stiffness, geometry stiffness, generalized eigenvalue analysis)
% [3] dynamic simulation

%EXAMPLE
clc; clear all; close all;
% Global variable
global M 
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% static analysis set
% substep=10;                                     %ºÉÔØ×Ó²½
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 

%dynamic analysis set
dt=1e-4;               % time step in dynamic simulation
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=1;                   % final time of dynamic simulation
out_dt=0.002;            % output data interval(approximately, not exatly)

amplitude=0;            % amplitude of external force of ground motion 
period=0.5;             %period of seismic

%% N C of the structure
% Manually specify node positions of double layer prism.
h=1;b=1;
N_1=b*[-1 1 0;1 1 0;1 -1 0;-1 -1 0]';
N_2=b*[-1 0 h;0 1 h;1 0 h;0 -1 h]'; 
N_3=b*[-1 -1 2*h;-1 1 2*h;1 1 2*h;1 -1 2*h]';
N=[N_1,N_2,N_3];

% Manually specify connectivity indices.
C_b_in = [1 8;4 7;3 6;2 5;7 10;6 9;5 12;8 11];  % Similarly, this is saying bar 1 connects node 1 to node 2,
C_s_in = [4 8;3 7;2 6;1 5;8 12;7 11;6 10;5 9;1 4;4 3;3 2;2 1;5 8;8 7;7 6;6 5;9 12;12 11;11 10;10 9];  % This is indicating that string connection

% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);%%
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);


%% Boundary constraints
pinned_X=[1:4]'; pinned_Y=[1:4]'; pinned_Z=(1:4)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% Group/Clustered information 
%generate group index
% gr={[9:16]};     % number of elements in one group
gr={[9,13],[10,14],[11,15],[12,16]};   
% gr={[3,4];[5,6]};     % number of elements in one group
% gr=[];                     %if no group is used
Gp=tenseg_str_gp(gr,C);    %generate group matrix
S=Gp';                      % clustering matrix
%% self-stress design
%Calculate equilibrium matrix and member length
[A_1,A_1c,A_1a,A_1ac,A_2,A_2c,A_2a,A_2ac,l,l_c]=tenseg_equilibrium_matrix_CTS(N,C,S,Ia);
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_2ac);
%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
% index_gp=[21,17:20]';                 % number of groups with designed force
index_gp=[9,13:16]';                 % number of groups with designed force
fd=[1e2;1e2*ones(4,1)];              % force in bar is given as -1000
[q_c,t_c,q,t]=tenseg_prestress_design_f(S,l,l_c,A_2ac,V2,w0a,index_gp,fd);

%% cross sectional design
index_b=find(t_c<0);              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings
[A_b,A_s,A_c,A,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
% A_c(9)=max(A);      %increase area of active member
A_c(17:20)=1e-8*ones(4,1);          % reduce the stiffness of middle string
E=S'*E_c;     %Young's modulus TTS
A=S'*A_c;     % Cross sectional area TTS
l0=(t+E.*A).\(E.*A.*l);
l0_c=S*l0;
mass=S'*rho.*A.*l0;
% % Plot the structure with radius
% R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.03,.1],r_b);
% R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.03,.1],r_s);
% R3Ddata.Nradius=0.1*ones(nn,1);
% tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);
tenseg_plot_CTS(N,C,index_b,S,[],[],[],[],[],t,[],[min(t),max(t)]);
tenseg_plot_CTS(N,C,index_b,S)
%% tangent stiffness matrix
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
% plot the mode shape of tangent stiffness matrix
num_plt=1:4;
plot_mode(K_mode,k,N,Ia,C_b,C_s,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.2,saveimg);

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0.1;     %damping coefficient
d_cri=damping_critical(rho,E_c,A_c);
d_c=d*d_cri;          %damping coefficient of all members
D=A_2c*diag(d_c)*A_2c';     %damping matrix

%% mode analysis
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
plot_mode(V_mode,omega,N,Ia,C_b,C_s,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg);

%% external force, forced motion of nodes, shrink of strings
% time step
if auto_dt
dt=pi/(8*max(omega)); 	% time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and 
ind_w=[];w=[];
% ind_dl0_c=[3]; dl0_c=[-1.0];
ind_dl0_c=[]; dl0_c=[];
[w_t,l0_ct]=tenseg_load_prestress_CTS(tspan,ind_w,w,'ramp',ind_dl0_c,dl0_c,l0_c,gravity,[0;0;0],C,mass);

% boundary node motion info
[~,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'vib_force',gravity,[0;0;9.8],C,mass,[1,2],amplitude,period);

% give initial speed of free coordinates
n0a_d=zeros(numel(a),1);                    %initial speed in X direction
    
%% Specify control objectives
ind_n_ct=[[9:12]'*3];n_ct1=[1*h*ones(4,1)];n_ct2=[1*h*ones(4,1)];
% ind_n_ct=[7:8];n_ct1=[1.8;0.3];n_ct2=[1.8;0.3];
Ic=transfer_matrix(ind_n_ct,a);         %transfer matrix for control coordinate
[n_ct_t,n_ct_dt,n_ct_ddt]=coord_vel_acc(tspan,n_ct1,n_ct2);     %nodal coordinate of control target
% num=numel(tspan);
% M_diff=triu(ones(num))-2*triu(ones(num),1)+triu(ones(num),2);   % matrix to calculate diff of vectors in time history
% M_diff(:,1)=0;                                      %(1st term 0)
% n_ct_t=n_ct1*ones(1,num)+(n_ct2-n_ct1)*[linspace(0,1,round(num/2)),linspace(1,1,num-round(num/2))];  % nodal coordinate of control target in time history
% n_ct_dt=n_ct_t*M_diff;
% n_ct_ddt=n_ct_dt*M_diff;

a_c=2*sqrt(100);    % coefficient in error dynamics(damping term)
b_c=100;            % coeffieient in error dynamics(stiffness term)

%% choose active, passive members
% ind_act=[9]; ind_pas=setdiff(1:size(S,1),ind_act);
ind_act=[9:12']; ind_pas=setdiff(1:size(S,1),ind_act);
I_act=transfer_matrix(ind_act,1:size(S,1));
I_pas=transfer_matrix(ind_pas,1:size(S,1));
% lower bound and upper bound of actuator
lb=zeros(numel(ind_act),1);    %lower bound of active string positive
% lb=-sigmas*I_act'*A_c;    %lower bound of active string positive
ub=sigmas*I_act'*A_c;             %upper bound is yielding force
%% dynamics calculation

% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia; data.Ib=Ib;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;           % external force
data.dnb_t=dnb_t; data.dnb_d_t=dnb_d_t;  data.dnb_dd_t=dnb_dd_t; % forced movement of pinned nodes
data.l0_t=l0_c;         % forced movement of pinned nodes
data.n0a_d=n0a_d;        %initial speed of free coordinates
data.M=M;data.D=D;
data.rho=rho_s;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;
data.Ic=Ic; data.n_ct_t=n_ct_t;data.n_ct_dt=n_ct_dt;data.n_ct_ddt=n_ct_ddt;
data.a_c=a_c;data.b_c=b_c;
data.I_act=I_act;data.I_pas=I_pas;
data.lb=lb;data.ub=ub;
%% shape control analysis
% solve dynamic equation with closed loop control
data_out=shape_control_CTS(data);
% time history of structure in shape control
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate 
l_t=data_out.l_t;   %time history of members' length 
l0c_t=data_out.l0_t; %time history of members' rest length 
nd_t=data_out.nd_t;   %time history of nodal coordinate
exitflag_t=data_out.exitflag_t;   % exitflag of lsqlin

%% plot member force 
tenseg_plot_result(out_tspan,t_t([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% Plot rest length l0_t
tenseg_plot_result(out_tspan,l0c_t([1,9,13,17,21],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Length (m)'},'plot_coordinate.png',saveimg);
%% plot exitflag
tenseg_plot_result(out_tspan,exitflag_t,{'Exitflag'},{'Time (s)','Length (m)'},'plot_coordinate.png',saveimg);

 
%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,90])
axis off 
%% save output data
if savedata==1
    save (['prism_CTS_1','.mat']);
end
%% make video of the dynamic
name=['CTS_prism_final'];
% tenseg_video(n_t,C_b,C_s,[],50,name,savevideo,material{2})
tenseg_video_CTS(n_t,C,index_b,S,[],[],[],[],[],[],t_t,[],min(numel(out_tspan),50),tf,name,savevideo)

%output data to tecplot
tenseg_tecplot(C,n_t,t_t,interp1([min(radius),max(radius)],[0.2,0.8],radius));
