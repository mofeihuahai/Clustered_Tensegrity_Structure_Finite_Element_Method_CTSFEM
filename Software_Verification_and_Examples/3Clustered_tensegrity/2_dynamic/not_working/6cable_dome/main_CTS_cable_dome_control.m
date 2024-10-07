%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%An double layer tensegrity tower with simplex%
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
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% dynamic analysis set
amplitude=50;            % amplitude of external force of ground motion 
period=0.5;             %period of seismic

dt=0.001;               % time step in dynamic simulation
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=2;                   % final time of dynamic simulation
out_dt=0.02;            % output data interval(approximately, not exatly)
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 

%% N C of the structure
% Manually specify node positions of double layer prism.
R=50;          %radius
p=12;          %complexity for cable dome
m=2;   %number of circle of the vertical bars
h=0.15*2*R;   %hight of the dome
beta=30*pi/180*ones(m,1);    %all angle of diagonal string
% [N,C_b,C_s,C] =generat_cable_dome(R,p,m,h,beta);
rate=0.3;
[N,C_b,C_s,C] =N_cable_dome(R,rate,p,m,h,beta);
[ne,nn]=size(C);% ne:No.of element;nn:No.of node

tenseg_plot(N,C_b,C_s);
title('Cable dome');
tenseg_plot(N,[],C);
%% Boundary constraints
pinned_X=([5:5:60])'; pinned_Y=([5:5:60])'; pinned_Z=([5:5:60])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% Group/Clustered information 
%generate group index
gr_whg=[1:13:13*(p-1)+1];    % whs
gr_nhg=[2:13:13*(p-1)+2];    % whs
gr_whs=[7:13:13*(p-1)+7];    % whs
gr_nhs=[12:13:13*(p-1)+12];     % nhs
gr_nhds=[13:13:13*(p-1)+13];     % nhds
gr_wxs=[kron(ones(1,p),[5,6])+13*kron(0:p-1,[1,1])];     % wxs
gr_nxs=[kron(ones(1,p),[10,11])+13*kron(0:p-1,[1,1])];     % nxs
gr_wjs=[kron(ones(1,p),[3,4])+13*kron(0:p-1,[1,1])];    % wjs
gr_njs=[kron(ones(1,p),[8,9])+13*kron(0:p-1,[1,1])];    %njs
gr={gr_whg;gr_nhg;gr_whs;gr_nhs;gr_nhds;gr_wxs;gr_nxs;gr_wjs;gr_njs};

% several clustered string in one group of ridge and diagonal strings
num_clu=3;
gr_whs2=reshape(gr_whs,numel(gr_whs)/num_clu,[])';
gr_nhs2=reshape(gr_nhs,numel(gr_nhs)/num_clu,[])';
gr_nhds2=reshape(gr_nhds,numel(gr_nhds)/num_clu,[])';
gr_wxs2=reshape(gr_wxs,numel(gr_wxs)/num_clu,[])';
gr_nxs2=reshape(gr_nxs,numel(gr_nxs)/num_clu,[])';
gr_wjs2=reshape(gr_wjs,numel(gr_wjs)/num_clu,[])';
gr_njs2=reshape(gr_njs,numel(gr_njs)/num_clu,[])';

gr2=[mat2cell(gr_whg',ones(1,p));
    mat2cell(gr_nhg',ones(1,p));
    mat2cell(gr_whs2,ones(1,num_clu));
    mat2cell(gr_nhs2,ones(1,num_clu));
    mat2cell(gr_nhds2,ones(1,num_clu));
    mat2cell(gr_wxs2,ones(1,num_clu));
    mat2cell(gr_nxs2,ones(1,num_clu));
    mat2cell(gr_wjs2,ones(1,num_clu));
    mat2cell(gr_njs2,ones(1,num_clu))];

Gp=tenseg_str_gp(gr,C);    %generate group matrix1
Gp2=tenseg_str_gp(gr2,C);    %generate group matrix2
S=Gp2';                      % clustering matrix
%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
A_1ac=A_1a*S';          %equilibrium matrix CTS
A_2ac=A_2a*S';          %equilibrium matrix CTS
l_c=S*l;
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
index_gp=[6];                   % number of groups with designed force
fd=20^2/3*1000;                        % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
t_c=pinv(S')*t;
q_c=pinv(S')*q;

%% cross sectional design
index_b=find(t_c<0);              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings
[A_b,A_s,A_c,A,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
E=S'*E_c;     %Young's modulus CTS
A=S'*A_c;     % Cross sectional area CTS
l0=(t+E.*A).\E.*A.*l;
mass=S'*rho.*A.*l0;
% % Plot the structure with radius
% R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.03,.1],r_b);
% R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.03,.1],r_s);
% R3Ddata.Nradius=0.1*ones(nn,1);
% tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);
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
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
plot_mode(V_mode,omega,N,Ia,C_b,C_s,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.8,saveimg);
%% control part
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
rate=0.8;
[N_1,C_b,C_s,C] =N_cable_dome(R,rate,p,m,h,beta);
rate=0.8;
[N_2,C_b,C_s,C] =N_cable_dome(R,rate,p,m,h,beta);

ind_n_ct=a;
n_ct1=N_1(:);n_ct1=n_ct1(ind_n_ct);     % control target start configuration
n_ct2=N_2(:);n_ct2=n_ct2(ind_n_ct);     % control target start configuration
% ind_n_ct=[7:8];n_ct1=[1.8;0.3];n_ct2=[1.8;0.3];
Ic=transfer_matrix(ind_n_ct,a);         %transfer matrix for control coordinate
[n_ct_t,n_ct_dt,n_ct_ddt]=coord_vel_acc(tspan,n_ct1,n_ct2);     %nodal coordinate of control target

a_c=2*sqrt(50);    % coefficient in error dynamics(damping term)
b_c=50;            % coeffieient in error dynamics(stiffness term)

%% choose active, passive members
ind_act=[3:6]; ind_pas=[1,2];
% ind_act=[1:6]; ind_pas=[];
I_act=transfer_matrix(ind_act,1:size(S,1));
I_pas=transfer_matrix(ind_pas,1:size(S,1));
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
%% shape control analysis
% solve dynamic equation with closed loop control
data_out=shape_control_CTS(data);
% time history of structure in shape control
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate 
l_t=data_out.l_t;   %time history of members' length 
l0_t=data_out.l0_t; %time history of members' rest length 
nd_t=data_out.nd_t;   %time history of nodal coordinate

%% plot member force 
tenseg_plot_result(out_tspan,t_t([1:6],:),{'1','2','3','4','5','6'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t([3*3-2,3*3-1],:),{'3X','3Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);

%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,90])

%% save output data
if savedata==1
    save (['Dbar_CTS_',material{1},'.mat']);
end
%% make video of the dynamic
name=['CTS_Tbar',material{1},'_slack_',num2str(material{2})];
% % tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t,C_b,C_s,[],50,name,savevideo,material{2})

%output data to tecplot
tenseg_tecplot(C,n_t,t_t,interp1([min(radius),max(radius)],[0.2,0.8],radius));
