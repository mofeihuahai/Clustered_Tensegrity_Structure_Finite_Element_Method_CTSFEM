%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%A Clustered Cable Dome Dynamics in deployment %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Aluminum');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% static analysis set
substep=1;                                     %�����Ӳ�
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 

%dynamic analysis set
dt=0.001;               % time step in dynamic simulation
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=1;                   % final time of dynamic simulation
% dt=tf/1e4;               % time step in dynamic simulation
out_dt=0.001;            % output data interval(approximately, not exatly)

amplitude=0;            % amplitude of external force of ground motion 
period=0.5;             %period of seismic
%% load data for dynamics
load('input_dyn.mat');
%% N C of the structure
% Manually specify node positions of double layer prism.
R=50;          %radius
p=12;          %complexity for cable dome
m=2;   %number of circle of the vertical bars
h=0.15*2*R;   %hight of the dome
beta=30*pi/180*ones(m,1);    %all angle of diagonal string
% [N,C_b,C_s,C] =generat_cable_dome(R,p,m,h,beta);

rate=0.2;
[N,C_b,C_s,C] =N_cable_dome(R,rate,p,m,h,beta);
[ne,nn]=size(C);% ne:No.of element;nn:No.of node

tenseg_plot(N,C_b,C_s);
axis off;
view(0,0)
 view(0,45)
 view(2)
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

gr2=[mat2cell(gr_whs2,ones(1,num_clu));
    mat2cell(gr_nhs2,ones(1,num_clu));
    mat2cell(gr_nhds2,ones(1,num_clu));
    mat2cell(gr_wxs2,ones(1,num_clu));
    mat2cell(gr_nxs2,ones(1,num_clu));
    mat2cell(gr_wjs2,ones(1,num_clu));
    mat2cell(gr_njs2,ones(1,num_clu))];

Gp=tenseg_str_gp(gr,C);    %generate group matrix1
Gp2=tenseg_str_gp(gr2,C);    %generate group matrix2
S=Gp2';                      % clustering matrix

tenseg_plot_CTS(N,C,[gr_whg,gr_nhg],S)
% view([0 90])
% view([0 30])
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
% index_gp=[6];                   % number of groups with designed force
% fd=20^2/3*1000;                        % force in bar is given as -1000
index_gp=[2];                   % nhg with designed force
fd=-5000;                        % force in bar is given as -5000N
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
t_c=pinv(S')*t;
q_c=pinv(S')*q;
%% cross sectional design
index_b=find(t_c<0);              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings
[A_b,A_s,~,~,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
E=S'*E_c;     %Young's modulus CTS
A=S'*A_c;     % Cross sectional area CTS
l0=(t+E.*A).\E.*A.*l;
mass=S'*rho.*A.*l0;
% % Plot the structure with radius
% R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.03,.1],r_b);
% R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.03,.1],r_s);
% R3Ddata.Nradius=0.1*ones(nn,1);
% tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);


tenseg_plot_CTS(N,C,[gr_whg,gr_nhg],S,[],[],[],[],[],t,[])

%% tangent stiffness matrix
num_plt=[];%1:4;
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
% plot the mode shape of tangent stiffness matrix
plot_mode(K_mode,k,N,Ia,C_b,C_s,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.8,saveimg);
%% input file of ANSYS
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0.01;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
plot_mode(V_mode,omega,N,Ia,C_b,C_s,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.8,saveimg);


%% dynamics:change rest length of strings
% time step
if auto_dt
dt=pi/(8*max(omega)); 	% time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and rest length
ind_w=[];w=[];
ind_dl0_c=[]; dl0_c=[];
% ind_dl0_c=[2]; dl0_c=[-120];
[w_t,l0_ct]=tenseg_load_prestress_CTS(tspan,ind_w,w,'ramp',ind_dl0_c,dl0_c,l0_c,gravity,[0;0;0],C,mass);

l0_c_t1=l0_c_t;     %recalculate the rest length
for i=1:numel(tspan)
l0_ct(:,i)=l0_c_t1(:,ceil(i*size(l0_c_t1,2)/numel(tspan)));       %replace rest length of clustered strings , from loading data
end
% boundary node motion info
[~,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'vib_force',gravity,[0;0;9.8],C,mass,[1,2],amplitude,period);

% give initial speed of free coordinates
n0a_d=zeros(numel(a),1);                    %initial speed in X direction
    %% dynamics calculation

% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia; data.Ib=Ib;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;           % external force
data.dnb_t=dnb_t; data.dnb_d_t=dnb_d_t;  data.dnb_dd_t=dnb_dd_t; % forced movement of pinned nodes
data.l0_t=l0_ct;         % forced movement of pinned nodes
data.n0a_d=n0a_d;        %initial speed of free coordinates
data.M=M;data.D=D;
data.rho=rho_s;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;

%% dynamic analysis
% solve dynamic equation
data_out=dynamic_solver_CTS(data);        %solve ODE of dynamic equation
% time history of structure
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate 
l_t=data_out.l_t;   %time history of members' length 
nd_t=data_out.nd_t;   %time history of nodal coordinate

%% load static result
output_sta=load('output_static.mat','n_t','t_gp_t');
n_t_sta=output_sta.n_t;%nodal coordinate
t_gp_sta=output_sta.t_gp_t;%member force


%% plot member force 
% tenseg_plot_result(out_tspan,t_t([1:3,2*p+1:2*p+3],:),{'1','2','3','4','5','6'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);

ratio=linspace(0.2,0.8,numel(out_tspan));
% ratio=linspace(0.2,0.8,100);
 plot_num=1:100;
t_gp_t=pinv(Gp)*t_t;    % members force in group 
tenseg_plot_result(ratio,t_gp_t,{'OB','IB','OHS','IHS','ITS','ODS','IDS','ORS','IRS'},{'{\itc}','{\itf} (N)'},'plot_member_force.png',saveimg);
% tenseg_plot_result(ratio,t_gp_t(:, 10*plot_num),{'OB','IB','OHS','IHS','ITS','ODS','IDS','ORS','IRS'},{'{\itc}','{\itf} (N)'},'plot_member_force.png',saveimg);
grid on;
columnlegend(3, {'OB','IB','OHS','IHS','ITS','ODS','IDS','ORS','IRS'}, 'location','southwest');
%% plot member force comparsion
tenseg_plot_result(out_tspan,[t_gp_t([1,2],:);n_t_sta([1,2],fix(linspace(1,size(n_t_sta,2),size(n_t,2))))]...
    ,{'OB-dynamic','IB-dymnamic','OB-static','IB-static'},{'Time (s)','Force (N)'},'plot_coordinate.png',saveimg);




%% Plot nodal coordinate curve X Y

tenseg_plot_result(out_tspan,[n_t([3*1-2,3*2-2],:);n_t_sta([3*1-2,3*2-2],fix(linspace(1,size(n_t_sta,2),size(n_t,2))))]...
    ,{'1X-dynamic','2X-dymnamic','1X-static','2X-static'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);

%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,90])

%% save output data
if savedata==1
    save (['CTS_cable_dome_','time_',num2str(tf),'s','.mat']);
end
%% make video of the dynamic
name=['CTS_cable_dome_CTS_mark_','time_',num2str(tf)];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t,C_b,C_s,[],min(numel(out_tspan),50),name,savevideo,[])
tenseg_video_CTS(n_t,C,[gr_whg,gr_nhg],S,[],[],[1,2],[],[],[],t_t,[],min(numel(out_tspan),50),tf,name,savevideo)
%output data to tecplot
tenseg_tecplot(C,n_t,t_t,interp1([min(radius),max(radius)],[0.2,0.8],radius));