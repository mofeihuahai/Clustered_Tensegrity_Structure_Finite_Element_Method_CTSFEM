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
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% static analysis set
p=1;
substep_1 = 20;             %荷载子步
substep=20*p;                                     %荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=1;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 

%dynamic analysis set
dt=2e-2;               % time step in dynamic simulation
auto_dt=1;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=2;                   % final time of dynamic simulation
out_dt=1e-3;            % output data interval(approximately, not exatly)

amplitude=0;            % amplitude of external force of ground motion 
period=0.5;             %period of seismic


%% %% N C of the structure
% Manually specify node positions
R=10; h=10; p=5; level=3;        % radius; height; number of edge, level
% AB=BC solve rotation angle
syms x
AB=2*R*sin(pi/p);
BC=norm([h,2*R*sin(x/2)]);
beta_x=eval(solve(BC==AB,x));
beta=beta_x(1);
% give rotation angle
% beta=15*pi/180; 	% rotation angle

angle=kron(ones(1,level+1),2*pi*((1:p)-1)./p)+kron(0:level,beta*ones(1,p));
N=R*[cos(angle); sin(angle); zeros(1,p*(level+1))]+kron(0:level,[0;0;h]*ones(1,p));

% Manually specify connectivity indices.
C_b_in_h=[(1:p*(level+1))',kron(ones(level+1,1),[2:p,1]')+kron((0:level)',p*ones(p,1))];% horizontal bar
C_b_in_v=[(1:p*level)',(p+1:(level+1)*p)'];% vertical bar
C_b_in_d=[(1:p*level)',kron(ones(level,1),[2:p,1]')+kron((1:level)',p*ones(p,1))];% diagonal bar
% strings
C_s_in=[(1:p*level)',kron(ones(level,1),[p,1:p-1]')+kron((1:level)',p*ones(p,1))]; %diagonal string

C_b_in = [C_b_in_h;C_b_in_v;C_b_in_d];   % This is indicating the bar connection
C_in=[C_b_in;C_s_in];
% Convert the above matrices into full connectivity matrices.
C = tenseg_ind2C(C_in,N);
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node

C_b=tenseg_ind2C(C_b_in,N);;C_s=tenseg_ind2C(C_s_in,N);
n_b=size(C_b,1);n_s=size(C_s,1);        % number of bars and string

% connectivity matrix for plot
C_bar_in=[];      % real bars
C_bar=tenseg_ind2C(C_bar_in,N);
C_rot_h_in=C_b_in;             %rotational hinges
C_rot_h=tenseg_ind2C(C_rot_h_in,N);

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% define hinge, rigid hinge
% C_in_h is the connectivity of higes, can be written in a function!!!!!!!!!
C_in_h=[C_b_in_h(p+1:(level)*p,:);C_b_in_v;C_b_in_d];
n_h=size(C_in_h,1);         % number of hinge

[~,index_h]=ismember(C_in_h,C_in,'rows');   % index_h is the index number of hinge
index_rh=[];index_rh_in_h=[];
% [~,index_rh]=ismember(C_in_2,C_in,'rows');   % index_h is the index number of rigid hinge
% [~,index_rh_in_h]=ismember(C_in_2,C_in_h,'rows');   % index_h is the index of rigid hinge in all hinge

C_h=tenseg_ind2C(C_in(setdiff([1:ne]',index_rh),:),N);     % connectivity matrix of all edges
C_rig_h=tenseg_ind2C(C_in(index_rh,:),N); % connectivity matrix of rigid edges
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% connectivity of triangle element Ca
% Ca can be written in a function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Ca=kron(ones(1,level),[[1:p;[2:p,1];[p+2:2*p,p+1]],[1:p;[p+2:2*p,p+1];p+1:2*p]])+kron(0:level-1,p*ones(3,p*2));
% Ca=generate_Ca(C_in,N);
% Ca=zeros(3,1)
[~,np]=size(Ca);        % ne:No.of element;np:No.of plate

% plot the origami configuration
tenseg_plot_ori(N,C_b,C_s,[],C_rig_h,[],[],[],[] ,[],Ca);
%%  不考虑节点顺序的连接关系矩阵
C_paper=tenseg_ind2C_paper(Ca,N);

%% 板的质量
rho_p=1;
A_p=tenseg_A_p(Ca,N);
mass_p=rho_p.*A_p.*thick;

M_p=1/12*kron((C_paper'*diag(mass_p)*C_paper+diag(diag(C_paper'*diag(mass_p)*C_paper))),eye(3));
%% transformation matrix from element to structure

E_n=cell(1,n_h);            %transformation matrix from element node to total node
node_in_hinge=zeros(n_h,4);
I=eye(3*nn);

for i=1:n_h
node2=C_in_h(i,1);  % start node of the hinge
node3=C_in_h(i,2);  % end node of the hinge
for j=1:np
    if (node2==Ca(1,j)&node3==Ca(2,j))|(node2==Ca(2,j)&node3==Ca(3,j))|(node2==Ca(3,j)&node3==Ca(1,j))
        node1=setdiff(Ca(:,j),[node2;node3]);
    elseif (node2==Ca(2,j)&node3==Ca(1,j))|(node2==Ca(3,j)&node3==Ca(2,j))|(node2==Ca(1,j)&node3==Ca(3,j))
        node4=setdiff(Ca(:,j),[node2;node3]);
    end
end
node_in_hinge(i,:)=[node1,node2,node3,node4];
E_n{i}=I(:,kron(node_in_hinge(i,:),3*ones(1,3))-kron(ones(1,4),[2,1,0]));
end
E_n_total=cell2mat(E_n);        % transformation matrix of the whole structure
%% Boundary constraints
 pinned_X=[1:p]'; pinned_Y=[1:p]'; pinned_Z=[1:p]';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%% Group/Clustered information 
%generate group index
% gr=[];
gr=[];     % number of elements in one group
Gp=tenseg_str_gp(gr,C);    %generate group matrix
% S=eye(ne);                  % no clustering matrix
S=Gp';                      % clustering matrix as group matrix

%tenseg_plot_CTS(N,C,[1,2],S);
%% equilibrium matrix


% equilibrium matrix of truss
[A_1,A_1c,A_1a,A_1ac,A_2,A_2c,A_2a,A_2ac,l,l_c]=tenseg_equilibrium_matrix_CTS(N,C,S,Ia);
% [A_1,A_1g,A_2,A_2g,l,l_gp]=tenseg_equilibrium_matrix2(N,C,Gp,Ia);

% equilibrium matrix of hinge
[phpn_e,phTpn,theta]=jacobian_ori(node_in_hinge,N,E_n_total);       % jacobian matrix
%A_o=[A_c,phTpn];
%A_o_a=Ia'*A_o;
A_co=[A_2c,phTpn];
A_co_a=Ia'*A_co;


%% SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_co_a);         % equilibrium of truss with hinge

%external force in equilibrium design
% w0=zeros(numel(N),1); w0a=Ia'*w0;
% 
%prestress design
% index_gp=[1];                   % number of groups with designed force
% fd=-1e2;                        % force in bar is given as -1000
% [q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
% [q_c,t_c,q,t]=tenseg_prestress_design_f_ORI(S,l,l_c,A_co_a,V2,w0a,index_gp, fd);
%% self-stress design (of truss)
t=zeros(ne,1);      %member force
t_c=S*t;                      
q=t./l;             % force density-
%% self-stress design (of hinge)
M=zeros(n_h,1);    

%% cross sectional design (of truss)
A=1e-6*ones(ne,1);
E=1e9*ones(ne,1);
A_c=pinv(S')*A;
E_c=pinv(S')*E;
l0=(t+E.*A).\E.*A.*l;
l0_c=S*l0;
index_b=[1:9]';              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings

%%  hinge section design  (of hinge)
k_h=1e-2*1/12*E_c(index_h).*l(index_h)*thick^3;
k_h(index_rh_in_h)=1/12*E_c(index_rh).*l(index_rh)*thick^3;      % increase stiffness of rigid hinge

%% rest length (of truss), initial angle (of hinge)
l0=l;                     %rest length of truss
theta_0=pi*ones(n_h,1);     % initial angle of hinge
% theta_0=[]; % initial angle of hinge
%% tangent stiffness matrix
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
[Kt_aa,Kg_aa,Ke_aa,K_mode_sort,k_sort]=tenseg_stiff_CTS3(Ia,C,S,t_c,A_2a,E_c,A_c,l0,l);

%% tangent stiffness matrix of hinge
%ph2px2 is the hessian matrix of theta to nodal coordinate
[ph2pn2_e,ph2pn2]=hessian_ori(node_in_hinge,N,E_n);         % calculate hessian matrix
G=cell2mat(ph2pn2);

%% tangent stiffness of the whole origami

K_t_oa=Kt_aa+Ia'*(phTpn*diag(k_h)*phTpn'+G*kron(M,eye(3*nn)))*Ia;

[K_mode,D1] = eig(K_t_oa);         % eigenvalue of tangent stiffness matrix
k=diag(D1); 
[k, ind] = sort(k);
K_mode = K_mode(:, ind);
% plot the mode shape of tangent stiffness matrix
num_plt=1:9;
plot_mode_ori(round(K_mode,12),k,N,Ia,C_bar,C_s,C_rot_h,C_rig_h,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.3,saveimg,3,Ca);

%% input file of ANSYS
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');

%% mass matrix and damping matrix
rho=1;
mass=rho.*A.*l;
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0.5;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l))*eye(3*nn);    %critical damping

%% vibration mode analysis
[V_mode,D1] = eig(K_t_oa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
if 0
plot_mode_ori(K_mode,k,N,Ia,C_b,C_s,[],C_rig_h,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg,3,Ca);
end

%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[];   %external force in Z 
ind_dnb=[]; dnb0=[];
%ind_dl0_c=[]; dl0_c=[];
ind_dl0_c=[n_b+1:ne]; dl0_c=-0.8*l0_c(ind_dl0_c)';
ind_theta_0=[]; dtheta_0=[];        % initial angel change with time
p=1;
% p=numel(dl0_c);  %加载次数
substep = substep*p;
% [w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
[w_t,dnb_t,l0_ct,theta_0_t,Ia_new,Ib_new]=tenseg_load_prestress_ori_p_1(substep,substep_1,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,ind_theta_0,dtheta_0,theta_0,l0_c,b,gravity,[0;9.8;0],C,M_p,p);
% figure
% plot(1:substep,l0_ct)
tenseg_plot_result(1:substep,l0_ct([51],:),{'string'},{'Substep','Rest length (m)'},'plot_coordinate.png',saveimg);
grid on;
%% Step1: statics: equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;           % forced movement of pinned nodes
data.l0_t=l0_ct;            % forced change of rest length
data.theta_0_t=theta_0_t;   % forced change of initial angle
data.k_h=k_h;               % stiffness of hinge
data.E_n=E_n;               % transfer matrix from matrix to structure
data.node_in_hinge=node_in_hinge;       % node in triangle element in hinge
data.substep=substep;    % substep
%data.InitialLoadFactor=0.001;
%data.MaxIcr=1000;
%.LoadType='Force'; % 'Force' or 'Displacement'
%data.StopCriterion=@(U)(norm(U)>0.5);

% nonlinear analysis
% data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
data_out1=static_solver_CTS_ORI(data);
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
M_t=data_out1.M_out;
N_out=data_out1.N_out;
%% plot member force 
tenseg_plot_result(1:substep,t_t([65],:),{'string'},{'Substep','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([58,59,60],:),{'x20','y20','z20'},{'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% plot hinge moment
 tenseg_plot_result(1:substep,M_t([1,11],:),{'hinge1','hinge2'},{'Substep','Moment (N·m)'},'plot_hinge_moment.png',saveimg);

%% Plot configuration
% for i=round(linspace(1,substep,3))
% tenseg_plot_CTS(reshape(n_t(:,i),3,[]),C,[1,2],S);
% axis off;
% end
% tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,90])


 j=linspace(1e-5,1,10);
for i=1:10
    num=ceil(j(i)*size(n_t,2));
tenseg_plot_ori(reshape(n_t(:,num),3,[]),C_bar,C_s,C_rot_h,C_rig_h,[],[],[30,30],[] ,[],Ca);
%  axis off;
end
%% save output data
if savedata==1
    save (['cable_net_CTS_static','.mat']);
end
%% make video of the dynamic
icrm=size(n_t,2); 
name=['krseling'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video_ori(n_t,C_b,C_s,C_h,C_rig_h,Ca,[],min(icrm,50),name,savevideo,[])

%% make video of the static
% name=['CTS_ORI'];
% % tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,material{2})
% tenseg_video_CTS(n_t,C,[1,2],S,[],[],[],[],[],[],t_t,[],min(substep,50),tf,name,savevideo)

%% Step 2: dynamics:change rest length of strings
% time step 
if auto_dt
dt=pi/(8*max(omega)); % time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;
tspan1=0:dt:tf/2;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and 
ind_w=[];w=[];
ind_dl0_c=[n_b+1:ne]; dl0_c=-0.8*l0_c(ind_dl0_c)';
ind_theta_0=[]; dtheta_0=[];        % initial angel change with time

[w_t,l0_ct,theta_0_t]=tenseg_load_prestress_CTS_ORI(tspan,ind_w,w,'ramp',ind_dl0_c,dl0_c,l0_c,ind_theta_0,dtheta_0,theta_0,gravity,[0;0;9.8],C,M_p);

% w_t=[w_t,w_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];   % second half no change of boundary info
% l0_ct=[l0_ct,l0_ct(:,end)*ones(1,numel(tspan)-numel(tspan1))];
% boundary node motion info
[~,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan1,a,b,'vib_force',gravity,[0;0;9.8],C,mass,[1,2],amplitude,period);
dnb_t=[dnb_t,dnb_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];
dnb_d_t=[dnb_d_t,dnb_d_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];
dnb_dd_t=[dnb_dd_t,dnb_dd_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];
% give initial speed of free coordinates
n0a_d=zeros(numel(a),1);                    %initial speed in X direction
    %% dynamics calculation

% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia; data.Ib=Ib;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.n_h=n_h;
data.w_t=w_t;           % external force
data.theta_0_t=theta_0_t;   % forced change of initial angle
data.k_h=k_h;               % stiffness of hinge
data.E_n=E_n;               % transfer matrix from matrix to structure
data.node_in_hinge=node_in_hinge;       % node in triangle element in hinge
data.dnb_t=dnb_t; data.dnb_d_t=dnb_d_t;  data.dnb_dd_t=dnb_dd_t; % forced movement of pinned nodes
data.l0_t=l0_ct;         % forced movement of pinned nodes
data.n0a_d=n0a_d;        %initial speed of free coordinates
data.M=M;data.D=D;
data.rho=rho_s;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;

%% dynamic analysis
% solve dynamic equation
data_out=dynamic_solver_CTS_ORI(data);        %solve ODE of dynamic equation
% time history of structure
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate 
l_t=data_out.l_t;   %time history of members' length 
nd_t=data_out.nd_t;   %time history of nodal coordinate

%% make video of the dynamic
icrm=size(n_t,2); 
name=['krseling_dynamic_1'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video_ori(n_t,C_b,C_s,C_h,C_rig_h,Ca,[],min(icrm,50),name,savevideo,[])

%% plot member force 
tenseg_plot_result(downsample(out_tspan,30),downsample(t_t([65],:),30),{'string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y

tenseg_plot_result(downsample(out_tspan,30),[downsample(n_t([58],:),30);downsample(n_t([59],:),30);downsample(n_t([60],:),30)],{'x20','y20','z20'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% plot member force 
tenseg_plot_result(out_tspan,t_t([65],:),{'string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([58,59,60],:),{'x20','y20','z20'},{'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,90])

%% save output data
if savedata==1
    save (['cable_net_CTS_dynamic',num2str(tf),'.mat']);
end

