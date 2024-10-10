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
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% static analysis set
substep_1=20;
substep=20;                                     %荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=1;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 

%dynamic analysis set
dt=1e-4;               % time step in dynamic simulation
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=1;                   % final time of dynamic simulation
out_dt=1e-4;            % output data interval(approximately, not exatly)

amplitude=0;            % amplitude of external force of ground motion 
period=0.5;             %period of seismic

%% N C of the structure
% Manually specify node positions of double layer prism.
% N=[0 0 0;1 1 0;2 0 0;1 -1 0]';
alpha=0;beta=0;gamma=0;
R=[cos(alpha)*cos(beta),cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma),cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma);
   sin(alpha)*cos(beta),sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma),sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma);
   -sin(beta),cos(beta)*sin(gamma),cos(beta)*cos(gamma)];
width=0.1;
% N=width*[2,1,3,0,2,4;0,sqrt(3),sqrt(3),0,0,0;0,0,0,0,0.05,0];       %nodal coordinate
N=width*[2,1,3,0,2,4;0,sqrt(3),sqrt(3),2*sqrt(3),2*sqrt(3),2*sqrt(3);0,0,0,0.05,0.01,0.05];       %nodal coordinate

% Manually specify connectivity indices.
C_t_in = [1,4;1,6;4,6];  % truss
C_l_in = [1,2;1,3;2,4;3,6;4,5;5,6;2,5;3,5;2,3];  % panel line

% Convert the above matrices into full connectivity matrices.
C_t = tenseg_ind2C(C_t_in,N);%%
C_l = tenseg_ind2C(C_l_in,N);
C=[C_l;C_t];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
[npl,~]=size(C_l);
[nt,~]=size(C_t);
% Plot the structure to make sure it looks right
tenseg_plot(N,C_t,C_l);

%% connectivity of triangle element Ca
% Ca can be written in a function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Ca=[1,2,2,3;3,5,3,6;2,4,5,5];
% Ca=generate_Ca(C_in,N);
% Ca=zeros(3,1)
[~,np]=size(Ca);        % ne:No.of element;np:No.of plate

% plot the origami configuration
tenseg_plot_ori(N,C_l,C_t,[],[],[],[],[],[] ,[],Ca);
%% Boundary constraints
 pinned_X=[1,2,3]'; pinned_Y=[1,2,3]'; pinned_Z=[1,2,3]';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%% C_pl_bar
 C_pl=zeros(np,3);  %Cpl
 C_pl=[1,2,9;3,7,5;9,8,7;8,4,6];
 [C_pl_bar,C_pl_bar_i]=tenseg_ind2C_bar(C_pl,C_l,Ca);

%%  不考虑节点顺序的连接关系矩阵
[C_pn_bar,n_pn_i,C_pn_i]=tenseg_ind2C_membrane(Ca,N);   %C_pn,i
% n=N(:)
% C_pn_bar=[];               
% for i=1:np
%   n_pn{i}=kron(cell2mat(C_pn(i)),eye(3))*n;  %n_pn,i
%   C_pn_bari=[1,1,1]*cell2mat(C_pn(i));       
%   C_pn_bar=[C_pn_bar;C_pn_bari];  %不考虑节点顺序的连接关系矩阵
% end

%% 板的质量
rho_p=0.1;
t_p=6e-3*ones(np,1);        % thickness of panel
A_p=tenseg_A_p(Ca,N);       % area of panel
mass_p=rho_p.*A_p.*t_p;
M_p=1/12*kron((C_pn_bar'*diag(mass_p)*C_pn_bar+diag(diag(C_pn_bar'*diag(mass_p)*C_pn_bar))),eye(3));

%% Group/Clustered information 
%generate group index
% gr=[];
gr={[10,11,12]};     % number of elements in one group
Gp=tenseg_str_gp(gr,C);    %generate group matrix
% S=eye(ne);                  % no clustering matrix
S=Gp';                      % clustering matrix as group matrix

gr_tc={[1,2,3]};     % number of elements in one group
Gp_tc=tenseg_str_gp(gr_tc,C_t);    %generate group matrix
% S=eye(ne);                  % no clustering matrix
S_tc=Gp_tc';                      % clustering matrix as group matrix
% tenseg_plot_CTS(N,C,[1,2],S);
%% self-stress design
%Calculate equilibrium matrix and member length
[A_2t,A_2ta,A_2tc,A_2tac,A_2l,A_2la,l_t,l_tc,l_l,l,l_c]=tenseg_equilibrium_matrix_truss_lines(N,C,C_t,C_l,S,S_tc,Ia);
for i=1:np
A_2p_i{i}=kron(R',eye(3))*kron(cell2mat(C_pn_i(i)),eye(3))*A_2l*cell2mat(C_pl_bar_i(i))';
B_lp_i{i}=A_2p_i{i}';

B_epsilon_i{i}=[cell2mat(B_lp_i(i))*kron(eye(3),[1,0,0;0,0,0;0,0,0])*cell2mat(n_pn_i(i)),cell2mat(B_lp_i(i))*kron(eye(3),[0,0,0;0,1,0;0,0,0])*cell2mat(n_pn_i(i)),cell2mat(B_lp_i(i))*kron(eye(3),[0,0,0;1,0,0;0,0,0])*cell2mat(n_pn_i(i))];
end
B_epsilon=blkdiag(B_epsilon_i{:});
E_p=1e6*ones(np,1);     % Young's modulus of panel
mu=0.5;                 
% D=zeros(3*np,3*np)
D=[diag(E_p/(1-mu^2)),diag(mu*E_p/(1-mu^2)),zeros(np,np);diag(mu*E_p/(1-mu^2)),diag(E_p/(1-mu^2)),zeros(np,np);zeros(np,np),zeros(np,np),diag(E_p/2*(1+mu))];
l0_l=0.19*ones(npl,1);
Delta_l_l=l_l-l0_l;
pVp_pn=A_2l*(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*Delta_l_l;  %partial_Vp/partial_n
t_l=(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*Delta_l_l;

E_t=1e6*ones(nt,1);    % Young's modulus of truss
A_t=1e-6*ones(nt,1)    % area of truss

l0_t=0.39*ones(nt,1);
% l0_tc=S_tc*l0_t;

l0=[l0_l;l0_t];
l0_c=S*l0;
Delta_l_t=l_t-l0_t;
% Delta_l_tc=l_tc-l0_tc;

% pVt_pn=A_2t*diag(E_t)*diag(A_t)*diag(l0_t.^-1)*Delta_l_t;   %partial_Vt/partial_n
pVt_pn=A_2tc*pinv(S_tc')*diag(E_t)*diag(A_t)*diag(l0_t.^-1)*Delta_l_t;   %partial_Vt/partial_n

t_t=diag(E_t)*diag(A_t)*diag(l0_t.^-1)*Delta_l_t;
t_tc=pinv(S_tc')*t_t;


% pV_pn=[A_2t,A_2l]*[t_t;t_l];  %  partial_V/partial_n=partial_Vt/partial_n+partial_Vp/partial_n
pV_pn=[A_2tc,A_2l]*[t_tc;t_l];  %  partial_V/partial_n=partial_Vt/partial_n+partial_Vp/partial_n


Ke_p=A_2l*(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*A_2l';
Kg_p=kron(C_l'*diag(l_l.^-1)*diag(t_l)*C_l,eye(3))-A_2l*diag(l_l.^-1)*diag(t_l)*A_2l';
Kt_p=Ke_p+Kg_p;   %partial_2Vp/partial_n*partial_n'
% p2Vp_pn2=A_2l*(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*A_2l'+kron(C_l'*diag(l_l.^-1)*diag(t_l)*C_l,eye(3))-A_2l*diag(l_l.^-1)*diag(t_l)*A_2l';

Ke_t=A_2tc*pinv(S_tc')*diag(E_t)*diag(A_t)*diag(l0_t.^-1)*pinv(S_tc)*A_2tc'
Kg_t=kron(C_t'*diag(l_t.^-1)*diag(S_tc*t_tc)*C_t,eye(3))-A_2tc*pinv(S_tc')*diag(l_t.^-1)*diag(S_tc'*t_tc)*pinv(S_tc)*A_2tc';
Kt_t=Ke_t+Kg_t;   %partial_2Vt/partial_n*partial_n'
% p2Vt_pn2=A_2t*diag(E_t)*diag(A_t)*diag(l0_t.^-1)*A_2t'+kron(C_t'*diag(l_t.^-1)*diag(t_t)*C_t,eye(3))-A_2t*diag(l_t.^-1)*diag(t_t)*A_2t';
% p2Vt_pn2=A_2tc*pinv(S_tc')*diag(E_t)*diag(A_t)*diag(l0_t.^-1)*pinv(S_tc)*A_2tc'+kron(C_t'*diag(l_t.^-1)*diag(S_tc*t_tc)*C_t,eye(3))-A_2tc*pinv(S_tc')*diag(l_t.^-1)*diag(t_t)*pinv(S_tc)*A_2tc';

K_t=Kt_p+Kt_t;   %partial_2V/partial_n*partial_n'
% p2V_pn2=p2Vp_pn2+p2Vt_pn2;

A_2tcl=[A_2tc,A_2l];
A_2tcla=Ia'*A_2tcl;
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_2tcla);

%external force in equilibrium design
% w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
% index_gp=[1];                   % number of groups with designed force
% fd=-1e2;                        % force in bar is given as -1000
% [q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
% [q_c,t_c,q,t]=tenseg_prestress_design_f(S,l,l_c,A_2ac,V2,w0a,index_gp,fd);

%% cross sectional design
% index_b=find(t_c<0);              % index of bar in compression
% index_s=setdiff(1:size(S,1),index_b);	% index of strings
% [A_b,A_s,A_c,A,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,t_p,hollow_solid);
% E=S'*E_c;     %Young's modulus CTS
% A=S'*A_c;     % Cross sectional area CTS
% l0=(t+E.*A).\E.*A.*l;
% l0_c=S*l0;
% mass=S'*rho.*A.*l0;
% % Plot the structure with radius
% R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.03,.1],r_b);
% R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.03,.1],r_s);
% R3Ddata.Nradius=0.1*ones(nn,1);
% tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);
index_l=(1:1:npl);              % index of bar in compression
index_t=setdiff(1:size(S,1),index_l);	% index of strings

%% tangent stiffness matrix
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS2(Ia,C,q,A_2ac,E_c,A_c,l0_c);
% plot the mode shape of tangent stiffness matrix

Kt_aa=Ia'*K_t*Ia;
[K_mode,D1] = eig(Kt_aa);         % eigenvalue of tangent stiffness matrix
k=diag(D1);  
[k_sort,I]=sort(k);
K_mode_sort=K_mode(:,I);
% plot the mode shape of tgent stiffness matrix
num_plt=1:9;
plot_mode_ori(K_mode_sort,k_sort,N,Ia,C_l,C_t,[],[],l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.1 ,saveimg,3,Ca);

%% input file of ANSYS
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');

%% mass matrix and damping matrix

% damping matrix
% d=1;     %damping coefficient
% D=d*2*max(sqrt(mass.*E.*A./l))*eye(3*nn);    %critical damping


%% mode analysis
[V_mode,D1] = eig(Kt_aa,Ia'*M_p*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
plot_mode_CTS(V_mode,omega,N,Ia,C,[1,2],S,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,0.1,saveimg);

%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[];   %external force in Z 
ind_dnb=[]; dnb0=[];
%ind_dl0_c=[]; dl0_c=[];
ind_dl0_c=[10];dl0_c=[-1.19];
% ind_theta_0=[]; dtheta_0=[];        % initial angel change with time
p=numel(dl0_c);  %加载次数
substep = substep*p;
% [w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress_membrane(substep,substep_1,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;0;9.8],C,M_p,p);

%% Step1: statics: equilibrium calculation
% input data
data.N=N; data.C=C; data.C_l=C_l;data.C_t=C_t;data.C_pl_bar=C_pl_bar;
data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;data.S_tc=S_tc;
data.E=E_p; data.A_p=A_p;data.t_p=t_p; data.index_l=index_l; data.index_t=index_t;data.B_epsilon=B_epsilon;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;           % forced movement of pinned nodes
data.l0_t=l0_ct;            % forced change of rest length
data.substep=substep;    % substep
%data.InitialLoadFactor=0.001;
%data.MaxIcr=1000;
%.LoadType='Force'; % 'Force' or 'Displacement'
%data.StopCriterion=@(U)(norm(U)>0.5);

% nonlinear analysis
% data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
data_out1=static_solver_CTS_membrane(data);
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
N_out=data_out1.N_out;
%% plot member force 
tenseg_plot_result(1:substep,t_t([1,2,3,5,6],:),{'1','2','3-4','5','6'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,0.5*(n_t([3*3-1],:)-n_t([3*4-1],:)-2),{'3Y'},{'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% Plot configuration
for i=round(linspace(1,substep,3))
tenseg_plot_CTS(reshape(n_t(:,i),3,[]),C,[1,2],S);
axis off;
end
tenseg_plot( reshape(n_t(:,end),3,[]),C_t,C_l,[],[],[0,90])

%% save output data
if savedata==1
    save (['cable_net_CTS_static','.mat']);
end
%% make video of the static
name=['Tbar_static2'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,material{2})
tenseg_video_CTS(n_t,C,[1,2],S,[],[],[],[],[],[],t_t,[],min(substep,50),tf,name,savevideo)

%% Step 2: dynamics:change rest length of strings
% time step
if auto_dt
dt=pi/(8*max(omega)); 	% time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;
tspan1=0:dt:tf/2;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and 
ind_w=[];w=[];
ind_dl0_c=[3,4,5];dl0_c=[-2,0.5,0.5];
% ind_dl0_c=[2]; dl0_c=[-120];
[w_t,l0_ct]=tenseg_load_prestress_CTS(tspan1,ind_w,w,'ramp',ind_dl0_c,dl0_c,l0_c,gravity,[0;0;0],C,mass);
w_t=[w_t,w_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];   % second half no change of boundary info
l0_ct=[l0_ct,l0_ct(:,end)*ones(1,numel(tspan)-numel(tspan1))];% boundary node motion info
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




%% plot member force 
tenseg_plot_result(out_tspan,t_t([1:5],:),{'1','2','3','4','5'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t([3*3-1],:)-n_t([3*4-1],:)-2,{'3Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
tenseg_plot( reshape(n_t(:,end),3,[]),C_t,C_l,[],[],[0,90])

%% save output data
if savedata==1
    save (['cable_net_CTS_dynamic',num2str(tf),'.mat']);
end
%% make video of the dynamic
name=['Tbar_dynamic_CTS',num2str(tf)];
% tenseg_video(n_t,C_b,C_s,[],min(numel(out_tspan),50),name,savevideo,material{2})
tenseg_video_CTS(n_t,C,[1,2],S,[],[],[],[],[],[],t_t,[],min(numel(out_tspan),50),tf,name,savevideo)

%output data to tecplot
tenseg_tecplot(C,n_t,t_t,interp1([min(radius),max(radius)],[0.2,0.8],radius));