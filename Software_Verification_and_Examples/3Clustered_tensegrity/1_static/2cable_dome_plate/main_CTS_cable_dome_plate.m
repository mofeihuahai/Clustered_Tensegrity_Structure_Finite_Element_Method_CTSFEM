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
material{2}=1; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

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
m=1;   %number of circle of the vertical bars
h=0.15*2*R;   %hight of the dome
beta=30*pi/180*ones(m,1);    %all angle of diagonal string
% [N,C_b,C_s,C] =generat_cable_dome(R,p,m,h,beta);

rate=0.8;
[N,C_b,C_s,C] =N_cable_dome_plate(R,rate,p,m,h,beta);
[ne,nn]=size(C);% ne:No.of element;nn:No.of node

tenseg_plot(N,C_b,C_s);
axis off
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

%% trajectory design
n_t=[];
substep=30;
rate_0=linspace(0.2,0.8,substep);
for i=1:substep
rate=rate_0(i);
[N,C_b,C_s,C] =N_cable_dome(R,rate,p,m,h,beta);

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

% plot_mode(K_mode,k,N,Ia,C_b,C_s,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.8,saveimg);
%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0.00;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz

%%record the result
n_t=[n_t,N(:)];
num_prestress_mode(i)=size(V2,2);
minimal_stiffness(i)=min(k);
minimal_frequency(i)=min(omega);
t_t(:,i)=t_c;
end
%% plot the prestress mode number and minimal stiffness
figure
plot(rate_0,num_prestress_mode,'-o','linewidth',2);
set(gca,'fontsize',18);grid on;
xlabel('deploying rate','fontsize',18,'Interpreter','latex');
ylabel('N0. of prestress mode','fontsize',18,'Interpreter','latex');
figure
plot(rate_0,minimal_stiffness,'-o','linewidth',2);
set(gca,'fontsize',18); grid on;
xlabel('deploying rate','fontsize',18,'Interpreter','latex');
ylabel('minimal stiffness eigenvalue','fontsize',18,'Interpreter','latex');
figure
plot(rate_0,minimal_frequency,'-o','linewidth',2);
set(gca,'fontsize',18); grid on;
xlabel('deploying rate','fontsize',18,'Interpreter','latex');
ylabel('minimal natural frequency','fontsize',18,'Interpreter','latex');



%% plot member force 
tenseg_plot_result(1:substep,t_t(:,:),{'whs','nhs','nhds','wxs','nxs','wjs','njs'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
% tenseg_plot_result(1:substep,t_t(:,:),{'IHS','nhs','nhds','wxs','nxs','wjs','njs'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y
tenseg_plot_result(rate_0,n_t([3*4-2],:),{'4X'},{'Deployment ratio','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% plot configuration
for i=round(linspace(1,substep,3))
tenseg_plot_CTS(reshape(n_t(:,i),3,[]),C,[gr_whg,gr_nhg],S,[],[],[45,30]);
% grid on;
axis off;
end
    %% make video of the trajectory
name=['CTS_cable_dome_trajectory'];
% % tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t,C_b,C_s,[],50,name,savevideo,material{2});