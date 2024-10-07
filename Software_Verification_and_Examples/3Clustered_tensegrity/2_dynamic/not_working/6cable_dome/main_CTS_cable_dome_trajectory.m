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
% [consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Aluminum');
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
savedata=0;             % save data or not (1) yes (0)no
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

Gp=tenseg_str_gp2(gr,C);    %generate group matrix1
Gp2=tenseg_str_gp2(gr2,C);    %generate group matrix2
S=Gp2';                      % clustering matrix

%% trajectory design
n_t=[];
dvd_num=20;    %number of dividing time
% ratio=linspace(0.2,0.8,dvd_num);
ratio=[linspace(0.2,0.8,dvd_num/2),0.8*ones(1,dvd_num/2)];
% ratio=linspace(0.01,0.99,100);      %deployment ratio
%% prestress design
for i=1:dvd_num
rate=ratio(i);
[N,C_b,C_s,C] =N_cable_dome(R,rate,p,m,h,beta);

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
l_c=S*l;
A_1ac=A_1a*diag(l.^-1)*S'*diag(l_c);          %equilibrium matrix CTS
A_2ac=A_2a*S';          %equilibrium matrix CTS

%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
% index_gp=[6];                   % number of groups with designed force
% fd=20^2/3*1000;                        % force in bar is given as -1000
index_gp=[2];                   % nhg with designed force
fd=-5000;                        % force in bar is given as -5000N
% index_gp=[6];                   % number of groups with designed force
% fd=20^2/3*1000;                        % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
t_c=pinv(S')*t;
q_c=pinv(S')*q;
%% cross sectional area design

index_b=find(t_c<0);              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings
[A_b,A_s,A_c,A,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
E=S'*E_c;     %Young's modulus CTS
A=S'*A_c;     % Cross sectional area CTS

%% record the result
num_prestress_mode(i)=size(V2,2);
t_t(:,i)=t_c;
A_t(:,i)=A;
A_c_t(:,i)=A_c;

end
%% cross sectional area redesign
A_c=max(A_c_t,[],2);  %choose the maximum area needed
A=max(A_t,[],2);  %choose the maximum area needed


%% evaluate the stiffness and vibration mode
for i=1:dvd_num
rate=ratio(i);
[N,C_b,C_s,C] =N_cable_dome(R,rate,p,m,h,beta);

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
l_c=S*l;
A_1ac=A_1a*diag(l.^-1)*S'*diag(l_c);          %equilibrium matrix CTS
A_2ac=A_2a*S';          %equilibrium matrix CTS

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
%% rest length and mass
l0=(t+E.*A).\E.*A.*l;
l0_c=S*l0;
mass=S'*rho.*A.*l0;

%% tangent stiffness matrix
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS2(Ia,C,q,A_2ac,E_c,A_c,l0_c);
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS3(Ia,C,S,t_c,A_2a,E_c,A_c,l0,l);
% plot the mode shape of tangent stiffness matrix
num_plt=1:4;
%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0.00;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
if i==0
plot_mode(V_mode,omega,N,Ia,C_b,C_s,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)',num_plt,1,saveimg,[0,50]);
end
%% record the result
n_t=[n_t,N(:)];
num_prestress_mode(i)=size(V2,2);
minimal_stiffness(i)=min(k);
minimal_frequency(i)=min(omega);
l0_t(:,i)=l0;
end
l0_c_t=S*l0_t;        %clusterd string length
t_gp_t=pinv(Gp)*S'*t_t;    % members force in group 
%% Output control target
ind_n_ct=a;
%% save data for dynamic analysis
if savedata
save('input_dyn.mat','l0_c_t','A','A_c');
save('output_static.mat');
end
%% plot the prestress mode number and minimal stiffness
plot_num=1:dvd_num;
% plot_num=1:100;
figure
plot(ratio(plot_num),num_prestress_mode(plot_num),'-o','linewidth',2);
set(gca,'fontsize',18);grid on;
xlabel('{\itc}','fontsize',18);   %number of cluster strings in one gruop
ylabel('{\itn_p}','fontsize',18);   %number of prestress mode 
figure
plot(ratio(plot_num),minimal_stiffness(plot_num),'-o','linewidth',2);
set(gca,'fontsize',18); grid on;
xlabel('{\itc}','fontsize',18);
ylabel('\lambda_{min}','fontsize',18); %minimal stiffness eigenvalue
figure
plot(ratio(plot_num),minimal_frequency(plot_num),'-o','linewidth',2);
set(gca,'fontsize',18); grid on;
xlabel('{\itc}','fontsize',18);
ylabel('\itf_{min}','fontsize',18); %minimal natural frequency

%% plot member force 

% tenseg_plot_result(1:90,t_t(:,1:90),{'whs','nhs','nhds','wxs','nxs','wjs','njs'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
tenseg_plot_result(ratio(plot_num),t_gp_t(:,(plot_num)),{'OB','IB','OHS','IHS','ITS','ODS','IDS','ORS','IRS'},{'Deployment ratio','Force(N)'},'plot_member_force.png',saveimg);
grid on;
columnlegend(3, {'OB','IB','OHS','IHS','ITS','ODS','IDS','ORS','IRS'}, 'location','southwest');
ylim([-5,6]*1e4);
%% plot rest length

tenseg_plot_result(ratio(plot_num),l0_c_t([3:3:21,22,23],plot_num),{'OHS','IHS','ITS','ODS','IDS','ORS','IRS','OB','IB'},{'Deployment ratio','Rest length (m)'},'plot_member_length.png',saveimg);
columnlegend(3, {'OHS','IHS','ITS','ODS','IDS','ORS','IRS','OB','IB'}, 'location','northeast');
grid on;
ylim([0,230]);
l0_gp=pinv(Gp)*l0_t;% string length in group
tenseg_plot_result(ratio(1:dvd_num),l0_gp(:,1:dvd_num),{'OB','IB','OHS','IHS','ITS','ODS','IDS','ORS','IRS'},{'Deployment ratio','Length (m)'},'plot_member_length.png',saveimg);
columnlegend(3, {'OB','IB','OHS','IHS','ITS','ODS','IDS','ORS','IRS'}, 'location','southwest');
grid on;
%% Plot nodal coordinate curve X Y
% tenseg_plot_result(1:1000,n_t([3*4-2,3*4],:),{'4X','4Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
tenseg_plot_result(1:dvd_num,n_t([3*1-2,3*2-2,3*1,3*2],:),{'1X','2X','1Z','2Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);

%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
for i=round(linspace(1,numel(ratio),4))
tenseg_plot(reshape(n_t(:,i),3,[]),C_b,C_s,[],[],[0,90]);
% axis off
% view([0,90]);
view([0,30]);
end
%% make video of the trajectory
name=['CTS_cable_dome_trajectory'];
% % tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t,C_b,C_s,[],dvd_num,name,savevideo,material{2});