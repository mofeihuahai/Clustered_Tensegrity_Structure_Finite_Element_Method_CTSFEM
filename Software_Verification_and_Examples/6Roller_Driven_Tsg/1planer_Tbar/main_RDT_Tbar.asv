%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%roller-driven planer T-bar%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
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
substep=1;                                     %substep
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
%% N C of the structure
% Manually specify node positions of double layer prism.
N=[1-0.5 0 0;1 1 0;1 1 0;1+0.5 0 0;1 -1 0]';    

% Manually specify connectivity indices.
 C_s_in = [1 3;3 4;4 5;5 1];  % This is indicating that string connection
 C_b_in = [1 4;2 5];  % Similarly, this is saying bar 1 connects node 1 to node 2,

 % Convert the above matrices into full connectivity matrices.
 C_b = tenseg_ind2C(C_b_in,N);%%
 C_s = tenseg_ind2C(C_s_in,N);
 C=[C_b;C_s];
 [ne,nn]=size(C);        % ne:No.of element;nn:No.of node

 % Plot the structure to make sure it looks right
 tenseg_plot(N,C_b,C_s);
 title('T-bar');
%% radius , side vector for plot
R=[0;0;0.05;0;0]; % radius vector
% mue=[0;0;-1;0;0]; %side vector

 %% unique node
 [N_unique,ia,ic]=unique(N','row','stable');
 I_u_temp=eye(numel(ic));
 I_unique0=zeros(numel(ic),numel(ia));
 for i=1:numel(ic)
     I_unique0(:,ic(i))=I_unique0(:,ic(i))+I_u_temp(:,i);
 end
 I_unique=kron(I_unique0,eye(3));
 %% Boundary constraints-transformation matrix of nodal coordinate
 pinned_X=[]; pinned_Y=[]; pinned_Z=(1:numel(ia))';
 [Ia_min,Ib_min,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,numel(ia));
 E_na=I_unique*Ia_min;
 E_nb=I_unique*Ib_min;
 %% sliding distance vector of pulley and transformation matrix
sld=zeros(nn,1);
% transformation matrix of sliding distance vector
free_sld=3; fixed_sld=setdiff(1:nn,free_sld);
E_s=eye(nn);
E_sa=E_s(:,free_sld);
E_sb=E_s(:,fixed_sld);
%% generalized coordinate q
% generalized coordinate 
q=[N(:);sld];
% transformation matrix of generalized coordinate
E_qa=blkdiag(E_na,E_sa);
E_qb=blkdiag(E_nb,E_sb);

%% Group/Clustered information 
%generate group index
gr={[3,4]};     % number of elements in one group
Gp=tenseg_str_gp(gr,C);    %generate group matrix
S=Gp';                      % clustering matrix

fig_handle=figure
tenseg_plot_CTS(N,C,[1,2],S,fig_handle);
%% equilibrium matrix and SVD
%Calculate equilibrium matrix and member length
H=N*C';                     % element's direction matrix
l=sqrt(sum(H.*H))';         % elements' length
Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H
A_2=kron(C',eye(3))*blkdiag(Cell_H{:})/diag(l);     % equilibrium matrix
A_RDT=E_qa'*[A_2;C'];

%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_RDT);

%%  cross section information
% cross sectional area
A=1e-5*[100*[1;1];ones(4,1)];
% young's modulus
E=[Eb*ones(2,1);Es*ones(4,1)];
% external force
f_ena=zeros(size(E_na,2),1);
f_esa=zeros(size(E_sa,2),1);

% member force design
index_gp=[1];
fd=-1e3;
e_nb=eye(ne);
e_d=e_nb(:,index_gp); % e_d is the matrix to select group of member with designed force
z=(e_d'*V2)\(e_d'*(fd-pinv(A_RDT)*[f_ena;f_esa]));   %self-stress coefficient
t=pinv(A_RDT)*[f_ena;f_esa]+V2*z;

index_b=find(t<0);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings

% member rest length
l0=E.*A.*l./(t+E.*A);
mass=rho_s.*A.*l0;
%% plot with pulley size
fig2=figure
tenseg_plot_RDT(N,C,R,index_b,S,fig2,[],[], [], [],[],[],[]);
% axis off;
%% tangent stiffness matrix
Kn=kron(C'*diag(l.\t)*C,eye(3));

K_T=[Kn+A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2',A_2*diag(E.*A./l0)*C;...
    C'*diag(E.*A./l0)*A_2',C'*diag(E.*A./l0)*C];
K_Tg=blkdiag(Kn-A_2*diag(l.\t)*A_2',zeros(numel(sld)));%geometry stiffness
K_Te=[A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2',A_2*diag(E.*A./l0)*C;...
    C'*diag(E.*A./l0)*A_2',C'*diag(E.*A./l0)*C];

K_T=0.5*(K_T+K_T');

[K_mode1,D1] = eig(E_qa'*K_T*E_qa);         % eigenvalue of tangent stiffness matrix
k=diag(D1);   
[k_sort1,e_nb]=sort(k);
K_mode_sort1=K_mode1(:,e_nb);
% plot mode
num_plt=4:9;
% plot_mode_CTS2(K_mode_sort(1:8,:),k_sort,N,E_na,C,1:2,S,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.2,saveimg,2);
plot_mode_RDT_sa(K_mode_sort1(1:8,:),K_mode_sort1(9,:),k_sort1,N,R,E_na,C,1:2,S,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.2,saveimg,2);

%% tangent stiffness of Clustered Tensegrity CTS
if 1
t_c=S'\t;E_c=S'\E;A_c=S'\A;l0_c=S*l0;A_2ac=E_na'*A_2*S';
[Kt_aa,Kg_aa,Ke_aa,K_mode2,k2]=tenseg_stiff_CTS3(E_na,C,S,t_c,E_na'*A_2,E_c,A_c,l0,l);
  if 0      % this is reduced order from RDT stiffness matrix
K_T11=E_na'*(Kn+A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2')*E_na;
K_T12=E_na'*(A_2*diag(E.*A./l0)*C)*E_sa;
K_T21=E_sa'*(C'*diag(E.*A./l0)*A_2')*E_na;
K_T22=E_sa'*(C'*diag(E.*A./l0)*C)*E_sa;
 
Kt_aa2=K_T11-K_T12/K_T22*K_T21;
[K_mode1,D1] = eig(Kt_aa2);         % eigenvalue of tangent stiffness matrix
k1=diag(D1);   
[k2,e_nb]=sort(k1);
K_mode2=K_mode1(:,e_nb);
  end

K_mode2(:,1)'*Kt_aa*K_mode2(:,1)      %eigenvalue
K_mode2(:,1)'*Kg_aa*K_mode2(:,1)      % eigenvalue of geometry stiffness
K_mode2(:,1)'*Ke_aa*K_mode2(:,1)      % eigenvalue of material stiffness

num_plt=1:8;
% plot_mode_CTS2(K_mode,k,N,Ia,C,1:2,S,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.2,saveimg,3);
% plot_mode(K_mode2,k2,N,E_na,C_b,C_s,l,'natrual vibration',...
%     'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg);
plot_mode_RDT(-K_mode2,k2,N,R,E_na,C,1:2,S,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.2,saveimg,2);


% compare stiffness eigen value
figure
semilogy(4:9,k_sort1(4:end),'-or',...
    4:8,k2(4:end),'-.ob','linewidth',2); %semilogy
set(gca,'fontsize',15,'LineWidth',2);
xlabel('Order of Stiffness','fontsize',18,'Interpreter','tex');
ylabel('Eigenvalue (N/m)','fontsize',18);
legend('PD-CTS','CTS','Location','southeast');
grid on;
fig=gcf;
% fig.Position(3:4)=[800,350];   %change fig size
fig.Position(1:4)=[100,100,650,550];   %change fig size
end

S*diag(l)*S'-diag(S*l)
S*inv(diag(l))*S'-inv(diag(S*l))
S'*diag(l0_c)*S-diag(S'*l0_c)
%% mass matrix 

M=tenseg_mass_matrix_RDT(mass,C,N,l); % generate mass matrix
%% mode analysis
[V_mode,D1] = eig(E_qa'*K_T*E_qa,E_qa'*M*E_qa);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
[omega_sort,e_nb]=sort(omega);
omega_sort(1:3)=0;
V_mode_sort=V_mode(:,e_nb);
num_plt=1:9;
plot_mode_RDT_sa(V_mode_sort(1:8,:),V_mode_sort(9,:),omega_sort,N,R,E_na,C,1:2,S,l,'tangent stiffness matrix',...
    'Order of Vibration Mode','Frequency (Hz)','Hz',num_plt,0.2,saveimg,2);
fig=gcf;
fig.Position(1:4)=[100,100,[650,350]];   %change fig size
% plot sliding distance
figure
bar(4:9,V_mode_sort(9,4:9))
set(gca,'fontsize',15,'LineWidth',2);
xlabel('Order of Stiffness','fontsize',18,'Interpreter','tex');
ylabel('Sliding distance (m)','fontsize',18);
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size
%% damping matrix
ksi=0.02;    %damping coefficient of steel
d_c=2/sqrt(3)*sqrt(rho_s)*A.*E.^0.5;                  % cricital damping 
D=[A_2;C']*diag(ksi.*d_c)*[A_2;C']';
%% statics analysis
substep=30;

ind_w=[];w=[];
ind_dqb=[[1;5]*3-2;[1;5]*3-1;18]; dqb0=[zeros(4,1);0.4];
ind_dl0=[]'; dl0=[]';
% ind_dl0_c=[1,2,3]'; dl0_c=[-40,-30,10]';
[w_t,dqb_t,l0_t,E_qa_new,E_qb_new]=tenseg_load_static_RDT(substep,ind_w,w,ind_dqb,dqb0,ind_dl0,dl0,l0,E_qa,E_qb,gravity,[0;9.8;0],C,mass);
[nq,nqa]=size(E_qa_new);
[~,nqb]=size(E_qb_new);
% modify external force(optional)
w_t(:,1:substep/2)=w_t(:,2:2:end); w_t(:,substep/2+1:end)=w_t(:,end)*ones(1,substep/2);
dqb_t(:,1:substep/2)=dqb_t(:,2:2:end); dqb_t(:,substep/2+1:end)=dqb_t(:,end)*ones(1,substep/2);
l0_t(:,1:substep/2)=l0_t(:,2:2:end); l0_t(:,substep/2+1:end)=l0_t(:,end)*ones(1,substep/2);

% input data
data.q=q; data.C=C; data.ne=ne; data.nn=nn; data.E_qa=E_qa_new; data.E_qb=E_qb_new;%data.S=S;
data.E=E; data.A=A; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dqb_t=dqb_t;% forced movement of pinned nodes
data.l0_t=l0_t;% forced movement of pinned nodes
data.substep=substep;    % substep


data_out=static_solver_RDT(data);
t_t=data_out.t_t;          %member force in every step
q_t=data_out.q_t;          %nodal coordinate in every step
l_t=data_out.l_t;          %member length in every step
K_t_t= data_out.Kt_t; %tangent stiffness of whole struct.

n_t=q_t(1:3*nn,:);          % nodal coordinate n
sld_t=q_t(3*nn+1:end,:);    % sliding distance
%% plot member force 
tenseg_plot_result2(1:substep,t_t,{'1', '2', '3','4','5','6'},{'Substep','Force (N)'} ...
    ,'plot_member_force.png',saveimg,{'-or','-+g','-xb','-*c','-^m','-vk'});
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size
%% Plot nodal coordinate curve X Y slid
% tenseg_plot_result(1:substep,[n_t([2*3-2,2*3-1],:);E_sa'*sld_t],{'2X','2Y','slide'},{'Substep','Coordinate /m)'},'plot_coordinate.png',saveimg);
%% Plot nodal coordinate curve X Y slid
tenseg_plot_result2(1:substep,[n_t([2*3-2,2*3-1],:);E_sa'*sld_t],{'$X_2$','$Y_2$','$s_3$'}, ...
    {'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg,{'-or',':+g','-xb'});
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size

%% plot stiffness in small deformation in XYZ direction
F_dir=zeros(nq,3);
F_dir([3*2-2,3*2-1],1:2)=eye(2);   % force with direction X Y 
F_dir(3*nn+3,3)=1;                      % force in sliding s_a

compliance_dir=zeros(3,substep);    % compliance with direction X Y s_a
stiff_dir=zeros(3,substep);         % stiff with direction X Y s_a
E_qa_temp=[E_qa_new,zeros(nq,1)];
E_qa_temp(3*nn+3,end)=1;              %add sliding freedom
for i=1:substep
    K_taa=E_qa_temp'*K_t_t{i}*E_qa_temp;

  F_dir=F_dir/diag(sqrt(diag(F_dir'*F_dir))); %normalized

    disp_dir=K_taa\(E_qa_temp'*F_dir);
compliance_dir(:,i)=diag(disp_dir'*K_taa'*disp_dir);     
stiff_dir(:,i)=1./compliance_dir(:,i);
end 

%plot stiffness
figure
semilogy(1:substep,stiff_dir(1,:),'-r',...
    1:substep,stiff_dir(2,:),'-.g',...
    1:substep,stiff_dir(3,:),'--b','linewidth',2); %semilogy
set(gca,'fontsize',18,'linewidth',2);
xlabel('Substep','fontsize',18,'Interpreter','tex');
ylabel('Stiffness (N/m)','fontsize',18);
lgd =legend('$X_2$','$Y_2$','$s_3$','location','best','fontsize',15,'Interpreter','latex');
grid on;
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size

%% Plot configurations
 j=linspace(0.01,1,4);
%  fig2=figure;% plot in one figure
for i=1:numel(j)
    num=ceil(j(i)*size(n_t,2));
%  tenseg_plot( reshape(n_t(:,num),3,[]),C_b,C_s,[],[],[]);
% tenseg_plot_CTS(reshape(n_t(:,num),3,[]),C,index_b,S);
% plot seperate figures
tenseg_plot_RDT(reshape(n_t(:,num),3,[]),C,R,index_b,eye(ne),[],[],[], [], [],t_t(:,num),[],[min(t_t),max(t_t)]);
% plot in one figure
% tenseg_plot_RDT(reshape(n_t(:,num),3,[]),C,R,index_b,eye(ne),fig2,[],[], [], [],t_t(:,num),[],[min(t_t),max(t_t)]);
 axis off;
end

%% make video of the statics
name=['Tbar_static_RDT'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,material{2})
% tenseg_video_CTS(n_t,C,[1,2],eye(ne),[],[],[],[],[],[],t_t,[],min(substep,50),5,name,savevideo)
tenseg_video_RDT(n_t,C,R,index_b,eye(ne),[],[],[],[],[],[],t_t,[],min(substep,50),5,name,savevideo)

%% dynamics analysis
% time step
dt=1e-4;
tf=0.05;
out_dt=1e-4;
tspan=0:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and 
ind_w=[];w=[];
ind_dl0=[]; dl0=[];
ind_dqb=[[1;5]*3-2;[1;5]*3-1;18]; dqb0=[zeros(4,1);0.4];
[w_t,dqb_t,dqb_d_t,dqb_dd_t,l0_t,E_qa_new,E_qb_new]=tenseg_load_dyna_RDT(tspan,ind_w,w,ind_dqb,dqb0,ind_dl0,dl0,l0,E_qa,E_qb,gravity,[0;9.8;0],C,mass);
% modify external force(optional)
step=numel(tspan)-1;
w_t(:,1:step/2)=w_t(:,2:2:end); w_t(:,step/2+1:end)=w_t(:,end)*ones(1,step/2+1);
dqb_t(:,1:step/2)=dqb_t(:,2:2:end); dqb_t(:,step/2+1:end)=dqb_t(:,end)*ones(1,step/2+1);
l0_t(:,1:step/2)=l0_t(:,2:2:end); l0_t(:,step/2+1:end)=l0_t(:,end)*ones(1,step/2+1);


% % boundary node motion info
% [~,dqb_t,dqb_d_t,dqb_dd_t,dz_a_t]=tenseg_ex_force_RDT(tspan,E_qa_new,E_qb_new,'vib_force',gravity,[0;0;0],C,mass,[1,2],amplitude,period);

% give initial speed of free coordinates
q0a_d=zeros(size(E_qa_new,2),1);                    %initial speed in X direction

   %% input data
% input data
data.q=q; data.C=C; data.ne=ne; data.nq=nq;data.nn=nn; data.E_qa=E_qa_new; data.E_qb=E_qb_new;%data.S=S;
data.E=E; data.A=A; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dqb_t=dqb_t; data.dqb_d_t=dqb_d_t;  data.dqb_dd_t=dqb_dd_t; % forced movement of pinned nodes
data.l0_t=l0_t;% forced movement of pinned nodes
data.q0a_d=q0a_d;        %initial speed of free coordinates
% data.M=M;data.D=D;
data.rho=rho_s; data.ksi=ksi;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;


%% solve ODE
% solve dynamic equation
data_out2=dynamic_solver_RDT(data);        %solve ODE of dynamic equation
% time history of structure
t_t2=data_out2.t_t;          %member force in every step
q_t2=data_out2.q_t;          %nodal coordinate in every step
l_t2=data_out2.l_t;          %member length in every step

n_t2=q_t2(1:3*nn,:);          % nodal coordinate n
sld_t2=q_t2(3*nn+1:end,:);    % sliding distance
if savedata==1
    save (['Tbr',num2str(tf),'.mat']);
end
%% plot member force 
tenseg_plot_result2(out_tspan,t_t2,{'1', '2', '3','4','5','6'},{'Time / s','Force / N'} ...
    ,'plot_member_force.png',saveimg,{'--r',':g','-b','-c','-.m','-.y'});
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size
%% Plot nodal coordinate curve X Y slid
f1=figure;
tenseg_plot_result2(out_tspan,[n_t2([2*3-2],:)],{'2X'}, ...
    {'Time / s','Coordinate / m'},'plot_coordinate.png',saveimg,{'-b','-+g'},f1);
% plot static result
tenseg_plot_result2([1:substep]*tf/substep,[n_t([2*3-2],:)],{'$X_2$'}, ...
    {'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg,{'--r','--g'},f1);
legend('$X_2$ dynamics','$X_2$ statics');
fig=gcf;
fig.Position(3:4)=[800,250];   %change fig size

%% Plot configurations
 j=linspace(0.01,1,5);
for i=1:numel(j)
    num=ceil(j(i)*size(n_t2,2));
%  tenseg_plot( reshape(n_t(:,num),3,[]),C_b,C_s,[],[],[]);
tenseg_plot_CTS(reshape(n_t2(:,num),3,[]),C,index_b,S);
%  axis off;
end
%% make video of the dynamic
name2=['Tbar_dynamic_RDT',num2str(tf)];
% tenseg_video(n_t,C_b,C_s,[],min(numel(out_tspan),50),name,savevideo,material{2})
% tenseg_video_CTS(n_t2,C,index_b,eye(ne),[],[],[],[],[],[],t_t2,[],min(numel(out_tspan),50),5,name,savevideo)
tenseg_video_RDT(n_t2,C,R,index_b,eye(ne),[],[],[],[],[],[],t_t2,[],min(substep,50),5,name2,savevideo);
