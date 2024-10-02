%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%A Clustered Cable Net full scale (deployable)%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% only plot the concept

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
substep=1;                                     %荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 
%% %% N C of the structure
% Manually specify node positions of double layer prism.
Ry=50;          %radius
Rx=1.2*Ry;
p=20;          %complexity for cable dome
% rate1=19.319/50;

% rate=linspace(0.2,0.8,4)
rate=linspace(0.2,0.3,1)
% rate2=0.25;
for j=1:numel(rate)

rate2=rate(j);
% rate1=2*rate2/(1+rate2);
rate1=(1+rate2)/2;
% generate node in one unit
beta1=2*pi/p;beta2=pi/p;
T1=[cos(beta1) -sin(beta1) 0
    sin(beta1) cos(beta1) 0
    0 0 1];
T2=[cos(beta2) -sin(beta2) 0
    sin(beta2) cos(beta2) 0
    0 0 1];
N00=Rx*[[1;0;0],rate1*T2*[1;0;0],rate2*[1;0;0]];      %initial N
N0=[];
for i=1:p    %rotate nodes
 N0=[N0,T1^(i-1)*N00];
end
%% Z coordinate for saddle shape
aaa=15;bbb=12;

N0(2,:)=Ry/Rx*N0(2,:);    %zoom Y

N0(3,:)=-(N0(1,:)/aaa).^2+(N0(2,:)/bbb).^2;

%% base node
N_base=diag([1,1,0])*N0(:,[1:3:3*p-2]);      %base node
N_base(3,:)=min(N0(3,:))*ones(1,p)-6;
N=[N0,N_base];
C_b_in1=[[1:3:3*p-2]',3*p+[1:p]'];
C_b1 = tenseg_ind2C(C_b_in1,N);%%     this is for vertical boundary
%% connectivity matrix 
% C_b_in=[];
C_b_in=[[1:3:3*p-2]',[4:3:3*p-2,1]'];
C_s_in=[[1:3:3*p-2]',[2:3:3*p-1]';[2:3:3*p-1]',[3:3:3*p]';[2:3:3*p-1]',[4:3:3*p-1,1]';[2:3:3*p-1]',[6:3:3*p,3]';[3:3:3*p]',[6:3:3*p,3]'];
% C_s_in=[[1:3:3*p-2]',[2:3:3*p-1]';[2:3:3*p-1]',[3:3:3*p]';[2:3:3*p-1]',[4:3:3*p-1,1]';[2:3:3*p-1]',[6:3:3*p,3]';[3:3:3*p]',[6:3:3*p,3]';[2:3:3*p-1]',[5:3:3*p-1,2]'];

C_b2 = tenseg_ind2C(C_b_in,N);%%   this is for circular boundary
C_s = tenseg_ind2C(C_s_in,N);
% index_b=[1:size(C_b,1)]';
C_b=[C_b1;C_b2];
C=C_s;                             % only consider string in calculation
C1=[C_b;C_s];                      % C1 is used in plot
nb=size(C_b,1);
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
%% Plot the structure to make sure it looks right
fig_handle=figure
% tenseg_plot(N1,C_b,C_s,fig_handle);
% tenseg_plot(N,C_b,C_s);
% title('Cable net');
%% %% Group/Clustered information 
%generate group index
% gr=[];
%gr={[1:p,2*p+1:3*p]',[p+1:2*p,3*p+1:4*p]',[4*p+1:5*p]',[5*p+1:6*p]};
gr={[1:p,2*p+1:3*p]';[p+1:2*p,3*p+1:4*p]';[4*p+1:5*p]'};  % outer diagonal, inner diagonal, inner hoop
Gp=tenseg_str_gp3(gr,C);    %generate group matrix
% S=eye(ne);                  % no clustering matrix
S=Gp';                      % clustering matrix is group matrix
gr1={nb+[1:p,2*p+1:3*p]';nb+[p+1:2*p,3*p+1:4*p]';nb+[4*p+1:5*p]'};  % outer diagonal, inner diagonal, inner hoop
Gp1=tenseg_str_gp(gr1,C1);    %generate group matrix
S_a=Gp1';     % this is used for plot

[nec,ne]=size(S);
tenseg_plot_CTS(N,C,[],S,fig_handle)

%% vertical brace
% N_base=diag([1,1,0])*N(:,[1:3:3*p-2]);      %base node
% N_base(3,:)=min(N(3,:))*ones(1,p)-6;
% N1=[N,N_base];
% C_b_in1=[[1:3:3*p-2]',3*p+[1:p]'];
% C_b1 = tenseg_ind2C(C_b_in1,N1);%%
% C_s1 = tenseg_ind2C(C_s_in,N1);
% % C_p=[C_b1;C_s1];
% % tenseg_plot(N1,C_b1,C_s1,fig_handle);% plot bars and string without group
% tenseg_plot(N1,C_b1,[],fig_handle);     %plot only bars (boundary)
% tenseg_plot(N,C_b,[],fig_handle);     %plot only bars (boundary)
% % axis off
%% vertical brace  method 2
% N_base=diag([1,1,0])*N(:,[1:3:3*p-2]);      %base node
% N_base(3,:)=min(N(3,:))*ones(1,p)-6;
% N1=[N(:,[1:3:3*p-2]),N_base];
% C_b_in1=[[1:p]',p+[1:p]'];
% C_b1 = tenseg_ind2C(C_b_in1,N1);%%

% C_p=[C_b1;C_s1];
% tenseg_plot(N1,C_b1,C_s1,fig_handle);% plot bars and string without group
tenseg_plot(N,C_b,[],fig_handle);     %plot only bars (vertical boundary)
% tenseg_plot(N,C_b,[],fig_handle);     %plot only bars (circle boundary)
axis off
%% plot hyperbolic paraboloid
if 0
xp=1.1*linspace(-Rx,Rx,40);
yp=1.1*linspace(-Ry,Ry,40);
[Xp,Yp]=meshgrid(xp,yp); %
Zp=-(Xp/aaa).^2+(Yp/bbb).^2;
ss=surf(Xp,Yp,Zp,'FaceAlpha',0.4);
ss.EdgeColor = 'none';

%  plot cylinder

[Xc,Yc,Zc] = cylinder;
ss=surf(Rx*Xc,Ry*Yc,Rx*(Zc)+min(N(3,:))-6,'FaceAlpha',0.3);
ss.EdgeColor = 'none';
end

%% view angle
% view([0,0]);
% view([90,0]);
% view([90,90]);
if saveimg==1
saveas(fig_handle,[num2str(j),'.png']);
end
end





%% %% Boundary constraints
pinned_X=([1:3:3*p-2,3*p+1:4*p])'; pinned_Y=([1:3:3*p-2,3*p+1:4*p])'; pinned_Z=([1:3:3*p-2,3*p+1:4*p])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);



%% %% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
A_1ac=A_1a*S';          %equilibrium matrix CTS
A_2ac=A_2a*S';          %equilibrium matrix CTS
l_c=S*l;                % length vector CTS
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_2ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
index_gp=[1]; % number of groups with designed force

fd=[1e4];                       % force in bar is given as -1000

% [q_gp,t_c,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V1(:,end),w0a,index_gp,fd);    %prestress design
[t_c,t]=tenseg_prestress_design2(Gp,l,l_gp,A_2ag,V1(:,end),w0a,index_gp,fd);    %prestress design

% 
% t_c=1e4*ones(nec,1);
% % t_c=1e7*[1;1;0.1];
% t=S'*t_c;




%% rest length design

index_b=find(t_c<0);              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings
[A_b,A_s,A_c,A,r_b,r_s,r_gp,radius,E_c,l0_c,rho,mass_c]=tenseg_minimass(t_c,l_c,eye(size(S,1)),sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
%% cross sectional design
% A_c=(6e-3)^2*ones(nec,1);
% E_c=Es*ones(nec,1);
E=S'*E_c;     %Young's modulus CTS
A=S'*A_c;     % Cross sectional area CTS

% l0=(t+E.*A).\E.*A.*l;
l0=0.9*l;
l0_c=S*l0;
mass=S'*rho.*A.*l0;
% % Plot the structure with radius
% R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.03,.1],r_b);
% R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.03,.1],r_s);
% R3Ddata.Nradius=0.1*ones(nn,1);
% tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);

%% tangent stiffness matrix
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS3(Ia,C,S,t_c,A_2a,E_c,A_c,l0,l);
% plot the mode shape of tangent stiffness matrix
num_plt=1:4;

% plot_mode2(K_mode,k,N,Ia,C_b,C_s,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.9,saveimg,3);

plot_mode_CTS2(K_mode,k,N,Ia,C1,1:nb,S_a,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.9,saveimg,3);

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
% sort the mode
[w_2_sort,I]=sort(w_2);
V_mode_sort=V_mode(:,I);

omega=real(sqrt(w_2_sort))/2/pi;                   % frequency in Hz

% plot_mode2(V_mode_sort,omega,N,Ia,C_b,C_s,l,'natrual vibration',...
%     'Order of Vibration Mode','Frequency (Hz)','Hz',num_plt,0.8,saveimg,3);

plot_mode_CTS2(V_mode_sort,omega,N,Ia,C1,1:nb,S_a,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)','Hz',num_plt,0.9,saveimg,3);

%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[];
ind_dnb=[]; dnb0=[];
ind_dl0_c=[]; dl0_c=[];
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);


%% Step1: equilibrium calculation
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
data_out1=static_solver_CTS(data);
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
N_out=data_out1.N_out;
tenseg_plot( N_out{:},C_b,C_s,[],[],[])
fig_handle=figure
tenseg_plot_CTS(N_out{:},C1,1:nb,S_a,fig_handle);
% tenseg_plot_CTS(N_out{:},[C_b],[],[],fig_handle);
% tenseg_plot(N,C_b1,[],fig_handle);     %plot only bars (vertical boundary)
% tenseg_plot(N_out{:},C_b,[],fig_handle);     %plot only bars (circle boundary)
tenseg_plot_CTS_dash(N,C,[],S,fig_handle);
% view([0,0]);
% view([90,0]);
% view([0,90]);

%% %% self-stress design

N2=N_out{end};         % final configuration 
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N2,C,Gp,Ia);
A_1ac=A_1a*S';          % equilibrium matrix CTS
A_2ac=A_2a*S';          % equilibrium matrix CTS
l_c=S*l;                % length vector CTS


%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_2ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
index_gp=[1]; % number of groups with designed force

fd=[1e4];                       % force in bar is given as -1000

% [q_gp,t_c,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V1(:,end),w0a,index_gp,fd);    %prestress design
[t_c,t]=tenseg_prestress_design2(Gp,l,l_gp,A_2ag,V2,w0a,index_gp,fd);    %prestress design

%% rest length design
l0=(t+E.*A).\E.*A.*l;
l0_c=S*l0;
% member mass
mass=S'*rho.*A.*l0;


%%   stiffness 
%% tangent stiffness matrix analysis
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS3(Ia,C,S,t_c,A_2a,E_c,A_c,l0,l);   %CTS
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS3(Ia,C,eye(ne),t,A_2a,E,A,l0,l);   % TTS
% plot the mode shape of tangent stiffness matrix
num_plt=1:4;
% plot_mode2(K_mode,k,N2,Ia,C_b,C_s,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.8,saveimg,3);
plot_mode_CTS2(K_mode,k,N2,Ia,C1,1:nb,S_a,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.9,saveimg,3);

figure
plot(1:10,k(1:10),'k-o','linewidth',1.5);
set(gca,'fontsize',18);
xlabel('Order of Eigenvalue','fontsize',18,'Interpreter','latex');
ylabel('Eigenvalue of Stiffness (N/m)','fontsize',18,'Interpreter','latex');
%% mode analysis
mass=S'*rho.*A.*l0;
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
[V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of 
[w_2_sort,I]=sort(w_2);
V_mode_sort=V_mode(:,I);

omega=real(sqrt(w_2_sort))/2/pi;                   % frequency in Hz
% plot_mode2(V_mode_sort,omega,N2,Ia,C_b,C_s,l,'natrual vibration',...
%     'Order of Vibration Mode','Frequency (Hz)','Hz',num_plt,0.8,saveimg,3);
plot_mode_CTS2(V_mode_sort,omega,N2,Ia,C1,1:nb,S_a,l,'natrual vibration',...
    'Order of Vibration Mode','Frequency (Hz)','Hz',num_plt,0.9,saveimg,3);

%% Step 2: External force
substep=8;
load=Rx*Ry*pi*100/(2*p);
ind_w=3*[[2:3:59]';[3:3:60]'];w=-load*ones(2*p,1);
ind_dnb=[]; dnb0=[];
%ind_dl0_c=[1,2,3,4]'; dl0_c=[-400,-300,200,100]';
ind_dl0_c=[]'; dl0_c=[]';
% ind_dl0_c=[1,2,3]'; dl0_c=[-40,-30,10]';
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_ct;% forced movement of pinned nodes
data.N=N_out{end};
data.substep=substep;    % substep

data_out3=static_solver_CTS(data);
t_t3=data_out3.t_out;          %member force in every step
n_t3=data_out3.n_out;          %nodal coordinate in every step
N_out3=data_out3.N_out;
t_c_t3=pinv(S')*t_t3;
%% plot member force 
tenseg_plot_result(1:substep,t_c_t3,{'ODC', 'IDC', 'HC'},{'Substep','Force / N'},'plot_member_force.png',saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t3([2*3,3*3],:),{'2Z','3Z'},{'Substep','Coordinate /m)'},'plot_coordinate.png',saveimg);
%% make video of the dynamic
name=['cable_net_CTS_load'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t3,C_b,C_s,[],min(substep,50),name,savevideo,material{2})



%% Step 3: change rest length of strings
substep=10;
ind_w=[];w=[];
ind_dnb=[]; dnb0=[];
%ind_dl0_c=[1,2,3,4]'; dl0_c=[-400,-300,200,100]';
ind_dl0_c=[1,2]'; dl0_c=[-700,-700]';
% ind_dl0_c=[1,2,3]'; dl0_c=[-40,-30,10]';
[w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_ct;% forced movement of pinned nodes
data.N=N_out{end};
data.substep=substep;    % substep

data_out2=static_solver_CTS(data);
t_t2=data_out2.t_out;          %member force in every step
n_t2=data_out2.n_out;          %nodal coordinate in every step
N_out2=data_out2.N_out;
t_c_t=pinv(S')*t_t2;
%% plot member force 
tenseg_plot_result(1:substep,t_c_t,{'ODC', 'IDC', 'HC'},{'Substep','Force / N'},'plot_member_force.png',saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t2([3*3-2],:),{'3X'},{'Substep','Coordinate /m)'},'plot_coordinate.png',saveimg);

%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
% tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[])
 j=linspace(0.01,1,3);
for i=1:3
    num=ceil(j(i)*size(n_t2,2));
%  tenseg_plot( reshape(n_t(:,num),3,[]),C_b,C_s,[],[],[]);
tenseg_plot_CTS(reshape(n_t2(:,num),3,[]),C1,1:nb,S_a);
 axis off;
end


%% %%%%%%%%%%%%%%%%% Redesign Prestress%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:substep
    N_new=reshape(n_t2(:,i),3,[]);
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N_new,C,Gp,Ia);
A_1ac=A_1a*S';          %equilibrium matrix CTS
A_2ac=A_2a*S';          %equilibrium matrix CTS
l_c=S*l;                % length vector CTS
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_2ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
index_gp=[1]; % number of groups with designed force

fd=[10000];                       % force in bar is given as -1000

[t_c_stp(:,i),t_stp(:,i)]=tenseg_prestress_design2(Gp,l,l_gp,A_2ag,V2,w0a,index_gp,fd);    %prestress design
%% rest length design
l0_stp(:,i)=(t_stp(:,i)+E.*A).\E.*A.*l;
l0_c_stp(:,i)=S*l0_stp(:,i);

end
%% plot member force 
% tenseg_plot_result(1:substep,t_t([1,2],:),{'外斜索','内环索'},{'加载步','力/ N'},'plot_member_force.png',saveimg);
tenseg_plot_result(1:substep,t_c_stp,{'ODC', 'IDC', 'HC'},{'Substep','Force / N'},'plot_member_force.png',saveimg);

%% plot rest length 
% tenseg_plot_result(1:substep,l0_c([1,2],:),{'外斜索','内环索'},{'加载步','原长/ N'},'plot_member_force.png',saveimg);
tenseg_plot_result(1:substep,l0_c_stp,{'ODC', 'IDC', 'HC'},{'Substep','Rest length /m'},'plot_member_force.png',saveimg);


%% save output data
if savedata==1
    save (['cable_net_CTS_',material{1},'.mat']);
end
%% make video of the dynamic
name=['cable_net_CTS'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t2,C_b,C_s,[],min(substep,50),name,savevideo,material{2})

%output data to tecplot
tenseg_tecplot(C,n_t,t_t,interp1([min(radius),max(radius)],[0.2,0.8],radius));