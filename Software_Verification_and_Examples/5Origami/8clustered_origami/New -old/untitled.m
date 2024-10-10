
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
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 

%dynamic analysis set
dt=1e-4;               % time step in dynamic simulation
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=4;                   % final time of dynamic simulation
out_dt=1e-4;            % output data interval(approximately, not exatly)

amplitude=0;            % amplitude of external force of ground motion 
period=0.5;             %period of seismic

%% N C of the structure
% Manually specify node positions of double layer prism.
% N=[0 0 0;1 1 0;2 0 0;1 -1 0]';
alpha=0;beta=0;gamma=0;         %对应x,y,z轴的角度
R=[cos(alpha)*cos(beta),cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma),cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma);
   sin(alpha)*cos(beta),sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma),sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma);
   -sin(beta),cos(beta)*sin(gamma),cos(beta)*cos(gamma)];   %转换矩阵

P_org=[1;1;1];  %局部坐标到整体左坐标的距离
width=0.1;
% N=width*[2,1,3,0,2,4;0,sqrt(3),sqrt(3),0,0,0;0,0,0,0,0.05,0];       %nodal coordinate
r=10; h=30; p=3; f=2;       % radius; height; number of edge;三角形边长划分数量

beta=180*(0.5-1/p); 	% rotation angle
for i=1:p               % nodal coordinate matrix N
    N_initial(:,i)=3*r*[cos(2*pi*(i-1)/p),sin(2*pi*(i-1)/p),0];
end

for i=p+1:2*p
    N_initial(:,i)=[3*r*cos(2*pi*(i-1)/p+beta*pi/180),3*r*sin(2*pi*(i-1)/p+beta*pi/180),h];
end

% for i=2*p+1:3*p
%     N(:,i)=[r*cos(2*pi*(i-1)/p+2*beta*pi/180),r*sin(2*pi*(i-1)/p+2*beta*pi/180),2*h];
% end

a=4;
b=5;    %需要划分三角形的三个节点位置
c=6;

N_f=N_initial([1 2 3],[a b c]); %需要划分三角形的三个节点坐标
N_new=membrane_f2(N_f,f);  %划分后行的节点坐标
N=[N_initial,N_new];

%自动选取对应节点连接panel lines      
C_l_in=[];
for i=1:(f-1)
    C_l_in_j=zeros(3*i,2); 
    for j=1:i 
        C_l_in_j(1+3*(j-1):(3*j),:)=[3+j+(1+(i-1))*(i-1)/2,3+j+(1+(i-1))*(i-1)/2+i;3+j+(1+(i-1))*(i-1)/2+i,3+j+(1+(i-1))*(i-1)/2+i+1;3+j+(1+(i-1))*(i-1)/2,3+j+(1+(i-1))*(i-1)/2+(i+1)];
    end 
    C_l_in = [C_l_in;C_l_in_j];    
end

C_l_in(C_l_in==4)=a;    %将4节点与划分节点的第一个位置调换

C_l_in(C_l_in==5)=-1;    %-1是一个不会出现的数，只是用来替换，无实际意义
C_l_in(C_l_in==(4+(1+(f-1))*(f-1)/2))=b;      %将底部第一个节点与划分节点的第二个位置调换
C_l_in(C_l_in==-1)=(4+(1+(f-1))*(f-1)/2);        
C_l_in(C_l_in==6)=-2;    %-2是一个不会出现的数，只是用来替换，无实际意义
C_l_in(C_l_in==(3+(1+f)*f/2))=c;              %将底部第二个节点与划分节点的第三个位置调换
C_l_in(C_l_in==-2)=(3+(1+f)*f/2);



C_b_in = [1,5;2,6;3,4];  %bar
C_s_in = [1,4;2,5;3,6];  %string
C_t_in = [C_b_in;C_s_in];  % truss

% Convert the above matrices into full connectivity matrices.
C_l = tenseg_ind2C(C_l_in,N);   %connectivity matrices of panel lines\
C_b = tenseg_ind2C(C_b_in,N);   %connectivity matrices of bar
C_s = tenseg_ind2C(C_s_in,N);   %connectivity matrices of string
C_t = tenseg_ind2C(C_t_in,N);   %connectivity matrices of truss

C=[C_l;C_b;C_s];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
[npl,~]=size(C_l);      % npl:No.of panel lines
[nt,~]=size(C_t);       % nt:No.of truss
% Plot the structure to make sure it looks right
tenseg_plot_membrane(N,C_b,C_s,C_l);

%% connectivity of triangle element Ca
% Ca can be written in a function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Ca=[];
for i=1:(f-1)
    Ca_i=zeros(3,1+2*(i-1));
    Ca_i(1,:)=floor(linspace(4+(1+(i-1))*(i-1)/2,4+(1+(i-1))*(i-1)/2+(i-1),1+2*(i-1))); %先确定两个端点然后中间的数向下取整
    Ca_i(2,:)=ceil(linspace(4+(1+i)*i/2,4+(1+i)*i/2+(i-1),1+2*(i-1)));                  %先确定两个端点然后中间的数向上取整
    Ca_i(3,1:2:end)=Ca_i(2,1:2:end)+1;      %选取第三行的奇数列，用第二行的奇数列+1赋值
    Ca_i(3,2:2:end)=Ca_i(2,2:2:end)-i;      %选取第三行的偶数列，用第二行的偶数列+i赋值
    Ca=[Ca,Ca_i];
end

Ca(Ca==4)=a;    %将4节点与划分节点的第一个位置调换

Ca(Ca==5)=-1;    %-1是一个不会出现的数，只是用来替换
Ca(Ca==(4+(1+(f-1))*(f-1)/2))=b;      %将底部第一个节点与划分节点的第二个位置调换
Ca(Ca==-1)=(4+(1+(f-1))*(f-1)/2);        
Ca(Ca==6)=-2;    %-2是一个不会出现的数，只是用来替换
Ca(Ca==(3+(1+f)*f/2))=c;              %将底部第二个节点与划分节点的第三个位置调换
Ca(Ca==-2)=(3+(1+f)*f/2);

% Ca=generate_Ca(C_in,N);
% Ca=zeros(3,1)
[~,np]=size(Ca);        % np:No.of plate

% plot the origami configuration
tenseg_plot_ori_membrane(N,C_b,C_s,C_l,[],[],[],[],[],[],[],Ca);
axis off;
%% Boundary constraints
 pinned_X=[1,2,3]'; pinned_Y=[1,2,3]'; pinned_Z=[1,2,3]';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%% C_pl_bar
 C_pl=zeros(np,3);  %Cpl
 C_pl=[];
for i=1:(f-1)
    C_pl_i=zeros(1+2*(i-1),3);
    C_pl_i(1,:)=[1+3*(1+(i-1))*(i-1)/2,2+3*(1+(i-1))*(i-1)/2,3+3*(1+(i-1))*(i-1)/2]; %先确定矩阵的第一行
    for j=1:(i-1)
    C_pl_i(2*j,1)=C_pl_i(2*j-1,1)-(2+3*(i-2));  %偶数行的第一项的等于奇数行减去(2+3*(i-2))
    C_pl_i(2*j,2:3)=C_pl_i(2*j-1,2:3)+1;        %偶数行的第二三项的等于奇数行加1
    C_pl_i(2*j+1,1)=C_pl_i(2*j,1)+(2+3*(i-1));  %奇数行的第一项的等于偶数行加上(2+3*(i-1))
    C_pl_i(2*j+1,2:3)=C_pl_i(2*j,2:3)+2;        %奇数行的第二三项的等于偶数行加2
    end
    C_pl=[C_pl;C_pl_i];
end
 [C_pl_bar,C_pl_bar_i]=tenseg_ind2C_bar(C_pl,C_l,Ca);

%%  不考虑节点顺序的连接关系矩阵
C_pn=Ca';
[ C_pn_bar,n_pn_i, n_pn_i_local,C_pn_i] = tenseg_ind2C_membrane( C_pn,N,R,P_org );
% [C_pn_bar,n_pn_i,C_pn_i]=tenseg_ind2C_membrane(Ca,N);  

%% 板的质量
rho_p=0.1;                  % density of panel
t_p=6e-3*ones(np,1);        % thickness of panel
A_p=tenseg_A_p(Ca,N);       % area of panel
mass_p=rho_p.*A_p.*t_p;     % mass of panel
M_p=1/12*kron((C_pn_bar'*diag(mass_p)*C_pn_bar+diag(diag(C_pn_bar'*diag(mass_p)*C_pn_bar))),eye(3)); % mass matrix of panel

%% Group/Clustered information 
%generate group index
% gr=[];
gr=[];     % number of elements in one group
Gp=tenseg_str_gp(gr,C);    %generate group matrix
% S=eye(ne);                  % no clustering matrix
S=Gp';                      % clustering matrix as group matrix

gr_tc=[];     % number of elements in one group(for truss)
Gp_tc=tenseg_str_gp(gr_tc,C_t);    %generate group matrix(for truss)
% S=eye(ne);                  % no clustering matrix(for truss)
S_tc=Gp_tc';                      % clustering matrix as group matrix(for truss)
% tenseg_plot_CTS(N,C,[1,2],S);
%% self-stress design
%Calculate equilibrium matrix and member length
[A_2t,A_2ta,A_2tc,A_2tac,A_2l,A_2la,l_t,l_tc,l_l,l,l_c]=tenseg_equilibrium_matrix_truss_lines(N,C,C_t,C_l,S,S_tc,Ia);
for i=1:np
A_2p_i{i}=kron(R',eye(3))*kron(cell2mat(C_pn_i(i)),eye(3))*A_2l*cell2mat(C_pl_bar_i(i))';
B_lp_i{i}=A_2p_i{i}';
B_epsilon_i{i}=[cell2mat(B_lp_i(i))*kron(eye(3),[1,0,0;0,0,0;0,0,0])*cell2mat(n_pn_i_local(i)),cell2mat(B_lp_i(i))*kron(eye(3),[0,0,0;0,1,0;0,0,0])*cell2mat(n_pn_i_local(i)),cell2mat(B_lp_i(i))*kron(eye(3),[0,0,0;1,0,0;0,0,0])*cell2mat(n_pn_i_local(i))];
end
B_epsilon=blkdiag(B_epsilon_i{:}); 
E_p=2.06e5*ones(np,1);     % Young's modulus of panel
mu=0.3;     %泊松比            
% D=zeros(3*np,3*np)
D=[diag(E_p/(1-mu^2)),diag(mu*E_p/(1-mu^2)),zeros(np,np);diag(mu*E_p/(1-mu^2)),diag(E_p/(1-mu^2)),zeros(np,np);zeros(np,np),zeros(np,np),diag(E_p/(2*(1+mu)))]; 

% l0_l=12.9*ones(npl,1);
% Delta_l_l=l_l-l0_l;

Delta_l_l=0.1*ones(npl,1);
l0_l=l_l-Delta_l_l;

pVp_pn=A_2l*(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*Delta_l_l;  %partial_Vp/partial_n
t_l=(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*Delta_l_l;


E_t=2.06e11*ones(nt,1);    % Young's modulus of truss
A_t=1e-6*ones(nt,1);    % area of truss

% l0_t=30*ones(nt,1);
l0_t=[65;65;65;33;33;33];
l0_tc=S_tc*l0_t;

E_tc=pinv(S_tc')*E_t;   % Young's modulus of truss
A_tc=pinv(S_tc')*A_t;    % area of truss


l0=[l0_l;l0_t];
l0_c=S*l0;
% Delta_l_t=l_t-l0_t;
Delta_l_tc=l_tc-l0_tc;

% pVt_pn=A_2t*diag(E_t)*diag(A_t)*diag(l0_t.^-1)*Delta_l_t;   %partial_Vt/partial_n
pVtc_pn=A_2tc*diag(E_tc)*diag(A_tc)*diag(l0_tc.^-1)*Delta_l_tc;   %partial_Vt/partial_n

t_tc=diag(E_tc)*diag(A_tc)*diag(l0_tc.^-1)*Delta_l_tc;



% pV_pn=[A_2t,A_2l]*[t_t;t_l];  %  partial_V/partial_n=partial_Vt/partial_n+partial_Vp/partial_n
pV_pn=[A_2tc,A_2l]*[t_tc;t_l];  %  partial_V/partial_n=partial_Vt/partial_n+partial_Vp/partial_n


Ke_p=A_2l*(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*A_2l';
Kg_p=kron(C_l'*diag(l_l.^-1)*diag(t_l)*C_l,eye(3))-A_2l*diag(l_l.^-1)*diag(t_l)*A_2l';
Kt_p=Ke_p+Kg_p;   %partial_2Vp/partial_n*partial_n'
% p2Vp_pn2=A_2l*(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*A_2l'+kron(C_l'*diag(l_l.^-1)*diag(t_l)*C_l,eye(3))-A_2l*diag(l_l.^-1)*diag(t_l)*A_2l';

Ke_tc=A_2tc*diag(E_tc)*diag(A_tc)*diag(l0_tc.^-1)*A_2tc';
Kg_tc=kron(C_t'*diag(l_t.^-1)*diag(S_tc'*t_tc)*C_t,eye(3))-A_2t*diag(l_t.^-1)*diag(S_tc'*t_tc)*A_2t';
Kt_tc=Ke_tc+Kg_tc;   %partial_2Vt/partial_n*partial_n'
% p2Vt_pn2=A_2t*diag(E_t)*diag(A_t)*diag(l0_t.^-1)*A_2t'+kron(C_t'*diag(l_t.^-1)*diag(t_t)*C_t,eye(3))-A_2t*diag(l_t.^-1)*diag(t_t)*A_2t';
% p2Vt_pn2=A_2tc*pinv(S_tc')*diag(E_t)*diag(A_t)*diag(l0_t.^-1)*pinv(S_tc)*A_2tc'+kron(C_t'*diag(l_t.^-1)*diag(S_tc*t_tc)*C_t,eye(3))-A_2tc*pinv(S_tc')*diag(l_t.^-1)*diag(t_t)*pinv(S_tc)*A_2tc';

K_t=Kt_p+Kt_tc;   %partial_2V/partial_n*partial_n'
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

%% bar和string的质量
A_b=5e-7*ones(size(C_b,1),1);
A_s=9e-10*ones(size(C_s,1),1);
A=[A_b;A_s];
rho=7870*ones(size(C_t,1),1);
mass=rho.*A.*l0_t;

M_b=1/6*kron((abs(C_t)'*diag(mass)*abs(C_t)+diag(diag(abs(C_t)'*diag(mass)*abs(C_t)))),eye(3));

M=M_p+M_b;  %总的质量矩阵
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
% % [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS2(Ia,C,q,A_2ac,E_c,A_c,l0_c);
% % plot the mode shape of tangent stiffness matrix
% 
% Kt_aa=Ia'*K_t*Ia;
% [K_mode,D1] = eig(Kt_aa);         % eigenvalue of tangent stiffness matrix
% k=diag(D1);  
% [k_sort,I]=sort(k);
% K_mode_sort=K_mode(:,I);
% % plot the mode shape of tgent stiffness matrix
% num_plt=1:9;
% plot_mode_ori_membrane(K_mode_sort,k_sort,N,Ia,C_b,C_s,C_l,[],[],l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.1 ,saveimg,3,Ca);

%% input file of ANSYS
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');

%% mass matrix and damping matrix

% damping matrix
% d=1;     %damping coefficient
% D=d*2*max(sqrt(mass.*E.*A./l))*eye(3*nn);    %critical damping

Cd=1e-5*eye(3*nn,3*nn);
%% mode analysis
% [V_mode,D1] = eig(Kt_aa,Ia'*M*Ia);         % calculate vibration mode
% w_2=diag(D1);                                    % eigen value of 
% omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
% % plot_mode_CTS(V_mode,omega,N,Ia,C,[1,2],S,l,'natrual vibration',...
% %     'Order of Vibration Mode','Frequency (Hz)',num_plt,0.1,saveimg);
% plot_mode_ori_membrane(V_mode,omega,N,Ia,C_b,C_s,C_l,[],[],l,'natrual vibration',...
%     'Order of Vibration Mode','Frequency (Hz)',num_plt,0.1 ,saveimg,3,Ca);
% 

%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[];   %external force in Z 
ind_dnb=[]; dnb0=[];
ind_dl0_tc=[4,5,6]; dl0_tc=[-32.9,-32.9,-32.9];
ind_dl0_l=[];dl0_l=[];
% ind_theta_0=[]; dtheta_0=[];        % initial angel change with time
p=1;  %加载次数
substep = substep*p;
% [w_t,dnb_t,l0_ct,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,l0_c,b,gravity,[0;9.8;0],C,mass);
[w_t,dnb_t,l0_tc_t,l0_l_t,Ia,Ib]=tenseg_load_prestress_membrane(substep,substep_1,ind_w,w,ind_dnb,dnb0,ind_dl0_tc,dl0_tc,ind_dl0_l,dl0_l,l0_tc,l0_l,b,gravity,[0;0;9.8],C,M,p);



U=zeros(size(N(:))); 
% 
% A=[1,1,1]';
% B=[2,1,2]';
% z=jacobian(A,B);
% ind_du=[10,11,12,13,14,15,16,17,18]; du=[0.1035,0.0602,-29.1305,-0.1039,0.0595,-29.1305,0.0004,-0.1197,-29.1305];

% ind_du=[10,11,12,13,14,15,16,17,18]; du=[0.1035,0.0602,-2.53e-74,-0.1039,0.0595,-2.53e-74,0.0004,-0.1197,-2.53e-74];
ind_du=[10,11,12,13,14,15,16,17,18]; du=[4.0634,-14.9741,-29.1305,10.9362,11.006,-29.1305,-14.9996,3.968,-29.1305];


    du_it=zeros(size(U));
    du_it(ind_du)=du;
    du_t=du_it*linspace(0,1,substep);
    U_t=du_t+U*linspace(1,1,substep);

    U=U_t(:,20);

    U_1=U(10:18,:);
   



%     X=N_exa+U;

    U_local=kron(R',eye(3))*U_1;    %局部节点位移

    n_pn_i_local_1=cell2mat(n_pn_i_local); %局部节点坐标

    xi=n_pn_i_local_1(1,:);
    yi=n_pn_i_local_1(2,:);
    xj=n_pn_i_local_1(4,:);
    yj=n_pn_i_local_1(5,:);     
    xk=n_pn_i_local_1(7,:);
    yk=n_pn_i_local_1(8,:);

    delta=0.5*(xj*yk-xk*yj-xi*yk+xk*yi+xi*yj-xj*yi);    %三角形面积

    P1=[yj-yk;yk-yi;yi-yj];
    P2=[xk-xj;xi-xk;xj-xi];     
    P3=[0;0;0];   

%     B=U(10:18,:);
    X1=U_local(1:3:end,:);
    X2=U_local(2:3:end,:);
    X3=U_local(3:3:end,:);

%     H=1/2*delta*[P1'*X1,P2'*X1;P1'*X2,P2'*X2];

    H=R*1/(2*delta)*[P1'*X1,P2'*X1,0;P1'*X2,P2'*X2,0;P1'*X3,P2'*X3,0]*R';
    
    E=0.5*(H+H'+H'*H);

%     E1=[E(1,1),E(2,2),2*E(1,2)]';
    E1=[E(1,1),E(2,2),E(3,3),2*E(1,2),2*E(2,3),2*E(3,1)]';

    lambda=E_p*mu/((1+mu)*(1-2*mu));
    G=E_p/(2*(1+mu));
    D1=[lambda+2*G,lambda,lambda,0,0,0;
        lambda,lambda+2*G,lambda,0,0,0;
        lambda,lambda,lambda+2*G,0,0,0;
        0,0,0,G,0,0;
        0,0,0,0,G,0;
        0,0,0,0,0,G];

%     T1=D*E1;
    T1=D1*E1;

% 
%     T=[T1(1,1),T1(3,1);
%        T1(3,1),T1(2,1)];

    T=[T1(1,1),T1(4,1),T1(6,1);
       T1(4,1),T1(2,1),T1(5,1);
       T1(6,1),T1(5,1),T1(3,1)];


    A=H+eye(3);

    sigma=det(A)^-1*A*T*A';


    Delta_l_l=0.2794*ones(3,1);

    sigma_l=D*inv(B_epsilon)*C_pl_bar*Delta_l_l;
%     rij=N(:,4)-N(:,5);  %b_ij
%     rjk=N(:,6)-N(:,5);  %b_jk
%     rik=N(:,6)-N(:,4);  %b_ik
% 
%     delta=rij(1,:)*(rjk(2,:)-rjk(3,:))+rij(2,:)*(rjk(3,:)-rjk(1,:))+rij(3,:)*(rjk(1,:)-rjk(2,:));
%     P1=1/delta*[rjk(2,:)-rjk(3,:);rjk(3,:)-rjk(1,:);rjk(1,:)-rjk(2,:)];
%     P2=1/delta*[rik(2,:)-rik(3,:);rik(3,:)-rik(1,:);rik(1,:)-rik(2,:)];
%     P3=1/delta*[rij(2,:)-rij(3,:);rij(3,:)-rij(1,:);rij(1,:)-rij(2,:)];
%     P=[P1,P2,P3];
%     
%     X=N(:,4:6);
%     X=X(:);
%     X1=X(1:3:end,:);
%     X2=X(2:3:end,:);
%     X3=X(3:3:end,:);
% 
% 
%     H=[(P*X1)';(P*X2)';(P*X3)'];
% 
%     A=H-eye(3);
% 
%     sigma=0.5*det(A).^-1*A*D*(H+H'+H'*H)*A';

%% Step1: statics: equilibrium calculation
% input data
 data.N=N; data.C=C; data.C_l=C_l;data.C_t=C_t;data.C_pl_bar=C_pl_bar;data.U_t=U_t;
data.ne=ne; data.nn=nn; data.Ia=Ia; data.Ib=Ib;data.S=S;data.S_tc=S_tc;data.D=D;
data.E_p=E_p; data.A_p=A_p;data.t_p=t_p;data.E_t=E_t; data.A_t=A_t; data.index_l=index_l; data.index_t=index_t;data.B_epsilon=B_epsilon;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;           % forced movement of pinned nodes
data.l0_tc_t=l0_tc_t;            % forced change of rest length
data.l0_l_t=l0_l_t;            % forced change of rest length
data.substep=substep;    % substep
%data.InitialLoadFactor=0.001;
%data.MaxIcr=1000;
%.LoadType='Force'; % 'Force' or 'Displacement'
%data.StopCriterion=@(U)(norm(U)>0.5);

% nonlinear analysis
% data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
data_out1=static_solver_CTS_membrane_tl(data);
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

% t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
sigma=data_out1.sigma_out;    %sigma_l
% t_tc=data_out1.t_tc_out;        %t_tc
N_out=data_out1.N_out;

%% plot member force 
% tenseg_plot_result(1:substep,t_tc([1,4],:),{'bar','string'},{'Substep','Force (N)'},'plot_member_force.png',saveimg);
% grid on;
%% Plot nodal coordinate curve X Y
% tenseg_plot_result(1:substep,n_t([10,11,12,22,23,24],:),{'x4','y4','z4','x8','y8','z8'},{'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg);
% grid on;
%% plot hinge moment
% tenseg_plot_result(1:substep,sigma_l([1,2,3],:),{'sigma_x','sigma_y','tau'},{'Substep','sigma (pa)'},'plot_hinge_moment.png',saveimg);

%% Plot configuration
% for i=round(linspace(1,substep,3))
% tenseg_plot_CTS(reshape(n_t(:,i),3,[]),C,[1,2],S);
% axis off;
% end
% tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,90])


 j=linspace(1e-5,1,10);
for i=1:10
    num=ceil(j(i)*size(n_t,2));
tenseg_plot_ori_membrane(reshape(n_t(:,num),3,[]),C_b,C_s,C_l,[],[],[],[],[30,30],[] ,[],Ca);
 axis off;
end
%% save output data
if savedata==1
    save (['cable_net_CTS_static','.mat']);
end
%% make video of the dynamic
icrm=size(n_t,2); 
name=['tower_m1'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video_ori_membrane(n_t,C_b,C_s,C_l,[],[],Ca,[],min(icrm,50),name,savevideo,[])

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
ind_dl0_tc=[4,5,6]; dl0_tc=[-32.9,-32.9,-32.9];
ind_dl0_l=[];dl0_l=[];


[w_t,l0_tc_t,l0_l_t]=tenseg_load_prestress_CTS_membrane(tspan,ind_w,w,'ramp',ind_dl0_tc,dl0_tc,l0_tc,ind_dl0_l,dl0_l,l0_l,gravity,[0;0;9.8],C,M);

% w_t=[w_t,w_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];   % second half no change of boundary info
% l0_ct=[l0_ct,l0_ct(:,end)*ones(1,numel(tspan)-numel(tspan1))];
% boundary node motion info
[~,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force_membrane(tspan1,a,b,'vib_force',gravity,[0;0;9.8],C,M,[1,2],amplitude,period);
dnb_t=[dnb_t,dnb_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];
dnb_d_t=[dnb_d_t,dnb_d_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];
dnb_dd_t=[dnb_dd_t,dnb_dd_t(:,end)*ones(1,numel(tspan)-numel(tspan1))];
% give initial speed of free coordinates
n0a_d=zeros(numel(a),1);                    %initial speed in X direction
    %% dynamics calculation

% input data
data.N=N; data.C_t=C_t;data.C_l=C_l; data.ne=ne; data.nn=nn;data.np=np; data.Ia=Ia; data.Ib=Ib;data.S_tc=S_tc;
data.E_tc=E_tc; data.A_tc=A_tc; data.A_p=A_p;data.t_p=t_p; data.B_epsilon=B_epsilon; data.C_pl_bar=C_pl_bar;
data.D=D; %转换矩阵  data.material=material; %constitue info
data.w_t=w_t;           % external force
data.dnb_t=dnb_t; data.dnb_d_t=dnb_d_t;  data.dnb_dd_t=dnb_dd_t; % forced movement of pinned nodes
data.l0_l_t=l0_l_t; data.l0_tc_t=l0_tc_t;         
data.n0a_d=n0a_d;        %initial speed of free coordinates
data.M=M;data.Cd=Cd;%阻尼矩阵
% data.rho=rho_s;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;

%% dynamic analysis
% solve dynamic equation
data_out=dynamic_solver_CTS_membrane(data);        %solve ODE of dynamic equation
% time history of structure
% t_t=data_out.t_t;   %time history of members' force
t_tc_t=data_out.t_tc_t;
sigma_l_t=data_out.sigma_l_t;
n_t=data_out.n_t;   %time history of nodal coordinate 
% l_tc_t=data_out.l_tc_t; 
% l_l_t=data_out.l_l_t;   %time history of members' length 
nd_t=data_out.nd_t;   %time history of nodal coordinate
%% Plot configuration
 j=linspace(1e-5,1,10);
for i=1:10
    num=ceil(j(i)*size(n_t,2));
tenseg_plot_ori_membrane(reshape(n_t(:,num),3,[]),C_b,C_s,C_l,[],[],[],[],[30,30],[] ,[],Ca);
%  axis off;
end
%% make video of the dynamic
icrm=size(n_t,2); 
name=['tower_m_dynamic'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video_ori_membrane(n_t,C_b,C_s,C_l,[],[],Ca,[],min(icrm,50),name,savevideo,[])


%% plot member force 
tenseg_plot_result(downsample(out_tspan,30),[downsample(t_tc_t([1],:),30);downsample(t_tc_t([4],:),30)],{'bar','string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y
tenseg_plot_result(downsample(out_tspan,30), [downsample(n_t([10],:),30);downsample(n_t([11],:),30);downsample(n_t([12],:),30);downsample(n_t([22],:),30);downsample(n_t([23],:),30);downsample(n_t([24],:),30)],{'x4','y4','z4','x8','y8','z8'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% plot hinge moment
tenseg_plot_result(downsample(out_tspan,30), [downsample(sigma_l_t([1],:),30);downsample(sigma_l_t([2],:),30);downsample(sigma_l_t([3],:),30)],{'sigma_x','sigma_y','tau'},{'Time (s)','sigma (pa)'},'plot_hinge_moment.png',saveimg);

%% plot member force 
% tenseg_plot_result(out_tspan,t_t([10],:),{'string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);
% grid on;
 tenseg_plot_result(downsample(out_tspan,30),downsample(t_t([10],:),30),{'string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);
grid on;
%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t([13,14,15],:),{'x5','y5','z5'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;

tenseg_plot_result(downsample(out_tspan,30),[downsample(n_t([13],:),30);downsample(n_t([14],:),30);downsample(n_t([15],:),30)],{'x6','y6','z6','x8','y8','z8'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% Plot final configuration
% tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_ct(index_s,end))
tenseg_plot( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,90])

%% save output data
if savedata==1
    save (['cable_net_CTS_dynamic',num2str(tf),'.mat']);
end