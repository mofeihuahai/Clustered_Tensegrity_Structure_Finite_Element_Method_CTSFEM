%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%CTS_consist_powe%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% 

% after running 'CTS_pulley_s', run this to recalculate member force and
% length
% Input: n_t
% Output  'l0_ct'


%EXAMPLE
clc; clear all; close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
% %写function
% [consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib_multi('material-1','material-2','material-3',''...);
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=1; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% static analysis set
substep=50;                                     %荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
z=0;
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 
pully =1;                      %pully influence 1 or 0
if pully 
  [consti_datap,Ep,sigmap,rho_p]=material_lib_p('Steel_Q345');  
end
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','nylon_rope');
%% input n_t 
 n_t=load('n_t.txt');
% w_t=load('w_t.txt');

 %% 
 l0=zeros(60,substep);
 t=zeros(26,substep);
 l0_c=zeros(26,substep);
 for i=1:substep
%% N C of the structure
N=reshape(n_t( :,i),3,[]);
p=12;
C_b_in=[[1:4:4*p-3]',[2:4:4*p-2]';[3:4:4*p-1]',[4:4:4*p]';];
C_s_in=[[2:4:4*p-2]',[3:4:4*p-1]';[3:4:4*p-1]',[6:4:4*p-2,2]';[4:4:4*p]',[8:4:4*p,4]'];

C_b = tenseg_ind2C(C_b_in,N);%%将上面输入c_b_in转化为关系矩阵
C_s = tenseg_ind2C(C_s_in,N);
C=[C_s;C_b];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node

% Plot the structure to make sure it looks right
% tenseg_plot(N,[],C);%绘图
% title('Cable dome');
%% Boundary constraints
pinned_X=([1:4:4*p-3])'; pinned_Y=([1:4:4*p-3])'; pinned_Z=([1:4:4*p-3])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%% Group/Clustered information 
%generate group index
% gr=[];
gr={[1:2*p]',[2*p+1:3*p]'};
Gp=tenseg_str_gp2(gr,C);    %generate group matrix
% S=eye(ne);                  % no clustering matrix
S=Gp';                      % clustering matrix is group matrix
%% cross sectional
E_s=Es*ones(3*p,1);
E_p=Ep*ones(2*p,1);
E=[E_s;E_p];
E_c=[Es*ones(2,1);Ep*ones(2*p,1)];
A_s=0.25e-6*ones(3*p,1);
A_p=1e-6*ones(2*p,1);
A=[A_s;A_p];
A_c=[0.25e-6*ones(2,1);1e-6*ones(2*p,1)];
rho=[rho_s*ones(3*p,1);rho_p*ones(2*p,1)];
%l0=(t+E.*A).\E.*A.*l;

%% external force
H=N*C';    
l=sqrt(diag(H'*H)); 
if i==1
g=[0;0;9.8];
mass=rho.*A.*l;
G=(gravity)*-0.5*kron(abs(C)'*mass,g);
w_t=zeros(size(G,1),1);
w_t=w_t+G;
w_ta=Ia'*w_t;
% w_ta=w_ta*linspace(0,1,substep);
end
% w_ta=Ia'*w_t;
%% equilibrium matrix and SVD of equilibrium matrix
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

l_c=S*l;                % length vector CTS  %同组构件的总长
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_2ag);%A_lag*V2=0   U2'*A_lag=0  V2即为自应力模态矩阵
%% pre-stress Calculate]
fd=75;                                     %constant force
t1=pinv(A_2ag)*w_ta;
I=eye(length(l_c));
index_h=2;                                %index of string with constant force
e_d=I(:,index_h);
%z=(t_t(2*p+1,1)-t1(2*p+1))/V2(2*p+1);
z(i)=(e_d'*V2)\(fd-e_d'*t1);        %self-stress coefficient
t2=z(i).*V2;
t_t(:,i)=t1+t2;
%external force in equilibrium design
%w0=zeros(numel(N),1); w0a=Ia'*w0;  

%% l0 and deltal0

l0_c(:,i)=(t_t(:,i)+E_c.*A_c).\E_c.*A_c.*l_c; 


 end
 save('l0_c.txt','l0_c','-ascii');
%  save('n_t.txt','n_t','-ascii');
 delta_l=zeros(size(l0_c,1),substep-1);
 for i=1:substep-1
 delta_l(:,i)=l0_c(:,i+1)-l0_c(:,i);
 end
 save('delta_l.txt','delta_l','-ascii');
 
%% plot member force 
tenseg_plot_result(1:substep,t_t([1,2],:),{'外斜索','内环索'},{'加载步','力/ N'},'plot_member_force.png',saveimg);
tenseg_plot_result(1:substep,t_t([1,2],:),{'Diagonal Cable','Hoop Cable'},{'Substep','Tension /N'},'plot_member_force.png',saveimg);

%% plot rest length 
tenseg_plot_result(1:substep,l0_c([1,2],:),{'外斜索','内环索'},{'加载步','原长/ N'},'plot_member_force.png',saveimg);
tenseg_plot_result(1:substep,l0_c([1,2],:),{'Diagonal Cable','Hoop Cable'},{'Substep','Rest length /m'},'plot_member_force.png',saveimg);



  %% plot 3 configurations
  syms h
h=solve(2.205^(1/2)*sin(pi/12)/sin(2*pi/3)*(h+0.16)==1.05*h,h);
h=eval(h);

Ng1=[(1.05^2/2)^(1/2) 0 -(1.05^2/2)^(1/2) 0 0 0];
Ng2=[0 (1.05^2/2)^(1/2) 0 -(1.05^2/2)^(1/2) (1.05^2/2)^(1/2) -(1.05^2/2)^(1/2)];
Ng3=[-h 0.32+h -h 0.32+h -h -h];
Ng=[Ng1;Ng2;Ng3];
C_z_in=[1 5; 1 2;2 3;2 5;3 4;3 5;3 6;4 1;4 6;6 1];
C_z = tenseg_ind2C(C_z_in,Ng);
Bz=Ng*C_z';
for k = 1:size(Bz,2)
        start_nodes(:,k) = Ng(:,C_z(k,:)==-1);
        end_nodes(:,k) = Ng(:,C_z(k,:)==1);
end

 j=linspace(0.01,1,3);
for i=1:3
    num=ceil(j(i)*size(n_t,2));
    
% tenseg_plot(reshape(n_t(:,num),3,[]),C_b,C_s); % plot strings and node
% 
% quiver3(start_nodes(1,:),start_nodes(2,:),start_nodes(3,:),Bz(1,:),Bz(2,:),Bz(3,:),'black.','Autoscale','off','LineWidth',3);

tenseg_plot([Ng,reshape(n_t(:,num),3,[])],blkdiag(C_z,C_b),[zeros(36,6),C_s]); %plot whole structure
axis off;
 view(-35,35);
% view(-15,60);
end




%% make video of the dynamic
name=['CTS_cable_net_animation2'];
n_t2=[Ng(:)*ones(1,size(n_t,2));n_t];  % add boundary node
tenseg_video(n_t2,blkdiag(C_z,C_b),[zeros(36,6),C_s],[],substep,name,savevideo,material{2});
% tenseg_video(flip(n_t2,2),blkdiag(C_z,C_b),[zeros(36,6),C_s],[],10,name,savevideo,material{2});
