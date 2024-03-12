%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%roller-driven 3D arm      %%%%%%%%%%%%%%%%%%%%
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
% Manually specify node positions
Height=1; b=1;h=0.5; level=3;        % shape parameters
N=zeros(3,5*(level+1)+4*(level+1));     %initialize N
% top nodes
N(:,1:5:5*level+1)=[0,0,1]'*Height*[[0:level]];
% side noses
N(:,kron([0:5:5*level],ones(1,4))+kron(ones(1,level+1),[2:5]))=...
    [kron(ones(1,level+1),0.5*[b,-b,-b,b;h,h,-h,-h]);kron(Height*[0:level]-Height/2,ones(1,4))];
% [kron(ones(1,level+1),0.5*[-b,b,b,-b;h,h,-h,-h]);zeros(1,4),kron(H*[1:level]-H/2,ones(1,4))];
% dulplicate side nodes
N(:,5*(level+1)+1:end)=[kron(ones(1,level+1),0.5*[b,-b,-b,b;h,h,-h,-h]);kron(Height*[0:level]-Height/2,ones(1,4))];
 
% Manually specify connectivity indices.
C_b_in=zeros(level*15+4,2);
% bottom bar
C_b_in(1:4*level,:) = [kron([1:5:5*(level-1)+1]',ones(4,1)),kron([5:5:5*level],ones(1,4))'+kron(ones(1,level),[2:5])'];  % Similarly, this is saying bar 1 connects node 1 to node 2,
% top bar
C_b_in(4*level+1:8*level+4,:) = [kron([1:5:5*(level)+1]',ones(4,1)),kron([0:5:5*level],ones(1,4))'+kron(ones(1,level+1),[2:5])'];  % Similarly, this is saying bar 1 connects node 1 to node 2,
% horizontal bar
C_b_in(8*level+5:12*level+4,:) =[kron([5:5:5*level],ones(1,4))'+kron(ones(1,level),[2:5])',kron([5:5:5*level],ones(1,4))'+kron(ones(1,level),[3,4,5,2])'];
% vertical bars
C_b_in(12*level+5:13*level+4,:) =[1:5:5*(level-1)+1;6:5:5*level+1]';
% cross bars
C_b_in(13*level+5:15*level+4,:)=[kron([5:5:5*level]',ones(2,1))+kron(ones(level,1),[2;3]),kron([5:5:5*level]',ones(2,1))+kron(ones(level,1),[4;5])];

C_s_in=zeros(level*12,2); 
% digaonal string
C_s_in(1:8*level,:) =[kron([5:5:5*level],ones(1,4))'+kron(ones(1,level),[2:5])',kron([0:5:5*(level-1)],ones(1,4))'+kron(ones(1,level),[3,4,5,2])';...
                     kron([5:5:5*level],ones(1,4))'+kron(ones(1,level),[2:5])',kron([0:5:5*(level-1)],ones(1,4))'+kron(ones(1,level),[5,2,3,4])'];
% vertical string
C_s_in(8*level+1:12*level,:) =[[5*(level+1)+1:5*(level+1)+4*(level)]',[5*(level+1)+5:5*(level+1)+4*(level+1)]'];

% Convert the above matrices into full connectivity matrices.
 C_b = tenseg_ind2C(C_b_in,N);%%
 C_s = tenseg_ind2C(C_s_in,N);
 C=[C_b;C_s];
 [ne,nn]=size(C);           % ne: No.of element;nn:No.of node
 nb=size(C_b,1);               % nb: No. of bars
 ns=size(C_s,1);             % ns: No. of strings
 % Plot the structure to make sure it looks right

 tenseg_plot(N,C_b,C_s);
 title('arm');
%% radius , side vector for plot
R=zeros(nn,1); % radius vector
% R(5*(level+1)+1:5*(level+1)+4*(level-1))=0.08;
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
 pinned_X=[1:5]; pinned_Y=[1:5]; pinned_Z=(1:5);
 [Ia_min,Ib_min,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,numel(ia));
 E_na=I_unique*Ia_min;
 E_nb=I_unique*Ib_min;
 %% sliding distance vector of pulley and transformation matrix
sld=zeros(nn,1);
% transformation matrix of sliding distance vector
free_sld=[5*(level+1)+1:5*(level+1)+4*(level)]; fixed_sld=setdiff(1:nn,free_sld);
E_s=eye(nn);
E_sa=E_s(:,free_sld);
E_sb=E_s(:,fixed_sld);
%% generalized coordinate q
% generalized coordinate 
q=[N(:);sld];
% transformation matrix of generalized coordinate
E_qa=blkdiag(E_na,E_sa);
E_qb=blkdiag(E_nb,E_sb);
[nq,nqa]=size(E_qa);
[~,nqb]=size(E_qb);
%% Group/Clustered information 
%generate group index
gr={[1:nb]';[nb+1:nb+8*level]';[nb+8*level+1:nb+12*level]'};     % number of elements in one group
Gp=tenseg_str_gp(gr,C);    %generate group matrix
S=Gp';                      % clustering matrix

fig_handle=figure;
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
A=1e-5*[100*ones(size(C_b,1),1);1e-2*ones(8*level,1);1*ones(4*level,1)];
% young's modulus
E=[Eb*ones(size(C_b,1),1);Es*ones(12*level,1)];
% external force
f_ena=zeros(size(E_na,2),1);
f_esa=zeros(size(E_sa,2),1);

% member force design
% index_gp=[1];
% fd=-1e3;
% e_nb=eye(ne);
% e_d=e_nb(:,index_gp); % e_d is the matrix to select group of member with designed force
% z=(e_d'*V2)\(e_d'*(fd-pinv(A_RDT)*[f_ena;f_esa]));   %self-stress coefficient
% t=pinv(A_RDT)*[f_ena;f_esa]+V2*z;
t=zeros(ne,1);

index_b=1:size(C_b,1);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings

% member rest length
l0=l.*[ones(size(C_b,1),1);0.5*ones(8*level,1);0.8*ones(4*level,1)];    % change rest length, change tension
mass=rho_s.*A.*l0;
%% plot with pulley size
fig2=figure
tenseg_plot_RDT(N,C,R,index_b,eye(ne),fig2,[],[], [], [],[],[],[])

%% tangent stiffness matrix
Kn=kron(C'*diag(l.\t)*C,eye(3));

K_T=[Kn+A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2',A_2*diag(E.*A./l0)*C;...
    C'*diag(E.*A./l0)*A_2',C'*diag(E.*A./l0)*C];
K_Tg=blkdiag(Kn-A_2*diag(l.\t)*A_2',zeros(numel(sld)));%geometry stiffness
K_Te=[A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2',A_2*diag(E.*A./l0)*C;...
    C'*diag(E.*A./l0)*A_2',C'*diag(E.*A./l0)*C];

K_T=0.5*(K_T+K_T');
K_Taa=E_qa'*K_T*E_qa;

[K_mode1,D1] = eig(K_Taa);         % eigenvalue of tangent stiffness matrix
k=diag(D1);   
[k_sort1,e_nb]=sort(k);
K_mode_sort1=K_mode1(:,e_nb);
% plot mode
num_plt=1:5;
% plot_mode_CTS2(K_mode_sort(1:8,:),k_sort,N,E_na,C,1:2,S,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.2,saveimg,2);
plot_mode_RDT(K_mode_sort1(1:size(E_na,2),:),k_sort1,N,R,E_na,C,index_b,eye(ne),l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.5,saveimg,3);
%% reduced order tangent stiffness of sliding distance vector
  if 1      % this is reduced order from RDT stiffness matrix
K_T11=E_na'*(Kn+A_2*diag(E.*A./l0)*A_2'-A_2*diag(l.\t)*A_2')*E_na;
K_T12=E_na'*(A_2*diag(E.*A./l0)*C)*E_sa;
K_T21=E_sa'*(C'*diag(E.*A./l0)*A_2')*E_na;
K_T22=E_sa'*(C'*diag(E.*A./l0)*C)*E_sa;
 
Kt_aa_nn=K_T11-K_T12/K_T22*K_T21;       % reduced order stiffness in nodal coordinate
Kt_aa_ss=-K_T21/K_T11*K_T12+K_T22;      % reduced order stiffness in sliding distance

[K_modess,D1] = eig(Kt_aa_ss);         % eigenvalue of tangent stiffness matrix
k1=diag(D1);   
[kss,e_nb]=sort(k1);
K_modess=K_modess(:,e_nb);          % mode vector in s_a
K_modens=-K_T11\K_T12*K_modess;     % mode vector in n_a

num_plt=1:12;
% plot_mode_CTS2(K_mode,k,N,Ia,C,1:2,S,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.2,saveimg,3);
% plot_mode(K_mode2,k2,N,E_na,C_b,C_s,l,'natrual vibration',...
%     'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg);
plot_mode_RDT(K_modens,kss,N,R,E_na,C,index_b,eye(ne),l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.5,saveimg,3);
  end
  %% reduced order tangent stiffness-general method: qa_a of end roller
  if 1                                    % reduced order-general method
I=eye(nqa);
index_qa_a=size(E_na,2)+1:size(E_na,2)+4; % active member index in qa
index_qa_p=setdiff(1:nqa,index_qa_a);       % passive members index in qa
E_qa_a=I(:,index_qa_a);                     % transform matrix for active member
E_qa_p=I(:,index_qa_p);                     % transform matrix for passive member

Kt_aa_11=E_qa_a'*K_Taa*E_qa_a;Kt_aa_12=E_qa_a'*K_Taa*E_qa_p;
Kt_aa_21=E_qa_p'*K_Taa*E_qa_a;Kt_aa_22=E_qa_p'*K_Taa*E_qa_p;
Kt_aa_a=Kt_aa_11-Kt_aa_12/Kt_aa_22*Kt_aa_21;    % reduced form stiffness matrix

[K_mode_a,D1] = eig(Kt_aa_a);         % eigenvalue of tangent stiffness matrix
k_a=diag(D1);   
[k_a_sort,e_nb]=sort(k_a);
K_mode_qaa_sort=K_mode_a(:,e_nb);          % mode vector in qa_a
K_mode_qa=(E_qa_a-E_qa_p*(Kt_aa_22\Kt_aa_21))*K_mode_qaa_sort;     % mode vector in qa
K_mode_na=K_mode_qa(1:size(E_na,2),:);          % mode shape in na
sqrt(sum(K_mode_na.^2))

num_plt=1:4;
% plot_mode_CTS2(K_mode,k,N,Ia,C,1:2,S,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.2,saveimg,3);
% plot_mode(K_mode2,k2,N,E_na,C_b,C_s,l,'natrual vibration',...
%     'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg);
plot_mode_RDT(K_mode_na,k_a_sort,N,R,E_na,C,index_b,S,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.5,saveimg,[45,30]);
% fig=gcf;
% fig.Position(1:4)=[100,100,[600,350]];   %change fig size
  end
    %% reduced order tangent stiffness-general method dqa_a=0
  if 1                                    % reduced order-general method
I=eye(nqa);
index_qa_a=size(E_na,2)+1:size(E_na,2)+4; % active member index in qa
index_qa_p=setdiff(1:nqa,index_qa_a);       % passive members index in qa
E_qa_a=I(:,index_qa_a);                     % transform matrix for active member
E_qa_p=I(:,index_qa_p);                     % transform matrix for passive member

Kt_aa_11=E_qa_a'*K_Taa*E_qa_a;Kt_aa_12=E_qa_a'*K_Taa*E_qa_p;
Kt_aa_21=E_qa_p'*K_Taa*E_qa_a;Kt_aa_22=E_qa_p'*K_Taa*E_qa_p;
Kt_aa_a=Kt_aa_11-Kt_aa_12/Kt_aa_22*Kt_aa_21;    % reduced form stiffness matrix

[K_mode_p,D1] = eig(Kt_aa_22);         % eigenvalue of tangent stiffness matrix
k_p=diag(D1);   
[k_p_sort,e_nb]=sort(k_p);
K_mode_qap_sort=K_mode_p(:,e_nb);          % mode vector in qa_a
K_mode_qa=E_qa_p*K_mode_qap_sort;     % mode vector in qa
K_mode_na=K_mode_qa(1:size(E_na,2),:);          % mode shape in na

num_plt=1:4;
% plot_mode_CTS2(K_mode,k,N,Ia,C,1:2,S,l,'tangent stiffness matrix',...
%     'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.2,saveimg,3);
% plot_mode(K_mode2,k2,N,E_na,C_b,C_s,l,'natrual vibration',...
%     'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg);
plot_mode_RDT(K_mode_na,k_p_sort,N,R,E_na,C,index_b,S,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.5,saveimg,[45,30]);
% fig=gcf;
% fig.Position(1:4)=[100,100,[650,350]];   %change fig size
end
%% tangent stiffness of Clustered Tensegrity CTS
if 0
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
plot_mode_RDT(-K_mode2,k2,N,R,E_na,C,1:2,eye(ne),l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)','N/m',num_plt,0.2,saveimg,2);


% compare stiffness eigen value
figure
semilogy(4:9,k_sort1(4:end),'-or',...
    4:8,k2(4:end),'-.ob','linewidth',2); %semilogy
set(gca,'fontsize',15,'LineWidth',2);
xlabel('Order of Stiffness','fontsize',18,'Interpreter','tex');
ylabel('Eigenvalue (N/m)','fontsize',18);
legend('RDT','CTS');
grid on;
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size
end

%% damping matrix
ksi=0.0;    %damping coefficient of steel
d_c=2/sqrt(3)*sqrt(rho_s)*A.*E.^0.5;                  % cricital damping 
D=[A_2;C']*diag(ksi.*d_c)*[A_2;C']';
%% statics analysis
substep=30;

ind_w=[];w=[];

ind_dqb=[find(sum(E_qa*E_qa_a,2))]; dqb0=[0.20*-K_mode_qaa_sort(:,1)];
ind_dl0=[]'; dl0=[]';
% ind_dl0_c=[1,2,3]'; dl0_c=[-40,-30,10]';
[w_t,dqb_t,l0_t,E_qa_new,E_qb_new]=tenseg_load_static_RDT(substep,ind_w,w,ind_dqb,dqb0,ind_dl0,dl0,l0,E_qa,E_qb,gravity,[0;9.8;0],C,mass);
[nq,nqa]=size(E_qa_new);
[~,nqb]=size(E_qb_new);
% modify external force(optional)
% w_t(:,1:substep/2)=w_t(:,2:2:end); w_t(:,substep/2+1:end)=w_t(:,end)*ones(1,substep/2);
% dqb_t(:,1:substep/2)=dqb_t(:,2:2:end); dqb_t(:,substep/2+1:end)=dqb_t(:,end)*ones(1,substep/2);
% l0_t(:,1:substep/2)=l0_t(:,2:2:end); l0_t(:,substep/2+1:end)=l0_t(:,end)*ones(1,substep/2);

% input data
rdm_q=[0.01*rand(size(E_na,1),1);zeros(size(E_sa,1),1)];  rdm_q([1:4,5*level+[1:4]])=0;% initial form with disturbance
q_disturb=q+rdm_q;
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
%% plot member force -
tenseg_plot_result2(1:substep,t_t(end-3:end,:),{'$f_{s_1}$','$f_{s_2}$','$f_{s_3}$','$f_{s_4}$'},{'Substep','Force (N)'} ...
    ,'plot_member_force.png',saveimg,{'-or','-+g','-vb','-^c','-^m','-vk'});
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size
%% Plot nodal coordinate curve X Y slid
% tenseg_plot_result(1:substep,[n_t([2*3-2,2*3-1],:);E_sa'*sld_t],{'2X','2Y','slide'},{'Substep','Coordinate /m)'},'plot_coordinate.png',saveimg);
%% Plot nodal coordinate curve X Y slid
tenseg_plot_result2(1:substep,[n_t([16*3-2:16*3],:);sld_t(21:24,:)],{'$X_{top}$','$Y_{top}$','$Z_{top}$','$s_1$','$s_2$','$s_3$','$s_4$'}, ...
    {'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg,{'-or','-+g','-xb','--<k','-->k','-.<k','-.>c'});
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size

%% plot stiffness in small deformation in XYZ direction
F_dir=zeros(nq,3);
F_dir([16*3-2:16*3],:)=kron(eye(3),ones(1,1)/norm(ones(1,1)));   % force with direction X Y Z


compliance_dir=zeros(3,substep);    % compliance with direction X Y Z
stiff_dir=zeros(3,substep);         % stiff with direction X Y Z
E_qa_temp=E_qa_new;
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
lgd =legend('$X_{top}$','$Y_{top}$','$Z_{top}$','location','best','fontsize',15,'Interpreter','latex');
grid on;
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size
% plot(stiff_dir(2,:))
%% Plot configurations
 j=linspace(0.01,1,4);
 fig2=figure;% plot in one figure
for i=1:numel(j)
    num=ceil(j(i)*size(n_t,2));
%  tenseg_plot( reshape(n_t(:,num),3,[]),C_b,C_s,[],[],[]);
% tenseg_plot_CTS(reshape(n_t(:,num),3,[]),C,index_b,S);
% plot seperate figures
% tenseg_plot_RDT(reshape(n_t(:,num),3,[]),C,R,index_b,eye(ne),[],[],[45,30], [], [],t_t(:,num),[],[min(t_t),max(t_t)]);
% plot in one figure
tenseg_plot_RDT(reshape(n_t(:,num),3,[]),C,R,index_b,eye(ne),fig2,[],[45,30], [], [],t_t(:,num),[],[min(t_t),max(t_t)]);
 axis off;
end

%% make video of the statics
name=['Tbar_static_RDT'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,material{2})
% tenseg_video_CTS(n_t,C,[1,2],eye(ne),[],[],[],[],[],[],t_t,[],min(substep,50),5,name,savevideo)
tenseg_video_RDT(n_t,C,R,index_b,eye(ne),[],[],[],[45,30],[],[],t_t,[],min(substep,50),5,name,savevideo)

%% dynamics analysis
% time step
dt=1e-5;
tf=0.1;
out_dt=1e-5;
tspan=0:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and 
ind_w=[];w=[];
ind_dl0=[]'; dl0=[]';
ind_dqb=[3*nn+find(sum(E_sa,2))]; dqb0=[2*K_modess(:,2)];    % 12 active roller
% ind_dqb=[find(sum(E_qa*E_qa_a,2))]; dqb0=[0.1*K_mode_qaa_sort(:,2)];    %4 active roller--unstable
[w_t,dqb_t,dqb_d_t,dqb_dd_t,l0_t,E_qa_new,E_qb_new]=tenseg_load_dyna_RDT(tspan,ind_w,w,ind_dqb,dqb0,ind_dl0,dl0,l0,E_qa,E_qb,gravity,[0;9.8;0],C,mass);
% modify external force(optional)
% step=numel(tspan)-1;
% w_t(:,1:step/2)=w_t(:,2:2:end); w_t(:,step/2+1:end)=w_t(:,end)*ones(1,step/2+1);
% dqb_t(:,1:step/2)=dqb_t(:,2:2:end); dqb_t(:,step/2+1:end)=dqb_t(:,end)*ones(1,step/2+1);
% l0_t(:,1:step/2)=l0_t(:,2:2:end); l0_t(:,step/2+1:end)=l0_t(:,end)*ones(1,step/2+1);
% figure
% plot(dqb_dd_t(end-9:end,:)')

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
tenseg_plot_result2(out_tspan,t_t2([size(C_b,1)+8*level+2:4:end-2,size(C_b,1)+8*level+1:4:end-3],:),{'$f_{s_1bottom}$','$f_{s_1middle}$','$f_{s_1top}$','$f_{s_2bottom}$','$f_{s_2middle}$','$f_{s_2top}$'},{'Time (s)','Force (N)'} ...
    ,'plot_member_force.png',saveimg,{'-r','-g','-b','-c','-m','-k'});
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size


%% Plot nodal coordinate curve X Y slid

tenseg_plot_result2(out_tspan,[n_t2([16*3-2:16*3],:)],{'$X_{top}$','$Y_{top}$','$Z_{top}$'}, ...
    {'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg,{'-r','--g','-.b'});
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size

%% Plot configurations
 j=linspace(0.01,1,5);
for i=1:numel(j)
    num=ceil(j(i)*size(n_t2,2));
%  tenseg_plot( reshape(n_t(:,num),3,[]),C_b,C_s,[],[],[]);
tenseg_plot_CTS(reshape(n_t2(:,num),3,[]),C,index_b,S);
%  axis off;
end
%% Plot configurations
 j=linspace(0.01,1,4);
%  fig2=figure;% plot in one figure
for i=1:numel(j)
    num=ceil(j(i)*size(n_t2,2));
%  tenseg_plot( reshape(n_t(:,num),3,[]),C_b,C_s,[],[],[]);
% tenseg_plot_CTS(reshape(n_t(:,num),3,[]),C,index_b,S);
% plot seperate figures
tenseg_plot_RDT(reshape(n_t2(:,num),3,[]),C,R,index_b,eye(ne),[],[],[45,30], [], [],t_t2(:,num),[],[min(t_t2),max(t_t2)]);
% plot in one figure
% tenseg_plot_RDT(reshape(n_t2(:,num),3,[]),C,R,index_b,eye(ne),fig2,[],[45,30], [], [],t_t2(:,num),[],[min(t_t2),max(t_t2)]);
 axis off;
end
%% make video of the dynamic
name2=['Tbar_dynamic_RDT',num2str(tf),'-3'];
% tenseg_video(n_t,C_b,C_s,[],min(numel(out_tspan),50),name,savevideo,material{2})
% tenseg_video_CTS(n_t2,C,index_b,eye(ne),[],[],[],[],[],[],t_t2,[],min(numel(out_tspan),50),5,name,savevideo)
tenseg_video_RDT(n_t2,C,R,index_b,eye(ne),[],[],[],[45,30],[],[],t_t2,[],min(substep,50),5,name2,savevideo);
