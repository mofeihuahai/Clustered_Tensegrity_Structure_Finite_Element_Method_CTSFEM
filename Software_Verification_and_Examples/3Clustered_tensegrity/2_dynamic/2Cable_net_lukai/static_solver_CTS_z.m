function data_out=static_solver_CTS_z(data)
%solve nonlinear equilibrium equations using modified Newton method
%converge to stable equilibrium, considering substep, for CTS( including
%TTS)

global E0 A l0 Ia Ib C S w ne Xb Xa dXa f_int l_int 
% minimize total energy? (1: use, 0: not use) it's time consuming
use_energy=1;
%% input data
p=12;
C=data.C;
ne=data.ne;
Ia=data.Ia;
Ib=data.Ib;
S=data.S;
index_b=data.index_b;
index_s=data.index_s;
substep=data.substep;
E0=data.E;
consti_data=data.consti_data;
material=data.material;
A=data.A;
% w0=data.w;
if  isfield(data,'w_t')
    if size(data.w_t,2)==substep
        w_t=data.w_t;
    elseif size(data.w_t,2)==1
        w_t=data.w_t*linspace(0,1,substep);
    end
else
    w_t=linspace(0,0,substep);
end
% if  isfield(data,'w_t')
%     w_t=data.w_t;
% end         %1130+

% dXb=data.dXb;
if  isfield(data,'dnb_t')
    if size(data.dnb_t,2)==substep
        dXb_t=data.dnb_t;
    elseif size(data.dnb_t,2)==1
        dXb_t=data.dnb_t;
    end
else
    dXb_t=linspace(0,0,substep);
end
% l0_0=data.l0;
if size(data.l0_t,2)==substep
    l0_t=data.l0_t;
elseif size(data.l0_t,2)==1
    l0_t=data.l0_t*linspace(1,1,substep);
end

if  isfield(data,'subsubstep')
    subsubstep=data.subsubstep;
else
    subsubstep=30;          %default ssubstep
end

X0=data.N(:);
data_out=data;     %initialize output data
data_out.E_out=E0*ones(1,substep);
%% calculate equilibrium
pinned_X1=([4:4:4*p])'; pinned_Y1=([4:4:4*p])'; pinned_Z1=([4:4:4*p])';
pinned_X2=([3:4:4*p-1])'; pinned_Y2=([3:4:4*p-1])'; pinned_Z2=([3:4:4*p-1])';
pinned_X3=([2:4:4*p-2])'; pinned_Y3=([2:4:4*p-2])'; pinned_Z3=([2:4:4*p-2])';
pinned_X4=([1:4:4*p-3])'; pinned_Y4=([1:4:4*p-3])'; pinned_Z4=([1:4:4*p-3])';
I=eye(3*48);
a1=sort([3*pinned_X1-2;3*pinned_Y1-1;3*pinned_Z1]);
a2=sort([3*pinned_X2-2;3*pinned_Y2-1;3*pinned_Z2]);
a3=sort([3*pinned_X3-2;3*pinned_Y3-1;3*pinned_Z3]);
b=sort([3*pinned_X4-2;3*pinned_Y4-1;3*pinned_Z4]);
I1=I(:,a1);
I2=I(:,a2);
I3=I(:,a3);
I4=I(:,b);
X=X0;               %initialize configuration
X10=I1'*X;
X20=I2'*X;
X30=I3'*X;
Xb0=I4'*X;           %pinned node

E=E0;
% lamda=linspace(0,1,substep);    %coefficient for substep
num_slack=ne*zeros(substep+1,1);    %num of string slack

cont=2;
for k=1:substep
  w=w_t(:,k); %external force
  
    X1=X10;
    X2=X20+0.1*dXb_t(:,k);
    X3=X30+0.9*dXb_t(:,k);
    Xb=Xb0+dXb_t(:,k);         %forced node displacement
    l0=l0_t(:,k);        %forced enlongation of string
    disp(k);
    u=1e-1;
   
    X=[I1';I2';I3';I4']\[X1;X2;X3;Xb];
    Xa=Ia'*X;
    l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length
    l_c=S*l;
    strain=(l_c-l0)./l0;        %strain of member 应变 
    sigma=strain.*E0;%++++
    f_c=sigma.*A;         %member force
%     if substep==1
%         e1=3*ones(1,1);e2=ones(1,1);e3=ones(12,1);e4=4*ones(12,1);ee=[e1;e2;e3;e4];
%         f_c=ee.*f_c;
%     end
    f=S'*f_c;
    q_c=f_c./l_c;
    q=f./l;      %reculate force density
    
    l_int=l;   f_int=f;
    for i=1:1e4
         X=[Ia';Ib']\[Xa;Xb];
        l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length
        l_c=S*l;
        %         q=E.*A.*(1./l0-1./l);      %force density
        strain=(l_c-l0)./l0;        %strain of member 
        sigma=strain.*E0;%++++
        f_c=sigma.*A;         %member force
        f=S'*f_c;
        q_c=f_c./l_c;
        q=f./l;      %reculate force density
        
        q_bar=diag(q);
        
        K=kron(C'*q_bar*C,eye(3));                      %stiffness matrix
        Fp=w-K*X;                                       %unbalanced force
        Fp_a=Ia'*Fp;                                 %see the norm of unbalanced force
        norm(Fp_a)
        if norm(Fp_a)<1e-4
            break
        end
        N=reshape(X,3,[]);
        H=N*C';
       Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H

A_1c=kron(C',eye(3))*blkdiag(Cell_H{:})*diag(l.^-1)*S'*diag(l_c);     % equilibrium matrix
        K_t=K+A_1c*diag(E0.*A./(l_c.^3))*A_1c';
%         for j=1:ne
%             Ki{j,1}=q_bar(j,j)*eye(3)+E(j)*A(j)*l(j)^(-3)*B(:,j)*B(:,j)';
%         end
%         K_t=kron(C',eye(3))*blkdiag(Ki{:})*kron(C,eye(3));
        K_taa=Ia'*K_t*Ia;
        
        %modify the stiffness matrix
        [V_mode,D]=eig(K_taa);                       %刚度矩阵特征根
        d=diag(D);                            %eigen value
        lmd=min(d);                     %刚度矩阵最小特征根
        if lmd>0
            Km=K_taa+u*eye(size(K_taa)); %修正的刚度矩阵
        else
            Km=K_taa+(abs(lmd)+u)*eye(size(K_taa));
        end
        dXa=Km\Fp_a;
        
        x=1;
        % line search
        if use_energy==1
            opt=optimset('TolX',1e-5);
            [x,V]=fminbnd(@energy_CTS,0,1,opt);
        end
        Xa=Xa+x*dXa;
    end
    %
    %     % change youngs mudulus if string slack
    %     strain=(l-l0)./l0;        %strain of member
    %     [E,stress]=stress_strain(consti_data,index_b,index_s,strain,material);
    % %     [E,sigma]=stress_strain(consti_data,index_b,index_s,strain,slack,plastic);
    %     f=stress.*A;         %member force
    %     q=f./l;      %reculate force density
    %     num_slack(k+1)=numel(find(E==0));
    %        % if string slack, recalculate with more steps
    %     if num_slack(k+1)>num_slack(k)
    %         p_s=k-1;
    %           p_e=k;
    %         [E,f,q] = nonlinear_solver(data,Xb0,w_t,dXb_t,l0_t,data_out.E_out(:,k-1),p_s,p_e,subsubstep,material);
    %     end
    %     num_slack(k+1)=numel(find(E==0));
    %
    %     if min(E)==0
    %         if cont<2
    %             [d_sort,idx]=sort(d);               %sorted eigenvalue
    %             D_sort=diag(d_sort);                  %sorted eigenvalue matrix
    %             V_mode_sort=V_mode(:,idx);              %sorted eigenvector
    %             index_bk=find(d_sort<1e-5);             %index for buckling mode
    %             cont=cont+1;
    %             Xa=Xa+0.0*min(l)*real(mean(V_mode_sort(:,index_bk),2));    %add unstable mode if needed
    %         end
    %     end
    
    
    %     if slack
    %         if sum(q_i(index_s)<1e-6)
    %             index_slack=find(q_i(index_s)<0);
    %             index_string_slack=index_s(index_slack);       %slack stings'number
    %             % change youngs mudulus of slack string E_ss=0
    %             E=E0;
    %             E(index_string_slack)=0;
    %             q=E.*A.*(1./l0-1./l);      %reculate force density
    %             q_bar=diag(q);
    %
    %             %give initial error in coordinate, prevent unstable solution
    %             if cont<3
    %             [d_sort,idx]=sort(d);               %sorted eigenvalue
    %             D_sort=diag(d_sort);                  %sorted eigenvalue matrix
    %             V_mode_sort=V_mode(:,idx);              %sorted eigenvector
    %             index_bk=find(d_sort<1e-5);             %index for buckling mode
    %             cont=cont+1;
    %             end
    %             Xa=Xa+0*min(l)*real(mean(V_mode_sort(:,index_bk),2));    %add unstable mode if needed
    %         else
    %             E=E0;              %use initial young's muldus
    %         end
    %     end
    
    
    %% output data
    
    data_out.N_out{k}=reshape(X,3,[]);
    data_out.n_out(:,k)=X;
    %     data_out.l_out(:,k)=l;
    %     data_out.q_out(:,k)=q;
    %     data_out.E_out(:,k)=E;
    data_out.t_out(:,k)=f;      %member force
    % data_out.V{k}=energy_cal(data_out);
    data_out.Fpn_out(k)=norm(Ia'*Fp);
end
data_out.E=E;
data_out.N=reshape(X,3,[]);