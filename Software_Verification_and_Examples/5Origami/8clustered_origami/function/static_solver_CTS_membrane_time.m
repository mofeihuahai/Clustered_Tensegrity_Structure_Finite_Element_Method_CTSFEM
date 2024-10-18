function data_out=static_solver_CTS_membrane_time(data)
%solve nonlinear equilibrium equations using modified Newton method
%converge to stable equilibrium, considering substep, for CTS( including
%TTS)

global E Ia Ib C S w ne Xb Xa dXa l_l_int l_tc_int C_t C_l  S_tc B_epsilon C_pl_bar A_p t_p E_tc A_tc l0_tc l0_l D
% minimize total energy? (1: use, 0: not use) it's time consuming
use_energy=0;

%% input data
C=data.C;
C_l=data.C_l;
C_t=data.C_t;
C_pl_bar=data.C_pl_bar;

ne=data.ne;
Ia=data.Ia;
Ib=data.Ib;
S=data.S;
S_tc=data.S_tc;
D=data.D;
% index_l=data.index_l;
% index_t=data.index_t;
B_epsilon=data.B_epsilon;
substep=data.substep;
E_p=data.E_p;
A_p=data.A_p;
t_p=data.t_p;
E_t=data.E_t;
A_t=data.A_t;
% consti_data=data.consti_data;
% material=data.material;
l0_tc_t=data.l0_tc_t;            
l0_l_t=data.l0_l_t;   
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

% dXb=data.dXb;
if  isfield(data,'dnb_t')
    if size(data.dnb_t,2)==substep
        dXb_t=data.dnb_t;
    elseif size(data.dnb_t,2)==1
        dXb_t=data.dnb_t*linspace(0,1,substep);
    end
else
    dXb_t=linspace(0,0,substep);
end
% l0_0=data.l0;
% if size(data.l0_t,2)==substep
%     l0_t=data.l0_t;
% elseif size(data.l0_t,2)==1
%     l0_t=data.l0_t*linspace(1,1,substep);
% end

if  isfield(data,'subsubstep')
    subsubstep=data.subsubstep;
else
    subsubstep=30;          %default ssubstep
end

X0=data.N(:);
data_out=data;     %initialize output data
data_out.E_out=E_p*ones(1,substep);


%% calculate equilibrium
X=X0;               %initialize configuration
Xb0=Ib'*X;           %pinned node
E=E_p;
% lamda=linspace(0,1,substep);    %coefficient for substep
num_slack=ne*zeros(substep+1,1);    %num of string slack
Xa=Ia'*X;
cont=2;
 u=1e-1;
E_tc=pinv(S_tc')*E_t;   % Young's modulus of truss
A_tc=pinv(S_tc')*A_t;    % area of truss
for k=1:substep
    w=w_t(:,k);               %external force
    Xb=Xb0+dXb_t(:,k);         %forced node displacement
    l0_tc=l0_tc_t(:,k);         %forced enlongation of string
    l0_l=l0_l_t(:,k);         %forced enlongation of string

%     disp(k);
%     X=[Ia';Ib']\[Xa;Xb];
% 
%     l_t=sqrt(sum((reshape(X,3,[])*C_t').^2))';  % elements' length of truss
%     l_tc=S_tc*l_t;
%     l_l=sqrt(sum((reshape(X,3,[])*C_l').^2))';  % elements' length of panel lines
% 
% 
%     l_l_int=l_l; 
%     l_tc_int=l_tc; 
% %     f_int=t;
%     
%     for i=1:1e3
%         X=[Ia';Ib']\[Xa;Xb];
%        
%         l_t=sqrt(sum((reshape(X,3,[])*C_t').^2))';  % elements' length of truss
%         l_tc=S_tc*l_t;
%         l_l=sqrt(sum((reshape(X,3,[])*C_l').^2))';  % elements' length of panel lines
%     
%         H_t=reshape(X,3,[])*C_t';
%         H_l=reshape(X,3,[])*C_l';
%         Delta_l_l=l_l-l0_l;
% 
%         Delta_l_tc=l_tc-l0_tc;
% 
%         Cell_H_t=mat2cell(H_t,3,ones(1,size(H_t,2)));          % transfer matrix H into a cell: Cell_H
%         Cell_H_l=mat2cell(H_l,3,ones(1,size(H_l,2)));          % transfer matrix H into a cell: Cell_H
% 
% 
%         A_2t=kron(C_t',eye(3))*blkdiag(Cell_H_t{:})*diag(l_t.^-1);
%         A_2tc=A_2t*S_tc';
%         A_2l=kron(C_l',eye(3))*blkdiag(Cell_H_l{:})*diag(l_l.^-1);
%         
% % %         t=cputime;
% %         tic;
% %         X=[Ia';Ib']\[Xa;Xb];
% %         l_l=sqrt(sum((reshape(X,3,[])*C_l').^2))';  % elements' length of panel lines
% %         Delta_l_l=l_l-l0_l;
% %         sigma_l=D/(B_epsilon)*C_pl_bar*Delta_l_l;
% %         toc
% % %         e=cputime-t;
% %           elapsedTime = toc;
% 
%         t_tc=diag(E_tc)*diag(A_tc)*diag(l0_tc.^-1)*Delta_l_tc;
%         t_l=(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*Delta_l_l;
%         
%         
%         Fp=w-A_2tc*t_tc-A_2l*t_l;                                       %unbalanced force
%         Fp_a=Ia'*Fp;                                 %see the norm of unbalanced force
%         norm(Fp_a)
%         if norm(Fp_a)<1e-5
%             break
%         end
%        
% Ke_p=A_2l*(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*A_2l';
% Kg_p=kron(C_l'*diag(l_l.^-1)*diag(t_l)*C_l,eye(3))-A_2l*diag(l_l.^-1)*diag(t_l)*A_2l';
% Kt_p=Ke_p+Kg_p;   %partial_2Vp/partial_n*partial_n'
% 
% Ke_tc=A_2tc*diag(E_tc)*diag(A_tc)*diag(l0_tc.^-1)*A_2tc';
% Kg_tc=kron(C_t'*diag(l_t.^-1)*diag(S_tc'*t_tc)*C_t,eye(3))-A_2t*diag(l_t.^-1)*diag(S_tc'*t_tc)*A_2t';
% Kt_tc=Ke_tc+Kg_tc;   %partial_2Vt/partial_n*partial_n'
% 
% K_t=Kt_p+Kt_tc;   %partial_2V/partial_n*partial_n'
% K_taa=Ia'*K_t*Ia;
% 
% % A_2c=kron(C',eye(3))*blkdiag(Cell_H{:})*diag(l.^-1)*S';     % equilibrium matrix
% % K_t=K+A_2c*diag(E.*A./(l0.^-1))*A_2c';
% %         K_taa=Ia'*0.5*(K_t+K_t')*Ia;
% 
% %         for j=1:ne
% %             Ki{j,1}=q_bar(j,j)*eye(3)+E(j)*A(j)*l(j)^(-3)*B(:,j)*B(:,j)';
% %         end
% %         K_t=kron(C',eye(3))*blkdiag(Ki{:})*kron(C,eye(3));
% 
%         
%         %modify the stiffness matrix
%         [V_mode,D1]=eig(K_taa);                       %刚度矩阵特征根
%         d=diag(D1);                            %eigen value
%         lmd=min(d);                     %刚度矩阵最小特征根
%         if lmd>0
%             Km=K_taa+u*eye(size(K_taa)); %修正的刚度矩阵
%         else
%             Km=K_taa+(abs(lmd)+u)*eye(size(K_taa));
%         end
%         dXa=Km\Fp_a;
% %          dXa=(lmd*eye(size(Km)))\Fp_a;
%         x=1;
%         % line search
%         if use_energy==1
%             opt=optimset('TolX',1e-5);
%             [x,V]=fminbnd(@energy_membrane,0,1e1,opt);
%         end
%         Xa=Xa+x*dXa;
%     end
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
    %         t=cputime;

    deltal2epsilon=(B_epsilon)\C_pl_bar;
       tic;
        X=[Ia';Ib']\[Xa;Xb];
        l_l=sqrt(sum((reshape(X,3,[])*C_l').^2))';  % elements' length of panel lines
        Delta_l_l=l_l-l0_l;
        epsilon=deltal2epsilon*Delta_l_l;
        
        toc
        sigma_l=D*epsilon;
%         e=cputime-t;

    data_out.N_out{k}=reshape(X,3,[]);
    data_out.n_out(:,k)=X;
    data_out.toc(:,k)=toc;
    %     data_out.l_out(:,k)=l;
    %     data_out.q_out(:,k)=q;
    %     data_out.E_out(:,k)=E;
    data_out.epsilon_out(:,k)=epsilon; 
    data_out.sigma_l_out(:,k)=sigma_l;     %sigma_l
%     data_out.t_tc_out(:,k)=t_tc;      %member force
    % data_out.V{k}=energy_cal(data_out);
%     data_out.Fpn_out(k)=norm(Ia'*Fp);
end
data_out.E=E;
data_out.N=reshape(X,3,[]);









