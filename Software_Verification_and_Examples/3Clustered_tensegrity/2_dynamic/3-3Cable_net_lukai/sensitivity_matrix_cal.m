function [K_loc,K_tloc,K_tw]=sensitivity_matrix_cal(C,E_c,A_c,l_oc,Ia,N,Gp,t)
%% 节点力对杆件原长的灵敏度矩阵K_loc
H=N*C'; 
l=sqrt(diag(H'*H));
l_c=Gp'*l;
t_c=pinv(Gp)*t;
Cell_H=mat2cell(H,3,ones(1,size(H,2)));  
A_2c=kron(C',eye(3))*blkdiag(Cell_H{:})*diag(l.^-1)*Gp;
K_loc=-1*A_2c*diag(E_c)*diag(A_c)/(diag(l_oc)^2);
%% 切线刚度矩阵K_taa
K_gaa=Ia'*kron(C'*diag(diag(l.^-1)*Gp*t_c)*C,eye(3))*Ia;
K_eaa=Ia'*A_2c*diag(E_c.*A_c.*(l_c.^-1))*A_2c'*Ia;
K_taa=K_gaa+K_eaa;
%% 索力对杆件原长的灵敏度矩阵k_tloc
K_tloc=-1*diag(E_c.*A_c.*(l_oc.^-1))*(A_2c'*Ia*pinv(K_taa)*Ia'*K_loc+diag(l_c.*(l_oc.^-1)));
%% 索力与节点外力的灵敏度矩阵K_tw
K_tw=diag(E_c.*A_c.*(l_oc.^-1))*A_2c'*Ia*pinv(K_taa)*Ia';
end
