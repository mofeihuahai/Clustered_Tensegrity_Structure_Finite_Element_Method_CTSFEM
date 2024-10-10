function V=energy_CTS_ORI(x)
% calculate the total energy, x is the cofficient in line search
global   A E l0 Ia Ib C Xa Xb dXa  w S theta_0 k_h node_in_hinge E_n_total

X=Ia*(Xa+x*dXa)+Ib*Xb;
l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length 
l_c=S*l;
[phpn_e,phTpn,theta]=jacobian_ori(node_in_hinge,reshape(X,3,[]),E_n_total);
V=0.5*(l_c-l0)'*diag(E.*A./l0)*(l_c-l0)+0.5*(theta-theta_0)'*diag(k_h)*(theta-theta_0)-w'*X;  %结构总势能
