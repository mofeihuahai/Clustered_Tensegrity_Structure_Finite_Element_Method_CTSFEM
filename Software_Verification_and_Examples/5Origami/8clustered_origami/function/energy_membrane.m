function V=energy_membrane(x)
% calculate the total energy, x is the cofficient in line search
global    Ia Ib C_t C_l Xa Xb dXa  w S_tc B_epsilon C_pl_bar A_p t_p E_tc A_tc l0_tc l0_l D

X=Ia*(Xa+x*dXa)+Ib*Xb;

l_t=sqrt(sum((reshape(X,3,[])*C_t').^2))';  % elements' length of truss

l_tc=S_tc*l_t;
l_l=sqrt(sum((reshape(X,3,[])*C_l').^2))';  % elements' length of panel lines
    
Delta_l_l=l_l-l0_l;
Delta_l_tc=l_tc-l0_tc;
Vp=0.5*Delta_l_l'*(inv(B_epsilon)*C_pl_bar)'*D*kron((diag(A_p)*diag(t_p)),eye(3))*inv(B_epsilon)*C_pl_bar*Delta_l_l;
Vtc=0.5*Delta_l_tc'*diag(E_tc)*diag(A_tc)*diag(l0_tc.^-1)*Delta_l_tc;
V=Vp+Vtc-w'*X;