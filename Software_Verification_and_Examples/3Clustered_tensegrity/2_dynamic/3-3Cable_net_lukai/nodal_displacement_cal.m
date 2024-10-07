function d_j =nodal_displacement_cal(ind_j,N_out1,N_out2)
l=length(ind_j);
d_j=zeros(4,l);
for i=1:l
    d_j(1,i)=N_out2(1,ind_j(i))-N_out1(1,ind_j(i));
    d_j(2,i)=N_out2(2,ind_j(i))-N_out1(2,ind_j(i));
    d_j(3,i)=N_out2(3,ind_j(i))-N_out1(3,ind_j(i));
    d_j(4,i)=sqrt(d_j(1,i)^2+d_j(2,i)^2+d_j(3,i)^2);
end