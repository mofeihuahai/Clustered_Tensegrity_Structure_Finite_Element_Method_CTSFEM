function [N0,C_b,n_qp] =N_plate(R,rate,p,H);
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
% generate node in one unit
alpha=pi/p;
N0=zeros(3,3*p);
N0(:,1:3)= R*  [rate 0 0; 
        sin(alpha)+rate -cos(alpha) 0;
        -sin(alpha)+rate -cos(alpha) 0]';
N0(3,:)=H;
T1=[cos(2*alpha) -sin(2*alpha) 0
    sin(2*alpha) cos(2*alpha) 0
    0 0 1];
for i=2:p
N0(:,3*i-2:3*i)=T1^(i-1)*N0(:,1:3);
end

C_b_in = kron(ones(p,1),[1 2;2 3;3 1])+kron(3*[0:p-1]',ones(3,2));  % Bar 
C_b = tenseg_ind2C(C_b_in,N0);

n_qp=mat2cell(N0(:),3*3*ones(1,p));
end