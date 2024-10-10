function [ C_pn_bar,n_pn_i, n_pn_i_local,C_pn_i] = tenseg_ind2C_membrane( C_pn,N,R,P_org )
% [ C_paper] = TENSEG_IND2C_PAPER( Ca, N ) creates a connectivity matrix
% from input index notation array and node matrix.
%
% Inputs:
%	Ca:  connectivity of triangle element 
%	N: node matrix (3 x n array for n nodes)
%
% Outputs:
%	C_paper: connectivity matrix of paper
%
% Example: If given four nodes (N is a 3x4 matrix)
%	C_ind would be: C_paper = [1,4;
%                              4,3;
%                              2,2]
% C_paper = tenseg_ind2C_paper([1 4;4 3;2 2], N);


nmembers = size(C_pn,1); % Number of members being created
nn = size(N,2); % Number of nodes in the structure
C_pni = zeros(3,nn);
[np,~]=size(C_pn);        % ne:No.of element;np:No.of plate

for i=1:nmembers
    side1=C_pn(i,1);
	side2=C_pn(i,2);        %选取节点数
    side3=C_pn(i,3);
	C_pni(1,side1) = 1;
	C_pni(2,side2) = 1;         %在对应位置赋值1
    C_pni(3,side3) = 1;
    C_pn_i{i}=C_pni;        %放入数组中
    C_pni = zeros(3,nn);
end
n=N(:);
C_pn_bar=[]; 

for i=1:np
  n_pn_i{i}=kron(cell2mat(C_pn_i(i)),eye(3))*n;  %n_pn,i  一个三角形单元的三个节点坐标
  n_pn_i_local{i}=kron(eye(3),R')*n_pn_i{i}-kron(R',eye(3))*kron(P_org,ones(3,1));
  C_pn_bari=[1,1,1]*cell2mat(C_pn_i(i));       
  C_pn_bar=[C_pn_bar;C_pn_bari];  %不考虑节点顺序的连接关系矩阵
end

