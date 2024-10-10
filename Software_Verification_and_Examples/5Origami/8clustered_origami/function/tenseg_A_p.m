function [ A_p ] = tenseg_A_p( Ca,N )
% [ C_paper] = TENSEG_IND2C_PAPER( Ca, N ) creates a connectivity matrix
% from input index notation array and node matrix.
%
% Inputs:
%	Ca:  connectivity of triangle element 
%	N: node matrix (3 x n array for n nodes)
%
% Outputs:
%	A_p: The area of the board


[~,np]=size(Ca);
A_p=zeros(np,1);
for i=1:np
    side1=Ca(1,i);
	side2=Ca(2,i);  %板的节点
    side3=Ca(3,i);
    V1=N(:,side1);
    V2=N(:,side2);  %板的坐标
    V3=N(:,side3);
    V4=V2-V1;       %板的向量
    V5=V3-V1;
    V6=cross(V4,V5);
    a=0.5*sqrt(sum(V6.*V6));  %板的面积
    % a=0.5*norm(V6);
    A_p(i)=a;

end 
A_p=A_p;
