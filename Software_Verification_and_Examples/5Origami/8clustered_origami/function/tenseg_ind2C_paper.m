function [ C_pn ] = tenseg_ind2C_paper( Ca,N )
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


nmembers = size(Ca,2); % Number of members being created
nn = size(N,2); % Number of nodes in the structure
C_pni = zeros(1,nn);
C_pn=[];
for i=1:nmembers
    side1=Ca(1,i);
	side2=Ca(2,i);
    side3=Ca(3,i);
	C_pni(1,side1) = 1;
	C_pni(1,side2) = 1;
    C_pni(1,side3) = 1;
    C_pn=[C_pn;C_pni];
    C_pni = zeros(1,nn);
end