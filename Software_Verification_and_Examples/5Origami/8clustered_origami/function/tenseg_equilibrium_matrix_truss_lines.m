function [A_2t,A_2ta,A_2tc,A_2tac,A_2l,A_2la,l_t,l_tc,l_l,l,l_c]=tenseg_equilibrium_matrix_truss_lines(N,C,C_t,C_l,S,S_tc,Ia)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function gives the equilibrium matrix of tensegrity structures
%
% Inputs:
%   N: Node coordinates
%	C: members' length
%	C_t: connectivity matrix of truss
%   C_l: connectivity matrix of panel lines
%   S: clustering matrix
%   S_tc: clustering matrix for truss
%   Ia: Tranform matrix to get free nodal coordinate: na=Ia'*n
%
% Outputs:
%	A_1: equilibrium matrix with no constraints, force density as variable
%	A_1g: Equilirium matrix withgroup constraints
%	density as variable.
%	A_2: equilibrium matrixno  with constraints, force as variable
%	A_2g: Equilirium matrix with group constraints, force
%   l: members' length vector
%   l_gp: members' length vector in group
%	as variable.
%%
H=N*C';
% element length
H_t=N*C_t';
H_l=N*C_l';
% element's direction matrix
l=sqrt(diag(H'*H));         % elements' length 
l_c=S*l;

l_t=sqrt(diag(H_t'*H_t));         % elements' length of truss
l_tc=S_tc*l_t;
l_l=sqrt(diag(H_l'*H_l));         % elements' length of panel lines
%l_c=pinv(S')*l;            % elements' length in group
Cell_H_t=mat2cell(H_t,3,ones(1,size(H_t,2)));          % transfer matrix H into a cell: Cell_H
Cell_H_l=mat2cell(H_l,3,ones(1,size(H_l,2)));          % transfer matrix H into a cell: Cell_H

A_2t=kron(C_t',eye(3))*blkdiag(Cell_H_t{:})*diag(l_t.^-1);
A_2l=kron(C_l',eye(3))*blkdiag(Cell_H_l{:})*diag(l_l.^-1);

A_2ta=Ia'*A_2t;
A_2la=Ia'*A_2l;
A_2tc=A_2t*S_tc';
A_2tac=A_2ta*S_tc';



end

