function [ C_pl_bar ,C_pl_bar_i] = tenseg_ind2C_bar( C_pl_in,C_l,Ca )

[npl,~]=size(C_l);

[~,np]=size(Ca);        % ne:No.of element;np:No.of plate

C_pl_bar=[]; %
 C_pli=zeros(3,npl);
 for i=1:np
      side1=C_pl_in(i,1);
      side2=C_pl_in(i,2);
      side3=C_pl_in(i,3);
      C_pli(1,side1)=1;
      C_pli(2,side2)=1;
      C_pli(3,side3)=1;
      C_pl_bar_i{i}=C_pli;        %C_pl,i
     C_pl_bar=[C_pl_bar;C_pli];
     C_pli=zeros(3,npl);
 end