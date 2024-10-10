clc;close all; clear;
p=3;
substep_1 = 5;
substep=substep_1*p;
l0=[2;2;2;2;3;4];
dl0_c=[-0.5,-1,2];
ind_dl0_c=[1,3,4];
dl0_i=zeros(size(l0));
% dl0_i(ind_dl0_c)=dl0_c;

[ne,nn]=size(l0);
num_1 = [0];
l0_t = zeros(ne,substep);
for i=1:p

dl0_i=zeros(size(l0))
        dl0_i(ind_dl0_c(i)) = dl0_c(:,i);
%     dl0_p1 = dl0_c(:,i);
%     B(i) = dl0_p1;
%     dl0_p{i}=B;
%    dl0_t{i}=cell2mat(dl0_p(i))'*linspace(0,1,substep/p);
   dl0_t{i}=dl0_i*linspace(0,1,substep/p);
   
   l0_t(:,[1+num_1:substep_1+num_1]) = cell2mat(dl0_t(i));

  num_1= num_1+substep_1;

end

 num_3 = [0];
  for i=1:p-1

      num_2 = l0_t(ind_dl0_c(i),substep_1+num_3)
      l0_t(ind_dl0_c(i),[substep_1+1+num_3:end])=num_2*ones(1,substep-substep_1-num_3)
     num_3= num_3+substep_1;
  end

  l0_t=l0_t+l0*linspace(1,1,substep);