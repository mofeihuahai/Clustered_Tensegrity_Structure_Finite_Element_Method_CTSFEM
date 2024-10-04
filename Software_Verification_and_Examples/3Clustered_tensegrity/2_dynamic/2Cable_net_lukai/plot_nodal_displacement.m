function [d_j]=plot_nodal_displacement(substep,data,ind_j,N_out2,N_out)
l=length(ind_j);
d_j=cell(1,l);
for i=1:l
    d_j0=zeros(4,substep-1);
    for j=1:substep-1
    d_j0(1,j)=N_out{j+1}(1,ind_j(i))-N_out2(1,ind_j(i));
    d_j0(2,j)=N_out{j+1}(2,ind_j(i))-N_out2(2,ind_j(i));
    d_j0(3,j)=N_out{j+1}(3,ind_j(i))-N_out2(3,ind_j(i));
    d_j0(4,j)=sqrt(d_j0(1,j)^2+d_j0(2,j)^2+d_j0(3,j)^2);
    end
    d_j{i}=d_j0;
end                                        %d_j0输出第一、二、三、四行分别为x、y、z轴和节点位移
%% plot
figure
count=1;
legend1=cell(1,l*size(data,2));
for i=1:l  
  for j=1:size(data,2)
      legend1{count}=['节点' num2str(ind_j(i)) '的' data{j} '方向位移' ];
      count=count+1;
  end
end
for i=1:l  
  for j=1:size(data,2)
      if strcmp(data{j},'dx')
          plot(1:substep-1,d_j{i}(1,:),'-x','linewidth',2);hold on
      elseif strcmp(data{j},'dy')
          plot(1:substep-1,d_j{i}(2,:),'-x','linewidth',2);hold on
      elseif strcmp(data{j},'dz')
          plot(1:substep-1,d_j{i}(3,:),'-x','linewidth',2); hold on
      else
          plot(1:substep-1,d_j{i}(4,:),'-x','linewidth',2); hold on
      end
  end
end
    
set(gca,'fontsize',18,'linewidth',1.15);
legend(legend1,'location','southwest');
ylabel('dz','fontsize',18);
xlabel('load step','fontsize',18);
