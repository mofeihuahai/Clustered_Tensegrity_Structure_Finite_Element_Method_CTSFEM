function [ N ] = membrane_f( N_initial,N_f,f )
% N_initial  结构初始节点坐标
% N_f   需要划分三角形的三个节点坐标
% f     划分次数  

 N_new=zeros(3,(1+f)*f/2);  %划分后三角形的节点坐标
 N_new(:,1)=N_f(:,1);
 N_new(:,(1+(1+(f-1))*(f-1)/2))=N_f(:,2);   %把三个顶点赋值到对应位置
 N_new(:,(1+f)*f/2)=N_f(:,3);

num_x1=linspace(N_new(1,1),N_new(1,(1+(1+(f-1))*(f-1)/2)),f);
num_x2=linspace(N_new(1,1),N_new(1,(1+f)*f/2),f);                     %把最外面的边长划分为f个（关于x）
num_x3=linspace(N_new(1,(1+(1+(f-1))*(f-1)/2)),N_new(1,(1+f)*f/2),f);

num_y1=linspace(N_new(2,1),N_new(2,(1+(1+(f-1))*(f-1)/2)),f);
num_y2=linspace(N_new(2,1),N_new(2,(1+f)*f/2),f);                     %把最外面的边长划分为f个（关于y）    
num_y3=linspace(N_new(2,(1+(1+(f-1))*(f-1)/2)),N_new(2,(1+f)*f/2),f);

num_z1=linspace(N_new(3,1),N_new(3,(1+(1+(f-1))*(f-1)/2)),f);
num_z2=linspace(N_new(3,1),N_new(3,(1+f)*f/2),f);                     %把最外面的边长划分为f个（关于z）    
num_z3=linspace(N_new(3,(1+(1+(f-1))*(f-1)/2)),N_new(3,(1+f)*f/2),f);


for i=1:f
    N_new(1,(1+(1+(i-1))*(i-1)/2))=num_x1(1,i);
    N_new(1,((1+i)*i/2))=num_x2(1,i);               %把上面划分的赋值到对应位置（关于x）
    N_new(1,((1+(f-1))*(f-1)/2+i))=num_x3(1,i);  

    N_new(2,(1+(1+(i-1))*(i-1)/2))=num_y1(1,i);
    N_new(2,((1+i)*i/2))=num_y2(1,i);               %把上面划分的赋值到对应位置（关于y）
    N_new(2,((1+(f-1))*(f-1)/2+i))=num_y3(1,i);

    N_new(3,(1+(1+(i-1))*(i-1)/2))=num_z1(1,i);
    N_new(3,((1+i)*i/2))=num_z2(1,i);               %把上面划分的赋值到对应位置（关于z）
    N_new(3,((1+(f-1))*(f-1)/2+i))=num_z3(1,i);

end

for i=2:(f-2)
    num_x=linspace(N_new(1,(1+(1+i)*i/2)),N_new(1,((1+(i+1))*(i+1)/2)),(i+1));
    num_y=linspace(N_new(2,(1+(1+i)*i/2)),N_new(2,((1+(i+1))*(i+1)/2)),(i+1));
    num_z=linspace(N_new(3,(1+(1+i)*i/2)),N_new(3,((1+(i+1))*(i+1)/2)),(i+1));
    N_new(1,((1+(1+i)*i/2):((1+(i+1))*(i+1)/2)))=num_x;                          %把中间的点按位置划分次数并赋值
    N_new(2,((1+(1+i)*i/2):((1+(i+1))*(i+1)/2)))=num_y;
    N_new(3,((1+(1+i)*i/2):((1+(i+1))*(i+1)/2)))=num_z;
end
N_new(:,[2,3,(1+(1+(f-1))*(f-1)/2),(1+f)*f/2])=N_new(:,[(1+(1+(f-1))*(f-1)/2),(1+f)*f/2,2,3]);  %将5,6列与底部两个调换

N_a=setdiff(N_initial',N_f','rows','stable');      
N_b=N_a';   %初始节点坐标去除三角形节点坐标后的节点坐标

N=[N_b,N_new];      %将划分后的节点坐标并入整体坐标内