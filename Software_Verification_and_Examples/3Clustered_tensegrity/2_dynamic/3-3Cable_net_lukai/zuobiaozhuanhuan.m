clear all
n_t=load('n_t.txt');
N=cell(1,10);
N2=[];
Ninew=zeros(3,48);
for i=1:50
   n_ti=n_t(:,i);
   Ni=reshape(n_ti,3,[]);
   Ninew(1,:)=Ni(1,:).*cos(pi/4)+Ni(2,:).*sin(pi/4)+1.05/2;
   Ninew(2,:)=Ni(2,:).*cos(pi/4)-Ni(1,:).*sin(pi/4)+1.05/2;
   Ninew(3,:)=Ni(3,:);
   N{i}=Ninew;
   N2=[N2;Ninew];
end
save('N2.txt','N2','-ascii');