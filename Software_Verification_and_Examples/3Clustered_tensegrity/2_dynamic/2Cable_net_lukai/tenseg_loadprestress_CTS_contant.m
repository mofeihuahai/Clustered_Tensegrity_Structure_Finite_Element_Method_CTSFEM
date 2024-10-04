function [w_t,l0_ct1]=tenseg_loadprestress_CTS_contant(tspan,ind_w,w,type,l0_c_ct,gravity,acc,C,mass)
substep=numel(tspan);
%% external force
%gravity force vector
G=(gravity)*-0.5*kron(abs(C)'*mass,acc);
%initialize force 
w0=zeros(size(G,1),1); %zero external force
w0(ind_w)=w;  %force exerted on bottom nodes
w_t=zeros(numel(G),substep);

switch type
    case 'impluse'
        w_t(:,find(tspan<0.05))=w0*20;        % impluse load in c_index
        w_t=w_t+G*ones(size(tspan));            % add gravity force
    case 'step'
        w_t=w0*ones(size(tspan));        % load in c_index
        w_t=w_t+G*ones(size(tspan));            % add gravity force
    case 'ramp'
        w_t=w0*linspace(0,1,numel(tspan));  % load in ind_w
        w_t=w_t+G*ones(size(tspan));            % add gravity force
end
%% rest length
s=size(l0_c_ct,2);
% c=round(numel(tspan)/(s-1));
c=numel(tspan)/(s-1);
% c_end=numel(tspan)-(s-2)*c;
delta_l=zeros(size(l0_c_ct,1),s-1);
l0_ct1=zeros(size(l0_c_ct,1),substep);
l0_ct1(:,1)=l0_c_ct(:,1);
 for i=1:s-1    
 delta_l(:,i)=l0_c_ct(:,i+1)-l0_c_ct(:,i);
 end
for i=1:substep-1
    j=ceil((i+1)/c);
    switch j
        case s-1
            l0_ct1(:,i+1)=l0_ct1(:,i)+delta_l(:,end)/c;
        case s
            l0_ct1(:,i+1)=l0_ct1(:,i)+delta_l(:,end)/c;
        otherwise
            l0_ct1(:,i+1)=l0_ct1(:,i)+delta_l(:,j)/c ;
    end

end
end