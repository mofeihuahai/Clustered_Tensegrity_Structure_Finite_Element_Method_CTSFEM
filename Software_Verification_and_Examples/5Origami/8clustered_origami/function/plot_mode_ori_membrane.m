function plot_mode_ori_membrane(Mode,value,N,Ia,C_b,C_s,C_l,C_h,C_rh,l,title,xlb,ylb,num_plt,ampli,saveimg,view_vec,Ca)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function plot the mode shape(eigenvalue) of tensegrity structure
% with eigenvalue

%Inputs:
%   C_h: hinge connectivity matrix
%   C_rh: rigid hinge connectivity matrix
%   Ca: connectivity matrix of triangle elements
switch nargin
    case 14+2
        view_vec=2;
end

figure
plot(1:numel(value),value,'k-o','linewidth',1.5); % original plot
% semilogy(1:numel(value),value,'k-o','linewidth',1.5); %semilogy
set(gca,'fontsize',18);
xlabel(xlb,'fontsize',18,'Interpreter','latex');
ylabel(ylb,'fontsize',18,'Interpreter','latex');
grid on;
if saveimg==1
    saveas(gcf,[title,'.png']);
end

for i=1:numel(num_plt)
    f1=figure;
     title2=({['Mode ',num2str(num_plt(i)),'  \lambda_{',num2str(num_plt(i)),'} = ',num2str(value(num_plt(i)),'%.1e')]});
% title2=[];
%plot buckling mode
    tenseg_plot_ori_membrane(N+ampli*max(l)*reshape(Ia*Mode(:,num_plt(i))/norm(Mode(:,num_plt(i))),3,[]),C_b,C_s,C_l,C_h,C_rh,f1,[],[],title2,[],Ca);
    tenseg_plot_ori_dash_membrane(N,C_b,C_s,C_l,C_h,C_rh,f1,[],[],title2,[],Ca);

     axis off;
        view(view_vec);
    if saveimg==1
        saveas(gcf,[title,' of ',num2str(num_plt(i)),'.png']);
    end
end

end

