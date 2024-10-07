%%%%% this file plot the dynamic and result of CTS prism in different
%%%%% material propeties
clc; clear all; close all;
%% static result
load prism_static_elasto_plastic.mat
t_t_sta=t_t;
n_t_sta=n_t;

% % plot member force 
% tenseg_plot_result(1:substep,t_t_sta([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
% grid on;
% % Plot nodal coordinate curve X Y
% tenseg_plot_result(1:substep,n_t_sta([[12]'*3-1],:),{'12Y'},{'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg);
% grid on;
%% dynamic rusult 50s case 1
saveimg=0
load prism_dynamic_1linear_elastic.mat
% t_t_dyn0=t_t;
n_t_dyn1=n_t;
out_tspan_1=0:1e-3:1;

% plot member force 
% tenseg_plot_result(out_tspan_0,t_t_dyn0([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
% tenseg_plot_result(out_tspan_1,n_t_dyn1([[12]'*3-1],:),{'12Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);

% plot(linspace(0,1,substep),n_t_sta([3*12-1],:),'--','linewidth',2);
figure
plot(out_tspan_1,n_t_dyn1([[12]'*3]-1,:),'-','linewidth',2);
hold on;
grid on


%% dynamic rusult 100s case 2
load prism_dynamic_1multielastic.mat
% t_t_dyn0=t_t;
n_t_dyn2=n_t;

plot(out_tspan_1,n_t_dyn2([[12]'*3]-1,:),'--','linewidth',2);

%% dynamic rusult 10s case 3
load prism_dynamic_1plastic.mat
% t_t_dyn0=t_t;
n_t_dyn3=n_t;


plot(out_tspan_1,n_t_dyn3([[12]'*3]-1,:),'-.','linewidth',2);

set(gca,'fontsize',18,'linewidth',1.15);
ylabel('Coordinate of 12Y (m) ','fontsize',18);
xlabel('Time (s)','fontsize',18);
legend('Linear elastic solution','Multi-linear elastic solution','Plastic solution','location','best','fontsize',15);

%% dynamic rusult 50s case 4
saveimg=0
load prism_dynamic_4linear_elastic.mat
% t_t_dyn0=t_t;
n_t_dyn4=n_t;
out_tspan_2=0:1e-3:4;

% plot member force 
% tenseg_plot_result(out_tspan_0,t_t_dyn0([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
% tenseg_plot_result(out_tspan_1,n_t_dyn1([[12]'*3-1],:),{'12Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);

% plot(linspace(0,1,substep),n_t_sta([3*12-1],:),'--','linewidth',2);
figure
plot(out_tspan_2,n_t_dyn4([[12]'*3]-1,:),'-.','linewidth',2);
hold on;
grid on


%% dynamic rusult 100s case 5
load prism_dynamic_4multielastic.mat
% t_t_dyn0=t_t;
n_t_dyn5=n_t;

plot(out_tspan_2,n_t_dyn5([[12]'*3]-1,:),'--','linewidth',2);

%% dynamic rusult 10s case 6
load prism_dynamic_4plastic.mat
% t_t_dyn0=t_t;
n_t_dyn6=n_t;


plot(out_tspan_2,n_t_dyn6([[12]'*3]-1,:),'-','linewidth',2);

set(gca,'fontsize',18,'linewidth',1.15);
ylabel('Y-Coordinate of node 12 (m) ','fontsize',18);
xlabel('Time (s)','fontsize',18);
legend('Linear elastic','Multi-linear elastic','Plastic','location','best','fontsize',15);