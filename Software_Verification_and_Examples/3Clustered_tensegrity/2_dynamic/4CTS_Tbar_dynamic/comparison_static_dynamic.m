%%%%% this file plot the dynamic and static comparsion of T bar 
clc; clear all; close all;
%% static result
load cable_net_CTS_static.mat
t_t_sta=t_t;
n_t_sta=n_t;

% plot member force 
tenseg_plot_result(1:substep,t_t_sta([1,2,3,5,6],:),{'1','2','3-4','5','6'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
grid on;
% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t_sta([3*3-1],:)-n_t_sta([3*4-1],:)-2,{'3Y'},{'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;

%% dynamic rusult 1s case1
load cable_net_CTS_dynamic1.mat
t_t_dyn1=t_t;
n_t_dyn1=n_t;
out_tspan_1=out_tspan;

% plot member force 
tenseg_plot_result(out_tspan_1,t_t_dyn1([1:5],:),{'1','2','3','4','5'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_1,0.5*(n_t_dyn1([3*3-1],:)-n_t_dyn1([3*4-1],:)-2),{'3Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,1,substep),0.5*(n_t_sta([3*3-1],:)-n_t_sta([3*4-1],:)-2),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on;
%% dynamic rusult 0.5s case 2
load cable_net_CTS_dynamic0.5.mat
t_t_dyn2=t_t;
n_t_dyn2=n_t;
out_tspan_2=out_tspan;

% plot member force 
tenseg_plot_result(out_tspan_2,t_t_dyn2([1:5],:),{'1','2','3','4','5'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_2,0.5*(n_t_dyn2([3*3-1],:)-n_t_dyn2([3*4-1],:)-2),{'3Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,0.5,substep),0.5*(n_t_sta([3*3-1],:)-n_t_sta([3*4-1],:)-2),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on;
%% dynamic rusult 0.1s case 3
load cable_net_CTS_dynamic0.1.mat
t_t_dyn3=t_t;
n_t_dyn3=n_t;
out_tspan_3=out_tspan;

% plot member force 
tenseg_plot_result(out_tspan_3,t_t_dyn3([1:5],:),{'1','2','3','4','5'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_3,0.5*(n_t_dyn3([3*3-1],:)-n_t_dyn3([3*4-1],:)-2),{'3Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,0.1,substep),0.5*(n_t_sta([3*3-1],:)-n_t_sta([3*4-1],:)-2),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on;
%% dynamic rusult 0.05s case 4
load cable_net_CTS_dynamic0.05.mat
t_t_dyn4=t_t;
n_t_dyn4=n_t;
out_tspan_4=out_tspan;

% plot member force 
tenseg_plot_result(out_tspan_4,t_t_dyn4([1:5],:),{'1','2','3','4','5'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_4,0.5*(n_t_dyn4([3*3-1],:)-n_t_dyn4([3*4-1],:)-2),{'3Y'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,0.05,substep),0.5*(n_t_sta([3*3-1],:)-n_t_sta([3*4-1],:)-2),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on;
