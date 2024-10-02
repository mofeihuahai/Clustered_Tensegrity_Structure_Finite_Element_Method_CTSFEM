%%%%% this file plot the dynamic and static comparsion of cable dome
clc; clear all; close all;

%% load static result
output_sta=load('output_static.mat','n_t','t_gp_t');
n_t_sta=output_sta.n_t;%nodal coordinate
t_gp_sta=output_sta.t_gp_t;%member force
%% load dynamic result 1s
output_dyn=load('CTS_cable_dome_time_1s.mat','n_t','out_tspan');
n_t=output_dyn.n_t;
out_tspan=output_dyn.out_tspan;
% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,[n_t([3*4-2],:)]...
    ,{'Dynamic solution'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',0);
hold on;
plot(out_tspan,[n_t_sta([3*4-2],fix(linspace(1,size(n_t_sta,2),size(n_t,2))))],'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on;

%% load dynamic result 2s
output_dyn=load('CTS_cable_dome_time_2s.mat','n_t','out_tspan');
n_t=output_dyn.n_t;
out_tspan=output_dyn.out_tspan;
% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,[n_t([3*4-2],:)]...
    ,{'Dynamic solution'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',0);
hold on;
plot(out_tspan,[n_t_sta([3*4-2],fix(linspace(1,size(n_t_sta,2),size(n_t,2))))],'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on;

%% load dynamic result 4s
output_dyn=load('CTS_cable_dome_time_4s.mat','n_t','out_tspan');
n_t=output_dyn.n_t;
out_tspan=output_dyn.out_tspan;
% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,[n_t([3*4-2],:)]...
    ,{'Dynamic solution'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',0);
hold on;
plot(out_tspan,[n_t_sta([3*4-2],fix(linspace(1,size(n_t_sta,2),size(n_t,2))))],'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on;
%% load dynamic result 8s
output_dyn=load('CTS_cable_dome_time_8s.mat','n_t','out_tspan');
n_t=output_dyn.n_t;
out_tspan=output_dyn.out_tspan;
% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,[n_t([3*4-2],:)]...
    ,{'Dynamic solution'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',0);
hold on;
plot(out_tspan,[n_t_sta([3*4-2],fix(linspace(1,size(n_t_sta,2),size(n_t,2))))],'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on;

