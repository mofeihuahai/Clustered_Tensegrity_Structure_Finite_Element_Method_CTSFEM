%%%%% this file plot the dynamic and static comparsion of CTS prism
clc; clear all; close all;
%% static result
load prism_static.mat
t_t_sta=t_t;
n_t_sta=n_t;

% plot member force 
tenseg_plot_result(1:substep,t_t_sta([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
grid on;
% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t_sta([[12]'*3],:),{'12Z'},{'Substep','Coordinate (m)'},'plot_coordinate.png',saveimg);
grid on;
%% dynamic rusult 50s case 7
load prism_dynamic_50.mat
% t_t_dyn0=t_t;
n_t_dyn7=n_t;
out_tspan_7=0:1e-3:50;

% plot member force 
% tenseg_plot_result(out_tspan_0,t_t_dyn0([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_7,n_t_dyn7([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,50,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on
%% dynamic rusult 100s case 8
load prism_dynamic_100.mat
% t_t_dyn0=t_t;
n_t_dyn8=n_t;
out_tspan_8=0:1e-3:100;
% plot member force 
% tenseg_plot_result(out_tspan_0,t_t_dyn0([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_8,n_t_dyn8([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,100,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on
%% dynamic rusult 10s case 0
load prism_dynamic_10.mat
% t_t_dyn0=t_t;
n_t_dyn0=n_t;
out_tspan_0=0:1e-3:10;

% plot member force 
% tenseg_plot_result(out_tspan_0,t_t_dyn0([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_0,n_t_dyn0([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,10,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on

%% dynamic rusult 5s case 5
load prism_dynamic_5.mat
% t_t_dyn5=t_t;
n_t_dyn5=n_t;
out_tspan_5=0:1e-3:5;

% plot member force 
% tenseg_plot_result(out_tspan_5,t_t_dyn5([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_5,n_t_dyn5([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,5,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on

%% dynamic rusult 1s case1
load prism_dynamic_1.mat
% t_t_dyn1=t_t;
n_t_dyn1=n_t;
out_tspan_1=0:1e-3:1;

% % plot member force 
% tenseg_plot_result(out_tspan_1,t_t_dyn1([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_1,n_t_dyn1([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,1,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on
%% dynamic rusult 0.5s case 2
load prism_dynamic_0.5.mat
% t_t_dyn2=t_t;
n_t_dyn2=n_t;
out_tspan_2=0:1e-3:0.5;

% % plot member force 
% tenseg_plot_result(out_tspan_2,t_t_dyn2([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_2,n_t_dyn2([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,0.5,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on
%% dynamic rusult 0.1s case 3
load prism_dynamic_0.1.mat
% t_t_dyn3=t_t;
n_t_dyn3=n_t;
out_tspan_3=0:1e-3:0.1;

% % plot member force 
% tenseg_plot_result(out_tspan_3,t_t_dyn3([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_3,n_t_dyn3([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,0.1,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on
%% dynamic rusult 0.05s case 4
load prism_dynamic_0.05.mat
% t_t_dyn4=t_t;
n_t_dyn4=n_t;
out_tspan_4=0:1e-4:0.05;

% % plot member force 
% tenseg_plot_result(out_tspan_4,t_t_dyn4([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_4,n_t_dyn4([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,0.05,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');
grid on
%% dynamic rusult 50s case 6
load prism_dynamic_50.mat
% t_t_dyn6=t_t;
n_t_dyn6=n_t;
out_tspan_6=linspace(0,50,size(n_t,2));

% % plot member force 
% tenseg_plot_result(out_tspan_6,t_t_dyn6([1,9,17,21,25],:),{'bar','diagonal string','bottom string','middle string','top string'},{'Time (s)','Force (N)'},'plot_member_force.png',saveimg);

% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan_6,n_t_dyn6([[12]'*3],:),{'12Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
hold on
plot(linspace(0,50,substep),n_t_sta([3*12],:),'--','linewidth',2);
legend('Dynamic solution','Quasi-static solution');