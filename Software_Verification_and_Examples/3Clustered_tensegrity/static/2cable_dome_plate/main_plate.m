%An double layer tensegrity tower with simplex%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% [1] structure design(calculate equilibrium matrix,
% group matrix,prestress mode, minimal mass design)
% [2] modal analysis(calculate tangent stiffness matrix, material
% stiffness, geometry stiffness, generalized eigenvalue analysis)
% [3] dynamic simulation

%EXAMPLE
clc; clear all; close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=1; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% dynamic analysis set
amplitude=50;            % amplitude of external force of ground motion 
period=0.5;             %period of seismic

dt=0.001;               % time step in dynamic simulation
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=2;                   % final time of dynamic simulation
out_dt=0.02;            % output data interval(approximately, not exatly)
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0) 

%% N C of the structure
% Manually specify node positions of double layer prism.
R=50;          %radius
p=12;          %complexity for cable dome

rate=0;
[N0,C_b,n_qp] =N_plate(R,rate,p);

tenseg_plot(N0,C_b,[]);
tenseg_plot_RBD(N0,C_b,[],[],[],[],[],[],n_qp);

% plot shell
if 1

x=[linspace(,R,20);

y=linspace(0,1,20);

%% y in 2D

z=y';

r=x;

th=linspace(0,2*pi,100);

[R,Th]=meshgrid(r,th);

[X,Y] = pol2cart(Th,R);
Z=(z*ones(size(th)))';

subplot(1,2,1)

mesh(X,Y,Z)

alpha(0.4)
%% y in 3D

[X2,Y2,Z2] = cylinder(y);

Z2=repmat(x',[1 length(X2(1,:))]);

subplot(1,2,2)

mesh(Y2,Z2,X2)
end
%% different deployment rate
num=4;
rate_n=linspace(0,0.8,num);
for i=1:num
    rate=rate_n(i);
[N0,C_b,n_qp] =N_plate(R,rate,p);

% tenseg_plot(N0,C_b,[]);
tenseg_plot_RBD(N0,C_b,[],[],[],[],[],[],n_qp);
axis(80*[-1 1 -1 1]);
end



%%
name=['plate_trajectory'];
% % tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video(n_t,C_b,C_s,[],50,name,savevideo,material{2});




