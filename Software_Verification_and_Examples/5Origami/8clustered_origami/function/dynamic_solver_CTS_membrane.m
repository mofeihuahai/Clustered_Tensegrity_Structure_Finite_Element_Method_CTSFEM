function data_out=dynamic_solver_CTS_membrane(data_in)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% This function solve dynamic equations using Runge_Kuta method
%
% Inputs:
%	data_in: data structure describing simulation task
%		[].N: initial node positions
%		[].Ia: transform matrix to get free nodal coordinate
%		[].tspan: time series 
% Outputs:
%	History: data structure containing simulation results
%		[].n_t: %time history of nodal coordinate
%		[].t_t: %time history of members' force
%		[].l_t: %time history of members' length
%% input data

Ia=data_in.Ia;
n0=data_in.N(:);
tspan=data_in.tspan;
n0a_d=data_in.n0a_d;    %initial speed of free coordinates
%% dynamic iteration

% initial value
n0a=Ia'*n0;
% n0a_d=zeros(size(n0a));
Y0a=[n0a;n0a_d];
% Perform simulation
data_out = ode4_CTS_membrane(@tenseg_dyn_x_xdot_CTS_membrane,tspan,Y0a,data_in);


