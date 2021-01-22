%% Clifford Algebra Geometric Library
% MSc. Edgar Macias Garcia (edgar.macias@cinvestav.mx)
% Dr. Julio Zamora Esquivel (julio.c.zamora.esquivel@intel.com )
% Prof. Eduardo José Bayro Corrochano (edb@gdl.cinvestav.mx)
% Centro de Investigación y Estudios Avanzados del Instituto Politécnico
% Nacional, Zapopan, México

%%
clc
clear

%Include libraries
addpath('include');
conformal;

%Define a conformal point directly
p1 = e1 + e2 + e3 + 1.5*ei + e0

%Convert euclidean vector to conformal point
p2 = blades.euc2confpoint([0, 3, -5])

%Clifford product 
p3 = p1*p2

%Wedge product
p4 = p1^p2

%Bivector operations
p5 = e12 + e23
p6 = ep + e23
p7 = p5*p6