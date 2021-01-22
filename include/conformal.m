%% Clifford Algebra Geometric Library
% MSc. Edgar Macias Garcia (edgar.macias@cinvestav.mx)
% Dr. Julio Zamora Esquivel (julio.c.zamora.esquivel@intel.com )
% Prof. Eduardo José Bayro Corrochano (edb@gdl.cinvestav.mx)
% Centro de Investigación y Estudios Avanzados del Instituto Politécnico
% Nacional, Zapopan, México

%% Configurate algebra

format short

%Setup Conformal Algebra basis
ep = blades(0);
e1 = blades(1); e2 = blades(2); e3 = blades(3); e4 = blades(4); e5 = blades(5);
e12 = blades([1,2]); e13 = blades([1,3]); e14 = blades([1,4]); e15 = blades([1,5]);
e23 = blades([2,3]); e24 = blades([2,4]); e25 = blades([2,5]);
e34 = blades([3,4]); e35 = blades([3,5]);
e45 = blades([4,5]); 
e0 = 0.5*(e5 - e4);
ei = e4 + e5;
