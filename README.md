# Conformal Geometric Algebra Library for MATLAB

Matlab library for Geometric Algebra calculations. 

**Key features** 

* Flexible and friendly user interface to develop operations using Conformal Algebra. 
* Includes most common geometric algebra operations (sum, subtraction, Clifford product, wedge product, dot product, norm). 

## Installation

In order to use the library, add the **/include** folder to your project folder, and include the library by: 

    %Include folder
    addpath('include');
    
    %Add library
    conformal;

## Use 

**Basic operations** 

By loading the conformal script, the elemental basis e1, e2, e3, e4, e5, e0 and ei are available to develop basic operations: 

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
    
## To Do List

* Make a full documentation.
* Currently improving operation speed performance.
* Enable graphic interface for points visualization.

## Citation

If you employ this library please cite as: 

    E.  Macias-Garcia,  J.  Zamora-Esquivel  and  E.  Bayro-Corrochano, Conformal Algebra library for MATLAB. Accessed on: Date, Year. [Online], Available:https://github.com/iqedgarmg/conformal_library, 2020.
