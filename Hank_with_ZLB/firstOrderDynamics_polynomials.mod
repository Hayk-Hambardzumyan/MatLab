// Dynare shell which declares model and solves for aggregate dynamics using
// first order approximation (when approximating conditional expectation with 
// polynomials)
//
// Thomas Winberry, July 26th, 2016

//% A simple new-keynesian model with government spending

warning off;

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parameters_polynomials.mod"

//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

@#include "variables_polynomials.mod"

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

@#include "equations_polynomials.mod"

end;

//----------------------------------------------------------------
// 4. Computation
//----------------------------------------------------------------

// Specify shock process

shocks;
//var aggregateTFPShock; stderr 0.1;
var someothershock ;   stderr 0.1;
//var eps_c ; stderr 0.01  ;
//var eps_g ; stderr 0.01  ;
//var eps_z ; stderr 0.01  ;

end;

options_.steadystate.nocheck = 1;

// Compute steady state (nocheck option ensures that Dynare runs even if steady
// state only computed approximately, i.e., with small numerical error)
//steady(nocheck);

// Check regularity conditions (turn on to check)
//check;
//model_diagnostics;
//model_info;

// Simulate
stoch_simul(order=1,hp_filter=100,irf=40) logAggregateOutput pi
	logAggregateConsumption logAggregateInvestment logWage r i;
stoch_simul(order=1,irf=0,periods=1000)   ;
