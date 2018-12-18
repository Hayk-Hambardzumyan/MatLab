% Shell which declares parameters and calls Dynare to solve model using
% first order approximation of aggregate dynamics
%
% Thomas Winberry, July 26th, 2016

clear all
close all
clc

cd('./Auxiliary Functions');

%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
computeGrids;

% Polynomials over grids (only if using polynomials to approximate conditional expectation)
if splineOpt == 0
	computePolynomials;
end

%----------------------------------------------------------------
% Save parameters in .mat files to import into Dynare 
%----------------------------------------------------------------

if splineOpt == 0	% if using polynomials to approximate individual decisions

	% Economic parameters
	save economicParameters.mat bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
		mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
		
	% Approximation parameters
	save approximationParameters.mat nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
		nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
		
	% Grids
	save grids.mat vAssetsGridZeros vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
		vAssetsGridFineZeros mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
		vAssetsGridQuadratureZeros vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
		
	% Polynomials
	save polynomials.mat vAssetsPoly vAssetsPolySquared vAssetsPolyFine vAssetsPolyQuadrature vAssetsPolyBC
	
else	% if using splines to approximate individual decisions

	% Economic parameters
	save economicParameters.mat bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
		mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
		
	% Approximation parameters
	save approximationParameters.mat nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
		nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
		
	% Grids
	save grids.mat vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
		mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
		vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
	
end

%----------------------------------------------------------------
% Run Dynare
%----------------------------------------------------------------

if splineOpt == 0	% if using polynomials to approximate individual decisions

	dynare firstOrderDynamics_polynomials
	
else	% if using splines to approximate individual decisions

	dynare firstOrderDynamics_splines
	
end

cd('../')

time = 1:40; 
z = zeros(1,40);

figure; 
subplot(5,2,1);
plot(time, oo_.irfs.logAggregateOutput_aggregateTFPShock, '-r',time, z, '--k');
title('Y')
subplot(5,2,2);
plot(time, oo_.irfs.logAggregateConsumption_aggregateTFPShock,'-r',time, z, '--k');
title('C')
subplot(5,2,3);
plot(time, oo_.irfs.logAggregateInvestment_aggregateTFPShock,'-r',time, z, '--k');
title('I')
subplot(5,2,4);
plot(time, oo_.irfs.logWage_aggregateTFPShock,'-r',time, z, '--k');
title('W')
subplot(5,2,5);
plot(time, oo_.irfs.r_aggregateTFPShock,'-r',time, z, '--k');
title('r')
subplot(5,2,6);
plot(time, oo_.irfs.pi_aggregateTFPShock,'-r',time, z, '--k');
title('pi')
subplot(5,2,7);
plot(time, oo_.irfs.i_aggregateTFPShock,'-r',time, z, '--k');
title('i')
subplot(5,2,8);
plot(time, oo_.irfs.mHat_1_aggregateTFPShock,'-r',time, z, '--k');
title('m1')
subplot(5,2,9);
plot(time, oo_.irfs.mHat_2_aggregateTFPShock,'-r',time, z, '--k');
title('m2')