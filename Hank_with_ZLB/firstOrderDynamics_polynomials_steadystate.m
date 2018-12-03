function [ys,check] = firstOrderDynamics_polynomials_steadyState(ys,exo)

% Computes stationary equilibrium of the model for Dynare; format is required
% to be called by Dynare (follows example of NK_baseline.mod in Dynare examples)
%
% Thomas Winberry, July 26th, 2016

tStart = tic;
fprintf('\nComputing steady state...\n')

%----------------------------------------------------------------
% Call parameters (the next set of commands will overwrite some)
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Read in parameters from Dynare declaration
%----------------------------------------------------------------

% Initialize indicator
check = 0;

% Read parameters from Dynare
global M_ 
% Read out parameters to access them with their name
for iParameter = 1:M_.param_nbr
  paramname = deblank(M_.param_names(iParameter,:));
  eval(['global ' paramname]);
  eval([ paramname ' = M_.params(' int2str(iParameter) ');']);
end
%----------------------------------------------------------------
% Steady State
%----------------------------------------------------------------
displayOpt = 'off';       % 'iter-detailed' or 'off'
coreSteadyState;

% Prices
r = aalpha * (aggregateCapital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = (aggregateCapital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

%----------------------------------------------------------------
% Save values of steady state variables for Dynare (must be exactly
% as declared in Dynare)
%----------------------------------------------------------------

% Coefficients on conditional expectation function
for iEpsilon = 1 : nEpsilon
	for iAsset = 1 : nAssets
		eval(sprintf('expectationCoefficient_%d_%d = mCoefficients(iEpsilon,iAsset);',...
			iEpsilon,iAsset));
	end
end

% Moments and parameters of density away from borrowing constraint
for iEpsilon = 1 : nEpsilon
	for iMoment = 1 : nMeasure
		eval(sprintf('moment_%d_%d = mMoments(iEpsilon,iMoment);',iEpsilon,iMoment));
		eval(sprintf('measureCoefficient_%d_%d = mParameters(iEpsilon,iMoment+1);',iEpsilon,iMoment));
	end
end

% Mass at borrowing constraint
for iEpsilon = 1 : nEpsilon
	eval(sprintf('mHat_%d = mHat(iEpsilon);',iEpsilon));
end

% Other variables
aggregateCapital = (1 - mHat(1,1)) * (1 - aggEmployment) * mMoments(1,1) + (1 - mHat(2,1)) * aggEmployment * mMoments(2,1);
aggregateTFP = 0;
logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (aggEmployment ^ (1 - aalpha)));
logAggregateInvestment = log(ddelta * aggregateCapital);
logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));
logWage = log(w);

% More Other variables that we specified in equations
beta = bbeta;
%r = 0;
i = 0.0001; %is small so that it's easy to get to ZLB by shocking beta
pApB = 1;
S = 0;
pi = 0;
beta=bbeta;
inot = 0.0001;
ma1 = 0; ma2 = ma1; ma3 = ma2; ma4 = ma3; ma5 = ma4;
ma6 = ma5; ma7 = ma6; ma8 = ma7; ma9 = ma8; ma10 = ma9;
ma11 = ma10; ma12 = ma11; ma13 = ma12; ma14 = ma13; ma15 = ma14;
ma16 = ma15; ma17 = ma16; ma18 = ma17; ma19 = ma18; ma20 = ma19;
ma21 = ma20; ma22 = ma21; ma23 = ma22; ma24 = ma23; ma25 = ma24;
ma26 = ma25; ma27 = ma26; ma28 = ma27; ma29 = ma28; ma30 = ma29;
ma31 = ma30; ma32 = ma31;

%another MA process for forward guidance (don't need for the usual case)
%m1 = 0;
%m2 = m1; m3 = m2; m4 = m3; m5 = m4;
%m6 = m5; m7 = m6; m8 = m7; m9 = m8; m10 = m9;
%m11 = m10; m12 = m11; m13 = m12; m14 = m13; m15 = m14;
%m16 = m15; m17 = m16; m18 = m17; m19 = m18; m20 = m19;
%m21 = m20; m22 = m21; m23 = m22; m24 = m23; 

% Save endogenous variables back into ys
for ii = 1 : M_.orig_endo_nbr
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']); 
end

fprintf('... Done!  Elapsed time: %2.2f seconds \n\n',toc(tStart))
