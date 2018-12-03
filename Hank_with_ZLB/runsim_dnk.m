%===================================================================
% Econ 611, U of M  
% 
% H. Hambardzumyan, J. Olson, A. Yadav 02/12/2018
% 
% Replication of "The Power of Forward Guidance Revisited", McKay, Nakamura, 
%   and Steinsson (2016) (adding capital to the model) (HANK+ZLB part)
%
% (Combined codes from Winberry (2016)(HANK part) and from Iacoviello
% (2013) (ZLB part))
%=====================================================================

% This part taken from Thomas Winberry, July 26th, 2016, HANK model

% Note that I didn't exclude the IRFs from Winberry part and this code is
% set up to run naive case (can be easilty modified for extended case in
% the model)


clear
close all
clc
set(0,'DefaultLineLineWidth',2)

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
    % Please note that spline method is not available yet

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


% This part is taken from Matteo Iacoviello (2013) RANK model with ZLB

%---------------------------------------------------------------------
% To compute multipliers at ZLB, solve two models
% simulation 1 is baseline  that takes us at the ZLB
% simulation 2 is baseline simulation that takes us at the ZLB plus G shock
%---------------------------------------------------------------------

nperiods=50;
maxiter=30;

% Pick color of charts
solution=1;

irfshock = char('someothershock'); % Shocks we look at: discount factor shock only (can add more shocks e.g. 'aggregateTFPShock')

% simulation 1 is a simulation with a zero baseline
baseline1=[ %0     0      0       0        0       0
            0.00  0.00   0.00    0.00     0.00    0.0 ]';

% In both models, there is a positive someother shock in period 6
scenario1=[%6  0.00   0.00    0.00     0.00    0.00
           6  0.0    0.0     0.0      0.0     0.0   ]';

if splineOpt == 0	% if using polynomials to approximate individual decisions

modnam = 'firstOrderDynamics_polynomials';
modnamstar = 'firstOrderDynamics_polynomials_zlb';

else	% if using splines to approximate individual decisions

modnam = 'firstOrderDynamics_splines';
modnamstar = 'firstOrderDynamics_splines_zlb'; %zlb

end    

constraint = '0<0';
constraint_relax ='0>-0';

% First time we solve simulation only with baseline shocks
[zdatabaseline_lin1 zdatabaseline_pie1 zdatass oobase_ Mbase_] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  baseline1,irfshock,nperiods,maxiter);

% Second time we solve simulation with baseline shocks and scenario
[zdatascenario_lin1 zdatascenario_pie1 zdatass oobase_ Mbase_ ] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  baseline1+scenario1,irfshock,nperiods,maxiter);


% constraint: the constraint (see notes 1 and 2 below). When the condition in constraint evaluates to true, the solution switches from the reference to the alternative regime.
% constraint relax: when the condition in constraint relax evaluates to true, the solution returns to the reference regime.
constraint = 'inot<-0'; 
constraint_relax ='inot>-0';  

% Pick color of charts
simulation=2;

irfshock = char('someothershock'); % Shocks we look at: discount factor shock only (can add more shocks, e.g. 'aggregateTFPShock')

baseline2=[%   0.00   0.00   0.00     0.00     0     0
              0.00   0.00   0.00     0.00     0.0   0.0  ]';

% In both simulations, there is a positive someother shock in period 1
scenario2=[% 6  0.00   0.00    0.00     0.00    0.00
            6  0.00   0.00    0.00     0.00    0.00   ]';


% First time we solve simulation only with baseline shocks
[zdatabaseline_lin2 zdatabaseline_pie2 zdatass oobase_ Mbase_] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  baseline2,irfshock,nperiods,maxiter);

% Second time we solve simulation with baseline shocks and scenario
[zdatascenario_lin2 zdatascenario_pie2 zdatass oobase_ Mbase_ ] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  baseline2+scenario2,irfshock,nperiods,maxiter);


% Note that we compute impulse responses in deviation from baseline
% In simulation=1, baseline1 has a no negative preference shock
% In simulation=2, baseline2 has a negative preference shock that takes economy to ZLB
for i=1:Mbase_.endo_nbr
  eval([deblank(Mbase_.endo_names(i,:)),'1 = zdatascenario_pie1(:,i)-zdatabaseline_pie1(:,i);']);
  eval([deblank(Mbase_.endo_names(i,:)),'2 = zdatascenario_pie2(:,i)-zdatabaseline_pie2(:,i);']);
end




titlelist = char('r','pi','Nom Interest rate','Notional Nom Interest rate','Output','Investment', 'Cons', 'Wage', 'beta');

ylabels = char('% deviation',...
  '% deviation',...
  '% deviation',...
  '% deviation',...
  '% deviation',...
  '% deviation',...
  '% deviation',...
  '% deviation',...
  '% deviation');

figtitle = '';
line1=[r1, pi1, i1, inot1, logAggregateOutput1, logAggregateInvestment1, logAggregateConsumption1, logWage1, beta1];
line2=[r2, pi2, i2, inot2, logAggregateOutput2, logAggregateInvestment2, logAggregateConsumption2, logWage2, beta2];

legendlist = cellstr(char('No ZLB','ZLB binds'));
figlabel = '';
makechart(titlelist,legendlist,figlabel,ylabels,line1,line2)

