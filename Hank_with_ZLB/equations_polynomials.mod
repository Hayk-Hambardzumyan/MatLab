// Specifies equations for firstOrderDynamics.mod
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Conditional expectation (#equations = nEpsilon * nAssets)
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon

	@#for iAssets in 1 : nAssets
	
		// Compute conditional expectation
		# expectation_@{iEpsilon}_@{iAssets} = exp(0
		@#for iPower in 1 : nAssets
			+ expectationCoefficient_@{iEpsilon}_@{iPower} * expectationPoly_@{iAssets}_@{iPower}
		@#endfor
		);
		
		// Compute savings policy
		# assetsPrime_@{iEpsilon}_@{iAssets} = max(w * (mmu * (1 - epsilonGrid_@{iEpsilon}) + (1 - ttau) * 
			epsilonGrid_@{iEpsilon}) + (1 + r) * assetsGrid_@{iAssets} - (expectation_@{iEpsilon}_@{iAssets} ^ 
			(-1 / ssigma)),aaBar);
			
		// Compute next period's consumption
		@#for iEpsilonPrime in 1 : nEpsilon
		
			# expectationPrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets} = exp(0
			@#for iPower in 1 : nAssets
				+ expectationCoefficient_@{iEpsilonPrime}_@{iPower}(+1) * cos((@{iPower} - 1) * acos(min(max(
					2 * ((assetsPrime_@{iEpsilon}_@{iAssets} - assetsMin) / (assetsMax - assetsMin)) - 1,-1),1)))
			@#endfor
			);
			
			# assetsPrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets} = max(w(+1) * (mmu * (1 - epsilonGrid_@{iEpsilonPrime}) + 
				(1 - ttau) * epsilonGrid_@{iEpsilonPrime}) + (1 + r(+1)) * assetsPrime_@{iEpsilon}_@{iAssets} - 
				(expectationPrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets} ^ (-1 / ssigma)),aaBar);
			
			# consumptionPrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets} = w(+1) * (mmu * (1 - epsilonGrid_@{iEpsilonPrime}) + 
				(1 - ttau) * epsilonGrid_@{iEpsilonPrime}) + (1 + r(+1)) * assetsPrime_@{iEpsilon}_@{iAssets} - 
				assetsPrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets};
				
		@#endfor
				
		// Functional equation
		log(expectation_@{iEpsilon}_@{iAssets}) = log(beta * (1 + r(+1)) * (0
		@#for iEpsilonPrime in 1 : nEpsilon
			+ epsilonTransition_@{iEpsilon}_@{iEpsilonPrime} * (consumptionPrime_@{iEpsilonPrime}_@{iEpsilon}_@{iAssets} ^ 
				(-ssigma))
		@#endfor
		));
		
	@#endfor

@#endfor

//----------------------------------------------------------------
// Compute various objects over quadrature grid for integrating distribution
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon

	@#for iAssets in 1 : nAssetsQuadrature
	
		// Compute conditional expectation
		# expectationQuadrature_@{iEpsilon}_@{iAssets} = exp(0
		@#for iPower in 1 : nAssets
			+ expectationCoefficient_@{iEpsilon}_@{iPower} * quadraturePoly_@{iAssets}_@{iPower}
		@#endfor
		);
		
		// Compute savings policy
		# assetsPrimeQuadrature_@{iEpsilon}_@{iAssets} = max(w * (mmu * (1 - epsilonGrid_@{iEpsilon}) + (1 - ttau) * 
			epsilonGrid_@{iEpsilon}) + (1 + r) * quadratureGrid_@{iAssets} - (expectationQuadrature_@{iEpsilon}_@{iAssets} ^ 
			(-1 / ssigma)),aaBar);
			
		// PDF of distribution
		# measurePDF_@{iEpsilon}_@{iAssets} = exp(0 + measureCoefficient_@{iEpsilon}_1 * (quadratureGrid_@{iAssets} - 
			moment_@{iEpsilon}_1(-1))
			@#for iMoment in 2 : nMeasure
				+ measureCoefficient_@{iEpsilon}_@{iMoment} * ((quadratureGrid_@{iAssets} - moment_@{iEpsilon}_1(-1)) ^ @{iMoment} - 
				moment_@{iEpsilon}_@{iMoment}(-1))
			@#endfor
			);
			
	@#endfor
	
	// Total mass of distribution
	# totalMass_@{iEpsilon} = 0
	@#for iAssets in 1 : nAssetsQuadrature
		+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilon}_@{iAssets}
	@#endfor
	;
	
@#endfor

//----------------------------------------------------------------
// Compute various objects at borrowing constraint for integrating distribution
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon

	// Compute conditional expectation
	# expectationBC_@{iEpsilon} = exp(0
	@#for iPower in 1 : nAssets
		+ expectationCoefficient_@{iEpsilon}_@{iPower} * bcPoly_@{iPower}
	@#endfor
	);
	
	// Compute savings policy
	# assetsPrimeBC_@{iEpsilon} = max(w * (mmu * (1 - epsilonGrid_@{iEpsilon}) + (1 - ttau) * 
		epsilonGrid_@{iEpsilon}) + (1 + r) * aaBar - (expectationBC_@{iEpsilon} ^ 
		(-1 / ssigma)),aaBar);
		
@#endfor

//----------------------------------------------------------------
// Relationship between moments of distribution and parameters
// (#equations = nEpsilon * nMeasure)
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon
	
	// First moments (uncentered)
	moment_@{iEpsilon}_1(-1) = (0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * quadratureGrid_@{iAssets} * measurePDF_@{iEpsilon}_@{iAssets}
		@#endfor
		) / totalMass_@{iEpsilon};
	
	// Higher order moments (centered)
	@#for iMoment in 2 : nMeasure
	moment_@{iEpsilon}_@{iMoment}(-1) = (0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilon}_@{iAssets} * ((quadratureGrid_@{iAssets} - 
				moment_@{iEpsilon}_1(-1)) ^ @{iMoment})
		@#endfor
		) / totalMass_@{iEpsilon};
		
	@#endfor
		
@#endfor

//----------------------------------------------------------------
// Law of motion for density away from borrowing constraint 
// (#equations = nEpsilon * nMeasure)
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon
	
	// First moment (uncentered)
	moment_@{iEpsilon}_1 = (0
	@#for iEpsilonTilde in 1 : nEpsilon
		+ ((1 - mHat_@{iEpsilonTilde}(-1)) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} * (0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilonTilde}_@{iAssets} *
				assetsPrimeQuadrature_@{iEpsilonTilde}_@{iAssets}
		@#endfor
		) / totalMass_@{iEpsilonTilde}) + mHat_@{iEpsilonTilde}(-1) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} 
			* assetsPrimeBC_@{iEpsilonTilde}
	@#endfor
	) / epsilonMass_@{iEpsilon};
	
	// Higher order moments (uncentered)
	@#for iMoment in 2 : nMeasure
		moment_@{iEpsilon}_@{iMoment} = (0
		@#for iEpsilonTilde in 1 : nEpsilon
			+ ((1 - mHat_@{iEpsilonTilde}(-1)) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} * (0
			@#for iAssets in 1 : nAssetsQuadrature
				+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilonTilde}_@{iAssets} * 
					(assetsPrimeQuadrature_@{iEpsilonTilde}_@{iAssets} - moment_@{iEpsilon}_1) ^ @{iMoment}
			@#endfor
			) / totalMass_@{iEpsilonTilde}) + mHat_@{iEpsilonTilde}(-1) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon}
				* (assetsPrimeBC_@{iEpsilonTilde} - moment_@{iEpsilon}_1) ^ @{iMoment}
		@#endfor
		) / epsilonMass_@{iEpsilon};
	@#endfor

@#endfor

//----------------------------------------------------------------
// Law of motion for mass at borrowing constraint
// (#equations = nEpsilon)
//----------------------------------------------------------------

@#for iEpsilon in 1 : nEpsilon

	mHat_@{iEpsilon} = (0
	@#for iEpsilonTilde in 1 : nEpsilon
		+ ((1 - mHat_@{iEpsilonTilde}(-1)) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} * (0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * measurePDF_@{iEpsilonTilde}_@{iAssets} *
				(assetsPrimeQuadrature_@{iEpsilonTilde}_@{iAssets} <= aaBar + 1e-8)
		@#endfor
		) / totalMass_@{iEpsilonTilde}) + mHat_@{iEpsilonTilde}(-1) * epsilonMass_@{iEpsilonTilde} * epsilonTransition_@{iEpsilonTilde}_@{iEpsilon} 
			* (assetsPrimeBC_@{iEpsilonTilde} <= aaBar + 1e-8)
	@#endfor
	) / epsilonMass_@{iEpsilon};

@#endfor

//----------------------------------------------------------------
// Factor prices and Capital (# equations = 2)
//----------------------------------------------------------------

# aggregateCapital = (1 - aggEmployment) * moment_1_1(-1) + aggEmployment * moment_2_1(-1);

//r = exp(aggregateTFP) * aalpha * (aggregateCapital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

//----------------------------------------------------------------
// Monetary policy and pricing (# equations = 5)
//----------------------------------------------------------------

// Inflation as a function of the relative prices selected by the firms that update their prices
1+pi=((1-0.15)/(1-0.15*(pApB)^(1/(1-1.2))))^(1-1.2);

//Fisher equation
(1+pi(+1))=(1+i)/(1+r);

//MA process for forward guidance (not used in this case)
//m1 = ssigmaTFP*aggregateTFPShock;
//m2 = m1(-1); m3 = m2(-1); m4 = m3(-1); m5 = m4(-1);
//m6 = m5(-1); m7 = m6(-1); m8 = m7(-1); m9 = m8(-1); m10 = m9(-1);
//m11 = m10(-1); m12 = m11(-1); m13 = m12(-1); m14 = m13(-1); m15 = m14(-1);
//m16 = m15(-1); m17 = m16(-1); m18 = m17(-1); m19 = m18(-1); m20 = m19(-1);
//m21 = m20(-1); m22 = m21(-1); m23 = m22(-1); m24 = m23(-1); 
//i = inot-0.014*m20-0.018*m21-0.017*m22;


//Simplified Taylor rule (assuming SS interest rate is 0.0001)
inot = 0.0001*exp(0)+1.5*logAggregateOutput ;

//Interest rate is 0 as we are at the ZLB
i = inot;

// Prices
pApB = ((1.2 * w * logAggregateOutput) + beta*(1-0.15) * (pi(+1))^(-1.2/(1-1.2)))/(logAggregateOutput + beta*(1-0.15) * (pi(+1))^(-1/(1-1.2)))*pApB(+1);

//----------------------------------------------------------------
// Law of motion for aggregate TFP (# equations = 1)
//----------------------------------------------------------------

//no TFP shock
aggregateTFP = rrhoTFP * aggregateTFP(-1); // + ssigmaTFP*aggregateTFPShock; 

//----------------------------------------------------------------
// Law of motion for beta (discount factor) (# equations = 25) 
//----------------------------------------------------------------

//once shock hits, people know it's gonna persist for 24 periods
log(beta) = log(bbeta)+ma24+ma23+ma22+ma21+ma20+ma19+ma18+ma17+ma16+ma15+ma14+ma13+ma12+ma11+ma10+ma9+ma8+ma7+ma6+ma5+ma4+ma3+ma2+ma1;

//MA process for the persistence of the shock to beta (MA for 24 periods so that we have i=0 exactly 20 periods)
ma1 = ssigmaTFP * someothershock;
ma2 = ma1(-1); ma3 = ma2(-1); ma4 = ma3(-1); ma5 = ma4(-1);
ma6 = ma5(-1); ma7 = ma6(-1); ma8 = ma7(-1); ma9 = ma8(-1); ma10 = ma9(-1);
ma11 = ma10(-1); ma12 = ma11(-1); ma13 = ma12(-1); ma14 = ma13(-1); ma15 = ma14(-1);
ma16 = ma15(-1); ma17 = ma16(-1); ma18 = ma17(-1); ma19 = ma18(-1); ma20 = ma19(-1);
ma21 = ma20(-1); ma22 = ma21(-1); ma23 = ma22(-1); ma24 = ma23(-1); 


//----------------------------------------------------------------
// Auxiliary variables of interest (# equations = 4)
//----------------------------------------------------------------

// Output
logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (aggEmployment ^ (1 - aalpha)));

// Investment
logAggregateInvestment = log(ddelta * aggregateCapital);

// Consumption
logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));

// Wage
logWage = log(w);