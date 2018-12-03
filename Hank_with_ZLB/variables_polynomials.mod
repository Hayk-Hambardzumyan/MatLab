// Declare variables for firstOrderDynamics.mod
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Conditional expectation coefficients
//----------------------------------------------------------------

@#for iPower in 1 : nAssets
	var expectationCoefficient_1_@{iPower} expectationCoefficient_2_@{iPower};
@#endfor

//----------------------------------------------------------------
// Density of households away from borrowing constraint
//----------------------------------------------------------------

// Moments of the distribution
@#for iMoment in 1 : nMeasure
    var moment_1_@{iMoment} moment_2_@{iMoment};
@#endfor

// Parameters of the distribution
@#for iParameter in 1 : nMeasure
    var measureCoefficient_1_@{iParameter} measureCoefficient_2_@{iParameter};  
@#endfor

//----------------------------------------------------------------
// Mass at borrowing constraint
//----------------------------------------------------------------

var mHat_1 mHat_2;

//----------------------------------------------------------------
// Prices and Monetary variables
//----------------------------------------------------------------

var r w i inot pi pApB;

//----------------------------------------------------------------
// Aggregate TFP and beta (+ MA variables for beta)
//----------------------------------------------------------------

var aggregateTFP beta;

var ma1 ma2 ma3 ma4 ma5 ma6 ma7 ma8 ma9 ma10 ma11 ma12 ma13 ma14 ma15; 
var ma16 ma17 ma18 ma19 ma20 ma21 ma22 ma23 ma24; 
//var ma25 ma26 ma27 ma28 ma29 ma30 ma31 ma32;

//----------------------------------------------------------------
// Auxiliary variables we're interested in
//----------------------------------------------------------------

var logAggregateOutput logAggregateInvestment logAggregateConsumption logWage;

//variables needed for forward guidance (extended case) (can add)
//var m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19 m20 m21 m22 m23 m24; 

//----------------------------------------------------------------
// Shocks
//----------------------------------------------------------------

varexo  someothershock; 
//aggregateTFPShock eps_c eps_g eps_z ; //can add more shocks
