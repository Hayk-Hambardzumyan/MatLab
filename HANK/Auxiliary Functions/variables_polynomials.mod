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
// Prices
//----------------------------------------------------------------

var r w i;

//----------------------------------------------------------------
// Aggregate TFP
//----------------------------------------------------------------

var aggregateTFP;
var ma1 ma2 ma3 ma4 ma5 pi pApB S;
var ma6 ma7 ma8 ma9 ma10; 
var ma11 ma12 ma13 ma14 ma15 ma16 ma17 ma18 ma19 ma20; 
//var ma21 ma22 ma23 ma24 ma25 ma26 ma27 ma28 ma29 ma30 ;

//----------------------------------------------------------------
// Auxiliary variables we're interested in
//----------------------------------------------------------------

var logAggregateOutput logAggregateInvestment logAggregateConsumption logWage;

//----------------------------------------------------------------
// Shocks
//----------------------------------------------------------------

varexo aggregateTFPShock;