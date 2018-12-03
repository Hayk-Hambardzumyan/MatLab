The main file to run is 'runsim_dnk.m', it is taken from Guerrieri, Iacoviello (2013) toolkit, and 
is combined with Winberry's (2016) algorithm.

The code in 'runsim_dnk.m' does the following:
1.Sets parameters and computes grids as needed in Winberry's algorithm and saves all these parameters and
grids in separate mat files.
2.Runs Dynare as in Guerrieri's toolkit, solves two models: simulation 1 - "No ZLB",
simulation 2 - "ZLB" (takes economy to ZLB because of the discount factor shock).
3.The two models now solve HANK type model (Winberry) in Dynare with 8 asset grids (can be increased) 
and two states - employed and unemployed. 
4.Simulation 1 has no constraints for Taylor rule and doesn't change when it hits the ZLB bound (i goes
negative), while simulation 2 has constraint that changes Taylor rule to i = 0 when discount factor shock
pushes i below zero (i=max{0,r_bar+1.5*pi}, i is the maximum of 0 and Taylor rule).
5.Two scenarios (No ZLB and ZLB) are ploted in the same graph.
6.The case mentioned above was for "naive" policy but it is easy to add the codes for forward guidance 
("extended" policy) (e.g. when nominal interest rate hits ZLB, the central bank sets the i = 0 for several
additional quarters beyond what is implied by the "naive" policy and then reverts back to the Taylor rule.  
