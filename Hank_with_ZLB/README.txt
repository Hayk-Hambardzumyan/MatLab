The main file to run is 'runsim_dnk.m', it is taken from (Guerrieri, Iacoviello)'s toolkit.

The modification is that instead of running RANK model with ZLB, it now runs a NK model incorporated in 
Winberry's Heterogeneous Agent model as a baseline (firstOrderDynamics_polinomials.mod) and a similar 
constrained model with ZLB (firstOrderDynamics_polinomials_zlb.mod) as the constrained model.

The only difference between these two files is that in the first one the equation for the interest rate
is derived from Taylor rule while in the second one that equation is replaced with i=0. 
Everything else in (Guerrieri, Iacoviello)'s toolkit is left as it was.


The code in 'runsim_dnk.m' does the following:
1.Sets parameters and computes grids as needed in Winberry's algorithm and saves all these parameters and
grids in separate mat files.
2.Runs Dynare as in Guerrieri's toolkit, solves two models: simulation 1 "No ZLB",
simulation 2 "ZLB" (takes economy to ZLB because of the discount factor shock).
3.The two models now solve HANK type model (Winberry) in Dynare with 8 asset grids (can be increased) 
and two states - employed and unemployed. 
4.Simulation 1 has no constraints for Taylor rule and doesn't change when it hits the ZLB bound (i goes
negative), while simulation 2 has constraint that changes Taylor rule to i = 0 when discount factor shock
pushes i below zero (i=max{0,r_bar+1.5*pi}, i is the maximum of 0 and Taylor rule).
5.Two scenarios (No ZLB and ZLB) are ploted in the same graph to make the effect of the ZLB visable.
6.That was the naive case but it is easy to add the codes for forward guidance (extended case) (e.g. when
i hits ZLB, the central bank sets the nominal rate to zero for several additional quarters beyond what is
implied by the naive policy and then reverts back to the Taylor rule.

  