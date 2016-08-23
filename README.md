# Cylinder-Temp-MATLAB-Sim
Simulate transient spatial and temporal dependence of a cylinder when exposed to a sudden 
change in environment. This task was accomplished through simulation via theconstruction 
of a MATLAB program. Temperature profiles for Biot numbers (Bi=0.01,0.1,1,10,100) were 
created for Fourier numbers (Fo=0.01,0.1,0.2,0.33,1,2) resulting in a total of 30 graphs 
available for analysis. 

Analysis of the results for the varying Biot numbers led to the general conclusion that the 
Lumped Capacitance Method (LCM) assumption was valid for Bi<0.1. For these small Biot values, 
the rate of convection between the cylinder’s surface and the surrounding fluid was determined 
to be negligible when compared to the internal rate of conduction within the object. These 
isothermal results can be visualized by the flat, planar temperature profiles obtained for Bi<0.1. 
Additionally, an increasing rate of cooling was found for increasing Bi values.

Similarly, the first term approximation was proved to be valid for Fo>0.2. Analyses of the results 
illustrate that for these large Fourier numbers, the discrepancy between the first term approximation 
and exact solution were non-existent. This can be attributed to two simple facts: the first of 
which deals with the direct proportionality between the time and Fourier number, while the
second is concerned with the relationship between the temperature profile obtained from the 
approximation and the midpoint temperature. For large Fo values corresponding to longer time 
intervals, the temperature of the centerpoint was determined to be an accurate representation 
of the entire cylinder.
