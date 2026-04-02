I started developing the "SolarSystem_2BP" code, but of course the outer planets' orbital paths and orbital periods are so large that results aren't that interesting.
Thus, I cut it down to the "SolarSystem_InnerPlanets" code to get a better, simpler look at 2BP propagation to simulate our solar system.

The "PlanetNine_Nept_3BP" code is pretty cool (if I do say so myself). 
It propagates Neptune's motion considering the effects of the Sun and the theorized "Planet Nine" using a 3BP ODE propagator, then compares these results to a simple Neptune-Sun 2BP ODE propagator (similar to the other two codes). 
While figure 1 shows that the difference between the two models is minimal within one of Neptune's orbital periods (and on the large scale of AU), figure 2 plots the difference between the two models' position predictions over time.
I recommend setting tspan to a large period of time (several of Nine's orbital periods [1 = 15,000 years]) to really see the difference. Currently it is set to 10 orbital periods of Planet Nine (150,000 years).
Focusing on smaller periods of time can expose interesting patterns as well.
Figure 3 plots numerical/truncation error, to show what amount of the drift in figure 2 is dominated by this error. Should be low, even for tspan = 150,000-year.

I created the "Mars Transfer and Ground Track" code for my final project in my orbital dynamics calss during my Master's program. Some cool graphs/animations in that one.
