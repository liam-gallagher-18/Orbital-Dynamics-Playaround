I started developing the "SolarSystem_2BP" code, but of course the outer planets' orbital paths and orbital periods are so large that results aren't that interesting.
Thus, I cut it down to the "SolarSystem_InnerPlanets" code to get a better, simpler look at 2BP propagation to simulate our solar system.

The "PlanetNine_Nept_3BP" code is pretty cool (if I do say so myself). 
It propagates Neptune's motion considering the effects of the Sun and the theorized "Planet Nine" using a 3BP ODE propagator, then compares these results to a simple Neptune-Sun 2BP ODE propagator (similar to the other two codes). 
While figure 1 shows that the difference between the two models is minimal within one of Neptune's orbital periods (and on the large scale of AU), figure 2 plots the difference between the two models' position predictions over time.
I recommend setting tspan to a large period of time (several of Nine's orbital periods [1 = 15,000 years]) to really see the difference. Currently it is set to 10 orbital periods of Planet Nine (150,000 years).
Focusing on smaller periods of time can expose interesting patterns as well.
Figure 3 plots numerical/truncation error, to show what amount of the drift in figure 2 is dominated by this error. Should be low, even for tspan = 150,000-year.'

The "Interplanetary Mission" code was developed in Python with help from Gemini, as this was not a problem I've tried to solve before. I used the contents of "SolarSystem_InnerPlanets" and added a spacecraft which can travel between two planets in the animation. If you scroll down the the "Spacecraft!" section, you can change the starting and ending planets for your mission and the time you want it to take (in Earth years). Be sure the time is <= the length of t_end - t_start (at the start of the code) so it fits in the animation time. The output is the total Delta v required to make that mission a reality... often it is a very large number! It was very cool to see how to calculate a Delta v when the only required inputs were starting and ending planets and mission time -- this is the math I needed an assist from Gemini to develop.

I created the "Mars_Transfer_and_Ground_Track" code for my final project in my orbital dynamics calss during my Master's program. It simulates a transfer orbit from a high Mars orbit into LMO, and generates some cool graphs/animations.
