# Homework 7
### due until 16.12., 8:30

Consider a light body which is moving in the gravitational field of two heavy objects. The motion of the light object will not influence the heavy objects.
The heavy bodies with mass ratio <img src="stuffy_stuff/f1.png" width="70"> circle in the (x,y) plane with frequency 1 around their common center of gravity, which we assume to be the origin.

The equations of motion of the light body are

<p align="center">
<img src="stuffy_stuff/f2.png" width="400">
</p>

Consider the motion of a light body in the field of earth and moon. The mass ratio is <img src="stuffy_stuff/f3.png" width="200">. We restrict ourself to planar motion (z=0) and may use the initial conditions

<p align="center">
<img src="stuffy_stuff/f4.png" width="250">
</p>

This (according to literature) should lead to a periodic motion with period T=17.065216560157.
These orbits are known as Arenstorf orbits.

***
Implement the Dormand-Prince 4/5 RK method with the coefficients

<p align="center">
<img src="stuffy_stuff/rk.png" width="500">
</p>

to solve the ODEs. The first line of *b*-coefficients corresponds to the fifth-order result and the second line to the fourth-order result.

Use the norm <img src="stuffy_stuff/f5.png" width="150">, where <img src="stuffy_stuff/f6.png" width="100">, to implement a step-size control. Limit
<img src="stuffy_stuff/f7.png" width="40"> to reasonable numbers *Tol*, i.e. try values in the range of *Tol* = 1e-3 to *Tol* = 1e-10.
After adjusting the time step continue the integration using the result from the fourth-order scheme (Check for yourself whether it matters which solution you use to continue the integration).


For *Tol* = 1e-5 supply a plot for the trajectory (x(t),y(t)) and plot dt versus t to see how the algorithm adjusts the step-size.
