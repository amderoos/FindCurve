FindCurve
=========

An R package for computing curves of solutions to non-linear systems of equations

This is a general purpose software package, written in C, that solves non-linear systems of equations of the form

*F*(*x*, *p*) = 0

where *F* is a vector-valued function of a vector *x* of unknowns and *p* represents a parameter. There is little limitation to the formulation of the system of equations, for example, *F*(*x*, *p*) may be defined as the result of an numerical integration of a system of ordinary differential equations (ODEs). For this purpose the package also includes a dedicated numerical solver for ODEs. When a solution point *x* to the system of equations has been found, the package will output all information about the solution point itself and all the values of additional functions *G*(*x*) that the user can define his/herself. Here again, *G*(*x*) is a vector-valued function.

The package is designed to compute curves of points that are solutions to the non-linear equation system as a function of the parameter *p*. In addition, the package has facilities to detect special points in these curves:

- **Branching points (BP):** These points represent special points where 2 different solution curves intersect for a particular parameter value. Although in generic systems of equations this happens only rarely, in systems of equations stemming from biological population dynamics this occurs frequently. These points then represent so-called invasion points, the intersection between a curve representing a zero-valued equilibrium density of a population and a curve representing equilibria with a positive equilibrium density. 

- **Limit points (LP):** These points represent special points where a curve reaches a maximum or minimum parameter value and bends back on its self. When the curve represents equilibrium points of a problem in biological population dynamics such a point marks the threshold where a stable equilibrium state and an unstable state (saddle point) merge and disappear.

- **Extreme points (EXT):** These points represent special points where a curve reaches a maximum or minimum in an objective function, which is one of the elements of the vector-valued function *G*(*x*) discussed above. More specifically these points represent values of the parameter *p* where the partial derivative of one of the objective functions with respect to any arbitrary parameter occurring in the system of equations becomes 0. This can be the parameter *p* but also any other parameter. 

The package also allows for the computation of the location of these special points (BP, LP, EXT) as a function of two parameters *p* and *q*.

The system of equations to be solved has to be implemented in C, but the package itself can be used from the R command-line, for example, using the Rstudio package. To start implementing a particular system of equations, use one of the example files (with an '.h' extension) provided with the package as a template (use the function `showsystem()` to get a list of all the examples). Unfortunately, at the current moment there is no manual available to describe the use of the package. Look at the examples that can be accessed with the function `demo()` for inspiration how to use it and for a demonstration of its capabilities. Furthermore, consult the help page via `?FindCurve` for the syntax to be used in R.

To install the package on your system use the command:

```
devtools::install_github("amderoos/FindCurve")
```
