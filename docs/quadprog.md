# quadprog

This module gives direct access to the Fortran routine in the [R quadprog package](https://cran.r-project.org/web/packages/quadprog/) underlying the function `solve.QP`.
It solves quadratic programming problems subject to linear constraints.
As it calls straight into the shared library, no R code is executed.

Its use is demonstrated by the program testing/testqp.d. That program
implements the solution to the example problem in the manual using D.

Compilation requires linking to quadprog.so. You can find the full path of
that library by running this command in R:

```
paste0(find.package("quadprog"), "/libs/quadprog.so")
```
