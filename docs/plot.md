# Plotting

Minimal plotting capabilities have been added to this point. The reason is that there is little motivation to do so. If you want to make a plot in R, you can do something like this:

```
evalRQ(`pdf("myplot.pdf")`);
evalRQ(`plot(` ~ x.name ~ `, main="A Plot of My Data")`);
evalRQ(`dev.off()`);
```

D doesn't bring a lot to the table when making R plots like this. Perhaps in the future I will decide that's not true, but for now I don't have too much motivation to add wrappers to existing R plotting functions.

## Functionality

The betterr.plot module provides a struct for simple plots. It holds the
following data, all of which are optional, other than `name`:

- string name: The R name of the variable to be plotted.
- string type: Plot type, for instance, `l` produces a line plot.
- string main: The main title
- string sub: Subtitle
- string xlab: Label on the x-axis
- string ylab: Label on the y-axis

## Example

This creates a file called myplot.pdf holding a plot of the data in x.

```
Vector x = [1.7, -0.2, 3.4, -0.2, 6.8, 1.2];
auto p = Plot(x.name);
p.main = "Plot of random points";
p.sub = "Just some test data, actually";
p.type = "b";
p.xlab = "x values";
p.ylab = "y values";
p.create("myplot");
```
