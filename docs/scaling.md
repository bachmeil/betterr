# Scaling to zero

With a language like D, there's a temptation to focus on performance.
This is done by thinking about overhead and speed as one or more
dimensions of the analysis goes to infinity:

- A very large number of observations
- A very large number of variables
- A very large number of models to fit
- A very large number of iterations
- A very large number of whatever

That's all well and good if it comes without a cost. In the real world,
much of the data analysis that gets done *does not* involve very large
values for any of those dimensions, and that's why languages like R,
Python, and SAS are so heavily used. What matters is convenience and
functionality. The cost of focusing on performance is the loss of
convenience and functionality.

One goal of betterr is to make D convenient and functional for those
situations where you have:

- A small number of observations
- A small number of variables
- A small number of models to fit
- A small number of iterations
- A small number of whatever

D should be a perfectly viable solution any time you would
reach for Python or R. I like writing programs in D because it's a good
language with static typing. I would use it even if performance was as good as Python. If you need
performance, you can drop down to a lower level to handle the
bottlenecks. For instance, if you want to do matrix algebra inside a
loop, you can use [dgretl](https://github.com/bachmeil/dgretl).

I view betterr as a library that scales to zero rather than infinity. 

