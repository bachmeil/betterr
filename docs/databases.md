# Working with databases

One of the big advantages of working with R is the extensive database
access. You have the usual relational databases (SQLite, Postgres, MariaDB, etc.),
plus non-relational databases like Mongo, and other storage formats such
as parquet. There are related projects such as DuckDB. This is a huge 
area with many options.

Basically, there's no reason to write any specific wrappers for these
packages. Take the RSQLite package as an example. If you do this query
from R

```
dbGetQuery(salesdb, "select customer, sum(amount) 
from sales group by customer order by sum(amount) desc")
```

the output is a data frame. You can do this from D:

```
auto df = DataFrame(`dbGetQuery(salesdb, "select customer, sum(amount) 
from sales group by customer order by sum(amount) desc")`);
```

and then you have the output you want as a data frame. It's possible that
this could be polished a bit with structs and methods, but it would be
a thin wrapper over R code however you do it, and at least at this stage,
it doesn't meet the value/time threshold for me to work on it. The important
thing is that you have access to all of the popular databases in D, with
tested, reliable interfaces using the functionality already supplied by
this library.
