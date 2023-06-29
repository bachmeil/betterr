# DataFrame

The DataFrame struct is, at this time, primarily used to store data read in from a file.

# DataFile

The DataFile struct is used to read in data from a file. At this time, it only reads in data from .csv files.

Example code:

```
auto df = DataFile("unrate.csv");
DataFrame unrate = df.read();
unrate.print("Unemployment rate, monthly, Jan 2022-Jan 2023");
```

# Comments

DataFrame and DataFile are currently works in progress. There is no point introducing a design that is just a thin wrapper over calling R directly. You can already do that with evalR.

DataFile is a higher priority than DataFrame, because R supports reading data from many formats, and then that data can be changed from DataFrame to other data structures.
