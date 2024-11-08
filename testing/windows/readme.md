About this directory

This directory is being used to develop the Windows version of betterr. There's not much of interest for anyone else in this directory. When it's ready for use, it'll be moved to the top directory. As of October 2024, the basics are working. It's a matter of getting everything working on Windows.

Requirements:

You have to copy all the R DLLs from your R installation into the directory that you're going to run the program. In addition, you have to install the RTools for your version of R and install RInside using devtools:

```
library(devtools)
install_github("https://github.com/cran/rinside")
```

You can find out where it was installed using the find.package function:

```
find.package("RInside")
```

From there you can dig around (or use search) and find libRInside.dll. On my computer, it's in `C:\Users\lanceb\AppData\Local\R\win-library\4.4\RInside\libs\x64`. Simply installing RInside from the usual channels has never worked for me - there have always been errors when running the program. The only way to get all the configuration set up properly is to compile the package yourself.

There is also a binding generator in this directory.

- makebindings.d contains the workhorse function makeBindings.
- creation.d is a sample program showing how to make the bindings and write the output to the terminal.
- rfunctions.d and rinsidefunctions.d contain structs that hold the functions we're binding. These modules are imported by makebindings.d.
- bindings.d is the finished product after polishing the formatting. The output of creation.d is copied into the file and edited.
- The compilation command is `ldmd2 -i creation.d`.
