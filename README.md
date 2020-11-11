The custom R-package *metabo.pwy* can be installed using the command
```
R CMD INSTALL metabo.pwy_0.1.0.tar.gz
```
Please note that this package has a couple of dependencies to other R-Packages, namely: data.table (>= 1.9.6), stringr (>= 1.0.0), bit64 (>= 4.0.2), tcltk (>= 3.6.1), CHNOSZ (>= 1.3.5), Rdisop (>= 1.46.0), rcdk (>= 3.5.0), lme4 (>= 1.1-21)


Once installed and within R, the metacyc database for metabolite annotation is build via
```
build.db("metacyc")
```

The script for the data analysis and visualisation is provided in `#all_stats.R`.
