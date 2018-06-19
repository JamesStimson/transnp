# **transnp**

### Introduction

This is **transnp**, which takes SNP alignments for clusters, together with sampling dates and some assumptions, and models transmission events and times.

There is markdown for an introductory vignette in the vignettes folder, which you can create on installation (see below).

This contains the work of Yuanwei Xu, James Stimson and Caroline Colijn.

### Installation

You can install **transnp** in **R** using the following command:
```{r}
devtools::install_github("JamesStimson/transnp", build_vignettes = TRUE)
```

There are some example input files that come with the installation. To find out where they are on your system, use system.file() like this:
```{r}
system.file("extdata", "demo_dates.csv", package = "transnp", mustWork = TRUE)
```

You will see something like this in response:
```{r}
[1] "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/transnp/extdata/demo_dates.csv"
```

### Getting help

To view the vignette once installed, run
```{r}
vignette("intro", package = "transnp")
```

Alternatively, you can run the R markdown *vignettes/intro.Rmd* yourself.

If you need further assistance using **transnp**, you can get in touch by emailing 
```{r}
james.stimson16@imperial.ac.uk
```
