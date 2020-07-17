# Seroprevalence

Code used in analysis of seroprevalence data in
> Disease-associated antibody phenotypes and probabilistic seroprevalence estimates during the emergence of SARS-CoV-2. De Castro et al,

Authors: [Chris Wallace](https://github.com/chr1swallace), [Stasia Grinberg](https://github.com/stas-g)

## Seroprevalence data

Clean and tidy input data

```{sh}
Rscript -e 'rmarkdown::render("load-and-explore-data.Rmd")'
```

Generate different ML probabilistic predictions

```{sh}

```

Compare different ML methods ("bake off")
```{sh}
Rscript -e 'rmarkdown::render("bake-off.Rmd")'
```

Estimate seroprevalence and plot
```{sh}
Rscript population-inference.R
Rscript figure-overlay-ben.R 
Rscript figure-probpositive.R
```


## Patient data


```{sh}
Rscript -e 'rmarkdown::render("explore-patient-data.Rmd")'
```

