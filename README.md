# RNA-seq Vignette: DESeq2 Analysis of the Airway Dataset

A tutorial walkthrough of RNA-seq differential expression analysis using DESeq2, built around the `airway` dataset from Himes et al. (2014).

## What this covers

This vignette goes through a full RNA-seq pipeline - from loading raw count data to visualizing results and comparing them to the original publication. The biological question is how dexamethasone (a glucocorticoid used for asthma) changes gene expression in human airway smooth muscle cells.

## Files

| File | What it is |
|------|-----------|
| `RNA-seq-Vignette_DESeq2-Analysis-of-the-Airway-Dataset.Rmd` | Source R Markdown file |
| `RNA-seq-Vignette_DESeq2-Analysis-of-the-Airway-Dataset.html` | Knitted HTML vignette - open this in a browser |
| `app.R` | Interactive Shiny app to explore results |

## View the vignette

👉 [Live vignette](https://yusraamena20.github.io/Vignette_DESeq2/RNA-seq-Vignette_DESeq2-Analysis-of-the-Airway-Dataset.html)

## Run the Shiny app locally

```r
# Install dependencies if needed
install.packages(c("shiny", "DT", "dplyr", "tibble"))
BiocManager::install(c("DESeq2", "airway", "EnhancedVolcano", "ComplexHeatmap"))

# Launch
shiny::runApp("app.R")
```

## Reference

Himes et al. (2014). RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid-Responsive Gene. *PLoS ONE* 9(6): e99625. https://doi.org/10.1371/journal.pone.0099625
