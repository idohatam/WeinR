
![Test](inst/images/Logo_Header.png)
---

Welcome to **WeinR** — an R workflow for evaluating and inspecting Oxford Nanopore long-read sequencing data.

This package provides a native R-based tool to calculate key QC metrics and generate visualizations for long-read datasets directly in **RStudio**.

**WeinR is currently in beta** (first release), and we greatly appreciate any feedback or suggestions.

Feedback form:  
https://forms.office.com/r/DZ9Nkrjspq

---
## Quick Start: 

### 1. Install and load the package

We will need to use dev tools. If your machine doesn't have it installed, install as shown below: 
```r
install.packages("devtools")

devtools::install_github("idohatam/WeinR", ref = "dev")
library(WeinR)
```

### Supported Input Files: 

WeinR accepts long-read sequencing files in the following formats:
- FASTQ (.fastq)
- Compressed FASTQ (.fastq.gz)
- BAM (.bam)
  *(BAM files must contain base quality scores in the QUAL field)*

### 2. Running QC with CreateReport()

This function runs QC on one or more sequencing files and optionally generates an HTML report.
