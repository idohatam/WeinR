
![Test](inst/images/Logo_Header.png)
---

Welcome to **WeinR** — an R workflow for evaluating and inspecting Oxford Nanopore long-read sequencing data.

It allows you to:

- Calculate key quality control (QC) metrics
- Generate HTML QC reports
- Filter and trim reads
- Export processed reads in FASTQ or BAM format

All directly within **R / RStudio.**

**WeinR is currently in beta** (first release), and we greatly appreciate any feedback or suggestions.

Feedback form:  
https://forms.office.com/r/DZ9Nkrjspq

---
## 1. Installation

We will need to use dev tools. If your machine doesn't have it installed, install as shown below: 
```r
install.packages("devtools")

devtools::install_github("idohatam/WeinR", ref = "main")
library(WeinR)
```

### Supported Input Files: 

WeinR accepts long-read sequencing files in the following formats:

- FASTQ (.fastq)
- Compressed FASTQ (.fastq.gz)
- BAM (.bam) *(BAM files must contain base quality scores in the QUAL field)*
  
You can provide:

- A vector of file paths
- OR a directory containing supported files

## 2. Run Quality Control with CreateReport()

This function runs QC on one or more sequencing files and optionally generates an HTML report.

### Example: Multiple Files in a Vector 
```r
obj <- CreateReport(
  files = c(
    "C:/path/to/sample_data_1.fastq",
    "C:/path/to/sample_data_2.fastq"
  ),
  report_name = "name_of_your_report",
  render_report = TRUE,
  force = TRUE
)
```
### Example: Using a Directory 

```r
obj <- CreateReport(
  files = "C:/path/to/sequencing_directory/",
  report_name = "name_of_your_report",
  render_report = TRUE,
  force = TRUE
)

```
### What This Produces

- QC metrics stored in a `LongReadQC` object
- An HTML report (if `render_report = TRUE`)
- A saved QC object:
    - WeinR_Outputs/qc.rds
    
## 3. Process Reads (Filtering / Trimming / Adapter Removal)

After QC, use `ProcessReads()` to filter, trim, or remove adapters.
```r
qc3 <- ProcessReads(
  obj,
  filter = TRUE,
  MinAvgQS = 20,
  MinLength = 100,
  MaxNumberNs = 2,
  AdapterSeq = "CCACGATAAATGCGAAAACTAG",
  MaxMismatchEnd = 3,
  MinOverlapEnd = 60,
  MinInternalDistance = 100,
  MinFragmentLength = 200,
  Start = 3,
  End = 6,
  OutFileType = c("fastq", "bam"),
  outpath = "_trimmed", 
  force = TRUE
)
```
You can customize this step depending on your needs:

- Only filtering
- Only trimming
- Only adapter removal
- Or a full workflow combining all steps

### Typical Workflow Summary

1. Install WeinR
2. Run `CreateReport()` to evaluate raw reads
3. Review QC metrics and report
4. Run `ProcessReads()` to filter/trim
5. Export cleaned FASTQ or BAM files

### Outuput Structre:

Within your working directory, you will see something similar to:

```r
WeinR_Outputs/run_YYYY-MM-DD_HHMMSS/
├── metrics/
|    └── sample_1.fastq
│        ├── meanprQscore.csv
│        ├── Ncount.csv
│        ├── ...
│   └── summary metrics.csv
│   └── additonal sample files...
├── plots/
|    └── sample_1
│       ├── length_hist.png (300 DPI)
│       ├── gc_hist.png
│   └──  additonal sample files...
├── reports/
│   └── weinr_qc_report.html
│   └── weinr_qc_report.rmd
├── processed/  (if ProcessReads run)
│   ├── sample1_processed_date_time.fastq
│   └── ...
└── objects/
    └── qc.rds
```

## Known limitations (Beta — please read!)

WeinR is currently in **beta**, so a few rough edges are expected. In particular:

- **Performance:** Implemented primarily in R, so it may run slower than C/C++ tools on large ONT datasets.
- **Memory usage:** Report generation and plotting can be resource-intensive for multi-GB FASTQ/BAM inputs.
- **Input assumptions:** BAM files must contain base quality scores in the QUAL field. Malformed or truncated FASTQ files may cause errors.

## Bioconductor Status

WeinR is currently in **beta**. We are actively working toward a future submission to **Bioconductor**.

Planned improvements for Bioconductor readiness include:

- **Performance optimization** (speeding up report generation and processing)
- **Improved robustness** for edge-case inputs (very small files, malformed FASTQ, BAM quality-field checks)
- **Enhanced documentation and examples** aligned with Bioconductor standards
- **Downstream compatibility** (output structures and objects that integrate smoothly with common Bioconductor workflows)

If you encounter bugs or have feature requests that would support Bioconductor submission, please share feedback:
https://forms.office.com/r/DZ9Nkrjspq


## License

This project is licensed under the MIT License.

