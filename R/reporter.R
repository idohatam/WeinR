CreateReport <- function(
    path          = "report.Rmd",
    title         = "My Report",
    author        = Sys.info()[["user"]],
    date          = format(Sys.Date(), "%B %d, %Y"),
    toc           = TRUE,
    code_folding  = c("none","show","hide"),
    theme         = "cosmo",
    highlight     = "tango",
    includes_css  = NULL,          # if NULL, we’ll write colorblind.css automatically
    overwrite     = FALSE,
    render_html   = TRUE,
    open_browser  = interactive(),
    mfa           = NULL           # LongReadQC object
) {
  stopifnot(!is.null(mfa))
  code_folding <- match.arg(code_folding)
  
  # Ensure directory
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  # Overwrite guard
  if (file.exists(path) && !overwrite) {
    stop("File already exists: ", normalizePath(path),
         "\nSet overwrite = TRUE to replace it.")
  }
  
  # If no CSS supplied, create an accessible “colorblind-friendly” CSS beside the Rmd
  if (is.null(includes_css)) {
    cb_css <- "
/* -------- Colorblind-friendly page theme -------- */
:root{
  --bg: #F7F8FA;         /* soft neutral page background */
  --card-bg: #FFFFFF;    /* white cards for contrast */
  --text: #0B0C0C;       /* near-black text */
  --muted: #4C4F52;      /* muted body text */
  --accent: #0072B2;     /* accessible blue (Okabe–Ito) */
  --accent-2: #E69F00;   /* orange (Okabe–Ito) */
  --focus: #FFBF47;      /* high-visibility focus */
  --border: #DADFE6;     /* light border */
}

html, body {
  background: var(--bg);
  color: var(--text);
}

a, a:visited { color: var(--accent); }
a:hover, a:focus { color: var(--accent-2); outline: none; }
:focus { outline: 3px solid var(--focus); outline-offset: 2px; }

h1, h2, h3, h4 {
  color: var(--text);
  letter-spacing: 0.2px;
}

h1.title { margin-bottom: 0.2rem; }
.subtitle, .author, .date { color: var(--muted); }

.section-card {
  background: var(--card-bg);
  border: 1px solid var(--border);
  border-left: 6px solid var(--accent);
  border-radius: 12px;
  padding: 1.0rem 1.25rem;
  margin: 1.2rem 0;
  box-shadow: 0 1px 3px rgba(0,0,0,0.06);
}

.section-card h2, .section-card h3 {
  margin-top: 0.2rem;
}

hr.section-sep {
  border: none;
  border-top: 1px solid var(--border);
  margin: 1.2rem 0;
}

.table { width: 100% !important; }
.dataTables_wrapper .dataTables_filter input {
  border: 1px solid var(--border);
  border-radius: 6px; padding: 4px 8px;
}
.dataTables_wrapper .dataTables_length select {
  border: 1px solid var(--border);
  border-radius: 6px; padding: 4px 8px;
}
"
css_path <- file.path(dirname(path), "colorblind.css")
writeLines(cb_css, css_path)
includes_css <- "colorblind.css"
  }

# YAML with params (pass S4 object directly)
yaml <- c(
  "---",
  sprintf('title: "%s"', title),
  sprintf('author: "%s"', author),
  sprintf('date: "%s"', date),
  "output:",
  "  html_document:",
  sprintf("    toc: %s", tolower(as.character(toc))),
  "    toc_depth: 3",
  sprintf("    theme: %s", theme),
  sprintf("    highlight: %s", highlight),
  sprintf("    code_folding: %s", code_folding),
  if (!is.null(includes_css)) sprintf("    css: %s", includes_css) else NULL,
  "df_print: paged",
  "self_contained: true",
  "params:",
  "  mfa: !r NULL",
  "---",
  ""
)

# Body (LongReadQC-aware) with section “cards”
body <- c(
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
  "qc <- params$mfa",
  "stopifnot(!is.null(qc))",
  "suppressPackageStartupMessages({",
  "  library(ggplot2)",
  "  library(DT)",
  "  library(gridExtra)",
  "})",
  "",
  "# Helper to compute layout up to 3 columns",
  "calc_layout <- function(n, max_cols = 3) {",
  "  ncol <- min(max_cols, max(1, n))",
  "  nrow <- ceiling(n / ncol)",
  "  list(ncol = ncol, nrow = nrow)",
  "}",
  "```",
  "",
  "<div class='section-card'>",
  "# MFA Analysis",
  "",
  "## Overview",
  "```{r}",
  "cat('Files:', length(qc@files), '\\n')",
  "cat('Metrics tables:', length(qc@metrics), '\\n')",
  "cat('Plot bundles:', length(qc@plots), '\\n')",
  "cat('Summary rows:', NROW(qc@summary_metrics), '\\n')",
  "```",
  "</div>",
  "",
  "<hr class='section-sep'/>",
  "",
  "<div class='section-card'>",
  "## Summary statistics by file",
  "```{r}",
  "summ <- qc@summary_metrics",
  "if (!is.null(summ) && NROW(summ) > 0) {",
  "  if ('file' %in% names(summ)) {",
  "    cols <- c('file', setdiff(names(summ), 'file'))",
  "    summ <- summ[, cols, drop = FALSE]",
  "  }",
  "  DT::datatable(summ, options = list(scrollX = TRUE, autoWidth = TRUE),",
  "                rownames = FALSE, caption = 'Per-file summary metrics')",
  "} else {",
  "  cat('No summary metrics available.')",
  "}",
  "```",
  "</div>",
  "",
  "<hr class='section-sep'/>",
  "",
  "<div class='section-card'>",
  "## Data previews",
  "```{r}",
  "files_data <- names(qc@metrics)",
  "if (is.null(files_data) || length(files_data) == 0) files_data <- as.character(seq_along(qc@metrics))",
  "for (i in seq_along(qc@metrics)) {",
  "  f <- files_data[i]",
  "  cat('--> ', f, '\\n\\n', sep = '')",
  "  df <- qc@metrics[[i]]",
  "  if (is.null(df)) {",
  "    cat('_No metrics table available._\\n\\n')",
  "  } else {",
  "    print(utils::head(df))",
  "    cat('\\n')",
  "  }",
  "}",
  "```",
  "</div>",
  "",
  "<hr class='section-sep'/>",
  "",
  "<div class='section-card'>",
  "## Plots",
  "",
  "### Histograms",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "hists <- lapply(qc@plots, function(ps) if (!is.null(ps$hist)) ps$hist else NULL)",
  "hists <- Filter(Negate(is.null), hists)",
  "if (length(hists)) {",
  "  lay <- calc_layout(length(hists), max_cols = 3)",
  "  gridExtra::grid.arrange(grobs = hists, ncol = lay$ncol, nrow = lay$nrow)",
  "} else {",
  "  cat('_No histograms available._')",
  "}",
  "```",
  "",
  "### Density plots",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "dens <- lapply(qc@plots, function(ps) if (!is.null(ps$density)) ps$density else NULL)",
  "dens <- Filter(Negate(is.null), dens)",
  "if (length(dens)) {",
  "  lay <- calc_layout(length(dens), max_cols = 3)",
  "  gridExtra::grid.arrange(grobs = dens, ncol = lay$ncol, nrow = lay$nrow)",
  "} else {",
  "  cat('_No density plots available._')",
  "}",
  "```",
  "",
  "### Boxplots",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "boxes <- lapply(qc@plots, function(ps) if (!is.null(ps$box)) ps$box else NULL)",
  "boxes <- Filter(Negate(is.null), boxes)",
  "if (length(boxes)) {",
  "  lay <- calc_layout(length(boxes), max_cols = 3)",
  "  gridExtra::grid.arrange(grobs = boxes, ncol = lay$ncol, nrow = lay$nrow)",
  "} else {",
  "  cat('_No boxplots available._')",
  "}",
  "```",
  "</div>",
  "",
  "<hr class='section-sep'/>",
  "",
  "<div class='section-card'>",
  "## Session info",
  "```{r}",
  "sessionInfo()",
  "```",
  "</div>"
)

# Write Rmd
writeLines(c(yaml, body), con = path)
message("Wrote: ", normalizePath(path, winslash = "/"))

# Render with S4 object passed via params
if (isTRUE(render_html)) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required to render HTML. Install it via install.packages('rmarkdown').")
  }
  html_out <- rmarkdown::render(
    path,
    output_format = "html_document",
    params = list(mfa = mfa),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  message("Rendered HTML: ", normalizePath(html_out, winslash = "/"))
  if (open_browser) utils::browseURL(html_out)
  return(invisible(list(rmd = normalizePath(path), html = normalizePath(html_out))))
}

invisible(list(rmd = normalizePath(path)))
}
