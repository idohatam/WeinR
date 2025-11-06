CreateReport <- function(
    path          = "report.Rmd",
    title         = "My Report",
    author        = Sys.info()[["user"]],
    date          = format(Sys.Date(), "%B %d, %Y"),
    toc           = TRUE,
    code_folding  = c("none","show","hide"),
    theme         = "cosmo",
    highlight     = "tango",
    includes_css  = NULL,
    overwrite     = FALSE,
    render_html   = TRUE,
    open_browser  = interactive(),
    mfa           = NULL
) {
  stopifnot(!is.null(mfa))
  code_folding <- match.arg(code_folding)
  
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(path) && !overwrite) {
    stop("File already exists: ", normalizePath(path),
         "\nSet overwrite = TRUE to replace it.")
  }
  
  if (is.null(includes_css)) {
    cb_css <- "
/* ==========================================================
   Accessible, colorblind-friendly stylesheet (Okabe–Ito)
   Applied to R Markdown html_document output
   ========================================================== */

/* -------- Color system (Okabe–Ito) -------- */
:root{
  --bg: #F8FAFB;
  --surface: #FFFFFF;
  --text: #0B0C0C;
  --muted: #4C4F52;
  --border: #E3E7ED;

  /* Okabe–Ito palette */
  --c-blue:   #0072B2;
  --c-orange: #E69F00;
  --c-verm:   #D55E00;
  --c-sky:    #56B4E9;
  --c-green:  #009E73;
  --c-yellow: #F0E442;
  --c-pink:   #CC79A7;
  --c-purple: #7F3C8D;

  --focus: #FFBF47;
}

/* -------- Global layout & typography -------- */
html, body {
  background: var(--bg);
  color: var(--text);
  font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, 'Apple Color Emoji','Segoe UI Emoji';
  line-height: 1.55;
  font-size: 16px;
}

body .main-container {
  max-width: 1100px;
  background: transparent;
  padding: 0 16px 32px 16px;
  margin: 0 auto;
}

/* Headings */
h1.title {
  font-weight: 800;
  letter-spacing: 0.2px;
  margin-bottom: 0.25rem;
}
.subtitle, .author, .date { color: var(--muted); }

/* -------- Header banner -------- */
.header-banner {
  background: linear-gradient(135deg, rgba(0,114,178,0.09) 0%, rgba(0,158,115,0.09) 100%);
  border: 1px solid var(--border);
  border-left: 6px solid var(--c-blue);
  border-radius: 14px;
  padding: 18px 20px;
  margin: 18px 0 8px;
  box-shadow: 0 2px 6px rgba(10, 31, 68, 0.06);
}
.header-banner .meta {
  display: flex;
  gap: 14px;
  flex-wrap: wrap;
  margin-top: 6px;
}
.badge {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  border: 1px solid var(--border);
  border-radius: 999px;
  padding: 4px 10px;
  font-size: 0.85rem;
  color: var(--muted);
  background: var(--surface);
}

/* -------- Cards -------- */
.section-card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-left: 6px solid var(--c-blue);
  border-radius: 14px;
  padding: 16px 18px;
  margin: 18px 0;
  box-shadow: 0 2px 6px rgba(10, 31, 68, 0.06);
}
.section-card h2, .section-card h3 { margin-top: 0.2rem; }

/* -------- TOC (sticky on wide screens) -------- */
#TOC {
  background: var(--surface) !important;
  border: 1px solid var(--border) !important;
  border-radius: 12px;
  padding: 12px !important;
  box-shadow: 0 1px 3px rgba(10,31,68,0.05);
}
@media (min-width: 1100px) {
  .row-fluid .span3 {
    position: sticky;
    top: 14px;
    height: calc(100vh - 20px);
    overflow: auto;
  }
}

/* -------- Tables (DT, base) -------- */
.table, table.dataTable { width: 100% !important; }
.dataTables_wrapper .dataTables_filter input,
.dataTables_wrapper .dataTables_length select {
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 6px 8px;
}
table.dataTable.stripe tbody tr.odd,
table.dataTable.display tbody tr.odd { background-color: rgba(86,180,233,0.06); }
table.dataTable.hover tbody tr:hover { background-color: rgba(0,158,115,0.08); }

/* -------- Figures -------- */
.figure, .rimage img, img {
  border-radius: 10px;
  border: 1px solid var(--border);
}
.caption, figcaption { color: var(--muted); }

/* -------- Code blocks -------- */
pre, code {
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, 'Liberation Mono', monospace;
}
pre code { border-radius: 10px !important; }

/* -------- Links & focus -------- */
a, a:visited { color: var(--c-blue); text-decoration: none; }
a:hover { color: var(--c-green); }
:focus { outline: 3px solid var(--focus); outline-offset: 2px; }

/* -------- Horizontal rule -------- */
hr.section-sep {
  border: none;
  border-top: 1px solid var(--border);
  margin: 1.2rem 0;
}

/* -------- Footer -------- */
.page-footer {
  color: var(--muted);
  text-align: center;
  font-size: 0.9rem;
  margin-top: 20px;
}
"
css_path <- file.path(dirname(path), "colorblind.css")
writeLines(cb_css, css_path)
includes_css <- "colorblind.css"
  }

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
  "    df_print: paged",
  "    self_contained: true",
  if (!is.null(includes_css)) sprintf("    css: %s", includes_css) else NULL,
  "params:",
  "  mfa: !r NULL",
  "---",
  ""
)

body <- c(
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
  "knitr::opts_knit$set(root.dir = dirname(knitr::current_input()))",
  "qc <- params$mfa",
  "stopifnot(methods::is(qc, 'LongReadQC'))",
  "suppressPackageStartupMessages({",
  "  library(ggplot2)",
  "  library(DT)",
  "  library(gridExtra)",
  "})",
  "",
  "# null-coalescing helper",
  "`%||%` <- function(a, b) if (!is.null(a)) a else b",
  "```",
  "",
  "<div class='section-card'>",
  "## Overview",
  "```{r}",
  "cat('Files:          ', length(qc@files), '\\n')",
  "cat('Metrics bundles:', length(qc@metrics), '\\n')",
  "cat('Plot bundles:   ', length(qc@plots), '\\n')",
  "cat('Summary rows:   ', NROW(qc@summary_metrics), '\\n')",
  "if (length(qc@metadata)) {",
  "  cat('Metadata keys:  ', paste(names(qc@metadata), collapse = ', '), '\\n')",
  "}",
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
  "                rownames = FALSE, caption = 'Per-file summary metrics',",
  "                class = 'stripe hover')",
  "} else {",
  "  cat('No summary metrics available.')",
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
  "hists <- lapply(qc@plots, function(ps) {",
  "  if (inherits(ps, 'ggplot')) return(ps)",
  "  if (is.list(ps)) ps[['hist']] %||% ps[['length_hist']] %||% NULL else NULL",
  "})",
  "hists <- Filter(Negate(is.null), hists)",
  "if (length(hists)) {",
  "  n <- length(hists); ncol <- min(3, n); nrow <- ceiling(n / ncol)",
  "  gridExtra::grid.arrange(grobs = hists, ncol = ncol, nrow = nrow)",
  "} else {",
  "  cat('No histograms available')",
  "}",
  "```",
  "",
  "### Density plots",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "dens <- lapply(qc@plots, function(ps) if (is.list(ps)) ps[['density']] else NULL)",
  "dens <- Filter(Negate(is.null), dens)",
  "if (length(dens)) {",
  "  n <- length(dens); ncol <- min(3, n); nrow <- ceiling(n / ncol)",
  "  gridExtra::grid.arrange(grobs = dens, ncol = ncol, nrow = nrow)",
  "} else {",
  "  cat('No density plots available')",
  "}",
  "```",
  "",
  "### Boxplots",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "boxes <- lapply(qc@plots, function(ps) if (is.list(ps)) ps[['box']] else NULL)",
  "boxes <- Filter(Negate(is.null), boxes)",
  "if (length(boxes)) {",
  "  n <- length(boxes); ncol <- min(3, n); nrow <- ceiling(n / ncol)",
  "  gridExtra::grid.arrange(grobs = boxes, ncol = ncol, nrow = nrow)",
  "} else {",
  "  cat('No boxplots available.')",
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
  "</div>",
  "",
  "<div class='page-footer'>",
  "  Generated with an accessible, colorblind-friendly theme (Okabe–Ito).",
  "</div>"
)

writeLines(c(yaml, body), con = path)
message("Wrote: ", normalizePath(path, winslash = "/"))

if (isTRUE(render_html)) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required to render HTML. Install it via install.packages('rmarkdown').")
  }
  html_out <- rmarkdown::render(
    path,
    output_format = "html_document",
    params        = list(mfa = mfa),
    knit_root_dir = dirname(path),
    envir         = new.env(parent = globalenv()),
    quiet         = TRUE
  )
  message("Rendered HTML: ", normalizePath(html_out, winslash = "/"))
  if (open_browser) utils::browseURL(html_out)
  return(invisible(list(rmd = normalizePath(path), html = normalizePath(html_out))))
}

invisible(list(rmd = normalizePath(path)))
} 