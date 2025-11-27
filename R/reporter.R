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
:root{
  --bg: #F8FAFB; --surface: #FFFFFF; --text: #0B0C0C; --muted: #4C4F52; --border: #E3E7ED;
  --c-blue:#0072B2; --c-orange:#E69F00; --c-verm:#D55E00; --c-sky:#56B4E9;
  --c-green:#009E73; --c-yellow:#F0E442; --c-pink:#CC79A7; --c-purple:#7F3C8D;
  --focus:#FFBF47;
}
html, body { background: var(--bg); color: var(--text); font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, 'Apple Color Emoji','Segoe UI Emoji'; line-height: 1.55; font-size: 16px; }
body .main-container { max-width: 1100px; background: transparent; padding: 0 16px 32px 16px; margin: 0 auto; }
h1.title { font-weight: 800; letter-spacing: 0.2px; margin-bottom: 0.25rem; }
.subtitle, .author, .date { color: var(--muted); }
.section-card { background: var(--surface); border: 1px solid var(--border); border-left: 6px solid var(--c-blue); border-radius: 14px; padding: 16px 18px; margin: 18px 0; box-shadow: 0 2px 6px rgba(10, 31, 68, 0.06); }
.section-card h2, .section-card h3 { margin-top: 0.2rem; }
#TOC { background: var(--surface) !important; border: 1px solid var(--border) !important; border-radius: 12px; padding: 12px !important; box-shadow: 0 1px 3px rgba(10,31,68,0.05); }
@media (min-width: 1100px) { .row-fluid .span3 { position: sticky; top: 14px; height: calc(100vh - 20px); overflow: auto; } }
.table, table.dataTable { width: 100% !important; }
.dataTables_wrapper .dataTables_filter input, .dataTables_wrapper .dataTables_length select { border: 1px solid var(--border); border-radius: 8px; padding: 6px 8px; }
table.dataTable.stripe tbody tr.odd, table.dataTable.display tbody tr.odd { background-color: rgba(86,180,233,0.06); }
table.dataTable.hover tbody tr:hover { background-color: rgba(0,158,115,0.08); }
.figure, .rimage img, img { border-radius: 10px; border: 1px solid var(--border); }
.caption, figcaption { color: var(--muted); }
pre, code { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, 'Liberation Mono', monospace; }
pre code { border-radius: 10px !important; }
a, a:visited { color: var(--c-blue); text-decoration: none; } a:hover { color: var(--c-green); }
:focus { outline: 3px solid var(--focus); outline-offset: 2px; }
hr.section-sep { border: none; border-top: 1px solid var(--border); margin: 1.2rem 0; }
.page-footer { color: var(--muted); text-align: center; font-size: 0.9rem; margin-top: 20px; }
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
  "                rownames = FALSE, caption = 'Per-file summary metrics', class = 'stripe hover')",
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
  "### Read length histograms",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "lens <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['length_hist']] else NULL))",
  "if (length(lens)) {",
  "  n <- length(lens); ncol <- min(3, n); nrow <- ceiling(n / ncol)",
  "  gridExtra::grid.arrange(grobs = lens, ncol = ncol, nrow = nrow)",
  "} else { cat('No length histograms available') }",
  "```",
  "",
  "### Average Q-score histograms",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "qhist <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['quality_hist']] else NULL))",
  "if (length(qhist)) {",
  "  n <- length(qhist); ncol <- min(3, n); nrow <- ceiling(n / ncol)",
  "  gridExtra::grid.arrange(grobs = qhist, ncol = ncol, nrow = nrow)",
  "} else { cat('No Q-score histograms available') }",
  "```",
  "",
  "### GC-content histograms",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "gch <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['gc_hist']] else NULL))",
  "if (length(gch)) {",
  "  n <- length(gch); ncol <- min(3, n); nrow <- ceiling(n / ncol)",
  "  gridExtra::grid.arrange(grobs = gch, ncol = ncol, nrow = nrow)",
  "} else { cat('No GC-content histograms available') }",
  "```",
  "",
  "### Read quality vs length",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "qvl <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['q_vs_length']] else NULL))",
  "if (length(qvl)) {",
  "  n <- length(qvl); ncol <- min(3, n); nrow <- ceiling(n / ncol)",
  "  gridExtra::grid.arrange(grobs = qvl, ncol = ncol, nrow = nrow)",
  "} else { cat('No quality-vs-length plots available') }",
  "```",
  "",
  "### Per-position quality",
  "```{r, fig.width=12, fig.height=4, fig.align='center'}",
  "ppq <- Filter(Negate(is.null), lapply(qc@plots, function(ps) if (is.list(ps)) ps[['per_pos_q']] else NULL))",
  "if (length(ppq)) {",
  "  n <- length(ppq); ncol <- min(3, n); nrow <- ceiling(n / ncol)",
  "  gridExtra::grid.arrange(grobs = ppq, ncol = ncol, nrow = nrow)",
  "} else { cat('No per-position quality plots available') }",
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
  if (open_browser) utils::browseURL(html_out)
  return(invisible(list(rmd = normalizePath(path), html = normalizePath(html_out))))
}

invisible(list(rmd = normalizePath(path)))
}
