# Path configuration only: this file does not load statistical implementations.
repro_load_config <- function(repro_dir, config_file = "") {
  repro_dir <- normalizePath(repro_dir, mustWork = TRUE)
  spec <- read.delim(file.path(repro_dir, "config", "paths.tsv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  user <- list()
  config_dir <- repro_dir
  if (nzchar(config_file)) {
    if (!file.exists(config_file)) stop("Configuration file not found: ", config_file)
    config_file <- normalizePath(config_file, mustWork = TRUE)
    config_dir <- dirname(config_file)
    user <- dget(config_file)
    if (!is.list(user)) stop("Configuration must return a named list().")
    if (length(user) && (is.null(names(user)) || any(!nzchar(names(user))) || anyDuplicated(names(user))))
      stop("Configuration fields must have unique non-empty names.")
  }
  allowed <- c("data_root", "output_root", "source_root", "paths")
  if (length(setdiff(names(user), allowed)))
    stop("Unknown configuration fields: ", paste(setdiff(names(user), allowed), collapse = ", "))
  absolute <- function(x) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
      stop("Each configured path must be one non-empty string.")
    x <- path.expand(x)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", x)) x <- file.path(config_dir, x)
    if (identical(x, "/")) x else sub("/+$", "", x)
  }
  roots <- list(data = file.path(repro_dir, "data"),
                output = file.path(repro_dir, "outputs"),
                source = file.path(repro_dir, "external"))
  for (kind in names(roots)) {
    field <- paste0(kind, "_root")
    if (!is.null(user[[field]])) roots[[kind]] <- absolute(user[[field]])
  }
  roots$repository <- repro_dir
  overrides <- user$paths
  if (is.null(overrides)) overrides <- list()
  if (!is.list(overrides) || (length(overrides) &&
      (is.null(names(overrides)) || any(!nzchar(names(overrides))) || anyDuplicated(names(overrides)))))
    stop("paths must be a list with unique non-empty path keys.")
  unknown <- setdiff(names(overrides), spec$key)
  if (length(unknown)) stop("Unknown path keys: ", paste(unknown, collapse = ", "))
  resolved <- setNames(vapply(seq_len(nrow(spec)), function(i) {
    key <- spec$key[i]
    if (!is.null(overrides[[key]])) absolute(overrides[[key]])
    else file.path(roots[[spec$root[i]]], spec$relative_path[i])
  }, character(1)), spec$key)
  list(directory = repro_dir, spec = spec, paths = resolved, file = config_file)
}

repro_path_function <- function(config) {
  force(config)
  function(key, relative = NULL) {
    if (length(key) != 1L || is.na(key) || !key %in% names(config$paths))
      stop("Unknown reproducibility path key: ", paste(key, collapse = ", "))
    if (is.null(relative)) return(config$paths[[key]])
    if (!is.character(relative) || length(relative) != 1L || is.na(relative) ||
        grepl("^/|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$)", relative))
      stop("The suffix must be a relative path without '..'.")
    # Preserve a final slash and filename prefixes used by historical paste0 calls.
    paste0(config$paths[[key]], "/", relative)
  }
}

repro_script_paths <- function(script) {
  result <- list()
  if (grepl("\\.py$", script)) {
    helper <- file.path(getOption("annokn.repro_dir"), "python", "run_analysis.py")
    lines <- system2(Sys.getenv("ANNOKN_PYTHON", "python3"),
                     shQuote(c(helper, "--references", script)), stdout = TRUE)
    if (!is.null(attr(lines, "status"))) stop("Could not inspect Python path references.")
    for (line in lines) {
      parts <- strsplit(paste0(line, "\tEND"), "\t", fixed = TRUE)[[1L]]
      result[[length(result) + 1L]] <- list(key = parts[1L], relative = parts[2L])
    }
    return(result)
  }
  walk <- function(x) {
    if (missing(x)) return(invisible(NULL))
    if (is.call(x) && identical(x[[1L]], as.name("repro_path"))) {
      if (length(x) < 2L || !is.character(x[[2L]])) stop("Non-literal path key in ", script)
      result[[length(result) + 1L]] <<- list(key = x[[2L]],
        relative = if (length(x) >= 3L) x[[3L]] else "")
    } else if (is.call(x) || is.expression(x) || is.pairlist(x)) {
      for (element in as.list(x)) walk(element)
    }
  }
  walk(parse(script, keep.source = FALSE))
  result
}

repro_preflight <- function(config, script) {
  refs <- repro_script_paths(script)
  keys <- unique(vapply(refs, `[[`, character(1), "key"))
  unknown <- setdiff(keys, names(config$paths))
  if (length(unknown)) stop("Unknown path keys in script: ", paste(unknown, collapse = ", "))
  spec <- config$spec[match(keys, config$spec$key), , drop = FALSE]
  missing <- character()
  for (i in seq_along(keys)) {
    key <- keys[i]
    path <- config$paths[[key]]
    is_output <- spec$root[i] == "output"
    present <- dir.exists(path)
    label <- if (present) "FOUND" else if (is_output) "WILL CREATE ON RUN" else "MISSING"
    cat(sprintf("[%s] %s = %s\n", label, key, path))
    if (!present && !is_output) missing <- c(missing, key)
  }
  # Validate exact source files, including the stem paths used by source_gitfile().
  path_fun <- repro_path_function(config)
  for (ref in refs) {
    row <- config$spec[match(ref$key, config$spec$key), , drop = FALSE]
    if (!row$root %in% c("source", "repository")) next
    p <- path_fun(ref$key, ref$relative)
    if (!grepl("\\.[Rr]$", p)) p <- paste0(p, ".R")
    if (!file.exists(p)) {
      cat("[MISSING SOURCE] ", p, "\n", sep = "")
      missing <- c(missing, ref$key)
    }
  }
  cat("Check scope: configured roots and directly referenced source files; not all data files, R packages, or manuscript-version compatibility.\n")
  list(ok = !length(missing), missing = unique(missing),
       output_keys = keys[spec$root == "output"])
}
