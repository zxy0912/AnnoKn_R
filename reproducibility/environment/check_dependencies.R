# Run from the repository root. This checks the installed library without changing it.
args <- commandArgs(trailingOnly = TRUE)
manifest <- if (length(args)) args[1] else "reproducibility/environment/r_dependencies.tsv"
spec <- read.delim(manifest, colClasses = "character", na.strings = "", check.names = FALSE)
rows <- lapply(seq_len(nrow(spec)), function(i) {
  p <- spec$package[i]
  message <- ""
  loaded <- tryCatch(requireNamespace(p, quietly = TRUE), error = function(e) {
    message <<- conditionMessage(e)
    FALSE
  })
  version <- if (loaded) as.character(packageVersion(p)) else NA_character_
  minimum <- spec$minimum_version[i]
  meets_minimum <- loaded && (is.na(minimum) || packageVersion(p) >= package_version(minimum))
  matches_checked <- loaded && packageVersion(p) == package_version(spec$version[i])
  data.frame(package = p, installed_version = version, checked_version = spec$version[i],
             loadable = loaded, meets_minimum = meets_minimum,
             matches_checked_version = matches_checked, message = message)
})
result <- do.call(rbind, rows)
write.table(result, stdout(), sep = "\t", quote = TRUE, row.names = FALSE)
if (!all(result$loadable & result$meets_minimum)) {
  cat("FAIL: install missing dependencies or resolve load/version errors.\n", file = stderr())
  quit(status = 1L)
}
cat("PASS:", nrow(result), "packages load and meet declared minimum versions.\n", file = stderr())
if (!all(result$matches_checked_version))
  cat("NOTE: some versions differ from the checked environment; see the table.\n", file = stderr())
