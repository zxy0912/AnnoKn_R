# Copy to config.local.R (git-ignored). Relative paths resolve beside that file.
# Active source dependencies are bundled; configure external datasets and outputs here.
list(
  data_root = "data",
  output_root = "outputs",
  paths = list(
    # Override a dataset independently when it is stored outside data_root:
    # gtex_expression = "/your/data/GTEx_Gene_expression_adjust_covariates",
    # gtex_plink = "/your/data/GTEx/ancestry",
    # ldsc_annotations = "/your/data/baselineLD_v2.2",
    # height_validation_workspace = "/your/writable/height/workspace",
    # scz_validation_workspace = "/your/writable/scz/workspace",
    # gtex_validation_workspace = "/your/writable/gtex/workspace"
  )
)
