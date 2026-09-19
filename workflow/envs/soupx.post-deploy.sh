#!env bash
set -o pipefail
# TODO:
# 1) Point to newer tagged version here, once any version contains: https://github.com/scverse/anndataR/pull/480
# 2) Build backported version 1.2.2 for bioconda, once it actually becomes available. See: https://github.com/scverse/anndataR/pull/503
Rscript --no-environ -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="false"); remotes::install_github("scverse/anndataR", ref ="243bfa78c2f6a5b48051f5dc76fdc3409264398b", upgrade = "never")'
