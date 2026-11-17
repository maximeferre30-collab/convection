#!/usr/bin/env bash
# Campagne du niveau 3 : trois répartitions des assemblages chauds (6 h),
# puis convergence de maillage sur le cas uniforme (2 h). ~3 h sur 4 cœurs.
set -euo pipefail
ICI=$(cd "$(dirname "$0")" && pwd)
for rep in uniforme bord coin; do
  "$ICI/lancer.sh" "$rep" 48 6
done
"$ICI/lancer.sh" uniforme 32 2
"$ICI/lancer.sh" uniforme 64 2
