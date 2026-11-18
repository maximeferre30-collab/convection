#!/usr/bin/env bash
# Niveau 3 : lance le cas Code_Saturne 3D de la piscine (sous Linux / WSL).
# Usage : ./lancer.sh <uniforme|bord|coin> [nx=48] [duree_h=6] [nprocs=cœurs physiques]
# Le cas tourne dans ~/cs/etudes/piscine (système de fichiers Linux, bien plus rapide que /mnt/c) ;
# séries temporelles et instantanés CSV sont recopiés dans cfd3d/resultats/<repartition>_nx<nx><SUFFIXE>/.
# Variables optionnelles : SUFFIXE (nom de variante), PISCINE_SORTIE_H, PISCINE_POST_H (voir SRC/).
set -euo pipefail

REP=${1:-uniforme}
NX=${2:-48}
DUREE=${3:-6}
NP=${4:-$(lscpu -p=core | grep -v '^#' | sort -u | wc -l)}   # cœurs physiques (OpenMPI refuse l'hyperthreading)
CS=${CS:-$HOME/cs/arch/bin/code_saturne}
ICI=$(cd "$(dirname "$0")" && pwd)
ETUDE=${ETUDE:-$HOME/cs/etudes/piscine}
NOM=CAS_${REP}_nx${NX}${SUFFIXE:-}

mkdir -p "$ETUDE"
if [ ! -d "$ETUDE/$NOM" ]; then
  (cd "$ETUDE" && "$CS" create -c "$NOM" --noref -q)
fi
cp "$ICI"/SRC/*.cpp "$ETUDE/$NOM/SRC/"

export PISCINE_REPARTITION=$REP PISCINE_NX=$NX PISCINE_DUREE_H=$DUREE
ID=$(date +%Y%m%d-%H%M%S)
(cd "$ETUDE/$NOM/DATA" && "$CS" run -n "$NP" --id "$ID")

RESU=$ETUDE/$NOM/RESU/$ID
DEST=$ICI/resultats/${REP}_nx${NX}${SUFFIXE:-}
rm -rf "$DEST" && mkdir -p "$DEST"
cp "$RESU"/*.csv "$DEST"/
cp "$RESU"/run_solver.log "$DEST"/ 2>/dev/null || true
echo "Résultats copiés dans $DEST"
echo "ParaView : $RESU/postprocessing/RESULTS.case"
