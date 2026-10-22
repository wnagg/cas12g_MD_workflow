#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG="${1:-$ROOT/env/Cas12g_3rep_400ns.env.example}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"

resolve_exe() {
  local value="$1"
  if [[ "$value" == */* ]]; then
    [[ -x "$value" ]] || { echo "ERROR: executable not found: $value" >&2; return 1; }
    printf '%s\n' "$value"
  else
    command -v "$value" || { echo "ERROR: command not found: $value" >&2; return 1; }
  fi
}
GMX_EXE="$(resolve_exe "${GMX:-gmx}")"
PYMOL_EXE="$(resolve_exe "${PYMOL:-pymol}")"
PYTHON_EXE="$(resolve_exe "${PYTHON:-python3}")"
WT_PATH="${WT_PDB:-}"
[[ -n "$WT_PATH" ]] || { echo "ERROR: set WT_PDB in the config." >&2; exit 1; }
[[ "$WT_PATH" = /* ]] || WT_PATH="$ROOT/$WT_PATH"
[[ -s "$WT_PATH" ]] || { echo "ERROR: WT PDB is missing or empty: $WT_PATH" >&2; exit 1; }

echo "=== Cas12g MD preflight ==="
echo "Repository: $ROOT"
echo "WT PDB: $WT_PATH"
echo "WORKDIR: ${WORKDIR:-unset}"
"$GMX_EXE" --version | head -n 20
"$PYMOL_EXE" -cq -d "print('PyMOL CLI OK')" -d quit 2>&1 | tail -n 10
"$PYTHON_EXE" - <<'PY'
import importlib
for name in ("numpy", "pandas", "matplotlib"):
    importlib.import_module(name)
    print("[OK] Python module:", name)
PY
"$PYTHON_EXE" - "$WT_PATH" "${CHAIN:-A}" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1])
chain = sys.argv[2]
ca = {}
hetatm = 0
for line in p.read_text(errors="replace").splitlines():
    if line.startswith("HETATM"):
        hetatm += 1
    if not line.startswith("ATOM  ") or line[12:16].strip() != "CA" or len(line) < 27:
        continue
    if (line[21].strip() or "_") != (chain.strip() or "_"):
        continue
    if line[26].strip():
        continue
    try:
        ca[int(line[22:26])] = line[17:20].strip().upper()
    except ValueError:
        pass
print(f"[OK] chain {chain}: {len(ca)} C-alpha residues; residue 623 = {ca.get(623)}")
if ca.get(623) != "VAL":
    raise SystemExit("ERROR: expected VAL at chain A residue 623. Verify numbering/chain before mutation.")
if len(ca) < 100:
    raise SystemExit("ERROR: unexpectedly few C-alpha atoms; check the input structure.")
if hetatm:
    print(f"[WARNING] found {hetatm} HETATM records; protein-only preparation may remove biologically relevant components.")
PY
mkdir -p "${WORKDIR:-$ROOT/work}"
df -h "${WORKDIR:-$ROOT/work}" || true
if command -v nvidia-smi >/dev/null 2>&1; then nvidia-smi; else echo "[WARNING] nvidia-smi not found; GPU availability not verified."; fi
echo "[OK] Preflight checks passed. Run --build-only, then --pilot-only before production."
