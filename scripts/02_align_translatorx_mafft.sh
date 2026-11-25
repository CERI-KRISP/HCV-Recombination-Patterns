#!/usr/bin/env bash
#
# Codon-aware alignment pipeline: MAFFT (amino-acid alignment) + back-translation
# via TranslatorX, for HCV query sequences + HCV-GLUE / LANL references.
#
# Replaces MACSE. Addresses Reviewer 1 comment #2 (align AA, translate back to nt)
# and avoids the internal-stop-codon problem caused by frame-shifted NT alignment.
#
# ---------------------------------------------------------------------------
# REQUIREMENTS (install once):
#   - MAFFT           (conda install -c bioconda mafft)  -- optional if you
#                        use TranslatorX's bundled copy instead
#   - TranslatorX      (https://github.com/nathanweeks/translatorx or the
#                        standalone perl script: translatorx_vLocal.pl)
#   - seqkit (optional, useful for trimming/QC)  conda install -c bioconda seqkit
#
# USAGE:
#   TRANSLATORX_PATH=/full/path/to/translatorx_vLocal.pl \
#     ./hcv_codon_aware_alignment.sh input_cds.fasta output_prefix
#
#   Or place this script in the SAME directory as translatorx_vLocal.pl
#   (and its bundled mafft-*-with-extensions folder) and just run:
#     ./hcv_codon_aware_alignment.sh input_cds.fasta output_prefix
#
#   input_cds.fasta  = HCV queries + reference sequences, ALREADY TRIMMED
#                       to start at the polyprotein ATG and end at the stop
#                       codon (no 5'/3' UTR). See 01_extract_filter_orfs.py.
#
# ---------------------------------------------------------------------------
# IMPORTANT -- TranslatorX invocation gotcha:
#   TranslatorX locates its own bundled MAFFT binary by resolving its own
#   install directory from how it was invoked (i.e. from $0). If you call
#   it with a bare filename and no directory component (e.g. just
#   "perl translatorx_vLocal.pl" from inside its own folder), that path
#   resolution silently breaks: it treats the string "translatorx_vLocal.pl"
#   itself as a directory, so its computed mafft path becomes something like
#     translatorx_vLocal.pl/../mafft-7.222-with-extensions/bin/mafft
#   which the shell rejects with "Not a directory". MAFFT then never runs,
#   the amino-acid alignment file comes out empty, and TranslatorX reports
#   EVERY sequence (queries and references alike) as "not found in amino
#   acid file", eventually crashing with "Illegal division by zero". This
#   looks like a data problem but isn't -- it's purely an invocation issue.
#   This script always calls TranslatorX via an absolute, realpath-resolved
#   path to avoid it.
# ---------------------------------------------------------------------------

set -euo pipefail

export PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:$PATH"

IN_FASTA="${1:?Usage: $0 input_cds.fasta output_prefix}"
PREFIX="${2:?Usage: $0 input_cds.fasta output_prefix}"

TRANSLATORX_PATH="${TRANSLATORX_PATH:-./translatorx_vLocal.pl}"
TRANSLATORX_PATH="$(realpath "$TRANSLATORX_PATH")"

if [ ! -f "$TRANSLATORX_PATH" ]; then
  echo "ERROR: translatorx_vLocal.pl not found at: $TRANSLATORX_PATH"
  echo "Set TRANSLATORX_PATH=/full/path/to/translatorx_vLocal.pl and re-run,"
  echo "or place this script in the same directory as translatorx_vLocal.pl."
  exit 1
fi

TX_DIR="$(dirname "$TRANSLATORX_PATH")"
if ! ls -d "$TX_DIR"/mafft-*-with-extensions >/dev/null 2>&1; then
  echo "WARNING: no bundled mafft-*-with-extensions folder found next to"
  echo "  $TRANSLATORX_PATH"
  echo "TranslatorX expects its bundled MAFFT alongside the script. If you"
  echo "unzipped TranslatorX and moved only the .pl file, re-unzip so the"
  echo "aligner subfolders stay next to translatorx_vLocal.pl."
fi

echo "Using TranslatorX at: $TRANSLATORX_PATH"

echo "== Step 0: sanity check input =="
# Every sequence length should be a multiple of 3 if UTRs are already trimmed.
# This just flags sequences to look at manually -- it does not remove anything.
if command -v seqkit >/dev/null 2>&1; then
  seqkit fx2tab -nl "$IN_FASTA" | awk '{if ($2 % 3 != 0) print $1"\tlength="$2" (not a multiple of 3 -- check UTR trimming/indels)"}' \
    > "${PREFIX}_length_warnings.txt" || true
  if [ -s "${PREFIX}_length_warnings.txt" ]; then
    echo "WARNING: some sequences are not multiples of 3 -- see ${PREFIX}_length_warnings.txt"
    echo "These likely still have UTR bases attached, or contain non-triplet indels."
  fi
else
  echo "(seqkit not found -- skipping length pre-check, not required)"
fi

# Duplicate-ID check -- a common silent cause of downstream weirdness
DUPES=$(grep ">" "$IN_FASTA" | sort | uniq -d || true)
if [ -n "$DUPES" ]; then
  echo "WARNING: duplicate sequence IDs found in $IN_FASTA:"
  echo "$DUPES"
  echo "Fix these before proceeding -- duplicate headers can silently break output parsing."
fi

echo "== Step 1: run TranslatorX (translates -> aligns AA in MAFFT -> back-translates to NT) =="
# -i  input nt fasta
# -o  output prefix
# -p  alignment program: F = MAFFT
# -c  genetic code: 1 = standard genetic code
#     (NOTE: some TranslatorX versions use -t for genetic code; the version
#     matching this script's -h output uses -c for genetic code and -t for
#     a separate F/T "guess reading frame" flag. Check `perl "$TRANSLATORX_PATH" -h`
#     if you're on a different version and this errors with an unrecognized option.)
perl "$TRANSLATORX_PATH" \
    -i "$IN_FASTA" \
    -o "$PREFIX" \
    -p F \
    -c 1

echo "== Step 2: outputs =="
echo "  ${PREFIX}.aaseqs.fasta        -> unaligned AA translations (check for real internal stops here)"
echo "  ${PREFIX}.aa_ali.fasta        -> AA alignment (MAFFT)"
echo "  ${PREFIX}.nt_ali.fasta        -> final codon-aware NT alignment -- USE THIS for RDP5 / recombination analysis"
echo "  ${PREFIX}.html                -> colour-coded alignment report, open in browser"
echo "  ${PREFIX}_summary.txt         -> lists any sequences excluded due to real internal stop codons/frameshifts"

echo "== Step 3: check for sequences TranslatorX flagged/excluded =="
if [ -f "${PREFIX}_summary.txt" ]; then
  grep -i -E "stop|exclud|remov|warn" "${PREFIX}_summary.txt" || echo "  No flagged sequences found."
fi

echo ""
echo "If a REFERENCE sequence gets flagged, do not delete it:"
echo "  1. Locate the offending codon position in ${PREFIX}.aaseqs.fasta"
echo "  2. Mask just that codon in the ORIGINAL nt fasta (replace 3 nt with NNN)"
echo "  3. Re-run this script -- the rest of that reference's sequence is preserved"
echo "     for parent/recombination typing, only the bad codon is neutralised."
