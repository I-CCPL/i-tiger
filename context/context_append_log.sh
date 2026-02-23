#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  context/context_append_log.sh [--plan|--decision] <title> <message>

Examples:
  context/context_append_log.sh "parser" "chk reader guard added"
  context/context_append_log.sh --plan "next" "start ham_q_builder implementation"
USAGE
}

if [ $# -lt 2 ]; then
  usage
  exit 1
fi

mode="decision"
if [[ "$1" == "--plan" ]]; then
  mode="plan"
  shift
elif [[ "$1" == "--decision" ]]; then
  mode="decision"
  shift
fi

if [ $# -lt 2 ]; then
  usage
  exit 1
fi

title="$1"
shift
message="$*"

timestamp="$(date '+%Y-%m-%d %H:%M:%S %z')"
line="- [${timestamp}] ${title}: ${message}"

if [[ "$mode" == "plan" ]]; then
  file="context/README.md"
  header="# i-tiger Context (Minimal)"
  section="## Plan Updates"
else
  file="context/decision-log.md"
  header="# Decision Log"
  section="## Session Updates"
fi

if [ ! -f "$file" ]; then
  printf "%s\n\n" "$header" > "$file"
fi

if ! grep -q "^${section}$" "$file"; then
  printf "\n%s\n" "$section" >> "$file"
fi

printf "%s\n" "$line" >> "$file"
echo "Appended to ${file}"
