#!/usr/bin/env bash
# Resolve the GeneTribe -s separator (default @, matching pair/submit defaults).
# Prints the separator string on stdout.
set -euo pipefail

sep="@"
if [[ "${GENETRIBE_EXTRA:-}" =~ (^|[[:space:]])-s[[:space:]]+([^[:space:]]+) ]]; then
  sep="${BASH_REMATCH[2]}"
fi
printf '%s\n' "${sep}"
