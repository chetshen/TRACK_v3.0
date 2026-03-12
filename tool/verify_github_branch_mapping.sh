#!/usr/bin/env bash
set -euo pipefail

# Verify whether a local commit/branch matches a GitHub remote branch tip,
# and print reconciliation commands to align local and remote.

if [[ $# -lt 2 || $# -gt 3 ]]; then
  cat <<'USAGE'
Usage:
  tool/verify_github_branch_mapping.sh <github_repo_url> <remote_branch> [local_ref]

Examples:
  tool/verify_github_branch_mapping.sh git@github.com:ORG/REPO.git codex/add-function-to-build-bridge-model
  tool/verify_github_branch_mapping.sh https://github.com/ORG/REPO.git codex/add-function-to-build-bridge-model work
USAGE
  exit 2
fi

repo_url="$1"
remote_branch="$2"
local_ref="${3:-HEAD}"

if ! git rev-parse --verify "$local_ref" >/dev/null 2>&1; then
  echo "[ERROR] Local ref '$local_ref' does not exist." >&2
  exit 1
fi

local_sha="$(git rev-parse "$local_ref")"
remote_line="$(git ls-remote --heads "$repo_url" "$remote_branch" || true)"

if [[ -z "$remote_line" ]]; then
  echo "[ERROR] Remote branch '$remote_branch' not found at '$repo_url'."
  echo
  echo "Check branch names with:"
  echo "  git ls-remote --heads $repo_url"
  exit 1
fi

remote_sha="$(awk '{print $1}' <<<"$remote_line")"

printf 'Local ref:    %s\n' "$local_ref"
printf 'Local SHA:    %s\n' "$local_sha"
printf 'Remote branch:%s\n' " $remote_branch"
printf 'Remote SHA:   %s\n' "$remote_sha"
printf '\n'

if [[ "$local_sha" == "$remote_sha" ]]; then
  echo "Status: MATCH (local ref equals remote branch tip)."
  exit 0
fi

if git merge-base --is-ancestor "$remote_sha" "$local_sha" 2>/dev/null; then
  echo "Status: LOCAL AHEAD (local contains remote tip)."
  echo
  echo "If desired, push local ref to update remote branch:"
  echo "  git push $repo_url $local_ref:refs/heads/$remote_branch"
  exit 0
fi

if git merge-base --is-ancestor "$local_sha" "$remote_sha" 2>/dev/null; then
  echo "Status: LOCAL BEHIND (remote contains local ref)."
  echo
  echo "To align local branch/ref to remote:"
  echo "  git fetch $repo_url $remote_branch"
  echo "  git checkout -B ${local_ref} FETCH_HEAD"
  exit 0
fi

echo "Status: DIVERGED (local and remote have different histories)."
echo
echo "Recommended reconciliation options:"
echo "  Option A (rebase local on remote):"
echo "    git fetch $repo_url $remote_branch"
echo "    git checkout ${local_ref}"
echo "    git rebase FETCH_HEAD"
echo
echo "  Option B (merge remote into local):"
echo "    git fetch $repo_url $remote_branch"
echo "    git checkout ${local_ref}"
echo "    git merge FETCH_HEAD"
echo
echo "  Option C (force remote to local, use with caution):"
echo "    git push --force-with-lease $repo_url $local_ref:refs/heads/$remote_branch"
