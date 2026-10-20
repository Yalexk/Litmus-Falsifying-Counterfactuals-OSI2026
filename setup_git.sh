#!/usr/bin/env bash
# Repo-local git identity. Every commit here is authored and committed by you.
# Nothing else writes commits; no Co-Authored-By trailers are added anywhere.
set -euo pipefail
read -rp "Your git name:  " NAME
read -rp "Your git email: " EMAIL
git init -q 2>/dev/null || true
git config user.name  "$NAME"
git config user.email "$EMAIL"
git config commit.gpgsign false
echo "Local identity set: $NAME <$EMAIL>"
echo "Verify with: git config --local --list | grep user"
