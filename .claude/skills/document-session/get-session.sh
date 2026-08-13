#!/bin/bash
# Get current session JSONL content for retrospective analysis.
#
# Usage:
#   get-session.sh [session-id]            # cat the session JSONL to stdout
#   get-session.sh --path [session-id]     # print the resolved file path instead
#
# Session resolution order:
#   1. explicit session-id argument
#   2. $CLAUDE_SESSION_ID env var
#   3. fallback: the most-recently-modified *.jsonl under ~/.claude/projects/
#      (when run mid-session this is the current session, since its transcript
#      is actively being written). Used by the document-session skill.

set -e

PATH_ONLY=0
if [ "$1" = "--path" ]; then
  PATH_ONLY=1
  shift
fi

SESSION_ID="${1:-$CLAUDE_SESSION_ID}"
PROJECTS_DIR="$HOME/.claude/projects"

if [ -n "$SESSION_ID" ]; then
  # Find the named session file (could be in any project subdirectory).
  SESSION_FILE=$(find "$PROJECTS_DIR" -name "${SESSION_ID}.jsonl" -type f 2>/dev/null | head -1)
else
  # No id and no env var: fall back to the most recently modified transcript.
  SESSION_FILE=$(find "$PROJECTS_DIR" -name '*.jsonl' -type f -printf '%T@ %p\n' 2>/dev/null \
    | sort -rn | head -1 | cut -d' ' -f2-)
fi

if [ -z "$SESSION_FILE" ] || [ ! -f "$SESSION_FILE" ]; then
  if [ -n "$SESSION_ID" ]; then
    echo "Error: Session file not found for ID: $SESSION_ID" >&2
  else
    echo "Error: No session id given and no session transcripts found under $PROJECTS_DIR" >&2
  fi
  exit 1
fi

if [ "$PATH_ONLY" = "1" ]; then
  echo "$SESSION_FILE"
else
  cat "$SESSION_FILE"
fi
