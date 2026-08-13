#!/usr/bin/env python3
"""Reduce a Claude Code session JSONL to a compact, line-addressable digest.

Used by the `document-session` skill. Does ALL the mechanical extraction so the
agent only has to do interpretation:

  * a compact EVENT STREAM (real user prompts verbatim, truncated assistant prose,
    one-line tool-call summaries) with the original JSONL line number on every line,
    so the agent can drill into the raw transcript at any point;
  * the mechanical WORK PRODUCTS index (files written/edited, commits, commands);
  * the measured FRICTION counts + session SIZE;
  * a ROUGH SPOTS list pointing at flagged lines (errors / interruptions / denials /
    retry loops) for targeted drill-down.

Design rules (from the brainstorm):
  * Reads the on-disk JSONL, so it is immune to in-context compaction.
  * Drops large *successful* tool outputs but KEEPS error outputs (the drill targets).
  * Only does cheap, unambiguous marks — no fuzzy "this went wrong" heuristics; that
    judgement is the agent's job.
  * Never silently truncates a category without saying so.

Usage:  session_digest.py <session.jsonl>
Output: plain text to stdout.
"""

import json
import re
import shlex
import sys
from datetime import datetime
from posixpath import basename

# --- tunables --------------------------------------------------------------
ASST_TEXT_CAP = 350  # chars of assistant prose kept per text block
ERR_OUTPUT_CAP = 600  # chars of an error output kept (tail)
CMD_PREVIEW_CAP = 160  # chars of a bash command shown inline
MAX_COMMANDS_LISTED = 50  # cap on the WORK PRODUCTS command list (excess noted)
MAX_ROUGH_SPOTS = 40  # cap on rough-spot list (excess noted)

# Tool-result content phrases that indicate a *permission denial* rather than a
# genuine runtime error. Best-effort; matched case-insensitively.
DENIAL_PHRASES = (
    "user doesn't want to proceed",
    "the user doesn't want to take this action",
    "tool use was rejected",
    "user has chosen not to",
    "requested permissions",
    "haven't granted",
)


def load_entries(path):
    entries = []
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            try:
                entries.append((lineno, json.loads(line)))
            except json.JSONDecodeError:
                # tolerate a partially-written trailing line
                continue
    return entries


def is_command_wrapper(text):
    """A slash-command invocation wrapper, e.g. <command-name>/brainstorm</…>."""
    t = text.lstrip()
    return t.startswith("<command-name>") or t.startswith("<command-message>") or t.startswith("<local-command-")


def extract_command(text):
    """Pull the slash-command name out of a command wrapper, if present."""
    marker = "<command-name>"
    if marker in text:
        start = text.index(marker) + len(marker)
        end = text.find("</command-name>", start)
        if end != -1:
            return text[start:end].strip()
    return None


def first_text(content_str, cap):
    s = " ".join(content_str.split())
    return s if len(s) <= cap else s[:cap] + "…"


# git's own global options, i.e. those valid *before* the subcommand. The ones
# listed here consume the following token as their value, so the subcommand is
# two positions on, not one.
GIT_GLOBAL_VALUE_FLAGS = {
    "-C",
    "-c",
    "--git-dir",
    "--work-tree",
    "--namespace",
    "--exec-path",
    "--super-prefix",
}


def is_git_commit(piece):
    """True only when `piece` actually invokes `git commit`.

    This used to be `re.search(r"\\bgit\\b.*\\bcommit\\b", piece)`, which matched
    any command containing both words — `git grep -n "do not commit"` was
    reported as a commit that never happened. Parse tokens instead: find the
    `git` invocation, step over git's global options, and require `commit` to
    be the subcommand.
    """
    try:
        tokens = shlex.split(piece)
    except ValueError:  # unbalanced quotes in a partial command
        tokens = piece.split()

    i = 0
    # step over leading VAR=value assignments, then an optional sudo
    while i < len(tokens) and "=" in tokens[i] and not tokens[i].startswith("-"):
        i += 1
    if i < len(tokens) and tokens[i] == "sudo":
        i += 1
    if i >= len(tokens) or basename(tokens[i]) != "git":
        return False

    i += 1
    while i < len(tokens):
        tok = tokens[i]
        if tok in GIT_GLOBAL_VALUE_FLAGS:
            i += 2
        elif tok.startswith("-"):
            i += 1
        else:
            break
    return i < len(tokens) and tokens[i] == "commit"


def tool_summary(name, inp):
    """One-line summary of a tool call, surfacing its key argument."""
    inp = inp or {}
    if name in ("Edit", "MultiEdit", "Read", "Write"):
        return f"{name}({inp.get('file_path', '?')})"
    if name == "NotebookEdit":
        return f"{name}({inp.get('notebook_path', '?')})"
    if name == "Bash":
        cmd = " ".join(str(inp.get("command", "")).split())
        if len(cmd) > CMD_PREVIEW_CAP:
            cmd = cmd[:CMD_PREVIEW_CAP] + "…"
        desc = inp.get("description")
        return f"Bash[{desc}]: {cmd}" if desc else f"Bash: {cmd}"
    if name in ("Grep", "Glob"):
        return f"{name}({inp.get('pattern', inp.get('query', '?'))})"
    if name == "Task":
        return f"Task({inp.get('subagent_type', '')}: {inp.get('description', '')})"
    if name == "Skill":
        return f"Skill({inp.get('skill', '?')})"
    # generic fallback: show a short key=val of the first scalar arg
    for k, v in inp.items():
        if isinstance(v, (str, int, float, bool)):
            sv = str(v)
            if len(sv) > 60:
                sv = sv[:60] + "…"
            return f"{name}({k}={sv})"
    return f"{name}()"


def result_text(content):
    """Flatten a tool_result content (str or list of blocks) to text."""
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts = []
        for b in content:
            if isinstance(b, dict):
                parts.append(b.get("text", ""))
            else:
                parts.append(str(b))
        return "\n".join(parts)
    return str(content)


def fmt_duration(start, end):
    if not start or not end:
        return "unknown"
    try:
        s = datetime.fromisoformat(start.replace("Z", "+00:00"))
        e = datetime.fromisoformat(end.replace("Z", "+00:00"))
        secs = int((e - s).total_seconds())
    except ValueError:
        return "unknown"
    if secs < 0:
        return "unknown"
    h, rem = divmod(secs, 3600)
    m = rem // 60
    return f"{h}h {m}m" if h else f"{m}m"


def main():
    if len(sys.argv) < 2:
        print("usage: session_digest.py <session.jsonl>", file=sys.stderr)
        sys.exit(2)

    entries = load_entries(sys.argv[1])

    # accumulators
    title = None
    session_id = None
    ts_first = ts_last = None
    ts_first_prompt = None  # timestamp of the first real typed prompt
    user_prompts = 0  # typed prompts (excludes meta / commands / tool results)
    asst_msgs = 0  # assistant JSONL entries containing text or tool_use
    tool_calls = 0
    tool_errors = 0
    denials = 0
    interruptions = 0
    files_written = []
    files_edited = []
    notebooks = []
    commands = []  # (lineno, preview)
    commits = []  # (lineno, message)
    rough_spots = []  # (lineno, kind, detail)
    events = []  # (lineno, "USER"/"ASST"/..., text)

    # retry-loop detection: a Bash command fired verbatim twice in a row.
    # Deliberately narrow — only the unambiguous case, to avoid flagging
    # legitimate consecutive edits / task calls as friction.
    prev_bash_cmd = None
    retry_loops = 0

    # map tool_use_id -> (lineno, name) so we can attribute errors to a tool
    toolname_by_id = {}

    for lineno, o in entries:
        t = o.get("type")
        ts = o.get("timestamp")
        if ts:
            ts_first = ts_first or ts
            ts_last = ts

        if t == "ai-title":
            title = o.get("aiTitle")
            continue
        if t in ("mode", "permission-mode"):
            session_id = session_id or o.get("sessionId")
            continue
        if t == "queue-operation":
            # a message typed while the agent was working == an interruption
            if o.get("operation") == "enqueue":
                interruptions += 1
                detail = first_text(str(o.get("content", "")), 80)
                rough_spots.append((lineno, "interruption", detail))
                events.append((lineno, "INTERRUPT", detail))
            continue
        session_id = session_id or o.get("sessionId")

        if t == "user":
            msg = o.get("message", {})
            content = msg.get("content")

            # string content == a typed prompt or a command wrapper
            if isinstance(content, str):
                if o.get("isMeta"):
                    continue
                cmd = extract_command(content) if is_command_wrapper(content) else None
                if cmd:
                    # slash-command: a user action, but not counted as a typed turn
                    events.append((lineno, "CMD", f"{cmd} (slash-command)"))
                elif is_command_wrapper(content):
                    continue  # other wrapper noise
                else:
                    events.append((lineno, "USER", first_text(content, 500)))
                    user_prompts += 1
                    if ts_first_prompt is None:
                        ts_first_prompt = ts
                continue

            # list content == tool results (and occasional injected text blocks)
            if isinstance(content, list):
                for b in content:
                    if not isinstance(b, dict):
                        continue
                    if b.get("type") == "tool_result":
                        txt = result_text(b.get("content"))
                        if b.get("is_error"):
                            low = txt.lower()
                            if any(p in low for p in DENIAL_PHRASES):
                                denials += 1
                                kind = "permission-denied"
                            else:
                                tool_errors += 1
                                kind = "tool-error"
                            tid = b.get("tool_use_id")
                            who = toolname_by_id.get(tid, "?")
                            tail = first_text(txt, ERR_OUTPUT_CAP)
                            rough_spots.append((lineno, kind, f"{who}: {tail}"))
                            events.append((lineno, "ERROR", f"{who}: {tail}"))
                    # injected text blocks in user turns are skill/system noise → skip
                continue

        if t == "assistant":
            msg = o.get("message", {})
            had_content = False
            saw_bash = False
            for b in msg.get("content", []):
                if not isinstance(b, dict):
                    continue
                bt = b.get("type")
                if bt == "text":
                    txt = b.get("text", "").strip()
                    if txt:
                        events.append((lineno, "ASST", first_text(txt, ASST_TEXT_CAP)))
                        had_content = True
                elif bt == "tool_use":
                    name = b.get("name", "?")
                    inp = b.get("input", {}) or {}
                    tool_calls += 1
                    had_content = True
                    toolname_by_id[b.get("id")] = name
                    events.append((lineno, "TOOL", tool_summary(name, inp)))

                    # work-product extraction
                    if name == "Write":
                        files_written.append(inp.get("file_path", "?"))
                    elif name in ("Edit", "MultiEdit"):
                        files_edited.append(inp.get("file_path", "?"))
                    elif name == "NotebookEdit":
                        notebooks.append(inp.get("notebook_path", "?"))
                    elif name == "Bash":
                        cmd = str(inp.get("command", ""))
                        commands.append((lineno, first_text(cmd, CMD_PREVIEW_CAP)))
                        saw_bash = True
                        # commit detection — handles `git -C <path> commit`, extra
                        # global flags, and multi-`-m` messages; skips --amend.
                        for piece in re.split(r"&&|;", cmd):
                            if is_git_commit(piece) and "--amend" not in piece:
                                m = ""
                                mo = re.search(r"-m\s+(['\"])(.*?)\1", piece, re.S)
                                if mo:
                                    m = " ".join(mo.group(2).split())[:80]
                                elif "-m" in piece:
                                    m = " ".join(piece.split("-m", 1)[1].split())[:80]
                                commits.append((lineno, m or "(message not parsed)"))
                        # retry loop: the exact same command fired twice in a row
                        norm = " ".join(cmd.split())
                        if norm and norm == prev_bash_cmd:
                            retry_loops += 1
                            rough_spots.append(
                                (lineno, "retry-loop", f"identical Bash repeated: {first_text(cmd, 80)}")
                            )
                        prev_bash_cmd = norm
            # a non-Bash assistant entry breaks a Bash retry streak
            if had_content and not saw_bash:
                prev_bash_cmd = None
            if had_content:
                asst_msgs += 1

    # ---- emit -------------------------------------------------------------
    out = []
    w = out.append

    start = ts_first_prompt or ts_first
    w("# SESSION DIGEST")
    w(f"session_id : {session_id}")
    w(f"title      : {title!r}")
    w(f"started    : {start} (first typed prompt)")
    w(f"ended      : {ts_last}")
    w(
        f"wall_clock : {fmt_duration(start, ts_last)} (first prompt → last entry; "
        "still includes mid-session idle — tool_calls is the truer activity signal)"
    )
    w(
        f"turns      : {user_prompts} typed prompts / {asst_msgs} assistant msgs "
        "(asst msgs are raw JSONL entries, not conversational turns)"
    )
    w(f"tool_calls : {tool_calls}")
    w("")
    w("## FRICTION (measured — exact counts, no judgement)")
    w(f"tool_errors        : {tool_errors}")
    w(f"interruptions      : {interruptions}")
    w(f"permission_denials : {denials}")
    w(f"retry_loops        : {retry_loops}")
    rate = f"{tool_errors}/{tool_calls}" if tool_calls else "0/0"
    pct = f"{(100 * tool_errors / tool_calls):.1f}%" if tool_calls else "n/a"
    w(f"error_rate         : {rate} ({pct} of tool calls)")
    w("")

    def dedup_counts(paths):
        """Ordered unique paths with occurrence counts."""
        order, counts = [], {}
        for p in paths:
            if p not in counts:
                order.append(p)
                counts[p] = 0
            counts[p] += 1
        return order, counts

    w("## WORK PRODUCTS (mechanical)")
    wo, wc = dedup_counts(files_written)
    w(f"files_written ({len(wo)} distinct):")
    for p in wo:
        w(f"  {p}" + (f"  (written ×{wc[p]})" if wc[p] > 1 else ""))
    eo, ec = dedup_counts(files_edited)
    w(f"files_edited ({len(eo)} distinct):")
    for p in eo:
        w(f"  {p}  (×{ec[p]} edit{'s' if ec[p] > 1 else ''})")
    if notebooks:
        no, nc = dedup_counts(notebooks)
        w(f"notebooks_edited ({len(no)} distinct):")
        for p in no:
            w(f"  {p}  (×{nc[p]})")
    w(f"commits ({len(commits)}):")
    for ln, m in commits:
        w(f"  L{ln}: {m}")
    w(f"commands_run ({len(commands)}):")
    shown = commands[:MAX_COMMANDS_LISTED]
    for ln, c in shown:
        w(f"  L{ln}: {c}")
    if len(commands) > MAX_COMMANDS_LISTED:
        w(
            f"  … {len(commands) - MAX_COMMANDS_LISTED} more commands not listed (cap "
            f"{MAX_COMMANDS_LISTED}) — read the raw transcript if needed"
        )
    w("")
    w("## ROUGH SPOTS (flagged for drill-down — open the raw JSONL at these lines)")
    if not rough_spots:
        w("  (none — no errors, interruptions, denials, or retry loops detected)")
    for ln, kind, detail in rough_spots[:MAX_ROUGH_SPOTS]:
        w(f"  L{ln} [{kind}] {detail}")
    if len(rough_spots) > MAX_ROUGH_SPOTS:
        w(f"  … {len(rough_spots) - MAX_ROUGH_SPOTS} more not listed (cap {MAX_ROUGH_SPOTS})")
    w("")
    w("## EVENT STREAM (LN = line in the raw JSONL; assistant prose truncated)")
    for ln, kind, text in events:
        w(f"L{ln:<5} {kind:<9} {text}")

    print("\n".join(out))


if __name__ == "__main__":
    main()
