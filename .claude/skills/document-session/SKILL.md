---
name: document-session
description: Capture the current Claude Code session — save its raw JSONL transcript and generate a critical, work-product-focused summary under documented_sessions/. Run before /clear, or when asked to document, summarise, or retrospect on the current session.
---

This skill documents the **current session** into two durable artefacts under a
`documented_sessions/` directory at the repository root (gitignored, but persistent
on disk):

- `documented_sessions/transcripts/<name>.jsonl` — the raw session transcript.
- `documented_sessions/summaries/<name>.md` — a readable summary serving two
  purposes: a **work-product index** (what the session changed) and a **critical
  retrospective** (what worked, what went sideways, what to do differently).

A rolling `documented_sessions/INDEX.md` gets one appended row per run.

## This skill is NON-INTERACTIVE

It is meant to run unattended right before `/clear`. **Do not ask the user
questions.** Resolve everything from the transcript and produce the artefacts. The
only acceptable reason to stop early is a hard failure (no transcript found).

## How the summary is made

A digest script does all the **mechanical** extraction (work products, measured
friction, an event stream with line numbers). You do the **interpretation** (goal,
category, TL;DR, retrospective), reading the on-disk transcript — which is the full
record even if your in-context memory was compacted. You only open the *raw*
transcript at the specific lines the digest flags as rough spots.

## Procedure

Run these steps in order. All paths are relative to the repository root (the
directory you launched Claude in). This skill is self-contained — it carries its own
`get-session.sh` and `session_digest.py` under `.claude/skills/document-session/`.

### 1. Resolve and copy the transcript

```bash
mkdir -p documented_sessions/transcripts documented_sessions/summaries
RAW="$(bash .claude/skills/document-session/get-session.sh --path)"   # most-recent jsonl = current session
echo "resolved transcript: $RAW"
```

If that errors, stop and report — there is nothing to document.

### 2. Generate the digest

```bash
python3 .claude/skills/document-session/session_digest.py "$RAW"
```

Read its output in full. It gives you: the header (session id, **started** time,
duration, tool-call count, typed-prompt count), the **measured friction** counts +
error-rate, the **work-products** index, a **rough-spots** list (lines to drill into),
and the **event stream**. The digest is ephemeral — do **not** save it.

### 3. Decide the filename

Build `NAME = YYMMDD-HHMM-<primary>-<slug>` where:
- `YYMMDD-HHMM` comes from the digest's **started** field (the first typed prompt).
- `<primary>` is the primary category (see Categorisation).
- `<slug>` is a 3–5 word lowercase kebab-case description of the session
  (e.g. `document-session-skill`). The digest's `title` is a decent starting point.

Both the transcript and the summary use this same `NAME`.

### 4. Trivial-session short-circuit

If the session did **near-nothing** — `files_written + files_edited == 0` **and**
`commits == 0` **and** `tool_errors + interruptions + denials == 0` **and**
`tool_calls <= 3` — don't manufacture a retrospective. Still save the transcript and
write a **stub** summary (header block + a one-line TL;DR + the line
`_Trivial session — no substantive work products; full retrospective skipped._`),
and still append the index row. Then skip to step 7.

### 5. Drill into rough spots

For each rough spot in the digest that you judge worth understanding (genuine
pushback, a real error, a retry loop — not every trivial blip), open the raw
transcript at that line to recover the detail the digest dropped, e.g.:

```bash
sed -n '60,62p' "$RAW" | python3 -m json.tool   # inspect entries around line 60
```

Interruptions and corrections are where the honest retrospective lives — read what
the user actually said and why.

### 6. Write the summary

Before writing, scan the **user turns** in the event stream and assess the
**course-corrections** count (see the template) — this is your judgement, not
something the digest computes. Then write `documented_sessions/summaries/<NAME>.md`
using the template below.

### 7. Copy the transcript and append the index

```bash
cp "$RAW" "documented_sessions/transcripts/<NAME>.jsonl"
```

Append a row to `documented_sessions/INDEX.md` (create it with the header first if it
doesn't exist — see Index format).

### 8. Report

Tell the user the two file paths and the one-line TL;DR. Keep it brief.

## Summary template

```markdown
# Session Summary: <title>

- **Date:** YYYY-MM-DD HH:MM
- **Session:** <session-id>
- **Category:** <primary>[ + <also>, …]
- **Size:** <N> tool calls · <M> typed prompts · <duration>
- **Mechanical friction (measured):** <X> tool errors · <Y> interruptions · <Z> denials · <W> retry loops — error rate <X>/<N> _(tool-level only; does not capture reasoning corrections; and `interruptions` is a raw count that over-reads as friction — decompose it per the reconcile note below)_
- **Course-corrections (agent-assessed):** <C> — times the user corrected a wrong assumption, rejected or redirected an output, or told you to change approach. This is the honest self-improvement signal; mechanical friction routinely misses it. Enumerate them in "what went sideways" so the count is auditable.

## Goal
<What the session was actually trying to achieve, reconstructed across ALL user
turns — not just the first prompt. One clear-cycle can be multi-phase (e.g. a
brainstorm that turned into implementation); note any pivots explicitly.>

## TL;DR
<2–3 sentences: what it set out to do and what came of it.>

## Work Products
**Files written** — <list, or "none">
**Files edited** — <list with (×N edits), or "none">
**Commits** — <list, or "none">
**Notable commands** — <only the ones worth remembering: migrations, scripts run,
  tests, git ops. Skip read-only probes. Or "none">

## Retrospective
**What worked** — <short, concrete bullets.>

**What went sideways** — <THE HEART of this section; never omit it. The landing
  zone for the drill-down findings: pushback, corrections, fumbles, dead-ends,
  errors, retry loops — each with a line on the *why*. **Enumerate each
  course-correction here** (so the header count is auditable). **Reconcile the two
  friction figures**: if mechanical friction is low but course-corrections are not,
  say so and explain — a clean tool record alongside repeated reasoning corrections
  is the common case, and the one most worth learning from. **Decompose
  `interruptions`**: the raw count lumps together (a) genuine redirects/failures,
  (b) background-task-completion notifications, and (c) solicited user-steering
  (new requests/context you invited) — only (a) is friction. Break the figure down
  and discount (b)/(c) explicitly rather than letting the raw number imply a rough
  session; for long-running / background-job sessions the count is routinely
  dominated by (b). If the session was genuinely clean on both, say exactly that in
  one line.>

**Would do differently** — <forward-looking, concrete lessons.>
```

## Tone & style

- **Flat, logical, declarative.** Plain prose. No hype, no marketing adjectives, no
  "successfully", no self-congratulation.
- **Lean critical.** Be blunt about the agent's own fumbles and where the user had
  to correct course. The **course-corrections count is the primary honesty anchor** —
  a session with 0 mechanical errors but several course-corrections did **not** "go
  smoothly," and saying so is the point. Never write "went smoothly" over a non-zero
  course-correction or mechanical-friction count without explaining it.
- Judge **delivered-vs-asked** against the reconstructed goal, not against the
  agent's own framing of what it set out to do.

## Categorisation

**Categorise by the session's output / purpose, not by the activities inside it.**
A session that reads a lot of code *in order to* produce a plan is `planning`, not
`analysis` or `investigation` — the reading served the plan. Tag what the session was
*for*, not what it did along the way.

Vocabulary (resist adding to it until a real session genuinely doesn't fit):

- `planning` — designing or scoping work to be done (brainstorms, scopes, designs).
- `implementation` — building or changing code/tooling as the deliverable.
- `investigation` — figuring out how something works or diagnosing, as the *primary
  aim* (not exploration done in service of a plan or fix).
- `fixes` — correcting defects in existing code.
- `analysis` — examining data or results to draw conclusions; the output is
  *findings*, not a plan or code.
- `other` — escape hatch.

**Default to a single primary tag.** Add a secondary ("also") tag only when a
*distinct, substantial* body of work with its own output happened alongside the
primary (e.g. a session that both implemented a feature *and* fixed unrelated bugs).
Do not add a secondary tag for activity that merely served the primary purpose. The
primary is the index grouping key.

## Index format

`documented_sessions/INDEX.md` is a markdown table, one row per documented session,
appended chronologically. Create it with this header if it doesn't exist:

```markdown
# Documented Sessions

| Date | Primary | Title | TL;DR | Size · friction | Links |
|---|---|---|---|---|---|
```

Each row shows size + friction compactly in the one column. Format
`<N> calls · ⚠<X>e/<Y>i · ✎<C>` where `⚠Xe/Yi` is mechanical (X tool errors /
Y interruptions) and `✎C` is the agent-assessed course-corrections count — the
trackable self-improvement signal. Example `44 calls · ⚠0e/3i · ✎1`:

```markdown
| 2026-06-23 22:16 | implementation | document-session skill | One-line gist. | 44 calls · ⚠0e/3i · ✎1 | [t](transcripts/<NAME>.jsonl) · [s](summaries/<NAME>.md) |
```
