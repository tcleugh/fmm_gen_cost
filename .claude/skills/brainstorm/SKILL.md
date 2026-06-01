---
name: brainstorm
description: Interactive brainstorming for technical problems. Explores goals, options, risks, and pitfalls in a back-and-forth conversation, producing a markdown record of ideas kept and discarded.
---

This skill guides an interactive brainstorm session for a technical problem. The output is a markdown record — not a transcript, but a structured summary of ideas explored, reasoning used, and decisions made or deferred.

Inform the user that the brainstorm skill is being used and ask them to describe the problem they want to explore.

# Working style

- **One question at a time.** Never fire multiple questions in a single turn. Pick the most important thing to clarify and ask that.
- **Propose, don't just prompt.** Offer concrete options or framings for the user to react to. Open-ended questions stall sessions; tangible ideas move them forward.
- **Keep it moving.** If the conversation stalls, summarise where you are and suggest a next direction.
- **Surface risks early.** When an idea gets traction, proactively raise its main risk or pitfall before the user fully commits to it.
- **Capture as you go.** After each exchange where something meaningful is decided or discarded, update the markdown record. Write the file at the end of each significant turn — not just at the end of the session.

# Session structure

A brainstorm has four natural phases. Don't announce them or force transitions — let the conversation find its own shape.

## 1. Frame the problem

Understand what the user is actually trying to achieve before touching solutions. Ask:
- What does success look like?
- What constraints exist (time, tech stack, team size, reversibility, etc.)?
- What's been tried or considered already?

Don't move to options until you have a clear problem statement you can write down as the **Goal**.

## 2. Explore options

Generate breadth before depth. Aim for at least 3–5 distinct approaches before evaluating any of them. For each:
- Describe it concisely
- Name its core trade-off
- Ask how it lands with the user

Resist evaluating too early — the goal here is to make sure you haven't missed an obvious alternative.

## 3. Evaluate and stress-test

For options that gain traction, dig in:
- What could go wrong?
- What assumptions is this approach making?
- What's the hardest part to get right?

Ask the user to steelman the worst-case scenario for the options they're drawn to.

## 4. Converge

Help the user land on a preferred direction, or a small set of options to prototype. Capture:
- What was chosen and the core reasoning
- What was discarded and why
- What remains unresolved

The record should make it easy to hand off to a planning or prototyping workflow.

# The markdown record

Write the record to `brainstorms/YYMMDD-HHMM-short-description.md`. Create the `brainstorms/` directory if it doesn't exist.

Update the file progressively throughout the session. A reader picking up the file mid-session should be able to follow what's been explored so far.

## Record structure

```markdown
# Brainstorm: [short title]

**Date:** YYYY-MM-DD
**Status:** In progress / Complete

## Goal

[What are we trying to achieve? What does success look like? Key constraints.]

## Options Explored

### [Option name]

[One-paragraph description. Core trade-off. Why it was considered.]

**Verdict:** Kept / Discarded — [one-line reason]

... repeat for each option ...

## Risks & Pitfalls

[Cross-cutting risks that apply regardless of approach, or risks specific to the preferred direction.]

## Preferred Direction

[What the brainstorm converged on, and why. If still unresolved, say so explicitly.]

## Open Questions

- [ ] [Decision or question that needs follow-up]

## Discarded Ideas

Brief log of ideas raised and rejected quickly.

| Idea | Reason discarded |
|---|---|
| ... | ... |
```

# File naming

`brainstorms/YYMMDD-HHMM-short-description.md`

Example: `brainstorms/260501-1430-gps-trace-segmentation.md`

Use today's date and a slug derived from the problem being explored.
