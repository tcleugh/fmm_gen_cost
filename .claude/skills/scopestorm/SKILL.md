---
name: scopestorm
description: Turns a brainstorm into a Jira-ready Epic + Task breakdown and a scope summary document suitable for generating a formal scope Word document.
---

This skill takes the output of a `/brainstorm` session and works with the user to produce two things:
1. A **Jira-ready breakdown** — one Epic with a set of Tasks underneath it
2. A **scope summary narrative** — a structured prose section suitable for handing to another tool to generate a formal scope document

Inform the user that the scopestorm skill is being used.

# Starting the session

1. **Find the brainstorm.** If the user hasn't specified a file, list what's in `brainstorms/` and ask which one to use.
2. **Read the brainstorm file** in full before saying anything else.
3. **Surface the preferred direction.** Quote or paraphrase it back to the user and confirm it's still the right starting point. If the brainstorm has no clear preferred direction, say so and resolve it before proceeding — don't guess.

Don't move forward until the preferred direction is confirmed.

# Working style

- **One question at a time.** Same rule as brainstorm — pick the most important thing and ask that.
- **Propose, don't just prompt.** Draft task titles or descriptions for the user to react to. Blank-slate questions stall the session.
- **Challenge scope creep.** If a task feels like it exceeds the brainstorm's preferred direction, flag it and ask whether it belongs here or in a future epic.
- **Be explicit about exclusions.** Every scope document is strengthened by a clear Out of Scope section — actively probe for things that might be assumed to be included but aren't.
- **Write progressively.** Update the output file after each significant decision. Don't accumulate everything and write at the end.
- **Translate internal shorthand into self-contained language.** Brainstorm sessions naturally generate option labels and code-names — things like `D-mesh` / `D-hub`, `Workstream A`, `V1 / V2 / V3`, `K-Refit`, `F-PSL + U-Hybrid`, `L-Light`. These are useful conversational shortcuts during the brainstorm, but they are meaningless to anyone reading the scope doc cold. When writing the scope doc, **translate every such label into what it actually is**:
    - `Workstream 0` → "the confident-sample task" or the task title itself
    - `D-mesh` → "the mesh-block-level anchor candidate set"
    - `K-Refit` → "a light end-to-end recalibration with no structural change"
    - `V1` → "the corridor checklist pass/fail check"
    - `F-PSL + U-Hybrid` → "a path-size logit with separate commute and leisure cost functions"
  The same principle applies to references to external documents — don't leave bare citations like `Tables 3-11/13/14` or `§3.5.2` of a report; describe what those sections actually contain in plain language.
  Standard project-management vocabulary (Must / Should / Could, MoSCoW tiers) and well-known technical terms (MNL, PSL, OD matrix) are fine, especially if introduced inline. The test: a stakeholder who has never seen the brainstorm should be able to read the scope doc end-to-end without needing to ask "what's a D-mesh?"

# Session structure

## 1. Confirm the direction

Read the brainstorm's **Preferred Direction** and **Open Questions** sections. If there are unresolved open questions that affect scope, resolve them now before defining the Epic.

## 2. Define the Epic

Work with the user to nail down:
- **Title** — short, action-oriented (e.g. "Implement GPS trace segmentation by transport mode")
- **Goal** — 2–3 sentences: what problem this solves and what success looks like
- **Acceptance criteria** — the conditions under which this epic is considered done

## 3. Break down into Tasks

Propose a draft task list based on the brainstorm. Each task should be:
- Independently actionable by one person
- Completable within a sprint (roughly)
- Named clearly enough that someone unfamiliar with the brainstorm understands what to do

For each task, work through:
- Title
- Description (what to do and why — enough for someone to pick it up cold)
- Acceptance criteria (how we know it's done)
- Any notes: dependencies, gotchas, or links to relevant brainstorm options

Iterate — add, remove, split, or merge tasks with the user until the list feels right.

## 4. Define Out of Scope

Explicitly list things that are adjacent to this epic but deliberately excluded. Pull from the brainstorm's **Discarded Ideas** and **Deferred** items as a starting point. Ask the user to confirm and add to this list.

## 5. Summarise risks

Pull relevant risks from the brainstorm. Trim to the ones that are still live given the chosen direction — don't carry over risks that were only relevant to discarded options.

## 6. Write the scope summary

Write a 4–6 paragraph narrative at the end of the document. This should be self-contained — a reader who has never seen the brainstorm should be able to understand:
- The problem being solved
- The chosen approach and why
- What will be delivered (the tasks, summarised)
- What is explicitly out of scope
- Key risks and mitigations

This section is the handoff to a formal scope document generator. Write it in clear, professional prose — not bullet points.

## 7. Self-containment review

Before declaring the scope doc done, read it end-to-end as if you'd never seen the brainstorm. Look specifically for:
- Option labels and code-names from the brainstorm (`D-mesh`, `V1`, `K-Refit`, `Workstream A`, etc.) — translate any that survived into plain-language descriptions.
- Bare citations to external documents (`Section 3.5`, `Tables 3-11/13/14`, `Eq. 5`) — replace with a brief description of what's actually being referenced.
- Acronyms used without introduction.
- Cross-task references that assume the reader knows what an earlier task produced.

If anything would prompt a reasonable stakeholder to ask "what does that mean?", rewrite it. Then ask the user for any final corrections.

# Output file

Write to `scopestorms/YYMMDD-HHMM-short-description.md`. Create the directory if it doesn't exist. Use the same slug as the source brainstorm file where possible.

## Output structure

```markdown
# Scope: [title]

**Date:** YYYY-MM-DD
**Source brainstorm:** [relative path to brainstorm file]
**Status:** Draft / Final

---

## Epic

**Title:** [Jira epic title]

**Goal:** [2–3 sentences: problem + success definition]

**Acceptance Criteria:**
- [ ] ...

---

## Tasks

### Task 1: [Title]

**Description:** [What to do and why — enough for a cold handoff]

**Acceptance Criteria:**
- [ ] ...

**Notes:** [Dependencies, gotchas, references to brainstorm options — omit if none]

---

... repeat per task ...

---

## Out of Scope

- [Explicit exclusion]
- ...

---

## Risks

| Risk | Impact | Mitigation |
|---|---|---|
| ... | High / Med / Low | ... |

---

## Scope Summary

[4–6 paragraphs of self-contained prose covering: problem, approach, deliverables, exclusions, risks. Written for a non-technical stakeholder who will use this to generate a formal scope document.]
```

# File naming

`scopestorms/YYMMDD-HHMM-short-description.md`

Example: `scopestorms/260501-1500-gps-trace-segmentation.md`

Match the slug to the source brainstorm file where possible.
