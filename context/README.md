# i-tiger Context (Minimal)

## Goal
- Implement Fortran MVP for Wannier-basis properties.
- Scope: Band + DOS + Berry curvature (non-SOC first).

## Constraints
- Input priority: `chk + eig + mmn`
- Default compute path: `chk + eig` (`mmn` optional unless advanced terms are requested)
- Formatted `chk` is out of MVP
- Output format: TXT + CSV
- Unit system: SI
- Parallel model: MPI-only

## Collaboration Rules
- User may directly edit any `*.f90` file at any time.
- AI must re-read target Fortran files before modifying.
- Never overwrite user edits blindly.

## AI Coding Guidelines
- Tradeoff: default to caution over speed, but use judgment for trivial tasks.

### 1) Think Before Coding
- State assumptions explicitly.
- If multiple interpretations exist, surface options instead of silently choosing one.
- If unclear, stop and ask instead of guessing.
- If a simpler path exists, mention it first.

### 2) Simplicity First
- Implement only what was requested.
- Avoid speculative abstractions/configurability for single-use code.
- Keep code minimal; remove avoidable complexity.

### 3) Surgical Changes
- Touch only lines directly needed for the request.
- Do not refactor unrelated areas.
- Match existing style.
- Remove only unused code created by your own change.

### 4) Goal-Driven Execution
- Define verifiable success criteria before implementation.
- For multi-step tasks, state concise step/check pairs.
- Verify with commands/tests whenever possible before claiming completion.

## Current Plan
1. Rebuild minimal scaffold after reset
2. Implement core types and I/O readers
3. Implement `H(q) -> H(R) -> H(k)` model
4. Implement Band/DOS/Berry modules
5. Add MPI reduction path
6. Add CLI/output/validation flow

## Workflow
- Update this file when goals/constraints/plan change.
- Append event-level history to `context/decision-log.md`.
- Use logger script:
  - `context/context_append_log.sh "title" "message"`
  - `context/context_append_log.sh --plan "title" "message"`

## Structure Policy
- Keep context docs minimal (single source here + change log).
- Avoid topic-split docs unless strictly necessary.

## Plan Updates
- [2026-02-23 13:45:27 +0900] structure: context/docs minimized and cases merged into tests
