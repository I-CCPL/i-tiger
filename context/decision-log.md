# Decision Log

## Architecture Decisions (2026-02-23)
1. MVP properties fixed to Band + DOS + Berry.
2. Non-SOC first.
3. MPI-only decomposition for MVP.
4. Input target: `chk+eig+mmn`, default compute path `chk+eig`.
5. Formatted `chk` excluded from MVP.
6. SI unit output selected.
7. User may directly modify `*.f90`; AI must re-read before edits.

## User Requests Applied
- Full MVP implementation with Fortran + MPI.
- Persistent context for cross-session restart.
- Explicit visibility that user edits Fortran sources directly.
- Context cleanup to remove duplicated topics across multiple files.

## Current Risks
1. `chk` binary layout/version variability.
2. Model/physics modules not completed yet.
3. Cluster-specific MPI/ifort runtime differences.
4. Concurrent edits between user and AI.

## Mitigation Notes
- Keep parser checks strict and explicit.
- Build in small compile-tested increments.
- Keep Makefile flags overridable.
- Re-open target files immediately before edit.

## Session Updates
- [2026-02-23 13:28:06 +0900] context-cleanup: script format updated for compact logging
- [2026-02-23 13:58:00 +0900] context-policy: added Karpathy-inspired AI coding guidelines into `context/README.md`
