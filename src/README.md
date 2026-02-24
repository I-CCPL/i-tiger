# src Layout

This structure mirrors the project workflow in `context/README.md`.

- core: common definitions and low-level utilities
- io: Wannier90 file readers (`chk`, `eig`, `mmn`) and TXT/CSV writers
- model: system container and Hamiltonian path (`H(q) -> H(R) -> H(k)`)
- physics: Band/DOS/Berry implementations and numerical integration helpers
- parallel: MPI decomposition and reduction
- cli: runtime orchestration and configuration parsing
