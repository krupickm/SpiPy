# CLAUDE.md
 
Read `insilico-screening.md` for full project specification.
 
## Scope boundary
 
All heavy computation (ORCA, CREST, Multiwfn) runs outside this package on an HPC cluster via PBS. This package manages job submission.
 
The Python code is responsible for:
- **Writing** structured input files (ORCA `.inp`, CREST input dirs, Multiwfn batch scripts)
- **Submitting** PBS jobs (`qsub`) and polling for completion (file-based, not PBS API)
- **Reading** and **parsing** output files (ORCA `.out`/`.fchk`, CREST `crest_conformers.xyz`, Multiwfn text dumps) into Python data structures for downstream ML
 
Do not scaffold installation or loading of external packages (ORCA, CREST, Multiwfn, Metacentrum modules). Assume they are available on the cluster and callable by the submitted PBS scripts.
 
