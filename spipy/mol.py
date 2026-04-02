"""Molecule handling for SpiPy: validation, 3D embedding, XYZ I/O, SP→MC transform.

The SP→MC transform uses programmatic RWMol bond manipulation (not Reaction SMARTS)
for debuggability. It targets the indolinospiropyran (BIPS) family specifically.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


class MolError(Exception):
    """Raised for molecule processing errors."""


# SMARTS for the spiropyran core: N-C(spiro) in a ring containing O, fused phenyl, and C=C.
# Atom maps: 1=N(indoline), 2=C(spiro), 3=O(pyran), 4=C3(alkene near spiro), 5=C4(alkene near phenyl)
#
# The pattern matches the 2H-chromene (benzopyran) fused ring with spiro junction:
#   N—C(spiro) is part of the indoline ring
#   C(spiro) is also in the pyran ring: O-c(ar)-c(ar)...-C=C-C(spiro)
_SPIRO_SMARTS = Chem.MolFromSmarts("[N:1][C:2]1([O:3]c2ccccc2[C:5]=[C:4]1)")


def validate_smiles(smiles: str) -> Chem.Mol:
    """Parse and sanitize a SMILES string.

    Returns a sanitized RDKit Mol object.

    Raises:
        MolError: If SMILES cannot be parsed.
    """
    if not smiles:
        raise MolError("Empty SMILES string")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise MolError(f"Cannot parse SMILES: {smiles!r}")
    return mol


def embed_3d(mol: Chem.Mol, random_seed: int = 42) -> Chem.Mol:
    """Generate 3D coordinates via ETKDG and optimize with MMFF94 (UFF fallback).

    Adds explicit hydrogens. The returned Mol includes Hs and a 3D conformer.

    Raises:
        MolError: If 3D embedding fails entirely.
    """
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    status = AllChem.EmbedMolecule(mol, params)

    if status == -1:
        # Retry with random initial coordinates for difficult molecules
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(mol, params)

    if status == -1:
        raise MolError("3D embedding failed — could not generate coordinates")

    # Optimize: MMFF first, UFF fallback
    mmff_status = AllChem.MMFFOptimizeMolecule(mol)
    if mmff_status == -1:
        # MMFF parameterization failed — fall back to UFF
        AllChem.UFFOptimizeMolecule(mol)

    return mol


def mol_to_xyz(mol: Chem.Mol, path: Path) -> Path:
    """Write molecule with 3D conformer to XYZ file format.

    Raises:
        MolError: If the molecule has no 3D conformer.
    """
    if mol.GetNumConformers() == 0:
        raise MolError("Molecule has no 3D conformer — call embed_3d first")

    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), Chem.MolToSmiles(Chem.RemoveHs(mol))]

    for i in range(mol.GetNumAtoms()):
        symbol = mol.GetAtomWithIdx(i).GetSymbol()
        pos = conf.GetAtomPosition(i)
        lines.append(f"{symbol:2s} {pos.x:14.8f} {pos.y:14.8f} {pos.z:14.8f}")

    path = Path(path)
    path.write_text("\n".join(lines) + "\n")
    return path


def sp_to_mc(mol: Chem.Mol, mc_form: str = "neutral") -> Chem.Mol:
    """Convert spiropyran (closed SP) to merocyanine (open MC) form.

    Uses SMARTS to identify the spiro centre, then programmatic RWMol
    manipulation for the ring-opening bond changes.

    Targets the indolinospiropyran (BIPS) family.

    Args:
        mol: RDKit Mol of a spiropyran.
        mc_form: ``"neutral"`` (phenol OH) or ``"zwitterionic"`` (O⁻ / N⁺).

    Returns:
        Mol of the merocyanine open form (no 3D coordinates).

    Raises:
        MolError: If the spiro centre is not found or the transform fails.
    """
    if mc_form not in ("neutral", "zwitterionic"):
        raise MolError(f"mc_form must be 'neutral' or 'zwitterionic', got {mc_form!r}")

    matches = mol.GetSubstructMatches(_SPIRO_SMARTS)
    if len(matches) == 0:
        raise MolError(
            "Spiro centre not found — is this an indolinospiropyran?"
        )
    if len(matches) > 1:
        raise MolError(
            f"Found {len(matches)} spiro centres — expected exactly 1 "
            f"(bis-spiropyrans are not supported)"
        )

    n_idx, spiro_idx, o_idx, c3_idx, c4_idx = matches[0]

    rwmol = Chem.RWMol(mol)

    # Kekulize so all bonds are explicit single/double — avoids conflicts
    # when modifying bonds adjacent to aromatic rings
    Chem.Kekulize(rwmol, clearAromaticFlags=True)

    # --- Bond surgery ---

    # 1. Break C(spiro)–O
    rwmol.RemoveBond(spiro_idx, o_idx)

    # 2. C(spiro)=C3: single → double
    rwmol.GetBondBetweenAtoms(spiro_idx, c3_idx).SetBondType(
        Chem.BondType.DOUBLE
    )

    # 3. C3–C4: double → single
    rwmol.GetBondBetweenAtoms(c3_idx, c4_idx).SetBondType(
        Chem.BondType.SINGLE
    )

    # 4. C4=C4a: single → double
    #    C4a is the neighbor of C4 that is not C3 and was in the aromatic ring
    c4a_idx = _find_c4a(rwmol, c4_idx, c3_idx)
    rwmol.GetBondBetweenAtoms(c4_idx, c4a_idx).SetBondType(
        Chem.BondType.DOUBLE
    )

    # --- Charges, protons, N–C bond order ---

    o_atom = rwmol.GetAtomWithIdx(o_idx)
    n_atom = rwmol.GetAtomWithIdx(n_idx)

    if mc_form == "neutral":
        # Phenol OH — O gets one H, stays uncharged
        o_atom.SetNumExplicitHs(1)
        o_atom.SetFormalCharge(0)
        # N–C(spiro) stays single (already is)
    else:
        # Zwitterionic: O⁻, N⁺, N=C(spiro) double bond
        o_atom.SetNumExplicitHs(0)
        o_atom.SetFormalCharge(-1)
        n_atom.SetFormalCharge(+1)
        rwmol.GetBondBetweenAtoms(n_idx, spiro_idx).SetBondType(
            Chem.BondType.DOUBLE
        )

    # Let RDKit recalculate implicit Hs and re-perceive aromaticity
    rwmol.UpdatePropertyCache(strict=False)

    try:
        Chem.SanitizeMol(rwmol)
    except Exception as e:
        raise MolError(f"Sanitization failed after ring opening: {e}") from e

    return rwmol.GetMol()


def _find_c4a(rwmol: Chem.RWMol, c4_idx: int, c3_idx: int) -> int:
    """Find C4a — the carbon neighbor of C4 that was part of the phenyl ring.

    C4a is the neighbor of C4 that is:
      - not C3
      - a carbon
      - bonded to at least 2 other carbons (part of the phenyl ring)
    """
    c4_atom = rwmol.GetAtomWithIdx(c4_idx)
    for neighbor in c4_atom.GetNeighbors():
        nidx = neighbor.GetIdx()
        if nidx == c3_idx:
            continue
        if neighbor.GetAtomicNum() == 6:
            # Should be in the phenyl ring — has multiple C neighbors
            c_neighbors = sum(
                1 for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 6
            )
            if c_neighbors >= 2:
                return nidx

    raise MolError("Could not identify C4a in the chromene ring")
