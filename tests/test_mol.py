"""Tests for spipy.mol — validation, 3D embedding, XYZ I/O, SP→MC transform."""

import pytest
from rdkit import Chem

from spipy.mol import MolError, embed_3d, mol_to_xyz, sp_to_mc, validate_smiles


class TestValidateSmiles:
    def test_valid_smiles_returns_mol(self):
        mol = validate_smiles("c1ccccc1")
        assert isinstance(mol, Chem.Mol)

    def test_invalid_smiles_raises(self):
        with pytest.raises(MolError):
            validate_smiles("this_is_not_smiles")

    def test_empty_string_raises(self):
        with pytest.raises(MolError):
            validate_smiles("")

    def test_roundtrip_canonical(self, bips_smiles):
        mol = validate_smiles(bips_smiles)
        out = Chem.MolToSmiles(mol)
        mol2 = Chem.MolFromSmiles(out)
        assert Chem.MolToSmiles(mol2) == out

    def test_spiropyran_parses(self, bips_smiles, nitro_bips_smiles, cl_bips_smiles):
        for smi in (bips_smiles, nitro_bips_smiles, cl_bips_smiles):
            mol = validate_smiles(smi)
            assert mol is not None


class TestEmbed3D:
    def test_produces_conformer(self, bips_mol):
        mol3d = embed_3d(bips_mol)
        assert mol3d.GetNumConformers() > 0

    def test_reasonable_bond_lengths(self, bips_mol):
        mol3d = embed_3d(bips_mol)
        conf = mol3d.GetConformer()
        for bond in mol3d.GetBonds():
            p1 = conf.GetAtomPosition(bond.GetBeginAtomIdx())
            p2 = conf.GetAtomPosition(bond.GetEndAtomIdx())
            dist = p1.Distance(p2)
            assert 0.8 < dist < 2.5, f"Unreasonable bond length: {dist:.3f} Å"

    def test_has_explicit_hydrogens(self, bips_mol):
        mol3d = embed_3d(bips_mol)
        h_count = sum(1 for a in mol3d.GetAtoms() if a.GetAtomicNum() == 1)
        assert h_count > 0

    def test_nitro_bips_embeds(self, nitro_bips_mol):
        """Nitro group can trip MMFF — verify UFF fallback works."""
        mol3d = embed_3d(nitro_bips_mol)
        assert mol3d.GetNumConformers() > 0


class TestMolToXyz:
    def test_writes_valid_xyz(self, bips_mol, tmp_path):
        mol3d = embed_3d(bips_mol)
        out = mol_to_xyz(mol3d, tmp_path / "test.xyz")
        assert out.exists()

        lines = out.read_text().splitlines()
        natoms = int(lines[0].strip())
        assert natoms == mol3d.GetNumAtoms()
        # Line 2 is comment, lines 3+ are coordinates
        for line in lines[2:]:
            parts = line.split()
            assert len(parts) == 4, f"Expected 4 fields, got {len(parts)}: {line}"

    def test_no_conformer_raises(self, tmp_path):
        mol = Chem.MolFromSmiles("C")
        with pytest.raises(MolError):
            mol_to_xyz(mol, tmp_path / "fail.xyz")


class TestSpToMc:
    def test_bips_produces_mc(self, bips_mol):
        mc = sp_to_mc(bips_mol)
        assert mc is not None
        # Heavy atom count should be preserved
        sp_heavy = sum(1 for a in bips_mol.GetAtoms() if a.GetAtomicNum() != 1)
        mc_heavy = sum(1 for a in mc.GetAtoms() if a.GetAtomicNum() != 1)
        assert sp_heavy == mc_heavy

    def test_nitro_bips_produces_mc(self, nitro_bips_mol):
        mc = sp_to_mc(nitro_bips_mol)
        assert mc is not None

    def test_cl_bips_produces_mc(self, cl_bips_mol):
        mc = sp_to_mc(cl_bips_mol)
        assert mc is not None

    def test_mc_no_longer_matches_spiro_smarts(self, bips_mol):
        mc = sp_to_mc(bips_mol)
        # The closed-form spiro pattern should not match in the open MC
        spiro_pat = Chem.MolFromSmarts("[N][C]1([O]c2ccccc2[C]=[C]1)")
        assert not mc.HasSubstructMatch(spiro_pat)

    def test_neutral_form_has_phenol_oh(self, bips_mol):
        mc = sp_to_mc(bips_mol, mc_form="neutral")
        # Look for oxygen with exactly 1 H (phenol)
        found_oh = False
        for atom in mc.GetAtoms():
            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:
                found_oh = True
                break
        assert found_oh, "Neutral MC should have a phenol OH"

    def test_zwitterionic_form_has_charges(self, bips_mol):
        mc = sp_to_mc(bips_mol, mc_form="zwitterionic")
        has_o_minus = any(
            a.GetFormalCharge() == -1 and a.GetAtomicNum() == 8
            for a in mc.GetAtoms()
        )
        has_n_plus = any(
            a.GetFormalCharge() == 1 and a.GetAtomicNum() == 7
            for a in mc.GetAtoms()
        )
        assert has_o_minus, "Zwitterionic MC should have O⁻"
        assert has_n_plus, "Zwitterionic MC should have N⁺"

    def test_non_spiropyran_raises(self):
        benzene = Chem.MolFromSmiles("c1ccccc1")
        with pytest.raises(MolError, match="Spiro centre not found"):
            sp_to_mc(benzene)

    def test_invalid_mc_form_raises(self, bips_mol):
        with pytest.raises(MolError, match="mc_form"):
            sp_to_mc(bips_mol, mc_form="invalid")

    def test_mc_passes_sanitization(self, bips_mol):
        mc = sp_to_mc(bips_mol)
        # Should not raise
        Chem.SanitizeMol(mc)

    def test_mc_zwitterionic_passes_sanitization(self, bips_mol):
        mc = sp_to_mc(bips_mol, mc_form="zwitterionic")
        Chem.SanitizeMol(mc)
