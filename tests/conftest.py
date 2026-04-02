"""Shared test fixtures for SpiPy.

Provides known spiropyran SMILES and Mol objects for use across test modules.
"""

import pytest
from rdkit import Chem

# 1',3',3'-trimethylindolino-2'-spirobenzopyran (unsubstituted BIPS)
BIPS_SMILES = "CN1c2ccccc2C(C)(C)C12Oc1ccccc1C=C2"

# 6-nitro-BIPS (the classic photochromic spiropyran)
NITRO_BIPS_SMILES = "CN1c2ccccc2C(C)(C)C12Oc1cc([N+](=O)[O-])ccc1C=C2"

# 5'-chloro-BIPS (halogen variant for robustness testing)
CL_BIPS_SMILES = "CN1c2ccc(Cl)cc2C(C)(C)C12Oc1ccccc1C=C2"


@pytest.fixture
def bips_smiles():
    return BIPS_SMILES


@pytest.fixture
def nitro_bips_smiles():
    return NITRO_BIPS_SMILES


@pytest.fixture
def cl_bips_smiles():
    return CL_BIPS_SMILES


@pytest.fixture
def bips_mol():
    return Chem.MolFromSmiles(BIPS_SMILES)


@pytest.fixture
def nitro_bips_mol():
    return Chem.MolFromSmiles(NITRO_BIPS_SMILES)


@pytest.fixture
def cl_bips_mol():
    return Chem.MolFromSmiles(CL_BIPS_SMILES)
