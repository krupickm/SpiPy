"""Configuration loading for SpiPy pipeline.

Reads a TOML config file and merges with sensible defaults.
All config sections are represented as dataclasses.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class MolConfig:
    mc_form: str = "neutral"  # "neutral" or "zwitterionic"


@dataclass
class IsomersConfig:
    clash_min_distance: float = 1.0  # Angstrom


@dataclass
class CrestConfig:
    wrapper_command: str = ""
    solvent: str = "chcl3"
    energy_window: float = 30.0  # kJ/mol


@dataclass
class OrcaTierConfig:
    method: str = "B3LYP"
    basis: str = "def2-SVP"
    solvent: str = "CPCM(CHCl3)"
    job_type: str = "opt"
    energy_window: float = 20.0  # kJ/mol, only meaningful for tier1


@dataclass
class OrcaConfig:
    wrapper_command: str = ""
    charge: int = 0
    multiplicity: int = 1
    tier1: OrcaTierConfig = field(default_factory=OrcaTierConfig)
    tier2: OrcaTierConfig = field(
        default_factory=lambda: OrcaTierConfig(
            basis="def2-TZVP",
            job_type="opt freq",
            energy_window=0.0,
        )
    )


@dataclass
class MultiwfnConfig:
    wrapper_command: str = ""
    analyses: list[str] = field(
        default_factory=lambda: ["esp_surface", "polarisability", "hirshfeld_charges"]
    )


@dataclass
class DescriptorsConfig:
    mordred: bool = True
    morfeus: bool = True
    soap: bool = True


@dataclass
class PipelineConfig:
    poll_interval: int = 60  # seconds


@dataclass
class SpiPyConfig:
    mol: MolConfig = field(default_factory=MolConfig)
    isomers: IsomersConfig = field(default_factory=IsomersConfig)
    crest: CrestConfig = field(default_factory=CrestConfig)
    orca: OrcaConfig = field(default_factory=OrcaConfig)
    multiwfn: MultiwfnConfig = field(default_factory=MultiwfnConfig)
    descriptors: DescriptorsConfig = field(default_factory=DescriptorsConfig)
    pipeline: PipelineConfig = field(default_factory=PipelineConfig)

    @classmethod
    def defaults(cls) -> SpiPyConfig:
        """Return config with all default values."""
        return cls()

    @classmethod
    def from_toml(cls, path: Path) -> SpiPyConfig:
        """Load config from a TOML file, merging with defaults."""
        with open(path, "rb") as f:
            raw = tomllib.load(f)

        cfg = cls.defaults()

        if "mol" in raw:
            cfg.mol = _merge(cfg.mol, raw["mol"])
        if "isomers" in raw:
            cfg.isomers = _merge(cfg.isomers, raw["isomers"])
        if "crest" in raw:
            cfg.crest = _merge(cfg.crest, raw["crest"])
        if "orca" in raw:
            orca_raw = raw["orca"]
            tier1_raw = orca_raw.pop("tier1", None)
            tier2_raw = orca_raw.pop("tier2", None)
            cfg.orca = _merge(cfg.orca, orca_raw)
            if tier1_raw:
                cfg.orca.tier1 = _merge(cfg.orca.tier1, tier1_raw)
            if tier2_raw:
                cfg.orca.tier2 = _merge(cfg.orca.tier2, tier2_raw)
        if "multiwfn" in raw:
            cfg.multiwfn = _merge(cfg.multiwfn, raw["multiwfn"])
        if "descriptors" in raw:
            cfg.descriptors = _merge(cfg.descriptors, raw["descriptors"])
        if "pipeline" in raw:
            cfg.pipeline = _merge(cfg.pipeline, raw["pipeline"])

        if cfg.mol.mc_form not in ("neutral", "zwitterionic"):
            raise ValueError(
                f"mol.mc_form must be 'neutral' or 'zwitterionic', "
                f"got {cfg.mol.mc_form!r}"
            )

        return cfg


def _merge(dc: object, overrides: dict) -> object:
    """Set attributes on a dataclass instance from a dict, ignoring unknown keys."""
    for key, value in overrides.items():
        if hasattr(dc, key):
            setattr(dc, key, value)
    return dc
