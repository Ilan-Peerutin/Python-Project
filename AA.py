"""
AA.py
Defines the AminoAcid class, which models the properties and behaviors of individual amino acids, including PTMs and molecular details.
"""
import math
import re

# AminoAcid class models a single amino acid and its properties
class AminoAcid:
    def __init__(self, name, three_letter, one_letter, polarity, charge, r_group, codon_list, pKa=None, volume=None, weight=None, PTM=False):
        # Initialize amino acid properties
        self.name = name
        self.three_letter = three_letter
        self.one_letter = one_letter
        self.polarity = polarity
        self.charge = charge
        self.r_group = r_group
        self.codon_list = codon_list
        self.pKa = pKa
        self.volume = volume
        self.PTM = PTM

        if weight is not None:
            self.weight = weight
        else:
            self.weight = self.calculate_weight()

    def __get_charge__(self, pH=None):
        # Calculate the charge of the amino acid at a given pH
        if pH is not None and self.pKa is not None:
            change = int(math.copysign(1, self.pKa - pH)) if self.pKa != pH else 1
            self.charge=self.charge-(change<0)
        return self.charge  
    
    def __has_codon__(self, codon):
        # Check if the codon is associated with this amino acid
        return codon.upper() in self.codon_list
    
    def is_identifier(self, identifier):
        # Check if the identifier matches the name, three-letter, or one-letter code

        identifier = identifier.upper()
        return identifier in {self.name.upper(), self.three_letter.upper(), self.one_letter.upper()}
    
    def calculate_weight(self):
        # Calculate the molecular weight of the amino acid

        atomic_weights = {
            'C': 12.01,
            'H': 1.008,
            'N': 14.01,
            'O': 16.00,
            'S': 32.06,
            'P': 30.97,
        }
        weight = 74.06
        if self.PTM == "Ubi":
            self.r_group=self.r_group[:-4]
            weight+=8565
        elif self.PTM == "GlcNAc":
            self.r_group=self.r_group[:-7]
            weight+=203
        cleaned = re.sub(r'[^A-Za-z0-9]', '', self.r_group)
        pattern = r'([A-Z][a-z]*)(\d*)'
        counts = {}
        
        for (elem, count) in re.findall(pattern, cleaned):
            count = int(count) if count else 1
            counts[elem] = counts.get(elem, 0) + count
        
        
        for atom, cnt in counts.items():
            w = atomic_weights.get(atom)
            if w is None:
                raise ValueError(f"Unknown element '{atom}' in formula")
            weight += w * cnt
        self.weight=weight
        return weight
    
    def add_PTM(self, modification):
        # Add a post-translational modification to the amino acid
        if self.one_letter in {"A","V","L","I","F","W"}:
            raise ("Amino acid does not undergo PTM")
        if modification in ("Phosphorylation", "p", "Phospo", "P"):
            if self.three_letter not in {"Ser", "Thr", "Tyr", "His", "Asp", "Glu", "Arg", "Lys", "Cys"}:
                raise Exception(f"{self.name} does not undergo Phosphorylation")
            self.r_group = self.r_group[:-1] + "-PO3"
            self.PTM = "Phospho"
            self.calculate_weight()
        elif modification in ("Acetylation", "a", "A", "Acetyl"):
            if self.three_letter not in {"Lys", "Met"}:
                raise Exception(f"{self.name} does not undergo Acetylation")
            self.r_group = self.r_group[:-1] + "-COCH3"
            self.PTM = "Acetyl"
            self.calculate_weight()
        elif modification in ("Methylation", "m", "M", "Methyl"):
            if self.three_letter not in {"Lys", "Arg", "His"}:
                raise Exception(f"{self.name} does not undergo Methylation")
            self.r_group = self.r_group[:-1] + "-CH3"
            self.PTM = "Methyl"
            self.calculate_weight()
        elif modification in ("Ubiquitination", "u", "U", "Ubi"):
            if self.three_letter not in {"Lys", "Met"}:
                raise Exception(f"{self.name} does not undergo Ubiquitination")
            self.r_group = self.r_group[:-1] + "-UBI"
            self.PTM = "Ubi"
            self.calculate_weight()
        elif modification in ("O-GlcNAcylation", "O-Glc", "GlcNAc"):
            if self.three_letter not in {"Ser", "Thr"}:
                raise Exception(f"{self.name} does not undergo O-GlcNAcylation")
            self.r_group = self.r_group[:-1] + "-GlcNAc"
            self.PTM = "GlcNAc"
            self.calculate_weight()
        else:
            raise Exception(f"Unknown modification: {modification}")
            
    def remove_PTM(self):
        if not hasattr(self, "PTM") or not self.PTM:
            raise ValueError(f"No PTM to remove from {self.name}")

        if self.PTM == "Phospho":
            self.r_group = self.r_group.replace("-PO3", "H")
        elif self.PTM == "Acetyl":
            self.r_group = self.r_group.replace("-COCH3", "H")
        elif self.PTM == "Methyl":
            self.r_group = self.r_group.replace("-CH3", "H")
        elif self.PTM == "Ubi":
            self.r_group = self.r_group.replace("-UBI", "H")
        elif self.PTM == "GlcNAc":
            self.r_group = self.r_group.replace("-GlcNAc", "H")
        else:
            raise ValueError(f"Unknown PTM type on {self.name}: {self.PTM}")

        self.PTM = None
        self.calculate_weight()


    def __str__(self):
        # String representation of the amino acid
        return (
            f"Amino Acid: {self.name} "
            f"({self.one_letter}, {self.three_letter})\n"
            f"  Polarity: {self.polarity}\n"
            f"  Charge: {self.charge}\n"
            f"  R-group: {self.r_group}\n"
            f"  Codons: {', '.join(self.codon_list)}\n"
            f"  pKa: {self.pKa if self.pKa is not None else 'N/A'}\n"
            f"  Volume: {self.volume if self.volume is not None else 'N/A'}\n"
        )
        
