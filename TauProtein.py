"""
TauProtein.py
Defines the TauProtein class, which models tau protein isoforms and their molecular transitions (phosphorylation, aggregation, truncation, etc.).
"""
import random
import math
import numpy as np
from typing import Optional
from Protein import Protein
from truncation import ProteinTruncator
from AA import AminoAcid
import re
from Environment import Environment as env
from collections import defaultdict

# TauProtein class models tau-specific logic, including isoforms, phosphorylation, aggregation, truncation, and simulation state
class TauProtein(Protein):
    def __init__(self, name="Tau", isoform="4R", sequence: Optional[np.array]=None,
                 weight=None, length=None, organism=None, location=None, expression_level=None):
        # Initialize tau protein properties, including isoform, phosphorylation, aggregation, and history
        super().__init__(name, sequence, weight, length, organism, location, expression_level)
        # Site-specific phosphorylation booleans
        self.S202_T205_phosphorylation = False
        self.T231_phosphorylation = False
        self.S396_S404_phosphorylation = False
        self.phosphorylation_rate = 0.1
        self.phosphatase_activity = 0.05
        self.temperature = 37
        self.aggregation_level = 0.0
        self.microtubule_binding = 1.0
        # Advanced/General attributes
        self.isoform = isoform
        self.phosphorylation_sites = {i: np.array([np.random.rand()]) for i in range(1, 80)}  # creating this for the 'update_state' func
        self.aggregation_state = "monomer"  # or "oligomer", "fibril"
        self.truncated_site = None
        self.soluble = True
        self.history = []
        self.age = 0
        self.is_truncated = False
        self.pathological = False
        self.sequence = []
        # If sequence is not provided, initialize as empty list
        if self.sequence is None:
            self.sequence = []

    def define_isoform(self, exon):
        # Define the tau isoform based on exon information
        if exon == "R2 at 10":
            self.isoform = "4R"
        if exon == "R1 at 10":
            self.isoform = "3R"
        return self.isoform

    # --- Site-specific phosphorylation logic ---
    def phosphorylate(self, site):
        # Site-specific phosphorylation logic for key tau residues
        if site == "S202_T205":
            self.S202_T205_phosphorylation = random.random() < self.phosphorylation_rate
            self.phosphorylation_sites[site] = self.S202_T205_phosphorylation
        elif site == "T231":
            if self.S202_T205_phosphorylation:
                self.T231_phosphorylation = random.random() < (self.phosphorylation_rate * 1.5)
            else:
                self.T231_phosphorylation = random.random() < self.phosphorylation_rate
            self.phosphorylation_sites[site] = self.T231_phosphorylation
        elif site == "S396_S404":
            if self.T231_phosphorylation:
                self.S396_S404_phosphorylation = random.random() < (self.phosphorylation_rate * 1.5)
            else:
                self.S396_S404_phosphorylation = random.random() < self.phosphorylation_rate
            self.phosphorylation_sites[site] = self.S396_S404_phosphorylation
        else:
            # General site logic
            self.phosphorylation_sites[site] = random.random() < self.phosphorylation_rate

    def dephosphorylate(self, site):
        # Site-specific dephosphorylation logic for key tau residues
        if site == "S202_T205":
            if random.random() < self.phosphatase_activity:
                self.S202_T205_phosphorylation = False
                self.phosphorylation_sites[site] = False
        elif site == "T231":
            if random.random() < self.phosphatase_activity:
                self.T231_phosphorylation = False
                self.phosphorylation_sites[site] = False
        elif site == "S396_S404":
            if random.random() < self.phosphatase_activity:
                self.S396_S404_phosphorylation = False
                self.phosphorylation_sites[site] = False
        else:
            if random.random() < self.phosphatase_activity:
                self.phosphorylation_sites[site] = False

    def simulate_phosphorylation(self, condition, temperature=None):
        # Simulate phosphorylation under different biological conditions
        if temperature is not None:
            self.temperature = temperature
        temperature_factor = 1 + (self.temperature - 37) * 0.01
        self.phosphorylation_rate *= temperature_factor
        if condition == "healthy":
            self.phosphorylation_rate = 0.1
            self.phosphatase_activity = 0.05
        elif condition == "pathological":
            self.phosphorylation_rate = 0.5
            self.phosphatase_activity = 0.01
        elif condition == "recovery":
            self.phosphorylation_rate = 0.05
            self.phosphatase_activity = 0.1
        self.phosphorylate("S202_T205")
        self.phosphorylate("T231")
        self.phosphorylate("S396_S404")
        self.update_aggregation_and_binding()

    def update_aggregation_and_binding(self):
        # Update aggregation and microtubule binding state based on phosphorylation
        total_phosphorylation = sum([self.S202_T205_phosphorylation, self.T231_phosphorylation, self.S396_S404_phosphorylation])
        if total_phosphorylation >= 0.7:
            self.aggregation_level = 1.0
            self.microtubule_binding = 0.0
        elif total_phosphorylation >= 0.5:
            self.aggregation_level = 0.7
            self.microtubule_binding = 0.3
        elif total_phosphorylation >= 0.3:
            self.aggregation_level = 0.3
            self.microtubule_binding = 0.6
        else:
            self.aggregation_level = 0.0
            self.microtubule_binding = 1.0

    def display_phosphorylation(self):
        # Print the current phosphorylation and aggregation state
        print(f"S202_T205 Phosphorylation: {self.S202_T205_phosphorylation}")
        print(f"T231 Phosphorylation: {self.T231_phosphorylation}")
        print(f"S396_S404 Phosphorylation: {self.S396_S404_phosphorylation}")
        print(f"Aggregation Level: {self.aggregation_level}")
        print(f"Microtubule Binding: {self.microtubule_binding}")

    # --- Advanced/General aggregation logic ---
    def count_phosphorylated_residues(self):
        # Count the number of phosphorylated sites (probability > 0.5)
        return sum(1 for p in self.phosphorylation_sites.values() if p > 0.5)
    
    def detect_aggregation_motifs(self):
        # Detect aggregation-prone motifs in the tau sequence
        sequence_str = ''.join(aa.one_letter for aa in self.sequence if hasattr(aa, 'one_letter'))
        motifs = ['VQIINK', 'VQIVYK']
        count = 0
        for motif in motifs:
            count += len(re.findall(motif, sequence_str))
        return count
    
    def compute_aggregation_score(self):
        # Compute an aggregation score based on phosphorylation and motifs
        phospho_count = self.count_phosphorylated_residues()
        motif_count = self.detect_aggregation_motifs()
        score = 0
        score += phospho_count * 1.5
        score += motif_count * 2
        if self.is_truncated:
            score += 2
        if self.isoform == '4R':
            score += 1
        return score
    
    def update_aggregation_state(self):
        # Update the aggregation state (monomer, oligomer, fibril) based on score
        score = self.compute_aggregation_score()
        if score >= 7:
            self.aggregation_state = 'fibril'
        elif score >= 4:
            self.aggregation_state = 'oligomer'
        else:
            self.aggregation_state = 'monomer'

    def update_state(self, environment, timepoints: np.array):
        """
        Main orchestrator: updates tau protein state over a series of timepoints by checking temperature, kinase, phosphatase, protease, and oxidative stress effects.
        """

        site_probabilities = defaultdict(list)
        temp_effect = self.check_temp(environment)
        kinase_effect = self.check_kinase(environment)
        phosphatase_effect = self.check_phosphatase(environment)
        protease_effect = self.check_protease(environment)
        oxidative_effect = self.check_oxidative_stress(environment)

        k_p = temp_effect * kinase_effect * oxidative_effect
        k_d = phosphatase_effect * protease_effect

        last_probs = self.phosphorylation_sites.copy()

        for site, prob_array in last_probs.items():
            prob_array = list(prob_array)
            for i in range(1, len(timepoints)):
                P_prev = prob_array[i - 1]
                P_new = P_prev + (k_p * (1 - P_prev) - k_d * P_prev)
                if P_new < 0:
                    P_new = P_new * -1
                
                prob_array.append(P_new)
                
            last_probs[site] = prob_array

        transposed = zip(*last_probs.values())
        site_probabilities = [sum(timepoint_probs) / len(timepoint_probs) for timepoint_probs in transposed]

        return site_probabilities
      
    def check_temp(self, environment):
        healthy_temp_range = (36, 38)
        if healthy_temp_range[0] <= environment.temperature <= healthy_temp_range[1]:
            return 1.0
        elif environment.temperature < healthy_temp_range[0]:
            return 0.7
        else:
            return 1.3

    def check_kinase(self, environment):
        healthy_kinase_range = (0.8, 1)
        if healthy_kinase_range[0] <= environment.temperature <= healthy_kinase_range[1]:
            return 1.0
        else:
            return 0.05 * environment.kinase_level

    def check_phosphatase(self, environment):
        phosphatase_range = (0.8, 1)
        if phosphatase_range[0] <= environment.kinase_level <= phosphatase_range[1]:
            return environment.kinase_level
        else:
         return 0.02 * environment.phosphatase_level

    def check_protease(self, environment):
        if 0.4 < environment.protease_level < 0.8:
            return 0.7
        elif environment.protease_level <= 0.4:
            return 0.2
        elif environment.protease_level <= 1:
            return 1.5
        return 1.0

    def check_oxidative_stress(self, environment):
        if 0.4 < environment.oxidative_stress < 0.8:
            return 0.7
        elif environment.oxidative_stress <= 0.4:
            return 0.2
        elif environment.oxidative_stress <= 1:
            return 1.5
        return 1.0

    def truncate(self, site):
        # Truncate the tau protein sequence at a given site
        if self.sequence is not None:
            truncated_seq, trunc_aa = ProteinTruncator.truncate(self.sequence, site)
            self.sequence = truncated_seq
            self.truncated_site = site
            self.is_truncated = True
            self.truncation_aa = trunc_aa
        else:
            raise ValueError("No sequence to truncate.") 
