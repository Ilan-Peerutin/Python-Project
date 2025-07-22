import random
import math
import numpy as np
from typing import Optional

class TauProtein(Protein):
    def __init__(self, name, sequence: Optional[np.array]=None, weight=None, length=None, organism=None, location=None, expression_level=None):
        super().__init__(name, sequence, weight, length, organism, location, expression_level)
        self.S202_T205_phosphorylation = False
        self.T231_phosphorylation = False
        self.S396_S404_phosphorylation = False
        self.phosphorylation_rate = 0.1
        self.phosphatase_activity = 0.05
        self.temperature = 37
        self.aggregation_level = 0.0
        self.microtubule_binding = 1.0

    def phosphorylate(self, site):
        if site == "S202_T205":
            self.S202_T205_phosphorylation = random.random() < self.phosphorylation_rate
        elif site == "T231":
            if self.S202_T205_phosphorylation:
                self.T231_phosphorylation = random.random() < (self.phosphorylation_rate * 1.5)
            else:
                self.T231_phosphorylation = random.random() < self.phosphorylation_rate
        elif site == "S396_S404":
            if self.T231_phosphorylation:
                self.S396_S404_phosphorylation = random.random() < (self.phosphorylation_rate * 1.5)
            else:
                self.S396_S404_phosphorylation = random.random() < self.phosphorylation_rate

    def dephosphorylate(self, site):
        if site == "S202_T205":
            if random.random() < self.phosphatase_activity:
                self.S202_T205_phosphorylation = False
        elif site == "T231":
            if random.random() < self.phosphatase_activity:
                self.T231_phosphorylation = False
        elif site == "S396_S404":
            if random.random() < self.phosphatase_activity:
                self.S396_S404_phosphorylation = False

    def simulate_phosphorylation(self, condition, temperature=None):
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
        print(f"S202_T205 Phosphorylation: {self.S202_T205_phosphorylation}")
        print(f"T231 Phosphorylation: {self.T231_phosphorylation}")
        print(f"S396_S404 Phosphorylation: {self.S396_S404_phosphorylation}")
        print(f"Aggregation Level: {self.aggregation_level}")
        print(f"Microtubule Binding: {self.microtubule_binding}")

tau = TauProtein(name="Tau Protein")

conditions = ["healthy", "pathological", "recovery"]
temperatures = [37, 40, 30]

for condition in conditions:
    for temperature in temperatures:
        print(f"\nSimulating under {condition} conditions at {temperature}°C:")
        tau.simulate_phosphorylation(condition, temperature)
        tau.display_phosphorylation()
