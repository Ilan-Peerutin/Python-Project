from AA import AminoAcid
from Protein import Protein
from typing import Optional, List, Tuple
import random
import numpy as np
import matplotlib.pyplot as plt

def phosphorylation_constants(protein: Protein, list_of_PTMs: Optional[List[Tuple[int, str]]] = None):
    initial_phospo_constant=0.005
    initial_dephospo_constant=0.02
    acetyl_constant_ranges=[2,5,0.8,1]
    methyl_constant_ranges=[0.5,0.9,1,1]
    ubi_constant_ranges=[0.5,1,1.2,2]
    glyco_constant_ranges=[0.1,0.3,1.5,2.5]
    phospho_residues=1
    possible_phopho=0
    acetyl_counts=0
    methyl_counts=0
    ubi_counts=0
    glyco_counts=0

    if list_of_PTMs is not None:
        for residue in list_of_PTMs:
            #print (protein.sequence[residue[0]-1])
            protein.sequence[residue[0]-1].add_PTM(residue[1])
    for aa in protein.sequence:
        match aa.PTM:
            case "Phospho":
                phospho_residues+=1
                possible_phopho+=1
            case "Acetyl":
                acetyl_counts+=1
            case "Methyl":
                methyl_counts+=1
            case "Ubi":
                ubi_counts+=1
            case "GlcNAc":
                glyco_counts+=1
            case __ if aa.three_letter in {"Ser", "Thr", "Tyr"}:
                possible_phopho+=1
    
    aK=random.uniform(acetyl_constant_ranges[0],acetyl_constant_ranges[1])
    adK=random.uniform(acetyl_constant_ranges[2],acetyl_constant_ranges[3])
    mK  = random.uniform(methyl_constant_ranges[0], methyl_constant_ranges[1])
    mdK = random.uniform(methyl_constant_ranges[2], methyl_constant_ranges[3])
    uK  = random.uniform(ubi_constant_ranges[0], ubi_constant_ranges[1])
    udK = random.uniform(ubi_constant_ranges[2], ubi_constant_ranges[3])
    gK  = random.uniform(glyco_constant_ranges[0], glyco_constant_ranges[1])
    gdK = random.uniform(glyco_constant_ranges[2], glyco_constant_ranges[3])
    
    pK = initial_phospo_constant * (aK ** acetyl_counts) * (mK ** methyl_counts) * (uK ** ubi_counts) * (gK ** glyco_counts)
    dpK = initial_dephospo_constant * (adK ** acetyl_counts) * (mdK ** methyl_counts) * (udK ** ubi_counts) * (gdK ** glyco_counts)

    return pK, dpK

def phosphorylation(protein: Protein, phospho_k, dephospho_k):
    phospho_residues=0
    possible_phopho=0
    for aa in protein.sequence:
        match aa.PTM:
            case "Phospho":
                phospho_residues+=1
                possible_phopho+=1
            case __ if aa.three_letter in {"Ser", "Thr", "Tyr"}:
                possible_phopho+=1
    dP=phospho_k*(1-(phospho_residues/possible_phopho))-dephospho_k*phospho_residues

    phospho_residues=0
    possible_phopho=0
    for aa in protein.sequence:
        match aa.PTM:
            case "Phospho":
                if random.random()<dephospho_k:
                    aa.remove_PTM()
                    possible_phopho+=1
                else:
                    phospho_residues+=1
                    possible_phopho+=1
            case __ if aa.three_letter in {"Ser", "Thr", "Tyr"}:
                if random.random()<dP:
                    aa.add_PTM("Phosphorylation")
                    phospho_residues+=1
                    possible_phopho+=1
                else:
                    possible_phopho+=1

    return phospho_residues/possible_phopho*100

def phospo_over_time(protein: Protein, time: int, list_of_PTMs: Optional[List[Tuple[int, str]]] = None):
    p_percentage=np.zeros(time)
    for i in range(time):
        pK, dpK=phosphorylation_constants(protein,list_of_PTMs)
        p_percentage[i]=phosphorylation(protein,pK,dpK)
    return p_percentage

def run_and_plot_disease_simulation():
    tau_seq="MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
    tau_prot=Protein("Tau",tau_seq)
    acetyl_sites = [
        (163, "Acetyl"),
        (174, "Acetyl"),
        (180, "Acetyl"),
    ]

    #phospho_data=phospo_over_time(tau_prot,180, acetyl_sites)
    phospho_data=phospo_over_time(tau_prot,180)



    plt.figure(figsize=(10, 6))
    plt.plot(phospho_data, label='Phosphorylation %', color='mediumblue', linewidth=2)
    plt.title('Tau Phosphorylation Simulation Over Time')
    plt.xlabel('Time Steps')
    plt.ylabel('Phosphorylated Residues (%)')
    plt.ylim(0, 5)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    print("Close the plot window to return to the menu.")
    plt.show()

# Remove or comment out any code at the top level that runs the simulation/plot!
# if __name__ == "__main__":
#     run_and_plot_disease_simulation()