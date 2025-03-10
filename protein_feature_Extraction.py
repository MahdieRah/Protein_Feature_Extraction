import os

# Set working directory
os.chdir(r"D:\biopython\protein_feature-Extraction")

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import requests
import pandas as pd

def fetch_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text.split("\n", 1)[-1].replace("\n", "")
        return fasta_data
    else:
        print(f"❌ Failed to fetch sequence for {uniprot_id}")
        return None

def compute_features(sequence):
    """Compute physicochemical properties of a protein sequence"""
    analysis = ProteinAnalysis(sequence)
    features = {
        "Molecular Weight": analysis.molecular_weight(),
        "Isoelectric Point": analysis.isoelectric_point(),
        "Aromaticity": analysis.aromaticity(),
        "Instability Index": analysis.instability_index(),
        "Helix Content": analysis.secondary_structure_fraction()[0],
        "Turn Content": analysis.secondary_structure_fraction()[1],
        "Sheet Content": analysis.secondary_structure_fraction()[2]
    }
    return features

def save_results(data, filename="protein_features.csv"):
    """Save extracted features to a CSV file"""
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"✅ Results saved in {filename}")

if __name__ == "__main__":
    protein_ids = ["P69905", "P68871"]  # Example UniProt IDs (Hemoglobin subunits)
    results = []
    
    for protein_id in protein_ids:
        sequence = fetch_sequence(protein_id)
        if sequence:
            features = compute_features(sequence)
            features["UniProt ID"] = protein_id
            results.append(features)
    
    save_results(results)
