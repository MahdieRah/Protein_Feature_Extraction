import os
import requests
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Function to fetch protein sequence from UniProt
def fetch_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text
        sequence = "".join(fasta_data.split('\n')[1:])  # Remove FASTA header
        return sequence
    else:
        print(f"❌ Error fetching sequence for {uniprot_id}")
        return None

# Function to predict disorder using IUPred3
def predict_disorder(sequence):
    url = "https://iupred3.elte.hu/iupred3"
    response = requests.post(url, data={"seq": sequence, "type": "long"})
    if response.status_code == 200:
        disorder_scores = [float(line.split()[2]) for line in response.text.split('\n') if line and not line.startswith('#')]
        return sum(disorder_scores) / len(disorder_scores)  # Return average disorder score
    else:
        print("❌ Error fetching disorder scores from IUPred3")
        return None

# Function to extract protein features
def extract_features(sequence):
    analysis = ProteinAnalysis(sequence)
    features = {
        "Length": len(sequence),
        "Molecular_Weight": analysis.molecular_weight(),
        "Isoelectric_Point": analysis.isoelectric_point(),
        "Instability_Index": analysis.instability_index(),
        "GRAVY": analysis.gravy(),
        "Disorder_Score": predict_disorder(sequence)
    }
    return features

# List of UniProt IDs to analyze
uniprot_ids = ["P69905", "P68871"]  # Example: Hemoglobin subunits

data = []
for uniprot_id in uniprot_ids:
    print(f"✅ Fetching sequence for {uniprot_id}...")
    sequence = fetch_sequence(uniprot_id)
    if sequence:
        print(f"✅ Extracting features for {uniprot_id}...")
        features = extract_features(sequence)
        features["UniProt_ID"] = uniprot_id
        data.append(features)

# Save results to CSV
if data:
    df = pd.DataFrame(data)
    df.to_csv("protein_features.csv", index=False)
    print("✅ Analysis completed! Results saved in protein_features.csv.")
else:
    print("❌ No valid protein sequences found.")

