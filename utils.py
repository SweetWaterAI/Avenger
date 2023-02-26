import requests
import os
from Bio import SeqIO
from Bio.Restriction import SmaI, RestrictionBatch
from Bio.SeqFeature import FeatureLocation,CompoundLocation


def load_gene(filename):
    """
    Load the input file as a Biopython SeqRecord object
    """
    seq_record = None
    file_ext = os.path.splitext(filename)[1].lower()
    if file_ext == ".gb" or file_ext == ".gbk" or file_ext == ".genbank":
        if os.path.exists(f"genes/{filename}"):
            seq_record = SeqIO.read(f"genes/{filename}", "genbank")
        else:
            seq_record = SeqIO.read(filename, "genbank")
    elif file_ext == ".fa" or file_ext == ".fasta":
        if os.path.exists(f"genes/{filename}"):
            seq_record = SeqIO.read(f"genes/{filename}", "fasta")
        else:        
            seq_record = SeqIO.read(filename, "fasta")
    elif os.path.exists(f"genes/{filename}.gb"):
        # be forgiving, maybe the user wants to load by name and not filepath:
        seq_record = SeqIO.read(f"genes/{filename}.gb", "genbank")        
    return seq_record

def get_mutation_data(gene_name):
    """
    Retrieves mutation data for a given gene from an external data source.

    Args:
        gene_name (str): Name of the gene to retrieve mutation data for.

    Returns:
        List of dictionaries, where each dictionary contains information about a single mutation.
    """
    url = f"https://magi.brown.edu/api/gene?hugo_symbol={gene_name}"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        mutations = []

        for mutation in data["mutations"]:
            mutation_data = {
                "amino_acid_change": mutation["amino_acid_change"],
                "chromosome": mutation["chromosome"],
                "start": mutation["start"],
                "stop": mutation["stop"],
                "reference": mutation["reference"],
                "variant": mutation["variant"],
                "classification": mutation["classification"],
            }
            mutations.append(mutation_data)

        return mutations
    else:
        return None