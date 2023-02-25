import os
import typer
import subprocess
import requests
from Bio import SeqIO, Entrez
from Bio.Restriction import SmaI, RestrictionBatch
from Bio.SeqFeature import FeatureLocation, CompoundLocation

app = typer.Typer()

def load_gene(filename):
    """
    Load the input file as a Biopython SeqRecord object
    """
    file_ext = os.path.splitext(filename)[1].lower()
    if file_ext == ".gb" or file_ext == ".gbk" or file_ext == ".genbank":
        seq_record = SeqIO.read(filename, "genbank")
    elif file_ext == ".fa" or file_ext == ".fasta":
        seq_record = SeqIO.read(filename, "fasta")
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    return seq_record

def generate_sgRNA_targets(genome_file: str, target_dir: str, output_file: str):
    # Load genome file
    genome_record = load_gene(genome_file)
    genome_seq = genome_record.seq
    
    # Load target files
    targets = []
    for filename in os.listdir(target_dir):
        if filename.endswith(".gb") or filename.endswith(".gbk"):
            target_record = load_gene(os.path.join(target_dir, filename))
            targets.append(target_record)
    
    # Generate potential sgRNA targets
    rb = RestrictionBatch([SmaI])
    potential_targets = []
    for target in targets:
        for feature in target.features:
            if feature.type == "CDS":
                feature_seq = feature.extract(genome_seq)
                for i in range(len(feature_seq)-20):
                    if rb.search(feature_seq[i:i+23]):
                        target_location = feature.location.parts[0].start + i
                        potential_targets.append(target_location)
                        
    # Write targets to output file
    with open(output_file, "w") as f:
        for target in potential_targets:
            f.write(f"{genome_record.id},{target}\n")

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

@app.command()
def download_gene_data(gene_name: str, email: str = typer.Option(None, help="Email address for NCBI queries")):
    """
    Downloads gene data for a given gene name from NCBI and saves all fields as annotations in a genbank file.
    """
    # If email is not provided, use a hard-coded email address
    if email is None:
        email = "your.email@example.com"
    
    # Set the email address for Entrez queries
    Entrez.email = email
    
    # Search for the gene name in NCBI's nucleotide database
    handle = Entrez.esearch(db="nucleotide", term=gene_name)
    record = Entrez.read(handle)
    handle.close()

    # If at least one record is found, download the first one and save all fields as annotations
    if record["IdList"]:
        handle = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="gb", retmode="text")
        gene_record = SeqIO.read(handle, "gb")
        handle.close()

        # Retrieve mutation data from MAGI - todo MAGI down.
        mutation_data = None #get_mutation_data(gene_name)

        # If there is mutation data, add it as an annotation to the gene record
        if mutation_data:
            for feature in gene_record.features:
                if feature.type == "gene":
                    feature.qualifiers["mutation_data"] = str(mutation_data)
                    break

        # Add annotations for gene name, accession, and taxonomy
        gene_record.annotations["gene_name"] = gene_name
        gene_record.annotations["source"] = "NCBI"
        check_accession = gene_record.id.split(":")
        if len(check_accession)>=2:
            gene_record.annotations["accession"] = check_accession[1]
        
        # Save the gene record as a genbank file with the gene name as the filename
        SeqIO.write(gene_record, f"{gene_name}.gb", "gb")
        typer.echo(f"Gene data for {gene_name} downloaded and saved to {gene_name}.gb")
    else:
        typer.echo(f"No records found for {gene_name} in NCBI's nucleotide database")

@app.command()
def find_mutations(reference_file: str = typer.Argument(..., help="Path to the reference genome FASTA file."),
                   gene_name: str = typer.Argument(..., help="Name of the gene of interest."),
                   output_file: str = typer.Argument(..., help="Path to the output VCF file.")):
    """
    Use GATK4 to find all mutations of a given gene within a genome.
    """
    # Use Conda to activate the GATK4 environment
    cmd = 'conda activate gatk4'
    subprocess.run(cmd, shell=True)

    # Use Biopython to get the gene coordinates
    reference_seq = SeqIO.read(reference_file, 'fasta')
    gene_seq = reference_seq.seq
    for feature in reference_seq.features:
        if feature.type == 'gene' and feature.qualifiers['gene'][0] == gene_name:
            start = feature.location.start.position
            end = feature.location.end.position
            gene_seq = feature.extract(reference_seq.seq)

    # Call variants using GATK4 HaplotypeCaller
    tmp_vcf = 'tmp.vcf'
    cmd = f'gatk HaplotypeCaller -R {reference_file} -I {reference_file} -L {gene_name} -O {tmp_vcf}'
    subprocess.run(cmd, shell=True)

    # Filter variants to retain only those overlapping with the specified gene
    cmd = f'bcftools view {tmp_vcf} -T {start}:{end} -O v -o {output_file}'
    subprocess.run(cmd, shell=True)

    # Clean up temporary files
    os.remove(tmp_vcf)


@app.command()
def modify_plasmid(
    backbone_file: str = typer.Option(..., "--backbone-file", "-b", help="The path to the GenBank file of the plasmid backbone to build from."),
    target_dir: str = typer.Option(..., "--target-dir", "-t", help="The path to the directory containing the list of target genes to replace."),
    replacement_dir: str = typer.Option(..., "--replacement-dir", "-r", help="The path to the directory containing the list of replacement genes to use."),
    crispr_file: str = typer.Option(..., "--crispr-file", "-c", help="The path to the GenBank file containing the desired multiplexed-CRISPR variant."),
    output_file: str = typer.Option(..., "--output-file", "-o", help="The path to the resulting plasmid.")
):
    """
    Modify a GenBank file of a plasmid backbone by replacing specified target genes with new sequences.

    Parameters:
        backbone_file (str): The path to the GenBank file of the plasmid backbone to build from.
        target_dir (str): The path to the directory containing the list of target genes to replace.
        replacement_dir (str): The path to the directory containing the list of replacement genes to use.
        crispr_file (str): The path to the GenBank file containing the desired multiplexed-CRISPR variant.
        output_file (str): The path to the resulting plasmid.
    """
    # Parse the CRISPR GenBank file to get the target sequence and PAM
    crispr_record = load_gene(crispr_file)
    target_sequence = crispr_record.features[0].qualifiers["target"][0]
    pam_sequence = crispr_record.features[0].qualifiers["PAM"][0]

    # Parse the GenBank file of the plasmid backbone
    backbone_record = load_gene(backbone_file)
    
    # Loop through each target gene and its replacement
    for target_file, replacement_file in zip(os.listdir(target_dir), os.listdir(replacement_dir)):
        # Parse the target gene and its replacement
        target_record = load_gene(os.path.join(target_dir, target_file))
        replacement_record = load_gene(os.path.join(replacement_dir, replacement_file))
        
        # Find the location of the target gene in the plasmid backbone
        target_feature = None
        for feature in backbone_record.features:
            if feature.type == "CDS" and feature.qualifiers.get("translation") == target_record.seq.translate():
                target_feature = feature
                break
        
        # Replace the target gene with the replacement sequence
        if target_feature is not None:
            # Find the location of the target sequence in the plasmid backbone
            target_index = backbone_record.seq.find(target_sequence)
            if target_index == -1:
                raise ValueError(f"Could not find target sequence {target_sequence} in plasmid backbone.")
            
            # Find the location of the PAM in the plasmid backbone
            pam_index = backbone_record.seq.find(pam_sequence, target_index+len(target_sequence))
            if pam_index == -1:
                raise ValueError(f"Could not find PAM sequence {pam_sequence} downstream of target sequence {target_sequence} in plasmid backbone.")
            
            # Replace the target sequence with the replacement sequence
            replacement_feature = SeqFeature(
                FeatureLocation(target_index, pam_index+len(pam_sequence)),
                type="CDS",
                qualifiers={
                    "translation": str(replacement_record.seq.translate()),
                    "product": target_feature.qualifiers["product"],
                    "note": target_feature.qualifiers["note"],
                },
            )
            backbone_record.features.remove(target_feature)
            backbone_record.features.append(replacement_feature)
        else:
            raise ValueError(f"Could not find target gene {target_record.id} in plasmid backbone.")
    
    # Write the modified plasmid GenBank file to disk
    with open(output_file, "w") as f:
        SeqIO.write(backbone_record, f, "genbank")

if __name__ == "__main__":
    app()
