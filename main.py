import os
import typer
import subprocess
from utils import load_gene 
import luaparser
import BIO
from Bio import SeqIO, Entrez, Seq, SeqRecord, SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation 
# todo DeepCRISPR needs to be ported to tf2
# from DeepCRISPR.deepcrispr import DCModelOntar
# todo wrap this into a better include:
# from sgrna_design import build_sgrna_library

app = typer.Typer( help="Build the future you choose." )

@app.command()
def get(gene_name: str, email: str = typer.Option(None, help="Email address for NCBI queries")):
    """
    Downloads gene data for a given gene name from NCBI and saves all fields as annotations in a genbank file.

    Parameters:
        gene_name (str): The name of the gene to download.
        email (str, optional): The email address to use for NCBI queries. Defaults to None.

    """
    if os.path.exists(f"genes/{gene_name}.gb"):
        typer.echo(f"found genes/{gene_name}.gb")
        return

    # If email is not provided, use a hard-coded email address
    if email is None:
        email = "contact@sweetwater.ai"
    
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
        SeqIO.write(gene_record, f"genes/{gene_name}.gb", "gb")
        typer.echo(f"Gene data for {gene_name} downloaded and saved to genes/{gene_name}.gb")
    else:
        typer.echo(f"No records found for {gene_name} in NCBI's nucleotide database")

@app.command()
def find(reference_file: str = typer.Argument(..., help="Path to the reference genome file."),
                   gene_name: str = typer.Argument(..., help="Name of the gene of interest."),
                   output_file: str = typer.Argument(..., help="Path to the output VCF file.")):
    """
    Use GATK4 to find all mutations of a given gene within a genome.

    Parameters:
        reference_file (str): Path to the reference genome file in FASTA format.
        gene_name (str): Name of the gene of interest.
        output_file (str): Path to the output VCF file.    
    """
    # Use Conda to activate the GATK4 environment
    cmd = 'conda activate gatk4'
    subprocess.run(cmd, shell=True)

    # Use Biopython to get the gene coordinates
    reference_seq = load_gene(reference_file)
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
def predict(
    gsRNA_candidates: str = typer.Argument(..., help="Path to a file containing gRNA sequences"),
    targets: str = typer.Argument(..., help="Path to a file containing target sequences"),
    genome: str = typer.Argument(..., help="Path to the full genome sequence"),
    output: str = typer.Option("ontar_results.tsv", help="Path to the output file"),
    model_path: str = typer.Option(..., help="Path to the DCModelOntar model file"),
    pam: str = typer.Option("NGG", help="PAM sequence used by the gRNA"),
):
    """
    Predict off-targets for a set of gRNA candidates and target sequences using the DCModelOntar model.

    Parameters:
        gsRNA_candidates (str): Path to a file containing gRNA sequences.
        targets (str): Path to a file containing target sequences.
        genome (str): Path to the full genome sequence.
        output (str): Path to the output file.
        model_path (str): Path to the DCModelOntar model file.
        pam (str): PAM sequence used by the gRNA. 
    """
    # Load the model
    model = DCModelOntar.load(model_path)

    # Read in the input files
    with open(gsRNA_candidates, "r") as f:
        gRNA_seqs = [line.strip() for line in f if line.strip()]

    with open(targets, "r") as f:
        target_seqs = [line.strip() for line in f if line.strip()]

    with open(genome, "r") as f:
        genome_seq = f.read().strip()

    # Run the prediction
    results = model.ontar_predict(gRNA_seqs, target_seqs, genome_seq, pam)

    # Write the results to a file
    with open(output, "w") as f:
        f.write("gRNA_sequence\ttarget_sequence\toff_target_count\toff_targets\n")
        for result in results:
            gRNA_seq, target_seq, off_target_count, off_targets = result
            off_target_str = ",".join(off_targets)
            f.write(f"{gRNA_seq}\t{target_seq}\t{off_target_count}\t{off_target_str}\n")

@app.command()
def design(
    target_genes: str = typer.Argument(..., help="List of target genes"),
    mutations: str = typer.Argument(..., help="List of mutations to introduce"),
    full_genome: str = typer.Argument(..., help="Full genome sequence"),
):
    """
    Design sgRNA candidates for given target genes and mutations using build_sgrna_library.

    Parameters:
        target_genes: List of target genes.
        mutations: List of mutations to introduce.
        full_genome: Full genome sequence.

    """
    # Generate sgRNA candidates for given target genes and mutations
    sgRNA_candidates = build_sgrna_library(target_genes, mutations, full_genome)

    # Write sgRNA candidates to sgRNA folder with unique names based on the target genes
    for gene, sgrna_list in sgRNA_candidates.items():
        with open(f"sgRNA/{gene}.txt", "w") as f:
            for sgrna in sgrna_list:
                f.write(sgrna + "\n")

    typer.echo("sgRNA candidates have been generated and written to the sgRNA folder.")

@app.command()
def build_dna(
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

@app.command(help="Build a gsRNA with a human cap 5' and tail as an annotated GenBank file")
def build_rna(
    pam: str = typer.Argument(..., help="The PAM sequence"),
    output_file: str = typer.Argument(..., help="The output file name"),
    poly_a_length: int = typer.Option(120, "-a", "--poly_a_length", help="The length of the poly-A' tail")
):
    """
    Build a gsRNA with a human cap 5' and tail as an annotated GenBank file.

    The PAM sequence must be specified as an argument, and the output file name
    must also be specified as an argument. The length of the poly-A' tail can be
    controlled with the --poly_a_length option (default is 120).

    Example:
    $ gsRNA_builder build NGG gsRNA.gb -a 20
    """

    # Transcribe the PAM sequence to RNA
    pam_rna = Seq(pam).reverse_complement().transcribe()

    # Define the cap and tail sequences
    cap_seq = "GAATTCGCGGCCGCTTCTAG"
    tail_seq = "A" * poly_a_length

    # Create the gsRNA sequence
    gsRNA_seq = str(cap_seq + str(pam_rna) + tail_seq)

    # Create a SeqRecord object for the gsRNA sequence
    gsRNA_record = SeqRecord(Seq(gsRNA_seq), id="gsRNA", name="gsRNA",
                             description="Guided RNA with human cap and tail")

    # Add an annotated feature for the PAM sequence
    pam_feature = SeqFeature(
        SeqFeature.FeatureLocation(
            start=len(cap_seq),
            end=len(cap_seq)+len(pam),
            strand=+1
        ),
        type="PAM",
        id="pam_feature",
        qualifiers={"note": f"PAM sequence: {pam}"}
    )
    gsRNA_record.features.append(pam_feature)

    # Add an annotated feature for the cap sequence
    cap_feature = SeqFeature(
        SeqFeature.FeatureLocation(
            start=0,
            end=len(cap_seq),
            strand=+1
        ),
        type="misc_feature",
        id="cap_feature",
        qualifiers={"note": "human cap sequence"}
    )
    gsRNA_record.features.append(cap_feature)

    # Add an annotated feature for the tail sequence
    tail_feature = SeqFeature(
        SeqFeature.FeatureLocation(
            start=len(cap_seq)+len(pam),
            end=len(gsRNA_seq),
            strand=+1
        ),
        type="misc_feature",
        id="tail_feature",
        qualifiers={"note": f"poly-A' tail length: {poly_a_length}"}
    )
    gsRNA_record.features.append(tail_feature)

    # Write the gsRNA sequence to a GenBank file
    with open(output_file, "w") as f:
        SeqIO.write(gsRNA_record, f, "genbank")

if __name__ == "__main__":
    app()
