# Import the necessary modules
from luaparser import ast
import lupa
import os
import Bio
from Bio import SeqIO, SeqFeature, Entrez
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import Seq


class Avenger:
    # Define Avenger's functionality for genetic engineering and gene editing operations
    def read(self, file_path):
        # Read genbank file and extract sequence
        genbank = SeqIO.read(file_path, 'genbank')
        sequence = str(genbank.seq)
        return genbank, sequence

    def delete(self, sequence, start, stop):
        # Delete sequence between start and stop positions
        sequence = sequence[:start] + sequence[stop:]
        return sequence

    def insert(self, sequence, position, to_insert):
        # Insert the to_insert sequence at specified position
        sequence = sequence[:position] + to_insert + sequence[position:]
        return sequence

    def replace(self, sequence, start, stop, replacement):
        # Replace sequence between start and stop positions with replacement sequence
        sequence = sequence[:start] + replacement + sequence[stop:]
        return sequence

    def create_sequence(self, seq, alphabet):
        # Create a new sequence object
        return Seq(seq, alphabet)

    def create_feature(self, feature_type, location, qualifiers):
        # Create a new feature object
        return SeqFeature(location, type=feature_type, qualifiers=qualifiers)

    def create_location(self, start, stop, strand=None):
        # Create a new location object
        return FeatureLocation(start, stop, strand=strand)

    def create_compound_location(self, locations):
        # Create a new compound location object
        return FeatureLocation(locations[0].start, locations[-1].end, parts=locations)

    def create_qualifier(self, qualifier_type, value):
        # Create a new qualifier object
        return {qualifier_type: [value]}

    def write(self, file_path, sequence, format):
        # Write sequence to file in specified format
        SeqIO.write(self.create_sequence(sequence, "generic"), file_path, format)

    
    def get(self, gene_name, email = None):
        """
        Downloads gene data for a given gene name from NCBI and returns the genetic code.

        Parameters:
            gene_name (str): The name of the gene to download.
            email (str, optional): The email address to use for NCBI queries. Defaults to None.

        Returns:
            int: The genetic code of the downloaded gene.
        """
        if os.path.exists(f"genes/{gene_name}.gb"):
            print(f"found genes/{gene_name}.gb")
        else:
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
                print(f"Gene data for {gene_name} downloaded and saved to genes/{gene_name}.gb")

        # Read the gene record from the genbank file
        gene_record = SeqIO.read(f"genes/{gene_name}.gb", "gb")

        # Extract the genetic code from the gene record
        genetic_code = gene_record.annotations["codon_start"]

        # Return the genetic code
        return genetic_code

def load_avenger(script_path):
    # Load the Lua script
    with open(script_path) as f:
        lua_script = f.read()

    # Create a Lua interpreter
    lua_engine = lupa.LuaRuntime()

    # Register the Avenger methods in the Lua environment
    avenger = Avenger()
    for name in dir(avenger):
        method = getattr(avenger, name)
        if callable(method):
            lua_engine.globals()[name] = method

    # Execute the Lua script in the Lua environment
    lua_engine.execute(lua_script)

# Use the parsed code to generate a final genbank file
if __name__ == '__main__':
    load_avenger("test.ave")
