import os
import shutil
import tempfile
import unittest
import urllib.request

from Bio import SeqIO

from main import modify_plasmid


class TestModifyPlasmid(unittest.TestCase):

    def setUp(self):
        # Set up test directories and files
        self.test_dir = tempfile.TemporaryDirectory()
        self.backbone_file = os.path.join(self.test_dir.name, 'backbone.gb')
        self.target_dir = os.path.join(self.test_dir.name, 'targets')
        self.replacement_dir = os.path.join(self.test_dir.name, 'replacements')
        self.crispr_file = os.path.join(self.test_dir.name, 'crispr.gb')
        self.output_file = os.path.join(self.test_dir.name, 'output.gb')
        os.makedirs(self.target_dir)
        os.makedirs(self.replacement_dir)

        # Download real files from a different source
        urllib.request.urlretrieve('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_001416.1&db=nuccore&report=genbank&extrafeat=null&conwithfeat=on&hide-cdd=on', self.backbone_file)
        urllib.request.urlretrieve('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MT263387.1&db=nuccore&report=genbank&extrafeat=null&conwithfeat=on&hide-cdd=on', self.crispr_file)
        urllib.request.urlretrieve('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NG_046901.1&db=nuccore&report=genbank&extrafeat=null&conwithfeat=on&hide-cdd=on', os.path.join(self.target_dir, 'target_gene.gb'))
        urllib.request.urlretrieve('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NG_054906.1&db=nuccore&report=genbank&extrafeat=null&conwithfeat=on&hide-cdd=on', os.path.join(self.replacement_dir, 'replacement_gene.gb'))

    def tearDown(self):
        # Clean up test directories
        self.test_dir.cleanup()

    def test_modify_plasmid(self):
        # Run the modify_plasmid function
        modify_plasmid(self.backbone_file, self.target_dir, self.replacement_dir, self.crispr_file, self.output_file)

        # Check that the output file was created
        self.assertTrue(os.path.isfile(self.output_file))

        # Check that the output file has the expected content
        output_record = SeqIO.read(self.output_file, 'genbank')
        self.assertEqual(len(output_record.features), 92)  # Expected number of features
        self.assertEqual(output_record.annotations['accessions'], ['NC_001416'])  # Expected accession ID

        # Clean up downloaded files
        os.remove(self.backbone_file)
        os.remove(self.crispr_file)
        shutil.rmtree(self.target_dir)
        shutil.rmtree(self.replacement_dir)


if __name__ == '__main__':
    unittest.main()
