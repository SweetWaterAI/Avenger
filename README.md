# PlasmidCompiler
This project simplifies the process of building plasmids and packaging them up into a single annoated genbank file that can be easily used on benchling or another workbench.  Although this task is commonly perfomred by hand in a browser or desktop application, performing senstive tasks like this with python allows us to be eaxct with our changes and to verify the changes.

### CLI Input
The script accepts the following input parameters:

- --backbone-file: Required. The path to the GenBank file of the plasmid backbone to build from.
- --target-dir: Required. The path to the directory containing the list of target genes to replace.
- --replacement-dir: Required. The path to the directory containing the list of replacement genes to use.
- --crispr-file: Required. The path to the GenBank file containing the desired multiplexed-CRISPR variant.
- --output-file: Required. The path to the final plasmid.


### Example Usage
To use the script, navigate to the directory containing the script and run the following command:

`python3 main.py --backbone-file ./genes/backbone.gb --target-dir ./genes/targets/ --replacement-dir ./genes/replacements/ --crispr-file ./genes/crispr.gb --output-file plasmid_1.gb`

This command assumes CRISPR is being used for a 1:1 replacment of a list of genes. To do this the code above will modify the plasmid GenBank file located at genes/backbone.gb and will add pair the replacements with targets based off of the alphabetical sort order of both directories.

### Reproduceablity / Flexablity
Commands can be saved within a .sh file and re-run to re-build a specific plasmid when needed. A team should be able to build any number of therapies or tests out of the same directory.  The code will take in fasta or genebank files.
