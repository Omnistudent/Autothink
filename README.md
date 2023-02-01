Autothink

Finds probable Insertion sequences (ISs) in genebank files
Requirements
Biopython
Standalone blast

In same folder as autothink script:
classholder.py
Genomeholder.py
helperobject.py
Running the Script

Takes about 10 minutes for example bmg5.1 file.

Running

At command line prompt on windows>

Python [scriptlocation]\autothink2_22jan.py --genome_directory A --outputdir B --files_to_check_for_is C --isblastdbfile D --blast_location E

Required parameters
--genome_directory	
Folder with genbank files to be treated. One sequence in one genbank file.

Example: 

C:\Users\Eris\Documents\autothinktestfolder\frankiatestgenomes\

A folder containing one genbank file, a shortened version of Frankia bmg51.gb 

--outputdir	
Directory for most of the output files and work files. A folder is made for each inputted genbank file here. The folder contains: 

(in folder X_curated) one genbank file, X.gbk_cleaned.gbk, with locations and classifications of insertion sequences.

One csv file, X.gb_firsthits_removed.csv with detailed info about the insertion sequences.

Example:
C:\Users\Eris\Documents\autothinktestfolder\outputdir\

--files_to_check_for_is	
The name of the file generated when the script looks at the genbank file.

Example:
C:\Users\Eris\Documents\autothinktestfolder\files_to_check_for_is.csv

--isblastdbfile 	
The location of the amino acid blast database containing known insertion sequences. This database is uploaded in the folder aadatabase


Example:
C:\Users\Eris\Documents\scripts\autothink\is_aa_30_nov2016.fa

This means three files named is_aa_30_nov2016.phr, .pin and .psq should be placed in the C:\Users\Eris\Documents\scripts\autothink folder.


--blast_location 	
Directory with standalone blast program 

Example:
C:\NCBI\blast-BLAST_VERSION+\bin\

Real line for running on my computer

C:\Users\Eris\Documents>python C:\Users\Eris\Documents\autothinktestfolder\autothink2_22jan.py --genome_directory C:\Users\Eris\Documents\autothinktestfolder\frankiatestgenomes\ --outputdir C:\Users\Eris\Documents\autothinktestfolder\outputdir\ --files_to_check_for_is C:\Users\Eris\Documents\autothinktestfolder\files_to_check_for_is.csv --isblastdbfile C:\Users\Eris\Documents\scripts\autothink\is_aa_30_nov2016.fa --blast_location C:\NCBI\blast-BLAST_VERSION+\bin\

As one line.


About the code
The script is very unorganized and contains parts that are in development, e.g. support for folders in input folder, building of databases.

In very short:

Goes through folder of genbank files and makes of list of files to treat
Blastxs file against a database of IS amino acids.
Curates the results, determining most likely IS sequence on nucleotide sequence.
Writes the information to a new genbank file.

The “summary” files should be double checked, headers and computed info may not be correct.
