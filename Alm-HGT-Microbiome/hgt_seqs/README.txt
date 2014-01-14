Description of horizontally transferred nucleotide sequence file (hgt.fst).

In order to determine whether a sequence of interest is present in our dataset,
a user can download and locally BLAST the sequence of interest against this
database.

hgt.fst follows the FASTA format, as in the following example: 
>12631 Yersinia_pestis_biovar_Orientalis_str_PEXU2 NZ_ACNS01000004
GTGTCAACGACGGATGAAAAGTGATCCACTTATATCTCCACCAACGGCC

The header line includes a unique numerical sequence identifier (12631),
followed by the genome name (Yersinia_pestis_biovar_Orientalis_str_PEXU2) to
which the sequence belongs and ending with the unique NCBI contig identifier
(NZ_ACNS01000004) on which the sequence was found (within the genome).

Each separate, horizontally transferred region in each genome is included
individually in this file. Thus an HGT between Bacteria_1, Bacteria_2 and
Bacteria_3 would appear three times. 