# **Tn3 Transposon/Toxin Finder (Tn3+TA_finder)**

Tn3 Transposon/Toxin Finder (Tn3+TA_finder) is a program for the automatic prediction of transposable elements of the Tn3 family associated with type II toxin and antitoxin pairs in bacteria and archaea. It compares bacterial and archaeal genome sequences to custom Tn3 transposase+resolvase and type II toxin+antitoxin databases. Results are writen in report files and pre-anotated GenBank files to help in subsequent manual curation.

---

## RELEASE

Version 1.0.0 - November 11, 2019.

Available from <https://github.com/danillo-alvarenga/tn3-ta_finder>.

---

## REQUIREMENTS

Tn3+TA_finder runs on Python 3.4 or newer and has been tested on Debian 9, Ubuntu 16.04 and CentOS 6.7. It should work on any modern GNU/Linux distro. Tn3+TA_finder depends on the Biopython library version 1.66 (for Python 3), the BLAST+ suite version 2.2.28 and Prodigal version 2.6.1. Newer versions of these packages are likely to work, but if you encounter problems you may try downgrading to the approved versions.

---

### DATABASES

Tn3+TA_finder depends on Tn3 transposase/resolvase and type II toxin/antitoxin databases. As an alternative to the databases included with this program, it is also possible to use custom databases. To do so, retrieve amino acid sequences for Tn3 transposases and resolvases and write them into a file called `Tn3+R.faa` using headers following the pattern `transposase_$NAME` or `resolvase_$NAME` (where `$NAME` is the name of the element). The same must be done for type II toxins and antitoxins, but using as headers `toxin_from_the_$FAMILY_family_$NAME` or `toxin_from_the_$FAMILY_family_$NAME` (where `$FAMILY` is the toxin family name and `$NAME` is the specific toxin name). Then move both files to a subdirectory called `db` in the directory where the `Tn3+TA_finder.py` script is stored. We recomend getting sequences from ISfinder (<https://www-is.biotoul.fr/>) for transposable elements and TADB (<http://202.120.12.135/TADB2/index.php>) for toxins and antitoxins.

---

## INSTALLATION

For Tn3+TA_finder to work the tblastn and prodigal executables need to be installed and available in your path. Follow installation instructions from distributors of these packages for your particular distro. In Debian and derivative distros, these programs are usually available in the official repositories and can be installed by issuing the following commands in a terminal in superuser mode:  
`apt install ncbi-blast+`  
`apt install prodigal`  
`apt install python3-biopython`  

Biopython can also be installed as superuser by pip3:  
`pip3 install biopython`  

For convenience, after downloading and extracting the latest Tn3+TA_finder release from GitHub you can move the Tn3+TA_finder directory to the desired destination and add it to your path:  
`mv Tn3+TA_finder-1.0.0/ bioinformatics/`  
`cd bioinformatics/Tn3+TA_finder-1.0.0/`  
`echo 'export PATH=$PATH'$(pwd) >> ~/.bashrc`  
`source ~/.bashrc`  

>**Note:** CentOS 6 only ships with the legacy Python 2 version. In case you are using a Red Hat-based distribution in which Python 3 is still unavailable, the latest Python version can be installed from source by issuing the following commands to the terminal as superuser:  
>1. Install the OpenSSL development package: `yum install openssl-devel`  
>2. Download the latest stable Python version (Python 3.5.1 at this moment): `wget https://www.python.org/ftp/python/3.5.1/Python-3.5.1.tgz`  
>3. Extract the package: `tar -xzf Python-3.5.1.tgz`  
>4. Change to the decompressed directory: `cd Python-3.5.1/`  
>5. Configure the build scripts: `./configure`  
>6. Build the executables: `make`  
>7. Install the built binaries: `make install`  

>**Note:** Older or newer versions of these programs as well as other operating systems might also work, but since they have not been tested this is not guaranteed.

---

## USAGE

`Tn3+TA_finder.py [-h] [-v] -f Sequences.fasta [Sequences.fasta ...] [-o Directory] [-g] [-t cores] [-e base pairs] [-p percentage] [-c percentage] [-d base pairs]`  

optional arguments:

`-h`, `--help` | show a help message and exit  
`-v`, `--version` | show version and exit  

`-f Sequences.fasta [Sequences.fasta ...]`, `--file Sequences.fasta [Sequences.fasta ...]` | target sequences  
`-o Directory`, `--out Directory` | output directory  
`-g`, `--gbk` | generate a genbank file with predictions  

`-t cores`, `--threads cores` | number of cores (for multifasta or multiple files)  

`-p percentage`, `--positives percentage` | minimum positives percentage  
`-c percentage`, `--coverage percentage` | minimum alignment coverage  
`-d base pairs`, `--distance base pairs` | maximum distance between ORFs  

`-m`, `--merge` | merge overlapping candidates  
`-e base pairs`, `--extend base pairs` | retrieve sequence with extended borders  

In order to run Tn3+TA_finder, indicate the target genome file in the nucleotide fasta format with the `-f`|`--file` argument and add flags adjusting the default parameters to your preferences. Optionally, indicate an output directory for the analyses results with the `-o`|`--output` argument. Tn3+TA_finder will output an annotated genbank file if you add the `-g`|`--genbank` flag to the command line.

Complete genomes can only be analyzed in a single processor core, but multiple files and diferent contigs from draft genomes in multifasta files can be analyzed concurrently in the number of cores you determine by adding the `-t`|`--threads` parameter followed by the number of threads to be used.

The `-p`|`--positives` parameter determines the minimum percentage of physicochemically similar amino acids between query and database sequences and defaults to 50 %. Likewise, the `-c`|`--coverage` parameter indicates the minimum allowed percentage of alignment size compared to the database sequence size, defaulting to 80 %. Finally, `-d`|`--distance` establishes the farthest the detected features can be located in the query sequence and defaults to 1000 bp. Borders for the candidates found can be extended by adding a number of base pairs with `-e`|`--extend`. Additionally, the `-m`|`--merge` option allows to merge candidates with overlapping coordinates into a single candidate.

**Examples**:  
`Tn3+TA_finder -f Escherichia\_coli\_A123\_complete\_genome.fasta`  
`Tn3+TA_finder -f Escherichia\_coli\_B456\_draft\_genome.fasta Escherichia\_coli\_C789\_draft\_genome.fasta -o results -g -t 4 -p 60 -c 70 -d 3000`  

---

## RESULTS

Tn3+TA_finder generates a report file indicating the query sequence followed by potential candidates, detailing length, order of features and general characteristics, partial sequences based on similarity to databases and predicted ORF coordinates. If you add the `-g`|`--genbank` flag, the program will also output a .gbk file with annotations corresponding to the analysis results which can be manually curated in programs such as the Artemis genome browser.

The report and gbk files include sequences for the predicted transposase, resolvase, toxin and antitoxin (reverse-complemented if necessary) and for the whole sequence, from first to last ORF. If the `-e`|`--extend` parameter was included, an additional sequence that is up and downstream longer than the whole sequence is provided.

Information about how the run was started is stored in the file `info.txt`, which is also written in the working directory.

>**Note:** The report file relies on BLAST alignments for providing feature sequences, which may be incomplete as a result of differences between extremities. Thus, the report file is only meant to allow a quick preview of these sequences. More accurate predictions are found in GenBank files, since they are based on Prodigal ORF detection.

>**Note:** There may be some discrepancy between tBLASTn and Prodigal sequence detection if stop codons are found in the middle of the region with similarity to the database. Consequently, annotation for one or more incomplete/divided features may be missing from the GenBank file.

---

## CITATION

If you find this program useful in your research, please cite the following references:
- Alvarenga DO, Varani AM (2019). Tn3+TA_finder: prokaryotic transposase element recognizer. Available from <https://github.com/danillo-alvarenga/tn3-ta_finder>.
- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421
- Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Freidberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25:1422-1423
- Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11:119
- Siguier P, Perochon J, Lestrade L, Mahillon J, Chandler M (2006) ISfinder: the reference centre for bacterial insertion sequences. Nucleic Acids Res 34:D32-D36
- Shao Y, Harrison EM, Bi D, Tai C, He X, Ou HY, Rajakumar K, Deng Z (2011) TADB: a web-based resource for type 2 toxin-antitoxin loci in bacteria and archaea. Nucleic Acids Res 39:D606-D611

---

## ISSUES AND REQUESTS

If you experience any issues or would like to see support for an additional feature implemented in Tn3+TA_finder, please file a request in the GitHub issue tracker or email it to the developer. Please feel free to contact the developer with any further questions regarding this software. You can reach him at <mailto:danillo.alvarenga@gmail.com>.

---

## LICENSE

Copyright © 2019 Danillo Oliveira Alvarenga

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/agpl-3.0.html>.
