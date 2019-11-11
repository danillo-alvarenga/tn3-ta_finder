#!/usr/bin/env python3
#
# Tn3+TA_finder: Tn3 Transposon/Toxin Finder
#
# Version 1.0.0 - Novemeber 11, 2019
#
# Copyright © 2019 Danillo Oliveira Alvarenga
#
# Tn3+TA_finder is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# Tn3+TA_finder is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
# License for more details.
#
# You should have received a copy of the GNU Affero GPL along with 
# Tn3+TA_finder. If not, see <http://www.gnu.org/licenses/agpl-3.0.html>.
#
import argparse
import os
import sys
import subprocess
from itertools import combinations
from multiprocessing import Pool
from operator import itemgetter
from shutil import rmtree
from time import strftime
from Bio import SeqFeature
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Get arguments from the command line.
parser = argparse.ArgumentParser(description="Tn3 transposon/toxin finder",
                                 formatter_class=lambda prog: 
                                 argparse.HelpFormatter(prog,
                                 max_help_position=10, width=100))

parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0", help="show version and exit")

parser.add_argument("-f", "--file", metavar="Sequences.fasta", required=True,
                    nargs='+', help="target sequences")
parser.add_argument("-o", "--out", metavar="Directory",
                    help="output directory")
parser.add_argument("-g", "--gbk", action="store_true",
                    help="generate a genbank file with predictions")

parser.add_argument("-t", "--threads", metavar="cores", type=int,
                    help="number of processor cores available for analyses")

parser.add_argument("-p", "--positives", metavar="percentage", type=float,
                    help = "minimum positives percentage")
parser.add_argument("-c", "--coverage", metavar="percentage", type=int,
                    help = "minimum alignment coverage")
parser.add_argument("-d", "--distance", metavar="base pairs", type=int,
                    help = "maximum distance between ORFs")

parser.add_argument("-m", "--merge", action="store_true",
                    help = "merge overlapping candidates")
parser.add_argument("-e", "--extend", metavar="base pairs", type=int,
                    help = "retrieve sequence with extended borders")

parser.set_defaults(out=os.getcwd(), threads=1, extend=0,
                    positives=40, coverage=60, distance=2000)

args = parser.parse_args()


# Run tBLASTn on databases according to input type and generate result lists.
def run_tblastn(fasta):

    # Skip this step if previous results are available.
    filename = fasta.split('/')[-1].rsplit('.', 1)[0]
    if "tblastn" in os.listdir(args.out) and \
        filename + ".tblastn" in os.listdir(args.out + "/tblastn/"):
            return

    db = os.path.dirname(os.path.realpath(__file__)) + "/db/"

    Tn3 = subprocess.check_output(
            ["tblastn", "-query", db+"Tn3+R.faa", "-subject", fasta,
             "-outfmt", "6 qseqid sstart send pident length qlen ppos",
             "-seg", "no", "-word_size", "7", "-evalue", "1e-5"],
             universal_newlines=True)
    TA = subprocess.check_output(
            ["tblastn", "-query", db+"T+A.faa", "-subject", fasta,
             "-outfmt", "6 qseqid sstart send pident length qlen ppos",
             "-seg", "no", "-word_size", "7", "-evalue", "1e-5"],
             universal_newlines=True)

    # Generate a list with both results.
    Tn3 = Tn3.split("\n")
    TA = TA.split("\n")
    results = [x for x in Tn3 + TA if x]

    return results


# Verify if any detected ORFs are within the parameters specified.
def find_candidates(results, distance, positives, length):

    transposases = []
    resolvases = []
    toxins = []
    antitoxins = []
    paired_transposases_resolvases = []
    paired_toxins_antitoxins = []
    candidates = []

    for item in results:
        item = item.rstrip().split("\t")
        category = item[0].split('_')[0]
        beginning = int(item[1])
        end = int(item[2])
        alignment_identity = float(item[3])
        alignment_coverage = round(int(item[4]) * 100 / int(item[5]), 1)
        positive_percentage = float(item[6])

        if positive_percentage < positives:
            continue
        elif alignment_coverage < length:
            continue

        if category == "transposase":
            info = item[0].split('_')[1]
            items = [(x[0], x[1]) for x in transposases]
            if (beginning, end) in items:
                continue
            else:
                transposases.append((beginning, end, info,
                                     positive_percentage, alignment_coverage))
        elif category == "resolvase":
            info = item[0].split('_')[1]
            items = [(x[0], x[1]) for x in resolvases]
            if (beginning, end) in items:
                continue
            else:
                resolvases.append((beginning, end, info,
                                   positive_percentage, alignment_coverage))
        elif category == "toxin":
            info = item[0].split('_')[3]
            items = [(x[0], x[1]) for x in toxins]
            if (beginning, end) in items:
                continue
            else:
                toxins.append((beginning, end, info,
                               positive_percentage, alignment_coverage))
        elif category == "antitoxin":
            info = item[0].split('_')[3]
            items = [(x[0], x[1]) for x in antitoxins]
            if (beginning, end) in items:
                continue
            else:
                antitoxins.append((beginning, end, info,
                                   positive_percentage, alignment_coverage))

    transposases = unduplicate(transposases)
    resolvaves = unduplicate(resolvases)
    toxins = unduplicate(toxins)
    antitoxins = unduplicate(antitoxins)

    for transposase in transposases:
        for resolvase in resolvases:
            spaces = (abs(transposase[0] - resolvase[0]),
                      abs(transposase[0] - resolvase[1]),
                      abs(transposase[1] - resolvase[0]),
                      abs(transposase[1] - resolvase[1]))
            if any(x <= distance for x in spaces):
                paired_transposases_resolvases.append((transposase, resolvase))

    for toxin in toxins:
        for antitoxin in antitoxins:
            spaces = (abs(toxin[0] - antitoxin[0]),
                      abs(toxin[0] - antitoxin[1]),
                      abs(toxin[1] - antitoxin[0]),
                      abs(toxin[1] - antitoxin[1]))
            if any(x <= distance for x in spaces):
                paired_toxins_antitoxins.append((toxin, antitoxin))

    for TR in paired_transposases_resolvases:
        transposase = TR[0]
        resolvase = TR[1]
        for TA in paired_toxins_antitoxins:
            toxin = TA[0]
            antitoxin = TA[1]
            spaces = (abs(transposase[0] - toxin[0]),
                      abs(transposase[0] - toxin[1]),
                      abs(transposase[1] - toxin[0]),
                      abs(transposase[1] - toxin[1]),
                      abs(transposase[0] - antitoxin[0]),
                      abs(transposase[0] - antitoxin[1]),
                      abs(transposase[1] - antitoxin[0]),
                      abs(transposase[1] - antitoxin[1]),
                      abs(resolvase[0] - toxin[0]),
                      abs(resolvase[0] - toxin[1]),
                      abs(resolvase[1] - toxin[0]),
                      abs(resolvase[1] - toxin[1]),
                      abs(resolvase[0] - antitoxin[0]),
                      abs(resolvase[0] - antitoxin[1]),
                      abs(resolvase[1] - antitoxin[0]),
                      abs(resolvase[1] - antitoxin[1]))
            if any(x <= distance for x in spaces):
                candidates.append((transposase, resolvase, toxin, antitoxin))

    candidates = eliminate_overlap(candidates)

    return candidates


# Remove overlapping features of lower identity.
def unduplicate(features):

    for x, y in combinations(features, 2):
        if x in features and (x[0] - 100 < y[0] < x[0] + 100 or \
                              x[1] - 100 < y[1] < x[1] + 100) and \
                              x[3] < y[3]:
            features.remove(x)
        if y in features and (y[0] - 100 < x[0] < y[0] + 100 or \
                              y[1] - 100 < x[1] < y[1] + 100) and \
                              y[3] < x[3]:
            features.remove(y)

    return features


# Ignore candidates accusing overlapping features.
def eliminate_overlap(candidates):

    non_overlapping_candidates = []

    for item in candidates:
        if ((item[0][0] - 50 < item[1][0] < item[0][0] + 50) or \
            (item[1][0] - 50 < item[0][0] < item[1][0] + 50)) and \
           ((item[0][1] - 50 < item[1][1] < item[0][1] + 50) or \
            (item[1][1] - 50 < item[0][1] < item[1][1] + 50)):
            continue
        elif ((item[0][0] < item[2][0] < item[0][1]) or \
              (item[0][1] < item[2][0] < item[0][0])) and \
             ((item[0][0] < item[2][1] < item[0][1]) or \
              (item[0][1] < item[2][1] < item[0][0])):
            continue
        elif ((item[0][0] < item[3][0] < item[0][1]) or \
              (item[0][1] < item[3][0] < item[0][0])) and \
             ((item[0][0] < item[3][1] < item[0][1]) or \
              (item[0][1] < item[3][1] < item[0][0])):
            continue
        elif ((item[1][0] < item[2][0] < item[1][1]) or \
              (item[1][1] < item[2][0] < item[1][0])) and \
             ((item[1][0] < item[2][1] < item[1][1]) or \
              (item[1][1] < item[2][1] < item[1][0])):
            continue
        elif ((item[1][0] < item[3][0] < item[1][1]) or \
              (item[1][1] < item[3][0] < item[1][0])) and \
             ((item[1][0] < item[3][1] < item[1][1]) or \
              (item[1][1] < item[3][1] < item[1][0])):
            continue
        else:
            non_overlapping_candidates.append(item)

    return non_overlapping_candidates


# Check in which order features are found.
def check_order(candidate, strand):

    transposase = (sorted((candidate[0], candidate[1]), key=int), strand[0])
    resolvase = (sorted((candidate[2], candidate[3]), key=int), strand[1])
    toxin = (sorted((candidate[4], candidate[5]), key=int), strand[2])
    antitoxin = (sorted((candidate[6], candidate[7]), key=int), strand[3])

    candidate = [transposase, resolvase, toxin, antitoxin]
    candidate = sorted(candidate, key=itemgetter(0))
    order = [x[1] for x in candidate]
    order = order[0] + " " + order[1] + " " + order[2] + " " + order[3]

    return (transposase[0], resolvase[0], toxin[0], antitoxin[0]), order


# Produce the reverse complement of a DNA sequence.
def reverse_complement(sequence):

    sequence = sequence[::-1]
    reverse_complemented_sequence = ''

    for character in sequence:
        if character is 'T':
            character = 'A'
        elif character is 'A':
            character = 'T'
        elif character is 'C':
            character = 'G'
        elif character is 'G':
            character = 'C'
        elif character is 'N':
            character = 'N'
        reverse_complemented_sequence += character

    return reverse_complemented_sequence


# Identify and retrieve candidate sequences.
def fetch_sequences(sequence, order, strands):

    transposase_sequence = sequence[order[0][0]:order[0][1]]
    if strands[0] is '-':
        transposase_sequence = reverse_complement(transposase_sequence)
    resolvase_sequence = sequence[order[1][0]:order[1][1]]
    if strands[1] is '-':
        resolvase_sequence = reverse_complement(resolvase_sequence)
    toxin_sequence = sequence[order[2][0]:order[2][1]]
    if strands[2] is '-':
        toxin_sequence = reverse_complement(toxin_sequence)
    antitoxin_sequence = sequence[order[3][0]:order[3][1]]
    if strands[3] is '-':
        antitoxin_sequence = reverse_complement(antitoxin_sequence)

    sequences = (transposase_sequence, resolvase_sequence,
                 toxin_sequence, antitoxin_sequence)

    return sequences


# Use Prodigal for predicting open reading frames in the detected sequence.
def predict_orfs(sequence, start):

    sequence_orfs = []
    sequence = ">sequence\n" + sequence
    prodigal = subprocess.check_output(
                            ["prodigal", "-f", "sco", "-p", "meta"],
                             input=sequence, universal_newlines=True,
                             stderr=subprocess.DEVNULL)
    prodigal = [x.lstrip(">") for x in prodigal.split("\n") if ">" in x]

    for item in prodigal:
        item = item.split("_")
        if item[3] is '+':
            beginning = int(item[1])
            end = int(item[2])
        elif item[3] is '-':
            beginning = int(item[2])
            end = int(item[1])
        sequence_orfs.append((beginning, end))

    complete_orfs = [(x+start, y+start) for (x, y) in sequence_orfs]

    return sequence_orfs, complete_orfs


# Combine BLAST and Prodigal outputs for detecting correct feature ORFs.
def match_candidate_to_orfs(candidate, predicted_orfs, beginning):

    matched = []
    for item in predicted_orfs:
        if abs((item[0]+beginning) - (candidate[0][0]+args.extend)) < 100 and\
           abs((item[1]+beginning) - (candidate[0][1]+args.extend)) < 100:
            matched.append((item[0], item[1], "transposase", candidate[0][2]))
        elif abs((item[0]+beginning) - (candidate[1][0]+args.extend)) < 100 and\
             abs((item[1]+beginning) - (candidate[1][1]+args.extend)) < 100:
            matched.append((item[0], item[1], "resolvase", candidate[1][2]))
        elif abs((item[0]+beginning) - (candidate[2][0]+args.extend)) < 100 and\
             abs((item[1]+beginning) - (candidate[2][1]+args.extend)) < 100:
            matched.append((item[0], item[1], "toxin", candidate[2][2]))
        elif abs((item[0]+beginning) - (candidate[3][0]+args.extend)) < 100 and\
             abs((item[1]+beginning) - (candidate[3][1]+args.extend)) < 100:
            matched.append((item[0], item[1], "antitoxin", candidate[3][2]))
        else:
            matched.append((item[0], item[1], False, None))

    return matched


# Verify features in common between candidates and merge them.
def merge_candidates(candidates):

    merged_candidates = []
    unique_candidates = []
    i = 0
    n = 1

    while i < len(candidates)-1:
        if i in merged_candidates:
            i += 1
            n = i + 1
            continue

        extra_features = []
        while n < len(candidates):
            if candidates[i][0] != candidates[n][0] and \
               candidates[i][1] != candidates[n][1] and \
               candidates[i][2] != candidates[n][2] and \
               candidates[i][3] != candidates[n][3]:
                candidate = list(candidates[i])
                n += 1
            else:
                if (candidates[i][0][0] != candidates[n][0][0]) and \
                   (candidates[i][0][1] != candidates[n][0][1]):
                    extra_feature = (candidates[n][0], "transposase")
                    if extra_feature not in extra_features:
                        extra_features.append(extra_feature)
                        merged_candidates.append(n)
                if (candidates[i][1][0] != candidates[n][1][0]) and \
                   (candidates[i][1][1] != candidates[n][1][1]):
                    extra_feature = (candidates[n][1], "resolvase")
                    if extra_feature not in extra_features:
                        extra_features.append(extra_feature)
                        merged_candidates.append(n)
                if (candidates[i][2][0] != candidates[n][2][0]) and \
                   (candidates[i][2][1] != candidates[n][2][1]):
                    extra_feature = (candidates[n][2], "toxin")
                    if extra_feature not in extra_features:
                        extra_features.append(extra_feature)
                        merged_candidates.append(n)
                if (candidates[i][3][0] != candidates[n][3][0]) and \
                   (candidates[i][3][1] != candidates[n][3][1]):
                    extra_feature = (candidates[n][3], "antitoxin")
                    if extra_feature not in extra_features:
                        extra_features.append(extra_feature)
                        merged_candidates.append(n)
            n +=1
        if extra_features:
            candidate = list(candidates[i]) + extra_features

        unique_candidates.append(candidate)
        i += 1
        n = i + 1

    return unique_candidates


# Use a sequence and matched ORFs to generate a genbank file.
def write_gbk(sequence, matched_orfs, filename, organism):

    date = strftime("%d-%b-%Y").upper()
    orfs = []
    features = []

    gbk_record = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna),
                        description=organism+" predicted Tn3 sequence",
                        annotations={"accession":'.', "version":'.',
                                     "organism":'.', "date":date,
                                     "data_file_division":"BCT"})

    for item in matched_orfs:
        if item[0] < item[1]:
            start = item[0]
            end = item[1]
            strand = 0
        else:
            start = item[1]
            end = item[0]
            strand = -1
        orfs.append((start, end, strand, item[2], item[3]))

    for item in orfs:
        if item[3] and item[4]:
            features.append(SeqFeature.SeqFeature(
                    SeqFeature.FeatureLocation(
                    item[0]-1, item[1], strand=item[2]),
                    type="misc_feature",
                    qualifiers={"note":item[4]+" family "+item[3]}))
        else:
            features.append(SeqFeature.SeqFeature(
                    SeqFeature.FeatureLocation(
                    item[0]-1, item[1], strand=item[2]),
                    type="CDS"))

    for item in features:
        gbk_record.features.append(item)

    SeqIO.write(gbk_record, filename, "gb")


# Verify gene orientation and write information for positive results.
def write_results(candidates, organism, filename, sequence, temp):

    if temp:
        filename = "tmp/" + filename

    with open(args.out + '/'  + filename + ".txt", "wt") as out:
        out.write("QUERY: " + organism + "\n")
        i = 1

        for item in candidates:
            if item[0][0] < item[0][1]:
                transposase_strand = '+'
                transposase_arrow = "Tn3---->"
            else:
                transposase_strand = '-'
                transposase_arrow = "<----Tn3"
            if item[2][0] < item[2][1]:
                toxin_strand = '+'
                toxin_arrow = "T-->"
            else:
                toxin_strand = '-'
                toxin_arrow = "<--T"
            if item[3][0] < item[3][1]:
                antitoxin_strand = '+'
                antitoxin_arrow = "A-->"
            else:
                antitoxin_strand = '-'
                antitoxin_arrow = "<--A"
            if item[1][0] < item[1][1]:
                resolvase_strand = '+'
                resolvase_arrow = "R--->"
            else:
                resolvase_strand = '-'
                resolvase_arrow = "<---R"

            candidate = "Candidate " + str(i)
            coordinates = (item[0][0], item[0][1], item[1][0], item[1][1],
                           item[2][0], item[2][1], item[3][0], item[3][1])
            arrows = (transposase_arrow, resolvase_arrow,
                      toxin_arrow, antitoxin_arrow)
            strands = (transposase_strand, resolvase_strand,
                      toxin_strand, antitoxin_strand)
            order = check_order(coordinates, arrows)
            candidate_length = max(coordinates) - min(coordinates)

            out.write("\n*" + candidate + "*\n\n")

            out.write("Length: " + str(candidate_length / 1000) + " kbp\n")
            out.write("Order: " + order[1] + "\n")
            out.write(("Feature\tPosition\tStrand\tLength(bp)\tType" +
                       "\tPositives(%)\tCoverage(%)").expandtabs(20) + "\n")
            out.write(("transposase\t" +
                       str(item[0][0]) + ".." + str(item[0][1]) + "\t" +
                       transposase_strand + "\t" +
                       str(abs(item[0][0] - item[0][1])) + "\t" +
                       item[0][2] + "\t" +
                       str(item[0][3]) + "\t" +
                       str(item[0][4])).expandtabs(20) + "\n")
            out.write(("resolvase\t" +
                       str(item[1][0]) + ".." + str(item[1][1]) + "\t" +
                       resolvase_strand + "\t" +
                       str(abs(item[1][0] - item[1][1])) + "\t" +
                       item[1][2] + "\t" +
                       str(item[1][3]) + "\t" +
                       str(item[1][4])).expandtabs(20) + "\n")
            out.write(("toxin\t" +
                       str(item[2][0]) + ".." + str(item[2][1]) + "\t" +
                       toxin_strand + "\t" +
                       str(abs(item[2][0] - item[2][1])) + "\t" +
                       item[2][2] + "\t" +
                       str(item[2][3]) + "\t" +
                       str(item[2][4])).expandtabs(20) + "\n")
            out.write(("antitoxin\t" +
                       str(item[3][0]) + ".." + str(item[3][1]) + "\t" +
                       antitoxin_strand + "\t" +
                       str(abs(item[3][0] - item[3][1])) + "\t" +
                       item[3][2] + "\t" +
                       str(item[3][3]) + "\t" +
                       str(item[3][4])).expandtabs(20) + "\n")

            sequences = fetch_sequences(sequence, order[0], strands)
            whole_sequence = sequence[min(coordinates):max(coordinates)]
            whole_seq_orfs, complete_sequence_orfs = \
                            predict_orfs(whole_sequence, min(coordinates))
            if args.extend:
                ext_sequence = sequence[max(0, min(coordinates)-args.extend):
                                             max(coordinates)+args.extend]
                extended_seq_orfs, complete_extended_orfs = \
                            predict_orfs(ext_sequence, min(coordinates))

            out.write("\n>transposase\n")
            out.write(sequences[0] + "\n")
            out.write(">resolvase\n")
            out.write(sequences[1] + "\n")
            out.write(">toxin\n")
            out.write(sequences[2] + "\n")
            out.write(">antitoxin\n")
            out.write(sequences[3] + "\n")
            out.write(">whole sequence\n")
            out.write(whole_sequence + "\n")
            if args.extend:
                out.write(">extended sequence\n")
                out.write(ext_sequence + "\n")

            out.write("\nPredicted ORFs in the whole sequence:\n")
            for orf in complete_sequence_orfs:
                out.write(str(orf[0]) + ".." + str(orf[1]) + "\n")

            if args.extend:
                out.write("\nPredicted ORFs in the extended sequence:\n")
                corrected_extended_orfs = [(x-args.extend, y-args.extend)
                                     for (x, y) in complete_extended_orfs]
                for orf in corrected_extended_orfs:
                    out.write(str(orf[0]) + ".." + str(orf[1]) + "\n")

            i += 1

            if args.gbk:
                candidate = candidate.replace(' ', '').lower() + ".gbk"
                filename = filename.replace("tmp/", '')

                if not args.extend:
                    name = args.out + '/' + filename + '_' + candidate
                    orfs = match_candidate_to_orfs(item,
                                                   whole_seq_orfs,
                                                   min(coordinates))
                    write_gbk(whole_sequence, orfs, name, organism)
                else:
                    name = args.out + '/' + filename + "_extended_" + candidate
                    orfs = match_candidate_to_orfs(item,
                                                   extended_seq_orfs,
                                                   min(coordinates))
                    write_gbk(ext_sequence, orfs, name, organism)


# Check for a multifasta file and split it if necessary.
def split_multifasta(fasta, filename):

    contigs = []
    contig_files = []
    report_files = []
    header = ''
    sequence = ''

    with open(fasta, "rt") as fasta:
        for line in fasta:
            if '>' in line[0] and not header:
                header = line
            elif '>' in line[0] and header:
                contigs.append((header, sequence))
                header = line
                sequence = ''
            elif line:
                sequence += line.rstrip()
        else:
            contigs.append((header, sequence))

    if len(contigs) > 1:
        n = 1
        os.makedirs(args.out + "/tmp", exist_ok=True)
        for item in contigs:
            contig_file = args.out + "/tmp/" + filename + "_contig" + str(n) + ".tmp"
            contig_files.append(contig_file)
            report_files.append(contig_file.replace(".tmp", ".txt"))
            with open(contig_file, "wt") as temp:
                temp.write(item[0])
                temp.write(item[1])
                n += 1

    return contig_files, report_files


# Remove temporary contig reports and concatenate them into a single report.
def merge_reports(report_files, filename):

    for item in report_files:
        basename = '_'.join(item.split('/')[-1].split('_')[:-1])
        if basename != filename:
            continue
        if item.split('/')[-1] not in os.listdir(args.out + "/tmp/"):
            continue
        else:
            with open(args.out + '/' + filename + ".txt", "wt") as report:
                with open(item, "rt") as contig_report:
                    for line in contig_report:
                        report.write(line)
                    else:
                        report.write("\n"+ 130*'-' + "\n")


# Check tBLASTn output, use it to search for candidates, and write results.
def analyze_fasta(fasta, temp):

    filename = fasta.split('/')[-1].rsplit('.', 1)[0]
    sequence = ''

    with open(fasta, "rt") as genome:
        organism = next(genome)
        organism = organism.split(',')[0].split('|')[-1].strip(" >.\n")
        for line in genome:
            sequence += line.rstrip()

    if "tblastn" not in os.listdir(args.out):
        os.mkdir(args.out + "/tblastn/")

    if not os.path.isfile(args.out + "/tblastn/" + filename + ".tblastn"):
        blast = run_tblastn(fasta)
        with open(args.out + "/tblastn/" +
                 filename + ".tblastn", "wt") as results:
            for item in blast:
                results.write(item + "\n")
    else:
        blast = []
        with open(args.out + "/tblastn/" +
                  filename + ".tblastn", "rt") as results:
            for line in results:
                blast.append(line.strip())

    candidates = find_candidates(blast, args.distance,
                                 args.positives, args.coverage)
    if len(candidates) > 1 and args.merge:
        candidates = merge_candidates(candidates)
    print(strftime("%c") + "\tFound " + str(len(candidates)) +
                                " candidate(s) in " + organism + ".")
    if candidates:
        write_results(candidates, organism, filename, sequence, temp)


# Verify input and call functions accordingly.
def main():

    try:
        print(strftime("%c") + "\tSearching for transposable elements.")
        os.makedirs(args.out, exist_ok=True)

        # File the command line to better allow reproducibility.
        with open(args.out + "/info.txt", "wt") as info:
            info.write("time:\t" + strftime("%c") + "\n")
            info.write("working directory:\t" + os.getcwd().replace(' ', "\ ") + "\n")
            info.write("command line:\t" + subprocess.list2cmdline(sys.argv[0:]) + "\n")

        # Run a multi-threaded contig analysis on a single multifasta genome.
        if len(args.file) == 1:
            fasta_file = args.file[0]
            filename = fasta_file.split('/')[-1].rsplit('.', 1)[0]
            contig_files, report_files = split_multifasta(fasta_file, filename)

            if contig_files:
                temp = [True] * len(contig_files)
                Pool(args.threads).starmap(analyze_fasta, zip(contig_files, temp))
                merge_reports(report_files, filename)
                if "tmp" in os.listdir(args.out):
                    rmtree(args.out + "/tmp")
            else:
                analyze_fasta(args.file[0], False)

        # Run parallel analyses if several files were provided.
        elif len(args.file) > 1:
            contigs_list = []
            reports_list = []
            temp = []
            split_reports = {}

            for item in args.file:
                filename = item.split('/')[-1].rsplit('.', 1)[0]
                contig_files, report_files = split_multifasta(item, filename)

                if contig_files:
                    for contig in contig_files:
                        contigs_list.append(contig)
                        temp.append(True)
                    for report in report_files:
                        reports_list.append(report)
                    split_reports.update({filename:reports_list})
                else:
                    contigs_list.append(item)
                    temp.append(False)

            Pool(args.threads).starmap(analyze_fasta, zip(contigs_list, temp))

            for multifasta, reports in split_reports.items():
                merge_reports(reports, multifasta)

            if "tmp" in os.listdir(args.out):
                rmtree(args.out + "/tmp")

    except KeyboardInterrupt:
        print(strftime("%c") + "\tProgram interrupted.")
    except:
        print(strftime("%c") + "\tError. Please check input files.")


if __name__ == "__main__":
    main()
