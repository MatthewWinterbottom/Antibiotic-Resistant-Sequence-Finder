#!/usr/bin/env python

from Bio import SeqIO
from subprocess import call
import os
from Bio.Blast import NCBIXML
import sys
import random
from glob import glob
from Bio import Entrez
import datetime

print(datetime.datetime.now().time())

def remove_stuff():

    path_to_assembly = os.path.expanduser('~/CarbFinder/temp_files/assembly/')

    path_to_forward_and_reverse = os.path.expanduser('~/CarbFinder/temp_files/forward_and_reverse/')

    path_to_output_files = os.path.expanduser('~/CarbFinder/OUTPUT/')

    path_to_ORFs = os.path.expanduser('~/CarbFinder/temp_files/ORFs/')

    path_to_blast_output = os.path.expanduser('~/CarbFinder/temp_files/blast_output/')

    path_to_HMM_output =  os.path.expanduser('~/CarbFinder/temp_files/HMM_output/')

    call(['rm', '-r', path_to_assembly])
    call(['rm', '-r', path_to_forward_and_reverse])
    call(['rm', '-r', path_to_output_files])
    call(['rm', '-r', path_to_ORFs])
    call(['rm', '-r', path_to_blast_output])
    call(['rm', '-r', path_to_HMM_output])

    call(['mkdir', path_to_ORFs])
    call(['mkdir', path_to_blast_output])
    call(['mkdir', path_to_HMM_output])
    call(['mkdir', path_to_output_files])

def fastqc():

    metagenome_path = sys.argv[1]

    fastqc_path = os.path.expanduser('~/CarbFinder/tools/fastq-dump')

    output_path = os.path.expanduser('~/CarbFinder/temp_files/forward_and_reverse')

    call([fastqc_path, '-I', '--split-files', metagenome_path, '-O', output_path])


def assembly():

    metagenome_path = sys.argv[1]

    metagenome_name = metagenome_path.split('/')[-1]

    metagenome_name = metagenome_name.split('.')[0]

    metaspades_path = os.path.expanduser('~/CarbFinder/tools/SPAdes-3.12.0-Darwin/bin/metaspades.py')

    forward_file = metagenome_name + '_1.fastq'

    reverse_file = metagenome_name + '_2.fastq'

    forward_file_path = os.path.expanduser('~/CarbFinder/temp_files/forward_and_reverse/' + forward_file)

    reverse_file_path = os.path.expanduser('~/CarbFinder/temp_files/forward_and_reverse/' + reverse_file)

    output_path_assembly = os.path.expanduser('~/CarbFinder/temp_files/assembly')

    call([metaspades_path, '-1', forward_file_path, '-2', reverse_file_path, '-o', output_path_assembly])

def same_frame_verification(start_codon, stop_codon):

    length_between_sequences = (stop_codon - (start_codon - 1))

    divide_by_three = length_between_sequences / 3

    if (divide_by_three).is_integer():

        return True

def write_out_ORFs(rec_id, seq, counter):

    counter = str(counter)

    output_path = os.path.expanduser('~/CarbFinder/temp_files/ORFs/')

    with open(output_path + 'ORF_' + counter + '_' + rec_id + '.fasta', 'w') as f:

        f.write('>Contig:|' + rec_id + '| uncharacterised_orf_' + counter + '\n' + seq\
                + '\n')

    f.close()

def sequence_extraction(org_log, record):

    for counter, item in enumerate(org_log):

        seq = str(record[item:org_log[item]].seq)

        rec_id = str(record.id)

        write_out_ORFs(rec_id, seq, counter)

def orf_log(start_codon, stop_codon, record):

    orf_log = {}

    for starts in start_codon:

        for stops in stop_codon:

            if stops > starts:

                if same_frame_verification(starts, stops) == True:

                    if (starts + 50) > stops:

                        break

                    elif stops in orf_log.values():

                        break

                    else:

                        orf_log[starts] = stops

                        break
    if len(orf_log) > 0:

        sequence_extraction(orf_log, record)


def start_and_stop_codon_finder():

    records = os.path.expanduser('~/CarbFinder/temp_files/assembly/contigs.fasta')

    for record in SeqIO.parse(records, 'fasta'):

        start_codon_positions = []

        stop_codon_positions = []

        list_stop_codon = ['TAG', 'TAA', 'TGA']

        range_of_record = range(len(record))

        for pos in range_of_record:

            codon = record[pos: pos + 3].seq

            if codon == 'ATG':

                start_codon_positions.append(pos)

            if codon in list_stop_codon:

                stop_codon_positions.append(pos)

        orf_log(start_codon_positions, stop_codon_positions, record)

def change_name_of_contigs():

    contigs = os.path.expanduser('~/CarbFinder/temp_files/assembly/contigs.fasta')

    seqs = []

    for record in SeqIO.parse(contigs, 'fasta'):

        seqs.append(str(record.seq))

    with open(contigs, 'w') as f:

        for counter, item in enumerate(seqs):

            f.write('>Contig_' + str(counter + 1) + '\n'  + item + '\n')

    f.close()

def parsing_blast_results_to_obtain_organism(upstream, downstream):

    Entrez.email = "M.Winterbottom2@ncl.ac.uk"

    try:

        result_handle = open(upstream)

        blast_record = NCBIXML.read(result_handle)

        description = str(blast_record.descriptions[0])

        accession = description.split('|')[3]

        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")

        records = Entrez.read(handle)

        organism_1 = (records[0]['GBSeq_organism'])

    except Exception:

        organism_1 = 'Organism from upstream sequence not found'

    try:


        result_handle = open(downstream)

        blast_record = NCBIXML.read(result_handle)

        description = str(blast_record.descriptions[0])

        accession = description.split('|')[3]

        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")

        records = Entrez.read(handle)

        organism_2 = (records[0]['GBSeq_organism'])

    except Exception:

        organism_2 = 'Organism from downstream sequence not found'

    return organism_1, organism_2


def blast2():

    blastn = os.path.expanduser('~/CarbFinder/tools/blastn')

    db = os.path.expanduser('~/CarbFinder/blast_db/nt.59')

    upstream = os.path.expanduser('~/CarbFinder/temp_files/upstream.fasta')

    downstream = os.path.expanduser('~/CarbFinder/temp_files/downstream.fasta')

    out_upstream = os.path.expanduser('~/CarbFinder/temp_files/blast_results_upstream')

    out_downstream = os.path.expanduser('~/CarbFinder/temp_files/blast_results_downstream')

    call([blastn, '-query', upstream, '-outfmt', '5' ,'-db', db, '-out', out_upstream])

    call([blastn, '-query', downstream, '-outfmt', '5', '-db', db, '-out', out_downstream])

    return parsing_blast_results_to_obtain_organism(out_upstream, out_downstream)

def write_files_for_blast(upstream, downstream):

    upstream2 = os.path.expanduser('~/CarbFinder/temp_files/upstream.fasta')

    downstream2 = os.path.expanduser('~/CarbFinder/temp_files/downstream.fasta')

    with open(upstream2, 'w') as f:

        f.write('>upstream\n' + str(upstream))

    f.close()

    with open(downstream2, 'w') as r:

        r.write('>downstream\n' + str(downstream))

    r.close()

    return blast2()

def extract_correct_contig_sequence(contig_name, ORF_sequence):

    path_to_contigs = os.path.expanduser('~/CarbFinder/temp_files/assembly/contigs.fasta')

    for record in SeqIO.parse(path_to_contigs, 'fasta'):

        if record.id == contig_name:

            start_of_orf = record.seq.find(ORF_sequence)

            end_of_orf = start_of_orf + len(ORF_sequence)

            upstream = record.seq[:start_of_orf]

            downstream = record.seq[end_of_orf:]

    return write_files_for_blast(upstream, downstream)

def organism_finder(query):

    record = SeqIO.read(query, 'fasta')

    orf_seq = str(record.seq)

    contig_name = record.id.split('|')[1]

    return extract_correct_contig_sequence(contig_name, orf_seq)

def write_blast_results(organism_one, organism_two, query, percentage, variant, e_values):

    path_to_results = os.path.expanduser('~/CarbFinder/Output/')

    record = SeqIO.read(query, 'fasta')

    seq = str(record.seq)

    rec_id = str(record.id)

    if organism_one == organism_two:

        with open(path_to_results + 'BLAST_RESULTS.txt', 'w') as f:

            f.write('\n\nBLAST RESULT\n'\
                    + '-----------------------------\n'\
                    + 'A POSITIVE RESULT WAS FOUND\n'\
                    + '\nRecord: ' + rec_id\
                    + '\nMatch: ' + variant\
                    + '\nPercentage Similiarity: ' + percentage\
                    + '\nE-value: ' + e_values\
                    + '\nOrganism: ' + organism_one\
                    + '\nSequence: ' + seq)

        f.close()

    if organism_one != organism_two:

        with open(path_to_results + 'BLAST_RESULTS.txt', 'a') as f:

            f.write('\n\nBLAST RESULT\n'\
                    + '-----------------------------\n'\
                    + 'A POSITIVE RESULT WAS FOUND\n'\
                    + '\nRecord: ' + rec_id\
                    + '\nMatch ' + variant\
                    + '\nPercentage Similiarity: ' + percentage\
                    + '\nE-value: ' + e_values\
                    + '\nOrganism one: ' + organism_one\
                    + '\nOrganism two: ' + organism_two\
                    + '\nSequence: ' + seq)

        f.close()

    with open(path_to_results + variant + '_new_variant_sequence.fasta', 'w') as f:

        f.write('>' + variant + '-like\n' + seq + '\n')

    f.close()

def add_to_file(query):

    record = SeqIO.read(query, 'fasta')

    random_variant_id = str(random.randint(0, 1000000000))

    rec_id = str(record.id) + random_variant_id

    rec_seq = str(record.seq)

    path_to_all_sequences_file = os.path.expanduser('~/CarbFinder/temp_files/all_sequences.fasta')

    with open(path_to_all_sequences_file, 'a') as f:

        f.write('>' + rec_id + '\n' + rec_seq)

    f.close()

def make_new_blast_db():

    makeblastdb = os.path.expanduser('~/CarbFinder/tools/makeblastdb')

    output = os.path.expanduser('~/CarbFinder/blast_db/carb_nucl_db')

    path_to_all_sequences_file = os.path.expanduser('~/CarbFinder/temp_files/all_sequences.fasta')

    call([makeblastdb, '-in', path_to_all_sequences_file, '-input_type', 'fasta', '-dbtype', 'prot', '-out', output])

def parse_output(query, output):

    length_of_record = len(SeqIO.read(query, 'fasta'))

    result_handle = open(output)

    blast_record = NCBIXML.read(result_handle)

    for alignment in blast_record.alignments:

        for hsp in alignment.hsps:

            if hsp.expect < 0.01:

                percentage = (hsp.identities/length_of_record) * 100

                percentage2 = str(percentage)

                variant = alignment.title.split(' ')[-1]

                e_value = str(hsp.expect)

                organism_one, organism_two = organism_finder(query)

                write_blast_results(organism_one, organism_two, query, percentage2, variant, e_value)

                if percentage != 100.0:

                    add_to_file(query)

                    make_new_blast_db()

        break

def blast(query, output):

    blastn = os.path.expanduser('~/CarbFinder/tools/blastn')

    blast_db = os.path.expanduser('~/CarbFinder/blast_db/carb_nucl_db')

    call([blastn, '-query', query, '-outfmt', '5', '-db', blast_db, '-out', output])


def run():

    path_to_ORFs = os.path.expanduser('~/CarbFinder/temp_files/ORFs/*')

    path_to_output = os.path.expanduser('~/CarbFinder/temp_files/blast_output/')

    for file in glob(path_to_ORFs):

        output_file = file.split('/')[-1]

        output_file = output_file.split('.')[0] + '_blast_results'

        output_file = path_to_output + output_file

        blast(file, output_file)

        parse_output(file, output_file)


def HMM_scan(test_file, orf_name):

    hmmscan = os.path.expanduser('~/CarbFinder/tools/hmmscan')

    hmms = os.path.expanduser('~/CarbFinder/HMMs/*.hmm')

    output_path = os.path.expanduser('~/CarbFinder/temp_files/HMM_output/')

    arg = '--domtblout'

    for hmm_path in glob(hmms):

        hmm = hmm_path.split('/')[-1]

        hmm = hmm.split('.')[0]

        output = output_path + 'OUTPUT_' + hmm + '_' + orf_name

        call([hmmscan, arg, output, hmm_path, test_file])

def create_new_HMM(HMM, query):

    hmm_build = os.path.expanduser('~/CarbFinder/tools/hmmbuild')

    hmm_press = os.path.expanduser('~/CarbFinder/tools/hmmpress')

    path_to_query = os.path.expanduser('~/CarbFinder/temp_files/ORFs/') + query + '*'

    input1 = glob(path_to_query)[0]

    random_no = str(random.randint(0, 1000000))

    output = os.path.expanduser('~/CarbFinder/HMMs/') + HMM + '_like_' + random_no + '.hmm'

    call([hmm_build, output, input1])

    call([hmm_press, output])

def write_HMM_result(HMM, query, e_value, sequence_score, organism_one, organism_two, orf):

    e_value = str(e_value)

    organism_one = str(organism_one)

    organism_two = str(organism_two)

    record = SeqIO.read(orf, 'fasta')

    seq = str(record.seq)

    path_to_results = os.path.expanduser('~/CarbFinder/Output/')

    with open(path_to_results + 'HMM_RESULTS.txt', 'a') as f:

        f.write('\n\nHMM RESULT\n'\
                + '-----------------------------\n'
                + 'A POSITIVE RESULT WAS FOUND\n'\
                + '\nRecord: ' + query\
                + '\nHMM: ' + HMM\
                + '\nSequence Score: ' + sequence_score\
                + '\nE-value: ' + e_value\
                + '\nOrganism one: ' + organism_one\
                + '\nOrganism two: ' + organism_two\
                + '\nSequence: ' + seq)

    f.close()

    with open(path_to_results + HMM + '_new_variant_sequence.fasta', 'w') as f:

        f.write('>' + HMM + '-like\n' + seq + '\n')

    f.close()

def gather_information(result_list, orf):

    orf2 = os.path.expanduser('~/CarbFinder/temp_files/ORFs/') + orf + '.fasta'

    e_values = [100.0]

    for result in result_list:

        result = result.split()

        e_value = float(result[6])

        if e_value < e_values[-1]:

            e_values.append(e_value)

            HMM = result[0]

            sequence_score = result[7]

            query = result[3]

    create_new_HMM(HMM, orf)

    organism_one, organism_two = organism_finder(orf2)

    write_HMM_result(HMM, query, e_value, sequence_score, organism_one, organism_one, orf2)

def parse_HMM_output(distinct_ORFs):

    path_to_HMM_output = os.path.expanduser('~/CarbFinder/temp_files/HMM_output/')

    save = False

    for orf2 in distinct_ORFs:

        orf = path_to_HMM_output + '*' + str(orf2)

        positive_result = []

        for file in glob(orf):

            with open(file, 'r') as f:

                for line in f.readlines():

                    if line.startswith('#\n'):

                        save = False

                    if save:

                        positive_result.append(line)

                    if line.startswith('#-------------------'):

                        save = True

            f.close()

        if len(positive_result) > 0:

            gather_information(positive_result, orf2)



def get_distinct_ORFs():

    path_to_HMM_output = os.path.expanduser('~/CarbFinder/temp_files/HMM_output/*')

    distinct_ORFs = []

    for file in glob(path_to_HMM_output):

        ORF = file.split('_')[-4:]

        ORF = '_'.join(ORF)

        if ORF not in distinct_ORFs:

            distinct_ORFs.append(ORF)

    parse_HMM_output(distinct_ORFs)

def HMM_run():

    path_to_ORFs = os.path.expanduser('~/CarbFinder/temp_files/ORFs/*')

    for file in glob(path_to_ORFs):

        ORF_name = file.split('/')[-1]

        ORF_name = ORF_name.split('.')[0]

        HMM_scan(file, ORF_name)

    get_distinct_ORFs()

remove_stuff()
fastqc()
assembly()
change_name_of_contigs()
start_and_stop_codon_finder()
run()
HMM_run()

print(datetime.datetime.now().time())
