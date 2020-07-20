import utils
from parsers import ClustalParser, MuscleParser
from matcher import AligmentDifferenceFinder
import os #to remove file in folder
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json

def runClustal(inputFile, reference_id, nseq = 3):
    """Call functions to parse a ClustalW alignment and produce a custom output structure"""
    parser = ClustalParser(nseq) #parser.py #(quante sequenze+file) --> sequenze lette
    alignment = parser.parse(inputFile, reference=reference_id)
    analyzer = AligmentDifferenceFinder()
    groups = analyzer.analyze(alignment)

    return utils.save(alignment, analyzer, reference_id=reference_id, tool='Muscle', path='output')


def runMuscle(inputFile, reference_id, nseq = 3):
    """Call functions to parse a MUSCLE alignment and produce a custom output structure"""
    muscle_parser = MuscleParser(nseq)  #parser.py #(quante sequenze+file) --> sequenze lette
    muscle_alignment = muscle_parser.parse(inputFile, reference=reference_id)
    analyzer = AligmentDifferenceFinder() # matcher.py
    groups = analyzer.analyze(muscle_alignment) # vettori indici mismatch rispetto al reference per sequenza

    return utils.save(muscle_alignment, analyzer, reference_id=reference_id, tool='Clustal', path='output')


def cleanUpOutputDir(folder=None):
    """clean folder output to make spaces for new output files"""
    dirname = os.path.dirname(__file__)
    path_to_dir = os.path.join(dirname, folder)  # path to directory you wish to remove
    files_in_dir = os.listdir(path_to_dir)     # get list of files in the directory

    for file in files_in_dir:# loop to delete each file in folder
        if os.name == 'nt': #if Windows
            os.remove(f'{path_to_dir}/{file}')     # delete file
        else: #if linux
            import sh
            sh.rm(f)

    #create new empty file for differences
    file = open(os.path.join(path_to_dir, "differences.txt"),"w+")
    file.write("Differences Clustal-Muscle alignment" + '\n\n')
    file.close()


def main():
    reference_id = open("input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2

    cleanUpOutputDir(folder='output')

    # # ClustalJ = file json from clustal alignment; MuscleJ = file json from muscle alignment
    # # Iran
    # ClustalJ = runClustal('analysis/iran-ref.txt', reference_id, 3)
    # MuscleJ = runMuscle('analysis/muscle-I20200523-084930-0610-44096621-p1m.clw', reference_id, 3)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt", "Iran", diff)

    # # Israel
    # ClustalJ = runClustal('analysis/israel-ref.txt', reference_id, 4)
    # MuscleJ = runMuscle('analysis/muscle-I20200523-085708-0753-28920419-p1m.clw', reference_id, 4)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt", "Israel ", diff)

    # # GISAID only
    # ClustalJ = runClustal('analysis/GISAID-all.txt', reference_id, 7)
    # MuscleJ = runMuscle('analysis/muscle-I20200523-090216-0023-41230765-p1m.clw', reference_id, 7)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt", "GISAID ", diff)

    # # NCBI only
    # ClustalJ = runClustal('analysis/all.txt', reference_id, 3)
    # MuscleJ = runMuscle('analysis/muscle-I20200512-170208-0225-69454386-p1m.clw', reference_id, 8)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt", "ncbi ", diff)

    # # All sequences
    # ClustalJ = runClustal('analysis/global.txt', reference_id, 14)
    # MuscleJ = runMuscle('analysis/muscle-I20200523-090837-0910-95910164-p1m.clw', reference_id, 15)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt","Global ", diff)

    # Compare global alignments
    print('[...]  Parsing ClustalW analysis...')
    ClustalJ = runClustal('analysis/clustal-global.clustal_num', reference_id, 14)
    print('[DONE] Parsing ClustalW analysis')

    print('[...]  Parsing Muscle analysis...')
    MuscleJ = runMuscle('analysis/muscle-global.clw', reference_id, 14)
    print('[DONE] Parsing Muscle analysis ')

    print('[...]  Comparing alignments...')
    diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ, path='output')
    print('[DONE] Comparing alignments')
    # utils.saveCompareFile("differences.txt","Global ", diff)

    analysis(ClustalJ)


def analysis(filename):
    with open(os.path.join('output', filename)) as f:
        variations = json.load(f)['unmatches'].items()

    ref_id, sequence = load_fasta(os.path.join('input', 'reference.fasta'))
    sequences = read_sequences(paths=[
        os.path.join('input', 'GISAID'),
        os.path.join('input', 'ncbi'),
    ])
    # sequences[ref_id] = sequence

    single_char_vars = []
    variation_rows = []
    for key, mismatch in variations:
        for sequence_id in mismatch['sequences']:
            variation_rows.append({
                'id': sequence_id,
                'alteration': mismatch['alt'],
            })

            for i, char in enumerate(mismatch['alt']):
                if mismatch['reference'][i] != '-':
                    alt_type = 'REPL' if char in ['A', 'C', 'G', 'T'] else 'DEL'
                else:
                    alt_type = 'INS'

                if char == '-':
                    change = '{} DEL'.format(mismatch['reference'][i])
                else:
                    change = '{} to {}'.format(mismatch['reference'][i], char)

                single_char_vars.append({
                    'id': sequence_id,
                    'alteration': char,
                    'change': change,
                    'type': alt_type,
                })

    # Convert the list to dataframe
    columns = ['id', 'alteration']
    df = pd.DataFrame(variation_rows, columns=columns)

    columns = ['id', 'alteration', 'change', 'type']
    df_chars = pd.DataFrame(single_char_vars, columns=columns)

    nas = []
    for seq_id, sequence in sequences.items():
        nas.append({
            'id': seq_id,
            'nas': sequence.count('N')
        })

    colors = ['mediumvioletred', 'darkorange', 'dodgerblue', 'green', 'coral', 'mediumpurple', 'gold', 'm', 'olivedrab', 'salmon','yellowgreen', 'c', 'lightgray']

    # Count NAs
    df_nas = pd.DataFrame(nas, columns=['id', 'nas']).set_index('id').sort_values(by='nas')
    ax = df_nas.plot.barh(legend=False)
    ax.set_title('Not Available Bases')
    ax.set_xlabel('Count')
    ax.set_xlim((0, 1450))
    ax.set_ylabel('Sequence ID')
    for i, (k, v) in enumerate(df_nas.iterrows()):
        v = v['nas']
        ax.text(v + 10, i - .15, str(v), color='black')
    plt.tight_layout()
    plt.savefig(os.path.join('relazione', 'images', 'plot-not-availables'))

    # Char Count per Variation Type
    chars_per_var = df_chars.groupby(['type']).count().sort_values(by='alteration')
    ax = chars_per_var.plot.pie(y='alteration', autopct='%.2f%%', labels=['Insertion', 'Replacement', 'Deletion'], colors=colors, legend=False)
    ax.set_title('Variation types')
    ax.set_ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join('relazione', 'images', 'plot-variations-per-type'))

    # Deletions per sequence
    chars_per_var = df_chars[df_chars['type']=='DEL'].groupby(['id']).count()
    ax = chars_per_var.plot.pie(y='alteration', autopct='%.2f%%', legend=False, colors=colors, explode=len(sequences.keys()) * [0.05])
    ax.set_title('Deletions per sequence')
    ax.set_ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join('relazione', 'images', 'plot-deletions-per-sequence'))

    # Replacements per sequence
    chars_per_var = df_chars[df_chars['type']=='REPL'].groupby(['id']).count()
    ax = chars_per_var.plot.pie(y='alteration', autopct='%.2f%%', legend=False, colors=colors, explode=(len(sequences.keys()) - 1) * [0.05])
    ax.set_title('Replacements per sequence')
    ax.set_ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join('relazione', 'images', 'plot-replacements-per-sequence'))

    # Alterations
    chars_per_var = df_chars.groupby(['change']).count().sort_values(by='type')
    ax = chars_per_var.plot.barh(y='type', legend=False)
    ax.set_title('Alterations')
    ax.set_xlabel('Count')
    ax.set_xlim((0, 470))
    ax.set_ylabel('Alteration')
    for i, (k, v) in enumerate(chars_per_var.iterrows()):
        v = v['type']
        ax.text(v + 4, i - .15, str(v), color='black')
    plt.tight_layout()
    plt.savefig(os.path.join('relazione', 'images', 'plot-alterations'))

    # Mismatches
    seq_per_var = df_chars.groupby(['id']).count().sort_values(by='alteration')
    ax = seq_per_var.plot.barh(y='alteration', legend=False)
    ax.set_title('Mismatches')
    ax.set_xlabel('Count')
    ax.set_xlim((0, 130))
    ax.set_ylabel('Sequence ID')
    for i, (k, v) in enumerate(seq_per_var.iterrows()):
        v = v['alteration']
        ax.text(v + 2, i - .15, str(v), color='black')
    plt.tight_layout()
    plt.savefig(os.path.join('relazione', 'images', 'plot-mismatches'))

    # plt.show()
    plt.close()


def read_sequences(paths=[]):
    sequences = {}

    for path in paths:
        for f in os.listdir(path):
            if f.endswith('.fasta'):
                seq_id, sequence = load_fasta(os.path.join(path, f))
                sequences[seq_id] = sequence

    return sequences


def load_fasta(path):
    with open(path, 'r') as f:
        lines = f.readlines()
        seq_id = lines[0].split('|')[1] if '|' in lines[0] else lines[0].split(' ')[0][1:]
        sequence = ''.join([line.rstrip() for line in lines[1:]])
        return seq_id, sequence

if __name__ == "__main__":
    main()
