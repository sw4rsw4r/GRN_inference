import pandas as pd
import csv

results_dir = 'data/'
geneannotfile = 'data/GRCh38_geneannot.txt'

enspTogenenameDict = dict()
geneannot = pd.read_csv(geneannotfile, sep='\t', header=0, dtype=str)
geneannot = geneannot.loc[geneannot.ensembl_peptide_id.notna()]
for index, row in geneannot.iterrows():
    if not (pd.isna(row['gene_name']) or pd.isna(row['ensembl_peptide_id'])):
        enspTogenenameDict[row['ensembl_peptide_id']] = row['gene_name']

rawnetworkfile = results_dir + 'v11.0_string_new.raw'
writenetworkfile = results_dir + 'human_string_v2.7_withdirection_new_noexpfilter'

score_dict = dict()
direction_dict = dict()

with open(rawnetworkfile) as csvfile:
    lines = csv.reader(csvfile, delimiter=' ')
    for line in lines:
        if line[0] == 'protein1':
            continue
        if int(line[2]) >= 500:
            try:
                # source = enspTogenenameDict[line[0].strip('^9606.|^10090.')]
                # target = enspTogenenameDict[line[1].strip('^9606.|^10090.')]
                source = enspTogenenameDict[line[0].replace(
                    '9606.', '').replace('10090.', '')]
                target = enspTogenenameDict[line[1].replace(
                    '9606.', '').replace('10090.', '')]
                if score_dict.get((source, target), 0) < int(line[2]):
                    score_dict[source, target] = int(line[2])
                    if (len(line) > 3):
                        direction_dict[source, target] = line[3]
                    else:
                        direction_dict[source, target] = 'NA'
            except KeyError:
                continue

fw = open(writenetworkfile, 'w')
for (tf, gene) in score_dict.keys():
    fw.writelines('\t'.join([tf.upper(), gene.upper(),
                             str(score_dict[(tf, gene)]), direction_dict[(tf, gene)]])+'\n')
fw.close()
