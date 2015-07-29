
# This script parses all the files created by Veksler-Lublinksys Target Prediction Algorithm from the target_prediction_auto.py
# Note that you need to know the list of all miRNA families created from the previous script (target_prediction_auto.py). This will be the folder names.
# Set them as an array of strings and equate them to families variable.
# Make sure to change os.chdir directory path which serves as the link to go through all the folders.


import glob
import os

targets = []
all_miRNA = []
families = ['miR-4270-4441']
for fam in families:

    os.chdir("/Users/Desktop/computational_work/target_prediction_III/hiv_linked_proteins/" + fam + "/")

    for file in glob.glob("*.aln"):
        f = open(file, 'r')
        text = f.read()
        rna = text.split('\n\n')
        fi = file.split(".")
        result = {}
        family= {}
        last = rna[-1].split('***************************************************************')
        last = filter(None, last)
        for lis in last:
            l = lis.split('\n')
            l = filter(None, l)
            if (l):
                first_line = l[0].split()
                las = first_line[0].replace("[", "")
                la = las.replace("]", "")
                miRNA = {}
                miRNA['miRNA'] = la
                l.pop(0)
                l = list(set(l))
                for i in l:
                    target_lines = i.split()
                    target = {}
                    target['target'] = target_lines[0].split('|')[1]
                    target['number_hits'] = target_lines[target_lines.index('hits')+1]
                    targets.append(target)

                miRNA['number_target_genes'] = len(targets)
                miRNA['targets'] = targets
                all_miRNA.append(miRNA)
                targets = []
        for i in all_miRNA:
            for a in i['targets']:
                if int(a['number_hits']) > 1:
                    print i['miRNA'] + '\t' + a['number_hits'] + '\t' + a['target']

        all_miRNA = []
