# This script automates Veksler-Lublinksys bi targeting prediction algorithm to allow to use large databases to run the algorithm with ease.
# For large lists of miRNA, the targeting prediction algorithm would require intensive computational work.
# This script simply groups miRNA into their respective families and run the algorithm altogether.
# More information about the algorith and all its requred files is available here http://www.cs.bgu.ac.il/~vaksler/TargetPrediction.htm

# Usage
# Simply group all the required into one folder. The only additional file reqired aside from those stated in the link above are the miRNA family list which can be downloaded from mirbase.org
# Change all the files names and paths as needed.
# miR_Family_Info.txt refers to the miR Families from miRbase.org.
# viral_miRNA_all_mature.txt refers to the miRNA list
# hiv_mrna.fasta refers to the mRNA list.
# Once all changes are made, simply run the script on terminal "python target_prediction_auto.py"

import os
import subprocess
import shutil

f = open('miR_Family_Info.txt', 'r')
m = open('viral_miRNA_all_mature.txt', 'r')

txt = f.read()
split = txt.split('\n')
testmiRNA = m.read()
splitmiRNA = testmiRNA.split('\n')

array = []
arraymiRNA = []
arrayFamily = []
members = []
go = []
y = 0
coun = 0
no_family = []


counter = 0
lengthmiRNA = len(splitmiRNA)

for i in split:
    array = i.split('\t')
    if array[2] == '9606':
        family = {}
        family['miR Family'] = array [0]
        family['Seed+m8'] = array[1]
        arrayFamily.append(family)


d = {}

for obj in arrayFamily:
    d[obj['miR Family']]= obj

desired_list = d.values()


for i in range(0, lengthmiRNA/2):
    miRNA={}
    miRNA['description'] = splitmiRNA[counter]
    miRNA['sequence'] = splitmiRNA[counter+1]

    counter=counter+2
    arraymiRNA.append(miRNA)


for a in desired_list:
    for i in arraymiRNA:
        str1 = i['sequence']
        str2 = a['Seed+m8']
        index = str1.find(str2)
        if index == 1:
            members.append(i)
    if not len(members) == 0:
        all_families = {}
        all_families['family'] = a['miR Family']
        y = y + len(members)
        all_families['members'] = members
        go.append(all_families)
        members = []


for i in arraymiRNA:
    for a in go:
        for b in a['members']:
            if i['description'] == b['description']:
                coun = 1
    if coun == 0:
        no_family.append(i)

    coun = 0




def write(title,text):
    replaced = title.replace('/','-')
    filename = "/Users/leolopelofranco/Desktop/computational_work/target_prediction_III/hiv_mrna/" + replaced + "/viral_miRNA_all_mature.fasta"
    config = "/Users/leolopelofranco/Desktop/computational_work/target_prediction_III/hiv_mrna/" + replaced + "/configparams.txt"
    mirsFile = 'mirsFile' + ' = ' + '/Users/leolopelofranco/Desktop/computational_work/target_prediction_III/hiv_mrna/' + replaced +'/viral_miRNA_all_mature.fasta' + '\n' + '\n'
    UTRsFile = 'UTRsFile' + ' = ' + '/Users/leolopelofranco/Desktop/computational_work/target_prediction_III/hiv_mrna/hiv_mrna.fasta' + '\n' + '\n'
    alignment = 'alignmentsFileLoc' + ' = ' + '/Users/leolopelofranco/Desktop/computational_work/target_prediction_III/hiv_mrna/' + replaced + '/' + '\n' + '\n'

    ### remains constant ###
    duplex = 'RNAduplexLoc = /Users/leolopelofranco/Desktop/computational_work/target_prediction/sample/ViennaRNA-2.1.8/Progs/'

    ####


    params = mirsFile + UTRsFile + alignment + duplex
    dir = os.path.dirname(filename)
    if not os.path.exists(dir):
        os.makedirs(dir)
    with open(filename, 'w') as x:
        x.write(text)
        x.close()
    with open(config, 'w') as i:
        i.write(params)
        i.close()
    dst = "/Users/leolopelofranco/Desktop/computational_work/target_prediction_III/hiv_mrna/" + replaced + "/"
    shutil.copy('/Users/leolopelofranco/Desktop/computational_work/target_prediction_III/targetPrediction.jar', dst)
    shutil.copy('/Users/leolopelofranco/Desktop/computational_work/target_prediction_III/parameters.txt', dst)
    if __name__ == "__main__":
       startingDir = os.getcwd()  # save our current directory
       testDir = dst
       os.chdir(testDir) # change to our test directory
       os.system("java -jar targetPrediction.jar configparams.txt parameters.txt")
       os.chdir(startingDir) # change back to where we started

text = ''

abc = []


for i in no_family:
    text = text + i['description'] + "\n" + i['sequence'] + "\n"

title = 'random'
abc.append(title);
stripped = text.rstrip()
write(title,stripped)
text = ''
title=''

for i in go:
    title = i['family']
    replaced = title.replace('/','-')
    for a in i['members']:
         text = text + a['description'] + "\n" + a['sequence'] + "\n"
    stripped = text.rstrip()
    abc.append(replaced)
    write(title,stripped)
    text = ''
    title=''

print abc
