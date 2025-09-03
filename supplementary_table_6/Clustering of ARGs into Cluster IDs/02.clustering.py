
import pandas as pd ##to read excel file
import csv

df = pd.read_csv('BlastN 95 thr.csv')

l =list(zip(df['qseqid'],df['sseqid']))


out = []
while len(l)>0:
    first, *rest = l
    first = set(first)

    lf = -1
    while len(first)>lf:
        lf = len(first)

        rest2 = []
        for r in rest:
            if len(first.intersection(set(r)))>0:
                first |= set(r)
            else:
                rest2.append(r)     
        rest = rest2

    out.append(first)
    l = rest


with open('Clustering of ARGs into Cluster IDs.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerows(out)
    

                    