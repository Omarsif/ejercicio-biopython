from Bio import SeqIO 

import os

'''
Archivo de prueba:
archivo = "/home/omar/ejercicio-biopython/orchid.gbk"
'''

def summarize_contents(archivo):

    nombreDeArchivo = os.path.basename(archivo)

    record = next(SeqIO.parse(archivo, "genbank"))
    print ("name:", nombreDeArchivo)
    print ("path:", archivo)

    records = list(SeqIO.parse(archivo, "genbank"))
    print("num_records: %i records" % len(records))

    for i, record in enumerate(SeqIO.parse(archivo, "genbank")):
        print("-id:", record.id)
        print("position:")
        for seq_feature in record.features :
            print('Start: %d, Stop: %d'%(int(seq_feature.location.start), int(seq_feature.location.end)))
summarize_contents(archivo)