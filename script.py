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
        print("name:", record.name)
        print("Description:", record.description)
summarize_contents(archivo)