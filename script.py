from Bio import SeqIO 

import os
import unittest

#Archivo de prueba:
archivo = os.path.abspath("data/AF323668.gbk")

def summarize_contents(archivo):

    nombreDeArchivo = os.path.basename(archivo)

    extension = (os.path.splitext(archivo))[1]

    if extension == ".gbk":
        tipo = "genbank"
    else:
        tipo = "fasta"

    record = next(SeqIO.parse(archivo, tipo))

    archiveName = "name:" + nombreDeArchivo
    
    pathName = "path:" + archivo

    records = list(SeqIO.parse(archivo, tipo))
    numRecords = ("num_records: %i records" % len(records))
    
    #lista para guardar los datos de cada iteración
    a = []
    for i, record in enumerate(SeqIO.parse(archivo, tipo)):
        id = "-id:" + record.id
        a.append(id)
        recordName = "name:" + record.name
        a.append(recordName)
        description = "Description:" + record.description
        a.append(description)
        
    return archiveName,pathName,numRecords,a
summarize_contents(archivo)

#print(summarize_contents(archivo))
