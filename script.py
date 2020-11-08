from Bio import SeqIO 
from Bio.Seq import Seq

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


#secuencias de ejemplo
DNA_sequence_1 = "GATCGATGGGCCTATATA"
DNA_sequence_2 = "GGATCGAAATCGC"

#Función que regresa el complemento inverso de la concatenación
#de dos secuencias de DNA
def concatenate_and_get_reverse_of_complement(DNA_sequence_1, DNA_sequence_2):
    sequence_1 = Seq(DNA_sequence_1)
    sequence_2 = Seq(DNA_sequence_2)
    sequence_3 = sequence_1 + sequence_2

    return sequence_3.reverse_complement()
concatenate_and_get_reverse_of_complement(DNA_sequence_1, DNA_sequence_2)
