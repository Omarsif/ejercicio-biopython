from Bio import SeqIO 
from Bio.Seq import Seq

import os
import unittest
import numpy

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



#Función que regresa un diccionario que incluye mRNA, secuencia de proteina y codones de paro usando código estándar
DNA_sequence = "GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA"

def print_protein_and_stop_codon_using_standard_table(DNA_sequence):
    DNA_sequence_1 = Seq(DNA_sequence)
    mRNA_sequence = DNA_sequence_1.transcribe()
    protein_sequence = mRNA_sequence.translate()
    
    dictionary = dict()
    dictionary['mRNA'] = mRNA_sequence.upper()
    DNA_sequence = DNA_sequence.upper()
    if len(DNA_sequence) %3 == 0 and ((DNA_sequence.find('ATG') >= 0) or (DNA_sequence.find('CTG') >= 0) or (DNA_sequence.find('TTG') >= 0)):
        dictionary['proteins'] = protein_sequence
        dictionary['stop_codons'] = []
        DNA_range = numpy.arange(0, len(DNA_sequence), 3)
        for i in DNA_range:
            codon = DNA_sequence[i:i+3]
            if codon in ('TAA', 'TAG', 'TGA'):
                dictionary['stop_codons'].append(codon)
    return dictionary
print_protein_and_stop_codon_using_standard_table(DNA_sequence)

#Función que regresa un diccionario que incluye mRNA, secuencia de proteina y codones de paro usando código mitocondrial de la levadura
def print_proteins_and_codons_using_mitocondrial_yeast_table(DNA_sequence):
    DNA_sequence_1 = Seq(DNA_sequence)
    mRNA_sequence = DNA_sequence_1.transcribe()
    protein_sequence = mRNA_sequence.translate()
    
    dictionary = dict()
    dictionary['mRNA'] = mRNA_sequence.upper()
    DNA_sequence = DNA_sequence.upper()
    if len(DNA_sequence) %3 == 0 and ((DNA_sequence.find('ATG') >= 0) or (DNA_sequence.find('ATA') >= 0) or (DNA_sequence.find('GTG') >= 0)):
        dictionary['proteins'] = protein_sequence
        dictionary['stop_codons'] = []
        DNA_range = numpy.arange(0, len(DNA_sequence), 3)
        for i in DNA_range:
            codon = DNA_sequence[i:i+3]
            if codon in ('TAA', 'TAG'):
                dictionary['stop_codons'].append(codon)
    print (dictionary)
    return dictionary
print_proteins_and_codons_using_mitocondrial_yeast_table(DNA_sequence)
