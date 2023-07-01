import csv
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import seaborn as sns
import matplotlib.pyplot as plt

def source():
    # Configurar las opciones de búsqueda en el NCBI
    Entrez.email = 'arnoldhurtado001@gmail.com'
    term = 'srcdb_refseq[PROP] AND biomol_genomic[PROP]'  # Término de búsqueda para obtener secuencias genómicas
    handle = Entrez.esearch(db='nucleotide', term=term, retmax=1000)  # Retmax limita el número de resultados a 1000
    record = Entrez.read(handle)
    ids = record['IdList']

    # Obtener información detallada de cada secuencia
    result = []
    for i, uid in enumerate(ids):
        handle = Entrez.efetch(db='nucleotide', id=uid, rettype='gb', retmode='text')
        record = SeqIO.read(handle, 'gb')
        organism = record.annotations.get('organism', 'N/A')
        result.append(organism)

    # Contar la frecuencia de cada especie
    species_count = {}
    for organism in result:
        species_count[organism] = species_count.get(organism, 0) + 1

    # Guardar los resultados en un archivo CSV
    with open('results/source.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Organismo Fuente', 'Frecuencia'])
        for organism, count in species_count.items():
            writer.writerow([organism, count])

    print('Extracción y conteo de especies completados. Los resultados se han guardado en "results/source.csv".')


def sequences(invalid_chars=['U']):
    # Configurar las opciones de búsqueda en el NCBI
    Entrez.email = 'arnoldhurtado001@gmail.com'  # Reemplaza con tu dirección de correo electrónico
    term = 'srcdb_refseq[PROP] AND biomol_genomic[PROP]'  # Término de búsqueda para obtener secuencias genómicas
    handle = Entrez.esearch(db='nucleotide', term=term, retmax=1000)  # Retmax limita el número de resultados a 1000
    record = Entrez.read(handle)
    ids = record['IdList']

    # Obtener información detallada de cada secuencia
    result = []
    for i, uid in enumerate(ids):
        handle = Entrez.efetch(db='nucleotide', id=uid, rettype='gb', retmode='text')
        record = SeqIO.read(handle, 'gb')
        organism = record.annotations.get('organism', 'N/A')
        for feature in record.features:
            if feature.type == 'CDS':
                translation = feature.qualifiers.get('translation', [''])[0]
                if translation.startswith('M') and not any(char in translation for char in invalid_chars):
                    analysis = ProteinAnalysis(translation)
                    molecular_weight = analysis.molecular_weight()
                    instability_index = analysis.instability_index()
                    result.append([organism, translation, molecular_weight, instability_index])

    # Guardar los resultados en un archivo CSV
    with open('results/glupeptides.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Organismo Fuente', 'Péptido', 'Peso Molecular', 'Índice de Inestabilidad'])
        for row in result:
            writer.writerow(row)

    # Generar un joinplot utilizando los valores de peso molecular e índice de estabilidad
    data = [row[2:] for row in result]
    sns.jointplot(x=[row[0] for row in data], y=[row[1] for row in data], kind='scatter', color='red', s=50)
    plt.xlabel('Peso Molecular')
    plt.ylabel('Índice de Inestabilidad')
    plt.title('Peso Molecular vs. Índice de Inestabilidad')
    plt.savefig('results/glupeptides.png')

    print('Extracción y análisis de péptidos completados. Los resultados han sido guardados en la careta results')