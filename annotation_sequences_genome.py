import sys
import getopt
import phylopandas as ph
import pandas as pd
import re
import csv

class Selecting_sequences():
    """This is a class to select a genome range from a coordinate file, in order to annotate

    Returns:
        [fasta] -- A fasta file with the gene sequences from the coordinates
    """

    def __init__(self):
        self.genome_fasta = None
        self.genes_coordinate = None

    def open_files(self, fasta_path, genome_path):
        self.genome_fasta = ph.read_fasta(fasta_path)
        self.genome_fasta.drop(columns=['label', 'uid', 'description'], inplace=True)

        self.genes_coordinate = pd.read_csv(genome_path, sep='\t',header=None, comment='#')
        self.genes_coordinate = self.genes_coordinate.rename(columns={0:'seqId_gene', 1:'source_gene', 2:'type_gene', 3:'start_gene', 4:'end_gene', 5:'score_gene', 6:'strand_gene', 7:'phase_gene', 8:'attributes_gene'})

    def selecting_and_saving(self, output):

        def annotations(df):
            scaffold = df.seqId_gene
            strand_gene = df.strand_gene
            start_gene = int(df.start_gene)
            end_gene = int(df.end_gene)
            sequencia_scaffold = self.genome_fasta [self.genome_fasta.id == scaffold].sequence.values[0]
            gene = sequencia_scaffold[start_gene:end_gene]
            id_sequencia = '>' + scaffold + ';' + str(start_gene) + ';' + str(end_gene) + ';' + strand_gene
            df['id_sequencia'] = id_sequencia
            df['gene_sequencia'] = gene
            return df

        coordinate_selected = self.genes_coordinate.apply(annotations, axis=1)
        coordinate_selected = coordinate_selected.loc[:,['id_sequencia', 'gene_sequencia']]

        coordinate_selected.to_csv(output, sep='\n', index=False, header=False, quoting=csv.QUOTE_NONE)

if __name__ == "__main__":

    try:
        OPTS, ARGS = getopt.getopt(sys.argv[1:], 'f:g:o:h', ['fasta_path', 'gff_path', 'output', 'help'])

    except getopt.GetoptError as err:
        print(err)
        sys.exit(1)

    SELECTING = Selecting_sequences()
    FASTA_PATH = None
    GFF_PATH = None
    OUTPUT = None

    for opt, arg in OPTS:

        if opt in ('-h', '--help'):
            print('''
            This program select a genome range (fasta file) from a some coordinates (.gff, .gtf or .txt file).
            annotation_sequences_genome.py -f --fasta_path -g --gff_path -o --output -h --help
            ''')
            sys.exit(2)

        elif opt in ('-f', '--fasta_path'):
            if re.match('.+fasta$|.+txt$', arg):
                FASTA_PATH = arg
            else:
                print("-f isn't a fasta or txt file")
                sys.exit(3)

        elif opt in ('-g', '--gff_path'):
            if re.match('.+gff$|.+gtf$|.+txt$', arg):
                GFF_PATH = arg
            else:
                print("-g isn't a gff, gtf or txt file")
                sys.exit(4)

        elif opt in ('-o', '--output'):
            if re.match('.+', arg):
                OUTPUT = arg
            else:
                print("-o missing the output name")
                sys.exit(5)

    if FASTA_PATH and GFF_PATH and OUTPUT:
        SELECTING.open_files(FASTA_PATH, GFF_PATH)
        SELECTING.selecting_and_saving(OUTPUT)

    else:
        print('missing arguments -f, -g or -o')
        sys.exit(6)