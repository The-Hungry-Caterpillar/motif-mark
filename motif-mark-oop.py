import cairo
import bioinfo

def get_args():
    import argparse
    parser = argparse.ArgumentParser(description = "this script outputs a map of motifs on a regions from a fasta file, exon sequence file, and motif marker files")
	#parser.add_argument('- command line variable', '-- python variable', description)
    parser.add_argument('-f', '--fasta_file', help='input fasta file to be mapped')
    # parser.add_argument('-e', '--exon_file', help='input exon file. exon file should be sequences of exons sorted by new lines')
    # parser.add_argument('-o', '--output_file', help='input output png file name')
    return parser.parse_args()
args=get_args() 

fasta_dictionary=bioinfo.fasta_reader(args.fasta_file) 



# intron1='atgtccacatgtagtcacgtttgacatcccagggccacctcagcaggccgtctctggggagaattttctctgatttcttccccttcccttgctggacccctgcacctgctggggaagatgtagctcactccgtctagcaagtgatgggagcgagtggtccagggtcaaagccagggtgcccttactcggacacatgtggcctccaagtgtcagagcccagtggtctgtctaatgaagttccctctgtcctcaaaggcgttggttttgtttccacag'
# exon1='AAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAG'
# intron2='gtatgactcacctgtgcgacccctggtgcctgctccgcgcagggccggcggcgtgccaggcagatgcctcggagaacccaggggtttctgtggctttttgcatgcggcgggcagctgtgctggagagcagatgcttcaccaattcagaaatccaatgccttcactctgaaatgaaatctgggcatgaatgtggggagaaaccttcactaacacactcttgctaaaacatagaatca'
# print(len(intron1+exon1+intron2))
# print(f'intron 1 start = 0, stop = {len(intron1)}')
# print(f'exon 1 start = {len(intron1)}, stop = {len(intron1)+len(exon1)}')
# print(f'intron 2 start = {len(intron1)+len(exon1)}, stop = {len(intron1)+len(exon1)+len(intron2)}')



class region:
    '''creates an object called region consisting of:
    - .contig: contig sequence
    - .length: length of contig
    - .exons: list of exons (start pos, stop pos) ordered by start pos
    input: contig sequence and exon file'''
    __slots__ = ['contig', 'length', 'exons','introns', 'surface_width', 'surface_height', 'level', 'header']

    def __init__(self, header, contig):
        self.header = header
        self.contig = contig
        self.length=len(contig)
        self.exons = self.find_exons()
        self.introns = self.find_exons()
        self.surface_width = 200 
        self.surface_height = 200
        self.level=100

    def find_exons(self):
        i=0
        exons = []
        while True:
            # update index of contig while the contig is lowercase (intron)
            # once the contig is no longer lowercase record index as start of exon
            try:
                while self.contig[i].islower():
                    i+=1
            # if the contig is over, break
            except IndexError:
                break
            start=i
            # update index of contig while the contig is uppercase (exon)
            # once the contig is no longer uppercase record index as finish of exon
            try:
                while self.contig[i].isupper():
                    i+=1
            # if the contig is over, record the final position as exon stop and break
            except IndexError:
                stop=i-1
                exons.append([start,stop])
                break
            stop=i-1
            exons.append([start,stop])
        return(exons)
          
    def draw_intron(self, start, stop, c):
        '''draws introns on cairo surface.
        input: intron start (i.e. previous exons stop), intron stop (i.e. next exon start), y coordinate of cairo surface
        output: a line representing intron, proportional to actual intron'''

        c.line_to(start, self.level)
        c.line_to(stop, self.level)
        c.set_source_rgb(1, 0, 0)
        c.set_line_width(1)
        c.stroke()

    def draw_exon(self, start, stop, c):
        '''draws exons on cairo surface.
        input: exon start, exon stop, y coordinate of cairo surface (same as draw_intron)
        output: a rectangle proportional to actual exon'''

        c.rectangle(start, self.level-5, stop-start, 10)
        c.set_source_rgb(0, .5, .5)
        c.fill()

    def normalize(self):
        '''function that transforms exon genome coordinates into cairo coordinates
        input: exon dictionary, cairo surface width
        output: exon dictionary with transformed coordinates'''

        # transform exon list
        for i in range(0,len(self.exons)):
            self.exons[i] = (
                round(self.surface_width * self.exons[i][0]/self.length, 0),
                round(self.surface_width * self.exons[i][1]/self.length, 0)
                )

    def draw_region(self):
        '''draws an entire region using 'draw_exon' and 'draw_intron' fucntions
        input: exon, likely from a 'region.exons' object'''

        with cairo.SVGSurface('example.svg', self.surface_width, self.surface_height) as surface:
            c = cairo.Context(surface)
            c.move_to(0, self.level)
            
            self.draw_intron(0, self.exons[0][0],c)

            for i in range( 0, len(self.exons) ): 
                try:
                    self.draw_exon(self.exons[i][0], self.exons[i][1], c)
                    self.draw_intron(self.exons[i][1], self.exons[i+1][0], c)
                
                except IndexError: # draws the final intron. If exons is final feature this intron will be length 0
                    self.draw_intron(self.exons[i][1], self.surface_width, c)
            
            surface.write_to_png(self.header + '.png')



# gene = region('>test', fasta_dictionary['>test'])
# ex = gene.exons
# gene.normalize()
# print(ex)
# gene.draw_region()

for key in fasta_dictionary:
    gene = region(key, fasta_dictionary[key])
    gene.normalize()
    gene.draw_region()

