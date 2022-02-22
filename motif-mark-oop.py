import cairo
import bioinfo
import re

def get_args():
    import argparse
    parser = argparse.ArgumentParser(description = "this script outputs a map of motifs on a regions from a fasta file, exon sequence file, and motif marker files")
    parser.add_argument('-f', '--fasta_file', help='input fasta file to be mapped')
    parser.add_argument('-m', '--motifs_file', help='input motifs file')
    return parser.parse_args()
args=get_args() 

fasta_dictionary=bioinfo.fasta_reader(args.fasta_file)
with open(args.motifs_file) as f:
    motifs=f.read().splitlines()
    motifs_list = [motif.upper() for motif in motifs]

IUPAC = bioinfo.IUPAC


class region:
    '''creates an object called region consisting of:
    - .contig: contig sequence
    - .length: length of contig
    - .exons: list of exons (start pos, stop pos) ordered by start pos
    input: contig sequence and exon file'''

    __slots__ = ['contig', 'length', 'exons','introns', 'header']

    def __init__(self, header, contig):
        self.header = header
        self.contig = contig
        self.length=len(contig)
        self.exons = self.find_exons()
        self.introns = self.find_exons()

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
          

class motif:
    __slots__=['motifs', 'motif_positions','contig']
    
    def __init__(self, motifs, contig):
        self.contig=contig
        self.motifs=motifs
        self.motif_positions=self.find_motifs()

    def find_motifs(self):
        motif_position=[]
        for motif in self.motifs: 
            reg_list = [IUPAC[letter.upper()] for letter in motif] # break the motif into list, and sub each letter for IUPAC replacement
            reg_pattern = ''.join(reg_list) # join the subbed list together into a string
            res = re.search(reg_pattern, self.contig, re.IGNORECASE) # search the above string in a contig
            try:
                res.span() # call re object
                start=res.start() # define start
                stop=res.end() # define end
                motif_position.append([int(start),int(stop)])
            except AttributeError:
                pass
        return(motif_position)


class doodle:
    '''still needs docstring'''

    __slots__ = ['level','surface_width','surface_height']

    def __init__(self):
        self.level=100
        self.surface_width=200
        self.surface_height=200
    
    def draw_motif(self, start, stop, c):
        c.rectangle(start, self.level-8, stop-start, 16)
        c.set_source_rgb(0.1, 0, .5)
        c.fill()

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

    def normalize(self, thing):
        '''need doc string'''

        # transform featurelist
        for i in range(0,len(thing)):
            thing[i] = (
                round(self.surface_width * thing[i][0]/gene.length, 0),
                round(self.surface_width * thing[i][1]/gene.length, 0)
                )

    def draw_region(self, gene_positions, motif_positions):
        '''draws an entire region using 'draw_exon' and 'draw_intron' fucntions
        input: exon, likely from a 'region.exons' object'''

        with cairo.SVGSurface('example.svg', self.surface_width, self.surface_height) as surface:
            c = cairo.Context(surface)
            c.move_to(0, self.level)
            
            self.draw_intron(0, gene_positions[0][0],c)

            for i in range( 0, len(gene_positions) ): 
                try:
                    self.draw_exon(gene_positions[i][0], gene_positions[i][1], c)
                    self.draw_intron(gene_positions[i][1], gene_positions[i+1][0], c)
                
                except IndexError: # draws the final intron. If exons is final feature this intron will be length 0
                    self.draw_intron(gene_positions[i][1], self.surface_width, c)

            for i in range( 0, len(motif_positions) ):
                self.draw_motif(motif_positions[i][0], motif_positions[i][1], c)
            
            surface.write_to_png(gene.header + '.png')


for key in fasta_dictionary:
    gene = region(key, fasta_dictionary[key])
    motifs = motif(motifs_list, gene.contig)
    draw=doodle()

    draw.normalize(gene.exons)
    draw.normalize(motifs.motif_positions)
    draw.draw_region(gene.exons, motifs.motif_positions)
 
