import cairo
import bioinfo
import re
import random


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
    motifs_list = [(motif.upper(), (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1))) for motif in motifs]
    max_motif_length=(len(max(motifs_list)[0]))

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
    __slots__=['motifs', 'motif_positions','contig', 'motif_name', 'color']
    
    def __init__(self, motif, contig, colors):
        self.contig=contig.upper()
        self.motifs=motif
        self.motif_name=motif
        self.motif_positions=self.find_positions() 
        self.color = colors

    def find_positions(self):
        motif_positions=[]
        reg_list = [IUPAC[letter] for letter in self.motifs.upper()] # break the motif into list, and sub each letter for IUPAC replacement
        reg_pattern = ''.join(reg_list) # join the subbed list together into a string
        for match in re.finditer(reg_pattern,self.contig):
            motif_positions.append([match.start(), match.end()])
        return(motif_positions)


class doodle:
    '''still needs docstring'''

    # __slots__ = ['level','surface_width','surface_height', 'font_size']

    def __init__(self, max_motif_length):
        self.level = 125
        # self.surface_width=400
        self.surface_width = 100 * len(motifs_list) + 200
        self.surface_height = 200
        self.font_size = self.surface_width/40
        
        # legend parameters
        self.legend_start=0
        self.legend_end = self.surface_width
        self.legend_top = self.surface_height/20
        self.legend_bottom = self.surface_height/5
        self.legend_middle = self.legend_top + (self.legend_bottom-self.legend_top)/2
        self.label_width=max_motif_length*self.font_size
        self.spacer=self.font_size/3
    
    def draw_motifs(self, motif_class, c):
        for i in range(0, len(motif_class.motif_positions)):
            start=motif_class.motif_positions[i][0]
            stop=motif_class.motif_positions[i][1]
            c.rectangle(start, self.level-8, stop-start, 16)
            c.set_source_rgb(motif_class.color[0], motif_class.color[1], motif_class.color[2])
            c.fill()

    def draw_intron(self, start, stop, c):
        '''draws introns on cairo surface.
        input: intron start (i.e. previous exons stop), intron stop (i.e. next exon start), y coordinate of cairo surface
        output: a line representing intron, proportional to actual intron'''

        c.line_to(start, self.level)
        c.line_to(stop, self.level)
        c.set_source_rgb(0, 0, 0)
        c.set_line_width(1)
        c.stroke()

    def draw_exon(self, start, stop, c):
        '''draws exons on cairo surface.
        input: exon start, exon stop, y coordinate of cairo surface (same as draw_intron)
        output: a rectangle proportional to actual exon'''

        c.rectangle(start, self.level-5, stop-start, 10)
  
        # setting color of the context for inside
        c.set_source_rgb(.8, .8, .8)
    
        # Preserving inside color of object
        c.fill_preserve()
    
        # setting color of the context for outline
        c.set_source_rgb(0, 0, 0.5)
    
        c.set_line_width(1)
    
        # stroke out the color and width property
        c.stroke()

    def normalize(self, thing):
        '''need doc string'''

        # transform featurelist
        for i in range(0,len(thing)):
            thing[i] = (
                round(self.surface_width * thing[i][0]/gene.length, 0),
                round(self.surface_width * thing[i][1]/gene.length, 0)
                )

    def draw_legend(self):
        c.rectangle(self.legend_start, self.legend_top, self.legend_end, self.legend_bottom)
        c.set_source_rgb(0, 0, 0)
        c.set_line_width(1)
        c.stroke()

    def draw_label(self, feature, i, color1, color2, color3,c):
        c.rectangle(self.legend_start + i*self.label_width + self.spacer, self.legend_top + .5*self.legend_middle, self.font_size, self.font_size)
        c.set_source_rgb(color1, color2, color3)
        c.fill_preserve()
        if feature=='Exon':
            c.set_source_rgb(0, 0, 0.5)
            c.set_line_width(1)
        x=c.get_current_point()[0]
        c.stroke()        
        c.move_to(x + self.spacer + self.font_size, self.legend_middle + .5*self.font_size)
        c.set_source_rgb(0,0,0)
        c.show_text(feature)
        c.stroke()

draw=doodle(max_motif_length)



for key in fasta_dictionary:
    gene = region(key, fasta_dictionary[key])
    draw.normalize(gene.exons)
    with cairo.SVGSurface('example.svg', draw.surface_width, draw.surface_height) as surface:
        c = cairo.Context(surface)
        c.move_to(0, draw.level)

        draw.draw_legend()

        # set font stuff
        c.set_source_rgb(0, 0, 0)
        c.set_font_size(draw.font_size)
        c.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        
        label_number=0
        draw.draw_label("Exon", label_number, .8,.8,.8, c)
        label_number+=1

        draw.draw_intron(0, gene.exons[0][0],c)

        for i in range( 0, len(gene.exons) ): 
            try:
                draw.draw_exon(gene.exons[i][0], gene.exons[i][1], c)
                draw.draw_intron(gene.exons[i][1], gene.exons[i+1][0], c)
            
            except IndexError: # draws the final intron. If exons is final feature this intron will be length 0
                draw.draw_intron(gene.exons[i][1], draw.surface_width, c)

        for motiv in motifs_list:
            moti = motif(motiv[0], gene.contig, motiv[1])
            if moti.motif_positions != []:
                draw.draw_label(motiv[0], label_number, motiv[1][0], motiv[1][1], motiv[1][2], c)
                label_number+=1
            draw.normalize(moti.motif_positions)
            draw.draw_motifs(moti, c)


        surface.write_to_png(gene.header + '.png')

 
