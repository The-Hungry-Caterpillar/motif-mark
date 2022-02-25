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

max_contig=(max([len(value) for value in fasta_dictionary.values()]))


with open(args.motifs_file) as f:
    motifs=f.read().splitlines()
    motifs_list = [(motif.upper(),
        (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1)), # assign random fill color for each motif
        (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1))) # assign random border color for each motif
        for motif in motifs
        ]
    max_motif_length=(len(max(motifs_list)[0])) # find the longest motif, used in creating legend

IUPAC = bioinfo.IUPAC # imports dictionary of all IUPAC symbols, e.g. 'H':'[ACT]'

class region:
    ''' insert docstring '''

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
    __slots__=[ 'motif_positions','contig', 'motif', 'color','border_color']
    
    def __init__(self, motif, contig, colors, border_color):
        self.contig=contig.upper()
        self.motif=motif
        self.motif_positions=self.find_positions() 
        self.color = colors
        self.border_color = border_color


    def find_positions(self):
        motif_positions=[]
        reg_list = [IUPAC[letter] for letter in self.motif.upper()] # break the motif into list, and sub each letter for IUPAC replacement
        reg_pattern = ''.join(reg_list) # join the subbed list together into a string
        # motif_positions = [[m.start(),m.end()] for m in re.finditer(reg_pattern, self.contig)]
        left = 0
        while True:
            match = re.search(reg_pattern, self.contig[left:])
            if not match:
                break
            left += match.start()+1
            motif_positions.append([left, left+len(self.motif)])

        return(motif_positions)

class doodle:
    '''still needs docstring'''

    # __slots__ = ['level','surface_width','surface_height', 'font_size']

    def __init__(self, max_motif_length, max_contig):
        # self.surface_width=400
        self.surface_width = 200 * len(motifs_list) + 200
        self.surface_height = 200 + 200*len(fasta_dictionary)
        self.font_size = 20
        self.label_width=max_motif_length*self.font_size
        self.max_contig=max_contig
        
        # legend parameters
        self.legend_start=0
        self.legend_end = self.surface_width
        self.legend_top = 5
        self.legend_bottom = 90
        self.legend_middle = 45
        self.spacer=self.font_size/3
    
    def draw_motif(self, motif_class, height, start, stop, c):
        c.rectangle(start, height-5, stop-start, 5)
        c.set_source_rgb(motif_class.color[0], motif_class.color[1], motif_class.color[2])
        c.fill_preserve()
        c.set_source_rgb(motif_class.border_color[0], motif_class.border_color[1], motif_class.border_color[2])
        c.set_line_width(1)
        c.stroke()

    def stack_motifs(self, motif_class, level, c):
        height = level
        for i in range(0, len(motif_class.motif_positions)):
            start=motif_class.motif_positions[i][0]
            stop=motif_class.motif_positions[i][1]

            height=random.randint(int(level-30), int(level+30))
            self.draw_motif(motif_class, height, start, stop, c)
                    
    def draw_intron(self, start, stop, level, c):
        '''draws introns on cairo surface.
        input: intron start (i.e. previous exons stop), intron stop (i.e. next exon start), y coordinate of cairo surface
        output: a line representing intron, proportional to actual intron'''

        c.line_to(start, level)
        c.line_to(stop, level)
        c.set_source_rgb(0, 0, 0)
        c.set_line_width(1)
        c.stroke()

    def draw_name(self, contig_name, level, c):
        c.move_to(0, level-60)
        c.set_source_rgb(0,0,0)
        c.show_text(contig_name)
        c.stroke()

    def draw_exon(self, start, stop, fill_colors, border_colors, level, c):
        '''draws exons on cairo surface.
        input: exon start, exon stop, y coordinate of cairo surface (same as draw_intron)
        output: a rectangle proportional to actual exon'''

        c.rectangle(start, level-30, stop-start, 60)
        c.set_source_rgb(fill_colors[0],fill_colors[1],fill_colors[2])
        c.fill_preserve()
        c.set_source_rgb(border_colors[0], border_colors[1], border_colors[2])
        c.set_line_width(1)
        c.stroke()

    def normalize(self, thing):
        '''need doc string'''

        # transform featurelist
        for i in range(0,len(thing)):
            thing[i] = (
                round(self.surface_width * thing[i][0]/self.max_contig, 0),
                round(self.surface_width * thing[i][1]/self.max_contig, 0)
                )

    def draw_legend(self):
        c.rectangle(self.legend_start, self.legend_top, self.legend_end, self.legend_bottom)
        c.set_source_rgb(0, 0, 0)
        c.set_line_width(1)
        c.stroke()

    def draw_label(self, feature, i, fill_colors, border_colors,c):
        c.rectangle(self.legend_start + i*self.label_width + self.spacer, #left
             self.legend_middle-7, #top
             self.font_size, #right
             20 #bottom
             )
        c.set_source_rgb(fill_colors[0],fill_colors[1],fill_colors[2])
        c.fill_preserve()
        c.set_source_rgb(border_colors[0], border_colors[1], border_colors[2])
        x=c.get_current_point()[0]
        c.stroke()        
        c.move_to(x + self.spacer + self.font_size, self.legend_middle + .5*self.font_size)
        c.set_source_rgb(0,0,0)
        c.show_text(feature)
        c.stroke()

draw = doodle(max_motif_length,max_contig)
with cairo.SVGSurface('example.svg', draw.surface_width, draw.surface_height) as surface:
    c = cairo.Context(surface)
    # c.move_to(0, draw.level)
    draw.draw_legend()
    # set font parameters
    c.set_source_rgb(0, 0, 0)
    c.set_font_size(draw.font_size)
    c.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    
    label_number=0
    draw.draw_label("Exon", label_number, (0.5,0.5,0.5), (0,0,0), c)
    label_number+=1
    for motiv in motifs_list:
        draw.draw_label(motiv[0], label_number, motiv[1], motiv[2], c)
        label_number+=1

    for h,key in enumerate(fasta_dictionary):
        
        level = h*200 + 250
        gene = region(key, fasta_dictionary[key])
        draw.normalize(gene.exons)
        length=round(draw.surface_width * gene.length/draw.max_contig, 0)

        draw.draw_name(key,level, c)
        draw.draw_intron(0, gene.exons[0][0], level, c)

        for i in range( 0, len(gene.exons) ): 
            try:
                draw.draw_exon(gene.exons[i][0], gene.exons[i][1], (0.5, 0.5, 0.5), (0, 0, 0), level, c)
                draw.draw_intron(gene.exons[i][1], gene.exons[i+1][0], level, c)
            
            except IndexError: # draws the final intron. If exons is final feature this intron will be length 0
                draw.draw_intron(gene.exons[i][1], length, level, c)

        for motiv in motifs_list:
            moti = motif(motiv[0], gene.contig, motiv[1], motiv[2])
            draw.normalize(moti.motif_positions)
            draw.stack_motifs(moti, level,c)


        surface.write_to_png(args.fasta_file + '.png')

 
