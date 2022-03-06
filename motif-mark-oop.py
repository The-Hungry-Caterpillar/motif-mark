import cairo
import bioinfo
import re
import random


def get_args():
    import argparse
    parser = argparse.ArgumentParser(description = "this script outputs a map of motifs on a regions from a fasta file, exon sequence file, and motif marker files")
    parser.add_argument('-f', '--fasta_file', help='input fasta file to be mapped')
    parser.add_argument('-m', '--motifs_file', help='input motifs file')
    parser.add_argument('-d', '--dark', action = "store_true", help = 'dark mode')
    return parser.parse_args()
args=get_args() 


# grabs the contigs from input fasta file
# header : contig
fasta_dictionary = bioinfo.fasta_reader(args.fasta_file)

# get max contig for normalizing lengths later on
max_contig=(max([len(value) for value in fasta_dictionary.values()]))

# get list of motifs and assign each of them random fill and border color
with open(args.motifs_file) as f:
    motifs=f.read().splitlines()
    motifs_list = [ # [ (motif 1 sequence, fill color 1, border color 1), ..., (motif N sequence, fill color N, border color N) ] 
        (motif.upper(),
        (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1)), # assign random fill color for each motif
        (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1))) # assign random border color for each motif
        for motif in motifs
        ]
    
    # find the longest motif, used in creating legend
    max_motif_length=(len(max(motifs_list)[0])) 


# import dictionary of all IUPAC symbols, e.g. 'H':'[ACT]'
IUPAC = bioinfo.IUPAC 

class region:
    ''' 
    Creates class used for plotting introns and exons. 
    Needs normalizng (doodle class function) before plotting.`

    Input: 
        - header: Contig header (used for labeling plots) 
        - contig: Contig sequence where introns denoted by lowercase and exons denoted by uppercase.
    Output:
        - length: Length of contig (used for drawing last intron)
        - exons: list of positions, e.g.
            [ [exon1 start, exon1 end], [exon2 start, exon2 end], ..., [exonN start, exonN end] ] 
    '''

    __slots__ = ['contig', 'length', 'exons','introns', 'header']

    def __init__(self, header, contig):
        self.header = header
        self.contig = contig
        self.length=len(contig)
        self.exons = self.find_exons()

    def find_exons(self):
        '''
        Looks through contig and returns start and stop positions of exons as a list
        '''
        i=0
        exons = []
        while True:
            try:
                # update index of contig while the contig is lowercase (intron)
                while self.contig[i].islower():
                    i+=1
            # if the contig is over, break
            except IndexError:
                break

            # now that the contig is no longer lowercase, mark the start with the current index
            start=i

            try:
                # update index of contig while the contig is uppercase (exon)
                while self.contig[i].isupper():
                    i+=1

            # if the contig is over, record the final position as exon stop and break
            except IndexError:
                stop=i-1
                exons.append([start,stop])
                break
            # now the contig is no longer uppercase, mark the stop with the current index
            stop=i-1

            # append the start and stop position of the exon
            exons.append([start,stop])

        return(exons)


class motif:
    '''
    Creates class used for plotting positions of a given motif across a contig.
    Positions need normalizing (a doodle class function) before plotting.

    Input:
        - motif: Motif sequence, must be IUPAC nucleic acid notation (https://en.wikipedia.org/wiki/Nucleic_acid_notation)
        - contig: Contig, most likely contig attribute from region class 
        - color: Motif fill color for plotting
        - border_color: Motif border color for plotting
    Output:
        - motif_positions: list of motif start and stop positions, including overlap, e.g.,
            [ [start1, stop1], [start2, stop2], ... , [startN, stopN] ] 
    '''

    __slots__=[ 'motif_positions','contig', 'motif', 'color','border_color']
    
    def __init__(self, motif, contig, colors, border_color):
        self.contig=contig.upper()
        self.motif=motif
        self.motif_positions=self.find_positions() 
        self.color = colors
        self.border_color = border_color


    def find_positions(self):
        motif_positions=[]

        # break the motif into list, and sub each letter for IUPAC replacement
        reg_list = [IUPAC[letter] for letter in self.motif.upper()] 
        
        # join the subbed list together into a string
        reg_pattern = ''.join(reg_list) 
        
        left = 0
        while True:
            # find the regex pattern in spliced string (splice by variable 'left')
            match = re.search(reg_pattern, self.contig[left:])
            
            # if there is no match, break
            if not match:
                break
            
            # update the splice position
            left += match.start()+1
            
            # append the [start, stop] to motif positions list
            motif_positions.append([left-1, left+len(self.motif)-1])

        return(motif_positions)


class doodle:
    '''
    Cairo drawing class which scales with number of items in legend (because legend is horizontal). 
    Plots multiple contigs and features on one surface.

    Input:
        - Longest motif, used for determining legend spacing
        - Max contig length, used for normalizing all feature lengths
    '''

    # __slots__ = ['level','surface_width','surface_height', 'font_size']

    def __init__(self, max_motif_length, max_contig):
        # since the legend is horizontal the width of the surface scales with number of motifs
        self.surface_width = 200 * len(motifs_list) + 200
        self.window_size = int(self.surface_width/5)
        self.surface_height = self.window_size + self.window_size*len(fasta_dictionary)
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
        '''
        draws a motif at specified location along contig with predifined fill and border colors
        '''
        c.rectangle(start, height-5, stop-start, 5)
        c.set_source_rgb(motif_class.color[0], motif_class.color[1], motif_class.color[2])
        c.fill_preserve()
        c.set_source_rgb(motif_class.border_color[0], motif_class.border_color[1], motif_class.border_color[2])
        c.set_line_width(1)
        c.stroke()

    def stack_motifs(self, motif_class, level, c):
        '''
        Uses self.draw motif function to draw all instances of a motif along a contig on a given cairo surface, c.
        This function draws the motifs at a random height around the contig as to allow for overlapping motifs.
        '''
        height = level
        for i in range(0, len(motif_class.motif_positions)):
            start=motif_class.motif_positions[i][0]
            stop=motif_class.motif_positions[i][1]

            height=random.randint(int(level - self.window_size/7), int(level + self.window_size/7))
            self.draw_motif(motif_class, height, start, stop, c)
                    
    def draw_intron(self, start, stop, level, c):
        '''
        This function draws the introns of a contig on a given cairo surface, c.
        '''
        c.line_to(start, level)
        c.line_to(stop, level)
        c.set_source_rgb(0, 0, 0) if not args.dark else c.set_source_rgb(1, 1, 1)
        c.set_line_width(1)
        c.stroke()



    def draw_exon(self, start, stop, fill_colors, border_colors, level, c):
        '''
        This fuction draws the exons of a contif onto a given cairo surface, c.
        '''
        c.rectangle(start, level-int(self.window_size/4), stop-start, 2*int(self.window_size/4))
        c.set_source_rgb(fill_colors[0],fill_colors[1],fill_colors[2])
        c.fill_preserve()
        c.set_source_rgb(border_colors[0], border_colors[1], border_colors[2])
        c.set_line_width(1)
        c.stroke()

    def normalize(self, thing):
        '''
        This function scales the locations of a given feature to the largest contig.
        '''

        # transform featurelist
        for i in range(0,len(thing)):
            thing[i] = (
                round(self.surface_width * thing[i][0]/self.max_contig, 0),
                round(self.surface_width * thing[i][1]/self.max_contig, 0)
                )

    def draw_legend(self):
        '''
        This function draws a legend box. Scales with number of features being plotted.
        '''
        c.rectangle(self.legend_start, self.legend_top, self.legend_end, self.legend_bottom)
        c.set_source_rgb(0, 0, 0) if not args.dark else c.set_source_rgb(1, 1, 1)
        c.set_line_width(1)
        c.stroke()

    def draw_label(self, feature, i, fill_colors, border_colors,c):
        '''
        This function draws the labels within the legend box onto a cairo surface, c.
        Uses predefined fill and border colors.
        '''
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
        c.set_source_rgb(0,0,0) if not args.dark else c.set_source_rgb(1, 1, 1)
        c.show_text(feature)
        c.stroke()

    def draw_name(self, contig_name, level, c):
        '''This function draws the name of the contig above contig on a given cairo surface, c.'''
        c.move_to(0, int(level-self.window_size/3))
        c.set_source_rgb(0,0,0) if not args.dark else c.set_source_rgb(1, 1, 1)
        c.show_text(contig_name)
        c.stroke()

draw = doodle(max_motif_length,max_contig)
with cairo.SVGSurface('example.svg', draw.surface_width, draw.surface_height) as surface:
    c = cairo.Context(surface)
    c.set_font_size(draw.font_size)

    # draw legend stuff
    # legend box
    draw.draw_legend() 
    label_number=0
    # legend exon
    draw.draw_label("Exon", label_number, (0.5,0.5,0.5), (0,0,0), c) if not args.dark else draw.draw_label("Exon", label_number, (0.8,0.8,0.8), (1,1,1), c) 
    label_number+=1
    
    # put motifs on legend
    for motiv in motifs_list:
        draw.draw_label(motiv[0], label_number, motiv[1], motiv[2], c)
        label_number+=1


    for h,key in enumerate(fasta_dictionary):
        
        # the h factor separates contigs vertically in cairo space
        level = h*draw.window_size + draw.window_size

        # initialize a region object based on the current contig
        gene = region(key, fasta_dictionary[key])

        # normalize all the exons and total length of current contig
        draw.normalize(gene.exons)
        length=round(draw.surface_width * gene.length/draw.max_contig, 0)

        # draw name of contig above contig
        draw.draw_name(key,level, c)

        # draw first intron 
        draw.draw_intron(0, gene.exons[0][0], level, c)

        # draws all introns and exons of a contig
        for i in range( 0, len(gene.exons) ): 
            try:
                draw.draw_exon(gene.exons[i][0], gene.exons[i][1], (0.5, 0.5, 0.5), (0, 0, 0), level, c) if not args.dark else draw.draw_exon(gene.exons[i][0], gene.exons[i][1], (0.8, 0.8, 0.8), (1, 1, 1), level, c)
                draw.draw_intron(gene.exons[i][1], gene.exons[i+1][0], level, c)
            
            except IndexError: # draws the final intron. If exons is final feature this intron will be length 0
                draw.draw_intron(gene.exons[i][1], length, level, c)

        # draws all motifs associated with current contig
        for motiv in motifs_list:
            moti = motif(motiv[0], gene.contig, motiv[1], motiv[2])
            draw.normalize(moti.motif_positions)
            draw.stack_motifs(moti, level,c)


        surface.write_to_png(args.fasta_file + '.png')