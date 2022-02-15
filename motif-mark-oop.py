import cairo
import bioinfo

# test_mode=False
# def get_args():
#     import argparse
#     parser = argparse.ArgumentParser(description = "still need description")
# 	#parser.add_argument('- command line variable', '-- python variable', description)
#     parser.add_argument('-f', '--file', help='sam filename to be deduped')
#     parser.add_argument('-u', '--umi', help='name of umi file')
#     return parser.parse_args()
# args=get_args()

fasta_dictionary=bioinfo.fasta_reader('test.fa') 
length_of_contig=len(fasta_dictionary['>test_contig'])


def exons_builder(file, contig):
    '''generates a dictionary of exons as described below
    input: a text file with exon seqeuences seperated by new lines
    return: a dictionary of exons as described below'''
    exons = {
        # exon number: (start position, stop position)
        }
    with open(file) as e:
        list_of_exons = e.read().splitlines()
        exons=[]
        for i in range(0,len(list_of_exons)):
            start = contig.find(list_of_exons[i])
            stop = len(list_of_exons[i]) + start
            exons.append([start,stop])
    exons.sort()
    return(exons)

class region:
    '''desc here'''
    __slots__ = ['contig', 'length', 'exons']
    def __init__(self, contig, exon_file):
        self.contig = contig
        self.length=len(contig)
        self.exons = exons_builder('exons.txt', contig)

region_dictionary ={
    # contig header: ( contig_length, exon_dictionary generated from exons builder )
}

for key in fasta_dictionary:
    x = region(fasta_dictionary[key],'exons.txt')
    region_dictionary[key]= (x.length, x.exons)


################################################
######## begin pycairo #########################
################################################


# set the surface dimensions
surface_width = 100 #NEED TO SET THIS TO MAX CONTIG LENGTH
scale = 20
surface_height = len(region_dictionary)*scale 


def normalize(exons, length_of_contig, surface_width):
    '''function that transforms exon genome coordinates into cairo coordinates
    input: exon dictionary, cairo surface width
    output: exon dictionary with transformed coordinates'''

    # transform exon dictionary
    for i in range(0,len(exons)):
        exons[i] = (
            round(surface_width * exons[i][0]/length_of_contig, 0),
            round(surface_width * exons[i][1]/length_of_contig, 0)
            )

def draw_intron(start, stop, level):
    '''draws introns on cairo surface.
    input: intron start (i.e. previous exons stop), intron stop (i.e. next exon start), y coordinate of cairo surface
    output: a line representing intron, proportional to actual intron'''
    c.line_to(start, level)
    c.line_to(stop,level)
    c.set_source_rgb(1, 0, 0)
    c.set_line_width(1)
    c.stroke()

def draw_exon(start, stop, level):
    '''draws exons on cairo surface.
    input: exon start, exon stop, y coordinate of cairo surface (same as draw_intron)
    output: a rectangle proportional to actual exon'''
    c.rectangle(start, level-scale/4, stop-start, 2*scale/4)
    c.set_source_rgb(0, .5, .5)
    c.fill()

def draw_region(exons, level):
    '''draws an entire region using 'draw_exon' and 'draw_intron' fucntions
    input: '''
    c.move_to(0, level)
    draw_intron(0, exons[0][0], level)
    for i in range( 0, len(exons) ): 
       
       # not all exons matched to the contig, this loop skips such exons
        if -1 in exons[i]:
            place = c.get_current_point()
            draw_intron(place[0], exons[i+1][0], level)

        else:
            try:
                draw_exon(exons[i][0], exons[i][1] , level)
                draw_intron(exons[i][1], exons[i+1][0], level)
            except IndexError: # draws the final intron. If exons is final feature this intron will be length 0
                draw_intron(exons[i][1], surface_width, level)



# normalize the exons within the region dictionary
for key in region_dictionary:
    contig_length = region_dictionary[key][0]
    exons= region_dictionary[key][1]
    normalize(exons, contig_length, surface_width)


with cairo.SVGSurface('example.svg', surface_width, surface_height) as surface:
    c = cairo.Context(surface)
    
    for i, key in enumerate(region_dictionary):
        print(region_dictionary[key])
        # set the y height for each region
        level = (scale*i) + (surface_height/(2 * len(region_dictionary)))

        # draw each region
        draw_region(region_dictionary[key][1], level)

    # draw_region(region_dictionary['>test_contig3'][1], 50)

    surface.write_to_png('test.png')
