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
print(fasta_dictionary)

length_of_contig=12000

exons = { # exon number: (start position, stop position)
    1: (0,5000),
    2: (7000,9000),
    3: (11000, 12000)
}


################################################
######## begin pycairo #########################
################################################


def normalize(exons, surface_width):
    '''function that transforms exon genome coordinates into cairo coordinates
    input: exon dictionary, cairo surface width
    output: exon dictionary with transformed coordinates'''

    # transform exon dictionary
    for key in exons:
        exons[key] = (
            round(surface_width * exons[key][0]/length_of_contig, 0),
            round(surface_width * exons[key][1]/length_of_contig, 0)
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
    c.rectangle(start, level-5, stop-start, 10)
    c.set_source_rgb(1, 1, 0)
    c.fill()

# set the surface dimensions
surface_width = 100
surface_height = 100

normalize(exons, surface_width)

with cairo.SVGSurface('example.svg', surface_width, surface_height) as surface:
    c = cairo.Context(surface)
    level=surface_height/2
    c.move_to(0, level)

    # draw the first intron. If starts with exon then this intron will be length 0
    draw_intron(0, exons[1][0], level)

    for i in range( 1, 1+len(exons) ): 
        try:
            draw_exon(exons[i][0], exons[i][1] , level)
            draw_intron(exons[i][1], exons[i+1][0], level)
        except KeyError: # draws the final intron. If exons is final feature this intron will be length 0
            draw_intron(exons[i][1], surface_width, level)

    surface.write_to_png('test.png')
