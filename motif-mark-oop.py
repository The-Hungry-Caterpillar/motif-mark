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



# def exons_builder(file, contig):
#     '''generates a list of exons sorted by start and stop relative to the input contif
#     input: a text file with exon seqeuences seperated by new lines
#     return: a sorted list of exons as described below'''

#     with open(file) as e:
#         list_of_exons = e.read().splitlines()
        
#         # generates a list: [ [exon start, exon stop], [ exon start, exon stop ], ... ]
#         exons=[]
#         for i in range(0,len(list_of_exons)):
#             start = contig.find(list_of_exons[i])
#             stop = len(list_of_exons[i]) + start
#             exons.append([start,stop])
#     exons.sort()
#     return(exons)












intron1='atgtccacatgtagtcacgtttgacatcccagggccacctcagcaggccgtctctggggagaattttctctgatttcttccccttcccttgctggacccctgcacctgctggggaagatgtagctcactccgtctagcaagtgatgggagcgagtggtccagggtcaaagccagggtgcccttactcggacacatgtggcctccaagtgtcagagcccagtggtctgtctaatgaagttccctctgtcctcaaaggcgttggttttgtttccacag'
exon1='AAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAG'
intron2='gtatgactcacctgtgcgacccctggtgcctgctccgcgcagggccggcggcgtgccaggcagatgcctcggagaacccaggggtttctgtggctttttgcatgcggcgggcagctgtgctggagagcagatgcttcaccaattcagaaatccaatgccttcactctgaaatgaaatctgggcatgaatgtggggagaaaccttcactaacacactcttgctaaaacatagaatca'
print(len(intron1+exon1+intron2))
print(f'intron 1 start = 0, stop = {len(intron1)}')
print(f'exon 1 start = {len(intron1)}, stop = {len(intron1)+len(exon1)}')
print(f'intron 2 start = {len(intron1)+len(exon1)}, stop = {len(intron1)+len(exon1)+len(intron2)}')



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
        self.surface_width = 100 
        self.surface_height = 100
        self.level=50

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

    def draw_region(self): #currently not working
        '''draws an entire region using 'draw_exon' and 'draw_intron' fucntions
        input: exon, likely from a 'region.exons' object'''

        with cairo.SVGSurface('example.svg', self.surface_width, self.surface_height) as surface:
            c = cairo.Context(surface)
            c.move_to(0, self.level)

            for i in range( 0, len(self.exons) ): 
                
                if self.exons[0][0] == 0:
                    pass
               
                else:
                    
                    self.draw_intron(0, self.exons[0][0], c)

                    # not all exons matched to the contig, this loop skips such exons
                    if -1 in self.exons[i]:
                        place = c.get_current_point()
                        self.draw_intron(place[0], self.exons[i+1][0], c)
                    
                    else:
                        try:
                            self.draw_exon(self.exons[i][0], self.exons[i][1], c)
                            self.draw_intron(self.exons[i][1], self.exons[i+1][0], c)
                        
                        except IndexError: # draws the final intron. If exons is final feature this intron will be length 0
                            self.draw_intron(self.exons[i][1], self.surface_width, c)
            
            surface.write_to_png(self.header + '.png')


class draw:
    def __init__(self):
        self.surface_height = 100
        self.surface_width = 100
        self.level = 50 

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

        c.rectangle(start, self.level-10, stop-start, 10)
        c.set_source_rgb(0, .5, .5)
        c.fill()

    def draw_region(self, header, exons): #currently not working
        '''draws an entire region using 'draw_exon' and 'draw_intron' fucntions
        input: exon, likely from a 'region.exons' object'''

        with cairo.SVGSurface('example.svg', self.surface_width, self.surface_height) as surface:
            c = cairo.Context(surface)
            c.move_to(0, self.level)

            for i in range( 0, len(exons) ): 
                
                if exons[0][0] == 0:
                    pass
               
                else:
                    
                    self.draw_intron(0, exons[0][0], c)

                    # not all exons matched to the contig, this loop skips such exons
                    if -1 in exons[i]:
                        place = c.get_current_point()
                        self.draw_intron(place[0], exons[i+1][0], c)
                    
                    else:
                        try:
                            self.draw_exon(exons[i][0], exons[i][1], c)
                            self.draw_intron(exons[i][1], exons[i+1][0], c)
                        
                        except IndexError: # draws the final intron. If exons is final feature this intron will be length 0
                            self.draw_intron(exons[i][1], self.surface_width, c)
            
            surface.write_to_png(header + '.png')


gene = region('>test', fasta_dictionary['>test'])
ex = gene.exons
print(ex)
gene.normalize()
print(ex)
gene.draw_region()


# class called draw doesn't work? whyyy tf?
try:
    draw.draw_region('>test', ex)
except TypeError:
    print("why tf doesn't this work")

# region_dictionary ={
#     # contig header: ( contig_length, exon_dictionary generated from exons builder )
# }

# for key in fasta_dictionary:
#     x = region(fasta_dictionary[key], args.exon_file)
#     region_dictionary[key]= (x.length, x.exons)


# ################################################
# ######## begin pycairo #########################
# ################################################


# # set the surface dimensions
# surface_width = 100 #NEED TO SET THIS TO MAX CONTIG LE`NGTH
# scale = 20
# surface_height = len(region_dictionary)*scale 









# # normalize the exons within the region dictionary
# for key in region_dictionary:
#     contig_length = region_dictionary[key][0]
#     exons= region_dictionary[key][1]
#     normalize(exons, contig_length, surface_width)



# with cairo.SVGSurface('example.svg', surface_width, surface_height) as surface:
#     c = cairo.Context(surface)
    
#     for i, key in enumerate(region_dictionary):

#         # set the y height for each region
#         level = (scale*i) + (surface_height/(2 * len(region_dictionary)))

#         # draw each region
#         draw_region(region_dictionary[key][1], level)


#     
