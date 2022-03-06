#!/usr/bin/env python

IUPAC={
    'A':"A",
    'C':"C",
    'G':"G",
    'T':"T",
    'U':'U',
    'W':'[AT]',
    'S':'[CG]',
    'M':'[AC]',
    'K':'[GT]',
    'R':'[AG]',
    'Y':"[CT]",
    'B':'[CGT]',
    'D':'[AGT]',
    'H':'[ACT]',
    'V':'[ACG]',
    'N':'[ACGT]'
}

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

if __name__ == '__bioinfo__':
	assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
	assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
	assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
	assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"




def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

if __name__ == '__bioinfo__':
	assert gc_content("GCGCGC") == 1
	assert gc_content("AATTATA") == 0
	assert gc_content("GCATGCAT") == 0.5




def convert_phred(letter):
    """Converts a single character into a phred score"""
    return(ord(letter)-33)

if __name__ == '__bioinfo__':
	assert convert_phred("I") == 40, "wrong phred score for 'I'"
	assert convert_phred("C") == 34, "wrong phred score for 'C'"
	assert convert_phred("2") == 17, "wrong phred score for '2'"
	assert convert_phred("@") == 31, "wrong phred score for '@'"
	assert convert_phred("$") == 3, "wrong phred score for '$'"



def average_qual_score(phred_score):
    '''Calculates the average quality score of an entire quality score string'''
    sum = 0
    
    #generates the sum
    for letter in phred_score:
        sum+=(convert_phred(letter))
    
    #calculates the average score
    average = sum/len(phred_score)
    
    return(average)

if __name__ == '__bioinfo__':
	assert average_qual_score('AABB') == 32.5, 'wrong average quality score'



def is_it_a_fasta_header (line):# determines if the input file line is a fasta header
	import re
	return True if len(re.findall('^>',line)) > 0 else False # if the line begins with > return True, else return False

if __name__ == '__bioinfo__':
	assert is_it_a_fasta_header('>adklfj;dlakjfd') == True, 'not a fasta header'
	assert is_it_a_fasta_header('ajfl;dsajf') == False, 'fasta header'




def fasta_reader (fasta_file):#this loops through an input fasta file line by line
# if the line is a header it populates the fasta_dictionary with a key = line
# if the line is not header it begins building a nucelotide string under previous header and continues until it runs into another header

    f=open(fasta_file,'r') #open the fasta file

    fasta_dictionary={}

    while True: 
        line=(f.readline()).strip()
        if line == '':
            break
        if is_it_a_fasta_header(line) == True: #this loop creates a key out of fasta header
            key=line
            fasta_dictionary[key]=''
			

        if is_it_a_fasta_header(line) == False: #this loop appends nucelotide sequences to the most recent fasta header key
            if validate_base_seq(line) == True:
                fasta_dictionary[key]+=line
            else:
                print("invalid sequence")
    f.close() #close fasta file
    return(fasta_dictionary)

if __name__ == '__bioinfo__':
	assert len(fasta_reader('fasta_reader_test.fa')) == 2; 'wrong fasta dictionary length'
	assert len(fasta_reader('fasta_reader_test.fa')['>NODE_2']) == 5, 'wrong nucleotide length'

