#from collections import defaultdict
import sys

def parseMapFile(mapFN, keyColumn=0):
    """
    Opens a QIIME mapping file and stores the contents in a dictionary
    keyed on SampleID (default) or a user-supplied one. The only 
    required fields are SampleID, BarcodeSequence, LinkerPrimerSequence 
    (in that order), and Description (which must be the final field).
    
    :@param mapFN: Full path to the map file
    :@param keyColumn: The column to key the returned dictionary on.
                       If keyColumn is larger than the number of columns,
                       it will be set to the default, 0.     
    
    Example data:
    #SampleID    BarcodeSequence    LinkerPrimerSequence    Treatment    SampleType    DOB    Description
    11.V13    ACGCTCGACA    GTTTGATCCTGGCTCAG    V13    Rat_Oral_Disease    111111    Rat_Oral
    """
    m = {}
    
    with open(mapFN) as mapfile:
        lines = mapfile.readlines()
        if keyColumn >= len(lines[0].split()): keyColumn = 0
        for line in lines[1:]:
            line = line.strip().split('\t')
            m[line[keyColumn]] = line
            
    return m


from Bio import SeqIO
from Bio.Seq import Seq
def barcodeReplace(mapFN, inFN, outFN):
    """
    In a mapping file of two primer sets, replace the barcodes from one with the other
    """
    orig_map = parseMapFile(mapFN)
    bmap = {}
    for key in orig_map:
        sample, v = key.split('.') 
        if v == 'v79':
            bmap[orig_map[key][1]] = orig_map[sample+'.V13'][1]
            
    
    # replace the v79 code with the v13 code
    outF = open(outFN,'w')
    out2 = open('/data/badcodes.fna', 'w')
    inF = open(inFN, 'rU')
    count = 0
    countGood = 0
    for record in SeqIO.parse(inF, 'fasta'):
        if str(record.seq)[:10] in bmap:
            record.seq = Seq(bmap[str(record.seq)[:10]]) + record.seq[10:]
            outF.write(record.format('fasta'))
            countGood += 1
        else:
            count += 1
            out2.write(record.format('fasta'))
        
    
    print str(countGood) + ' good records.'
    print str(count) + ' bad records.'
    inF.close()
    outF.close()
    out2.close()

        
def fastaSizeFilter(inFN, outFN, minLen):
    out = open(outFN, 'w')
    
    with open(inFN) as fasta:
        for line in fasta:
            if line[0] == '>':
                descr = line
            else:
                if len(line) >= minLen:
                    out.write(descr)
                    out.write(line)
    
    out.close()
    


def fnaConvertTo(inFN):
    counter = 1
    printSeq = False
    with open(inFN) as fasta:
        for line in fasta:
            if line[0] == '>':
                try:
                    if int(line[1]) <= 4:
                        info = line[1:].split()[0].split('_')
                        sampleID = info[0]
                        hID = info[1]
                        print '>HNS%s_%i %s orig_bc=AAAAAAAA new_bc=AAAAAAAA bc_diffs=0\n' % (sampleID, counter, hID)
                        counter  += 1
                        printSeq = True
                        
                except ValueError:
                    pass
            else:
                if printSeq:
                    print(line.strip())
                printSeq = False



def fnaRenumber(inFN):
    counter = 1

    with open(inFN) as fasta:
        for line in fasta:
            if line[0] == '>':
                line = line[1:].split()
                sampleID = line[0].split('_')[0]
                print '>%s_%i %s\n' % (sampleID, counter, ' '.join(line[1:]))
                counter  += 1
            else:
                print(line.strip())








if __name__ == '__main__':
#    filterFN = sys.argv[1]
#    toFilterFN = sys.argv[2]
#    outFN = sys.argv[3]
#    
#    fastaFilter(filterFN, toFilterFN, outFN)
#    fastaSizeFilter(sys.argv[1], sys.argv[2], int(sys.argv[3]))
#    fnaConvertTo(sys.argv[1])
#    fnaRenumber(sys.argv[1])
    barcodeReplace(sys.argv[1], sys.argv[2], sys.argv[3])