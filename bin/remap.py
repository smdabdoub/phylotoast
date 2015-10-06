import sys

def parseMaps(goodMapFile, badMapFile):
    with open(goodMapFile) as gmf:
        goodMap = gmf.readlines()
        del goodMap[0]

    with open(badMapFile) as bmf:
        badMap = bmf.readlines()
        del badMap[0]

    bcMap = {}

# create a map from the sample ID to the good barcode
    for line in goodMap:
        line = line.split('\t')
        bcMap[line[0]] = line[1]

    for line in badMap:
        line = line.split('\t')
        if line[0] in bcMap:
            bcMap[line[0]] = (line[1], bcMap[line[0]])

return dict([t for t in bcMap.values() if isinstance(t, tuple)])


if __name__ == '__main__':
    goodMapFile = sys.argv[1]
    badMapFile = sys.argv[2]
    fastaFile = sys.argv[3]

    bcMap = parseMaps(goodMapFile, badMapFile)

    newFasta = open('new.fasta', 'w')

    with open(fastaFile) as f:
        for line in f:
            if line[0] == '>':
                header = line
                newFasta.write(header)
            else:
                # sequence info
                seq = line
                badCode = seq[:8]
                if badCode in bcMap:
                    newFasta.write(bcMap[badCode] + seq[8:])
                else:
                    print 'Mapping error: ' + header

    newFasta.close()
