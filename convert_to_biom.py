 '''
Created on Nov 15, 2012

@author: Shareef Dabdoub
'''
from collections import OrderedDict
from datetime import datetime
import sys

def generate_default_fields():
    return OrderedDict({'id': 'ethnicity_data',
            'format': 'Biological Observation Matrix 1.0',
            'format_url': 'http://biom-format.org/documentation/format_versions/biom-1.0.html',
            'type': 'OTU table',
            'generated_by': 'Mjolnir 0.1',
            'date': datetime.isoformat()})
    
def process_input(inF, fields):
    headerLines = 2
    columns = []
    rows = []
    groups = []
    
    for line in inF:
        if headerLines == 2:
            columns = line.strip().split()[1:]
            headerLines = 1
        if headerLines == 1:
            groups = list(set(line.strip().split()[1:]))
        

def write_BIOM_file(outF, fields):
    outF.write('{\n')
    for f in fields:
        outF.write('%s: %s,\n' % fields[f])

    outF.write('}\n')



if __name__ == '__main__':
    fields = generate_default_fields()
    with open(sys.argv[1]) as inF:
        process_input(inF,fields)
    with open('out.biom', 'w') as outF:
        write_BIOM_file(fields)