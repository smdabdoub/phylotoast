import sys
import json
from collections import defaultdict

def biom_summary(biom):
    tamtcounts = defaultdict(int)
    tot_seqs = 0.0
    
    for row, col, amt in biom['data']:
        tot_seqs += amt
        rtax = biom['rows'][row]['metadata']['taxonomy']
        for i,t in enumerate(rtax):
            t = t.strip()
            if i == len(rtax)-1 and len(t) > 3 and len(rtax[-1]) > 3:
                t = 's__'+rtax[i-1].strip().split('_')[-1]+'_'+t.split('_')[-1]
            tamtcounts[t] += amt        
    
    print 'Total sequence count: %i\n' % tot_seqs
    
    def levelData(counts, seqCount, level='g'):
        lvlData = sorted([(k, counts[k]) for k in counts if level+'_' in k], key=lambda x: x[1], reverse=True)
        return ['{}: {:4.2f}%'.format(p[0], p[1]/seqCount*100) for p in lvlData]
    
    def lookup(s, d):
        for idx, item in enumerate(d):
            if s in item:
                return d[idx]
        return s+': None'
    
    print 'Top families:\n'
    lvlDataF = levelData(tamtcounts, tot_seqs, 'f')
    print lvlDataF[:20]
            
    print '\nTop genera:\n'
    lvlDataG = levelData(tamtcounts, tot_seqs)
    print lvlDataG[:20]
            
    print '\nTop species:\n'
    lvlData = levelData(tamtcounts, tot_seqs, 's')
    print lvlData[:20]
    
    
    print '\nRed Complex:'
    print lookup('s__Porphyromonas_gingivalis', lvlData)
    print lookup('s__Treponema_denticola', lvlData)
    print lookup('s__Tannerella_forsythia', lvlData)
    
    print '\nOrange Complex: '
    print lookup('s__Fusobacterium_nucleatum', lvlData)
    print lookup('s__Prevotella_intermedia', lvlData)
    print lookup('s__Prevotella_nigrescens', lvlData)
    print lookup('s__Peptostreptococcus_micros', lvlData)
    
    print '\nYellow Complex: '
    print lookup('s__Streptococcus_sanguis', lvlData)
    print lookup('s__Streptococcus_oralis', lvlData)
    print lookup('s__Streptococcus_mitis', lvlData)
    print lookup('s__Streptococcus_gordonii', lvlData)
    print lookup('s__Streptococcus_intermedius', lvlData)
    
    print '\nGreen Complex: '
    print lookup('s__Aggregatibacter_actinomycetemcomitans', lvlData)
    print lookup('s__Campylobacter_concisus', lvlData)
    print lookup('s__Eikenella_corrodens', lvlData)
    print lookup('g__Capnocytophaga', lvlDataG)
    
    print '\nPurple Complex: '
    print lookup('s__Veillonella_parvula', lvlData)
    print lookup('s__Actinomyces_odontolyticus', lvlData)
    print lookup('s__Selenomonas_noxia', lvlData)
    print lookup('s__Actinomyces_naeslundii', lvlData)

def main():
    with open(sys.argv[1], 'rU') as bF:
        biom = json.loads(bF.readline())
        biom_summary(biom)
    
if __name__ == '__main__':
    main()