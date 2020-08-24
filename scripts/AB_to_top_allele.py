__author__='mlin'
from optparse import OptionParser
import pandas as pd

USAGE = """
AB_to_top_allele.py --bim 
                    --strand <Illumina strand report>
                    --pos <want chr, pos be updated from the report or not, default y >
                    --out 
"""

parser = OptionParser(USAGE)
parser.add_option('--bim')
parser.add_option('--strand')
parser.add_option('--pos', default='y')
parser.add_option('--out')

(options,args)=parser.parse_args()

print('Options in effect:\n--bim', options.bim, '\n--strand', options.strand, '\n--pos', options.pos, '\n--out', options.out, '\n')

outfile = open(options.out, 'w')

bim = pd.read_csv(options.bim, header=None, dtype=str, sep='\s+')
strand = pd.read_csv(options.strand, header=5, names = ["SNP_Name","Ilmn_ID", "Build","Chr","Coord","Forward_Seq","Forward_Allele1","Forward_Allele2"],dtype=str, sep='\s+', comment='#')
A_dict = dict(zip(strand['SNP_Name'], strand['Forward_Allele1']))
B_dict = dict(zip(strand['SNP_Name'], strand['Forward_Allele2']))

for key , value in A_dict.items():
    print(key, value)

if options.pos=='y':
   chr_dict = dict(zip(strand['SNP_Name'], strand['Chr']))
   pos_dict = dict(zip(strand['SNP_Name'], strand['Coord']))


for i in range(0, bim[0].size):
    if i%1000==0:
        print ('On line: ', i)
    A = A_dict.get(bim[1][i], 'NA')
    B = B_dict.get(bim[1][i], 'NA')
    if bim[4][i]=="A":
        bim[4][i]= A
    elif bim[4][i]=="B":
        bim[4][i]= B
    if bim[5][i]=="A":
        bim[5][i]= A
    elif bim[5][i]=="B":
        bim[5][i]= B
    if options.pos=='y':
        bim[0][i] = chr_dict.get(bim[1][i], bim[0][i])
        bim[3][i] = pos_dict.get(bim[1][i], bim[3][i])
    
    
    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (bim[0][i], bim[1][i], bim[2][i], bim[3][i], bim[4][i],bim[5][i]))            


