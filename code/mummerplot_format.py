#！mummerplot_format.py
# -*- coding：utf-8 -*-
"""
Created on 2/7/23 1:24 PM

@author: cuixue

IDE: PyCharm 
"""

import re
import argparse

###对>30 depth的ref.fa and qry.fa alignmented file, bam --> coverage do filtering
### condition：depth >2 and juge if have diff chr
def main(argparse):
    parser = argparse.ArgumentParser(description='Duplication block similar and oringnal finding')
    parser.add_argument('--cover', '-c', help='inputfile(coverage file)', required=True, type=str)

    parser.add_argument('--pos', '-p', help='coverage allfilter position file',\
                        type=str, required=True)
    args = parser.parse_args()
    return args.cover, args.pos


if __name__ == '__main__':
    covername, posname = main(argparse)
    coverfile = open(covername)
    posfile = open(posname, 'w')
    allchr = set()
    lastchr = ''
    lastpos = ''
    while 1:
        lines = coverfile.readlines(10000)
        if not lines:
            break
        for line1 in lines:
            line1 = line1.rstrip()
            cut1 = line1.strip().split('\t')
            if re.search('^>', line1):
                posfile.write(line1+'\n')
            else:
                posfile.write(' '+line1+'\n')