#!/usr/bin/env python3

'''
extracts desired ultraconserved region using regular expressions
'''

import re, sys

args = sys.argv
file = open(args[1]).read()
uc_number = args[2]

uce = f'> uc\.{uc_number}\+'

match = re.search(r'(%s\n[ATCG]+)\n>' % uce, file)
final = match.group(1)

print(final)


