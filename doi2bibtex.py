#!/usr/bin/env python2.7
"""
Take one argument -- the doi, and convert it to bibtex using an API
call to dx.doi.org.
"""
from sys import argv
import os


if argv[0].find('doi') != -1:  # run from commandline as executable
    doi = argv[1]
else:  # run from python
    doi = argv[2]

cmd = (
    'curl -sLH "Accept: text/bibliography; style=bibtex" ' +
    'http://dx.doi.org/' + doi)
bib_oneliner = os.popen(cmd).read()

# convert bib_oneliner to formatted (multiline) bibtex
bib = ''

# extract type
entry_type = bib_oneliner[
    bib_oneliner.find('@') + 1: bib_oneliner.find('{')]

bib += '@' + entry_type + '{' + doi + ',\n'  # use doi as cite key
# parse body
body = bib_oneliner[bib_oneliner.find(',') + 2:-2] + ','
while body:
    # match curly braces
    left_minus_right = 0
    i = 0
    while True:
        if body[i] == '{':
            left_minus_right += 1
        if body[i] == '}':
            left_minus_right -= 1
            if left_minus_right == 0:
                # outermost level matched up, one entry finished
                # advance one char for the trailing comma
                i += 1
                break
        i += 1

    bib += '  ' + body[:i + 1] + '\n'
    body = body[i + 1:].strip()

bib += '}'
print(bib)
