import utils
from sys import argv


filename = argv[1]
outname = filename[:-3] + 'sorted.bed'


bedfile = utils.parseTable(filename, '\t')
out = []
for line in bedfile:

    coords = [int(line[1]), int(line[2])]
    start = min(coords)
    end = max(coords)

    newline = [line[0], start, end] + line[3:]
    out.append(newline)

utils.unParseTable(out, outname, '\t')
    
