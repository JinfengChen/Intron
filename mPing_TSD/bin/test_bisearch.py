from collections import defaultdict
import re

def binarySearch(data, val):
    highIndex = len(data)-1
    lowIndex = 0
    while highIndex > lowIndex:
            index = (highIndex + lowIndex) / 2
            sub = int(data[index])
            #print highIndex, index, lowIndex, sub, val
            if data[lowIndex] == val:
                    return [lowIndex, lowIndex]
            elif sub == val:
                    return [index, index]
            elif data[highIndex] == val:
                    return [highIndex, highIndex]
            elif sub > val:
                    if highIndex == index:
                            return sorted([highIndex, lowIndex])
                    highIndex = index
            else:
                    if lowIndex == index:
                            return sorted([highIndex, lowIndex])
                    lowIndex = index
    return sorted([highIndex, lowIndex])

#Chr1    DHS     DHsites 10537   10645   .       +       .       ID=1;
#Chr1    DHS     DHsites 10785   11022   .       +       .       ID=2;
def read_gff(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                idx  = '%s_%s_%s' %(unit[0], unit[3], unit[4])
                data[unit[3]] = idx
                data[unit[4]] = idx
    return data


snp  = 1010000
data = read_gff('../input/GSM655033_Rice_Seedling_DHsites.MSU7.Corrected.1chr.gff')
positions = sorted(data.keys(), key=int)
interval  = binarySearch(positions, snp)

dhs_s = positions[interval[0]]
dhs_e = positions[interval[1]]
print data[dhs_s]
print data[dhs_e]
