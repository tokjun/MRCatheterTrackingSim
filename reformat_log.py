#! /usr/bin/python

#
# reformat_log.py
#  Reformat log file to properly insert separators
#

import sys, csv
import operator


inputFileName = None
outputFileName = None
delimiter = ' '

for arg in sys.argv:
    param = arg.split('=', 2)
    if param[0] == 'IN':  # Input file
        inputFileName = param[1]
    if param[0] == 'OUT': # Output file
        outputFileName = param[1]
    if param[0] == 'DEL': # Deliminter (default = ' ')
        delimiter = param[1]

if inputFileName==None or outputFileName==None:
    print('Usage: %s IN=<input CSV file> OUT=<outputCSV file> KEY=<key> DIR=<direction>' % (sys.argv[0]))
    print('      <input CSV file>  : Input CSV file.')
    print('      <output CSV file> : Output CSV file.')
    sys.exit()


with open(inputFileName) as inputFile:
    csvReader = csv.reader(inputFile, delimiter=delimiter)

    with open(outputFileName, "wb") as outputFile:
        fileWriter = csv.writer(outputFile, delimiter=delimiter)
        for row in csvReader:
            outrow = []
            for col in row:
                if col.count('.') > 1:
                    digit = 0
                    str1 = ''
                    str2 = ''
                    for c in col:
                        if digit < 6:  # TODO: the number of digits are not necessarily 6
                            if c.isdigit():
                                digit = digit + 1
                            str1 = str1 + c
                        else:
                            str2 = str2 + c
                    outrow.append(str1)
                    outrow.append(str2)
                else:
                    outrow.append(col)
            fileWriter.writerow(outrow)



    
