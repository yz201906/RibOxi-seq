#!/usr/bin/python
import sys
try:
    input_bed = sys.argv[1]
    output_file_name = sys.argv[2]
except:
    print("Usage: ./")
output = open(output_file_name, 'w+')
count_dict={}
for x in range (109078, 110947):
    count_dict.update({str(x):'0'})
for y in range (113348, 118418):
    count_dict.update({str(y):'0'})

with open (input_bed) as data:
    for line in data:
        new_count = 0
        line_list=(line.strip("\n")).split("\t")
        if line_list[0] == 'chrUn_gl000220':
            if (109078<int(line_list[1])<=110946 or 113348<int(line_list[1])<=118417):
                new_count= 1+int(count_dict[line_list[1]])
                count_dict[line_list[1]]=str(new_count)
for z in range (109078, 110947):
    output.write('18S_'+str(z-109077) + '\t' + count_dict[str(z)] + '\n')
for i in range (113348, 118418):
    output.write('28S_'+str(i-113348) + '\t' + count_dict[str(i)] + '\n')

output.close
