#!/usr/bin/python
import sys
try:
    input_bed = sys.argv[1]
    output_file_name = sys.argv[2]
except:
    print("Usage: ./")
output = open(output_file_name, 'w+')
count_dict={}
for x in range (0, 9075):
    count_dict.update({str(x+1):'0'})

with open (input_bed) as data:
    for line in data:
        new_count = 0
        line_list=(line.strip("\n")).split("\t")
        if (0<int(line_list[1])<=9075) and (line_list[5]=="+"):
            new_count= 1+int(count_dict[line_list[1]])
            count_dict[line_list[1]]=str(new_count)
for z in range (0, 9075):
    output.write('base_'+str(z+1) + '\t' + count_dict[str(z+1)] + '\n')

output.close
