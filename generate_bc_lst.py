#!/usr/bin/env/python

from sys import argv

def create_barcodes_fasta(location, file_name):
    import os
    files_list = os.listdir(path=location)
    list1 = []
    for i in files_list:
        i = ''.join((item for item in i if not item.isdigit()))
        i = i[5:22]
        list1.append(i)
    list1.pop()
    list1.pop()
    set1 = set(list1)
    list2 = []
    for i in set1:
        i = i.split('_')
        list2.append(i)
    bc_list = []
    for i in list2:
        for j in i:
            bc_list.append(j)
    final_list = set(bc_list)
    os.chdir('/home/stud9/NAS/data')
    with open (file_name, 'w') as f:
        count = 0
        for i in final_list:
            count += 1
            f.write(f'>Adapter {count}: \n{i}\n')
            

path_to = argv[1]
file = argv[2]

create_barcodes_fasta(path_to, file)