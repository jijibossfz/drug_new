# coding=utf-8
import re
import os

# files = '1mol2'
# reg = '@<TRIPOS>MOLECULE[\s|\S]*?ROOT'
# re_reg = re.compile(reg)
#
# with open(files, 'r') as f:
#     date = f.read()
#
# con = re_reg.findall(date)
#
# os.chdir('1mol2')
# for mol in con:
#     file_name = mol.split('\n')[1] + '.mol2'
#     with open(file_name, 'w') as w:
#         w.write(mol)

# file_lst = os.listdir(files)
# res = '/home/jianping/pywork/drug/media/autoduck_db/1'
# for n in file_lst:
#     path = '/home/jianping/pywork/drug/media/autoduck_db/1mol2/%s' % n
#     os.system("python /home/jianping/pywork/drug/extra_apps/vina/prepare_ligand4.py -l %s -v"
#                   % path)
#     os.system("mv %s %s" % (path.split('/')[-1].split('.')[0]+'.pdbqt', res))

