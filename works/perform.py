# coding=utf-8
from __future__ import unicode_literals
import os
import cPickle
import re
import math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from drug.settings import BASE_DIR, TARGET_FOLDER_BASE, drugdb
username = 'wibrow'
work_name = 'iiooi'
center_x = 5
center_y = 6
center_z = 6
size_x = 6
size_y = 9
size_z = 10
pdb_file = '1hsg.pdb'
user_db_name = '3.sdf'
mol_db = '2'
lig_file = 'CHEMBL16326.mol2'
pdb_path = os.path.join(BASE_DIR, 'media', 'screen', username, work_name, pdb_file)  # (username, work_name)
lig_path = os.path.join(BASE_DIR, 'media', 'screen', username, work_name, lig_file)
res_path = os.path.join(BASE_DIR, 'media', 'screen', username, work_name)

FP_PARAM = {
    'topological_hashed': {
        "mean": 0.000026936641120031898,
        "sd": 0.00398734520007183,
        "sd_exp": 0.5145149536573593,
        "tc": 0.66,
        "fp_func": lambda m: AllChem.GetHashedTopologicalTorsionFingerprint(m)
    },
    'atompair_hashed': {
        "mean": 0.00002541280101763038,
        "sd": 0.003864080986890242,
        "sd_exp": 0.5118760617288712,
        "tc": 0.61,
        "fp_func": lambda m: AllChem.GetHashedAtomPairFingerprint(m)
    },
    'maccs': {
        "mean": 0.00003686608558399457,
        "sd": 0.00484043909504593,
        "sd_exp": 0.5199608769322853,
        "tc": 0.88,
        "fp_func": lambda m: AllChem.GetMACCSKeysFingerprint(m)
    },
    'morgan_hashed': {
        "mean": 0.00002318732496154415,
        "sd": 0.003829515743309366,
        "sd_exp": 0.5101931872278405,
        "tc": 0.65,
        "fp_func": lambda m: AllChem.GetHashedMorganFingerprint(m, 2)
    }
}


def raw_score(target_mol_pkl, mol_fp, cutoff):
    # try:
    sim_list = list()
    for el_fp in target_mol_pkl:
        sim = TanimotoSimilarity(mol_fp, el_fp)
        if sim >= cutoff:
            sim_list.append(sim)
    rawscore = sum(sim_list)
    if rawscore > 0:
        return round(rawscore, 3)
    return None


def z_score(rs, size, mean, sd, sd_exp):
    return (rs - size * mean) / (sd * size ** sd_exp)


def p_value(z):
    x = -math.exp(-z * math.pi / math.sqrt(6) - 0.577215665)
    if z > 28:
        return -(x + x ** 2 / 2 + x ** 3 / 6)
    else:
        return 1 - math.exp(x)


def pred(smiles, target_list):
    result = list()
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        # todo: invalid rdkit molecule
        print 'invalid smiles'
        return None
    # for fp_name, fp_parm in FP_PARAM.iteritems():
    #     mol_fp = fp_parm['fp_func'](mol)
    #     print fp_name
    #     fp_result = dict()
    fp_dict = dict()
    for fp_name, fp_param in FP_PARAM.iteritems():
       fp_dict[fp_name] = fp_param['fp_func'](mol)

    for idx, chembl_id in enumerate(target_list):
        # print idx
        target_result = dict()
        for fp_name, fp_param in FP_PARAM.iteritems():
            mol_fp = fp_dict[fp_name]
            target_mol_pkl = cPickle.load(open(os.path.join(TARGET_FOLDER_BASE, fp_name, chembl_id), 'r'))
            rs = raw_score(target_mol_pkl, mol_fp, fp_param['tc'])
            if rs:
                zscore = z_score(rs, len(target_mol_pkl), fp_param['mean'], fp_param['sd'], fp_param['sd_exp'])
                pvalue = p_value(zscore)
                target_result[fp_name] = pvalue

        if target_result:
            target_result['chembl_id'] = chembl_id
            result.append(target_result)
    return result


def perform_dock(work_name, center_x, center_y, center_z, size_x, size_y, size_z, pdb_path, lig_path, res_path):
    # dock_status(work_name=work_name, status='正在处理')
    res = '%s/res' % res_path
    if not os.path.exists(res):
        os.mkdir(res)
    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r %s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    os.system("mv %s %s" % (pdbqt, res_path))
    os.system("python %s/extra_apps/vina/prepare_ligand4.py -l %s -v" % (BASE_DIR, lig_path))
    ligqt = lig_path.split('/')[-1].split('.')[0] + '.pdbqt'
    os.system("mv %s %s" % (ligqt, res_path))
    os.system("%s --receptor %s/%s --ligand %s/%s --center_x %s --center_y %s --center_z %s --size_x %s --size_y %s"
              " --size_z %s" % (vina_path, res_path, pdbqt, res_path, ligqt, center_x, center_y, center_z,
                                size_x, size_y, size_z))
    outqt = ligqt.split('.')[0]+'_out.pdbqt'
    outqt_path = os.path.join(res_path, outqt)
    if os.path.exists(outqt_path):
        reg = 'REMARK VINA RESULT:(.*?)\n'
        re_reg = re.compile(reg)
        with open(outqt_path, 'r') as f:
            data = f.read()
        out_lst = re_reg.findall(data)
        if out_lst:
            med = []
            for model in out_lst:
                med.append(float(model.split()[0]))
            affinity = min(med)
        os.system("mv %s %s" % (outqt_path, res))
        os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res, outqt)))
        out_path = os.path.join(res.split('media/')[1], outqt[:-2])


def gen_user_db_qt_smiles(input_file, output_path):
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    sdf = Chem.SDMolSupplier(input_file.encode('utf-8'))
    sdf = [n for n in sdf if n is not None]
    print len(sdf)
    smiles = []
    for mol in sdf:
        name = mol.GetProp('_Name'.encode('utf-8'))
        # print type(name)
        name = name + '.pdb'.encode('utf-8')
        print type(name)
        mol_path = os.path.join(output_path, name)

        Chem.MolToPDBFile(mol=mol, filename=mol_path.encode('utf-8'))
        smile = Chem.MolToSmiles(mol)
        smiles.append([name, smile])
    df = pd.DataFrame(smiles)
    smiles_file = os.path.join(output_path, 'smiles.csv')
    df.to_csv(smiles_file, index=False, header=False, encoding='utf-8')

    pdblst = os.listdir(output_path)
    pdblst = [pdb for pdb in pdblst if pdb.endswith('.pdb')]
    for pdb in pdblst:
        pdb_path = os.path.join(output_path, pdb)
        os.system("python %s/extra_apps/vina/prepare_ligand4.py -l %s -v" % (BASE_DIR, pdb_path))
        pdbqt = pdb.split('.')[0] + '.pdbqt'
        if os.path.exists(pdbqt):
            os.system("mv %s %s" % (pdbqt, output_path))


def perform_screen_user(work_name, center_x, center_y, center_z, size_x, size_y, size_z, user_db_name, pdb_path, res_path):
    """
    用户提供数据库以及中心坐标和盒子大小进行筛选
    :param work_name:
    :param center_x:
    :param center_y:
    :param center_z:
    :param size_x:
    :param size_y:
    :param size_z:
    :param user_db_name:
    :param pdb_path:
    :param res_path:
    :return:
    """

    screen_out = 'screen_out.csv'
    user_db = os.path.join(res_path, 'userdb')
    res = os.path.join(res_path, 'res')

    if not os.path.exists(user_db):
        os.mkdir(user_db)

    if not os.path.exists(res):
        os.mkdir(res)

    input_file = os.path.join(res_path, user_db_name)

    gen_user_db_qt_smiles(input_file, user_db)

    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r /%s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    os.system("mv %s %s" % (pdbqt, res_path))

    target_list = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))
    smile_file = os.path.join(user_db, 'smiles.csv')
    df = pd.read_csv(smile_file, header=None, encoding='utf-8')
    smile_data = df.values.tolist()

    for ligand in smile_data:
        smiles = ligand[1]
        sea_res = pred(smiles, target_list)
        if sea_res:
            ligand = ligand[0].split('.')[0]+'.pdbqt'
            ligand_path = os.path.join(user_db, ligand)
            if os.path.exists(ligand_path):
                os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                          " --center_z %s --size_x %s --size_y %s --size_z %s" % (vina_path, res_path, pdbqt,
                                                                                  ligand_path, center_x, center_y, center_z,
                                                                                  size_x, size_y, size_z))
                os.system("mv %s %s" % (ligand_path.split('.')[0]+'_out.pdbqt', res))

    res_lst = os.listdir(res)
    reg = 'REMARK VINA RESULT:(.*?)\n'
    re_reg = re.compile(reg)
    screen_res = []
    for out in res_lst[:]:
        out_path = os.path.join(res, out)
        with open(out_path, 'r') as f:
            data = f.read()
        out_lst = re_reg.findall(data)
        if out_lst:
            med = []
            for model in out_lst:
                med.append(float(model.split()[0]))
            file_name = out.split('_out')[0]
            screen_res.append([file_name, min(med)])
    arr = np.array(screen_res)
    df = pd.DataFrame(arr, columns=['id', 'Affinity (kcal/mol)'])
    df = df.sort_values("Affinity (kcal/mol)")
    df.to_csv(screen_out, index=False)
    os.system("mv %s %s" % (screen_out, res))

# smi = os.path.join(BASE_DIR, 'media', 'sdf', 'smiles.csv')
# sdf = pd.read_csv(smi, header=None, encoding='utf-8')
# print sdf.values.tolist()


#
#
# perform_screen_user(work_name, center_x, center_y, center_z, size_x, size_y, size_z, user_db_name, pdb_path, res_path)

# os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f compound_sdf-18_2_53_39_out.pdbqt -v" % BASE_DIR)

# perform_dock(work_name, center_x, center_y, center_z, size_x, size_y, size_z, pdb_path, lig_path, res_path)

def perform_screen(work_name, center_x, center_y, center_z, size_x, size_y, size_z, mol_db, pdb_path, res_path):
    """
    用户指定筛选中心以及盒子大小
    :param work_name:
    :param center_x:
    :param center_y:
    :param center_z:
    :param size_x:
    :param size_y:
    :param size_z:
    :param mol_db:
    :param pdb_path:
    :param res_path:
    :return:
    """
    # screen_status(work_name, status='正在处理')

    screen_out = 'screen_out.csv'
    res = os.path.join(res_path, 'res')

    if not os.path.exists(res):
        os.mkdir(res)

    ligand_db = os.path.join(drugdb, mol_db)

    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r %s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    os.system("mv %s %s" % (pdbqt, res_path))

    target_list = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))
    smile_file = os.path.join(ligand_db, 'smiles.csv')
    df = pd.read_csv(smile_file, header=None, encoding='utf-8')
    smile_data = df.values.tolist()

    for ligand in smile_data:
        smiles = ligand[1]
        sea_res = 1
        sea_res = pred(smiles, target_list)
        if sea_res:
            ligand = ligand[0].split('.')[0] + '.pdbqt'
            ligand_path = ligand_db + '/' + ligand
            os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                      " --center_z %s --size_x %s --size_y %s --size_z %s" % (vina_path, res_path, pdbqt,
                                                                              ligand_path, center_x, center_y, center_z,
                                                                              size_x, size_y, size_z))
            os.system("mv %s %s" % (ligand_path.split('.')[0]+'_out.pdbqt', res))
            os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                    ligand_path.split('/')[-1].split('.')[0]+'_out.pdbqt')))

    res_lst = os.listdir(res)
    reg = 'REMARK VINA RESULT:(.*?)\n'
    re_reg = re.compile(reg)
    screen_res = []
    insert_lst = []
    for out in res_lst[:]:
        out_path = os.path.join(res, out)
        with open(out_path, 'r') as f:
            data = f.read()
        out_lst = re_reg.findall(data)
        if out_lst:
            med = []
            for model in out_lst:
                med.append(float(model.split()[0]))
            file_name = out.split('_out')[0]
            screen_res.append([file_name, min(med)])
            # insert_lst.append(Screen(work_name=work_name, affinity=min(med), path=os.path.join(res.split('media/')[1],
                                                                                               # out[:-2])))
    # Screen.objects.bulk_create(insert_lst)
    arr = np.array(screen_res)
    df = pd.DataFrame(arr, columns=['id', 'Affinity (kcal/mol)'])
    df = df.sort_values("Affinity (kcal/mol)", ascending=False)
    df.to_csv(screen_out, index=False)
    os.system("mv %s %s" % (screen_out, res))


# perform_screen(work_name, center_x, center_y, center_z, size_x, size_y, size_z, mol_db, pdb_path, res_path)
