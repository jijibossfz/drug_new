# coding=utf-8
from __future__ import division
import re
import os
import math
import cPickle
import multiprocessing as mp
import shutil

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from celery import shared_task
from works.models import Screen, Dock, VirScreen, SeaTarget, SeaVirScreen, Gbsa, UserProfile, AutoDock, AutoDock2, \
    VirtualScreen, VirtualScreen2
from drug.settings import BASE_DIR, TARGET_FOLDER_BASE, drugdb
from extra_apps.utils.email_send import email_status, email_status_sea

# 获取SEA算法预测的靶点列表
TargetList = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))
# from perform import gen_user_db_qt_smiles


# SEA建模后四种指纹分别计算出来的参数
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


# SEA打分函数
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


# 计算z_score公式
def z_score(rs, size, mean, sd, sd_exp):
    return (rs - size * mean) / (sd * size ** sd_exp)


# 计算 p_value公式
def p_value(z):
    x = -math.exp(-z * math.pi / math.sqrt(6) - 0.577215665)
    if z > 28:
        return -(x + x ** 2 / 2 + x ** 3 / 6)
    else:
        return 1 - math.exp(x)


# SEA算法预测靶点函数
def pred(smiles, target_list):
    """
    :param smiles: 输入需要预测分子的Smiles
    :param target_list: 输入SEA算法中能够预测的全部的靶点列表
    :return: 返回预测分子的靶点,字典形式返回
    """
    result = list()
    try:
        # 读取输入的分子Smiles,转化为rdkit内部可以识别别的分子
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
    # 获取所有建模的指纹以及计算出来的值
    for fp_name, fp_param in FP_PARAM.iteritems():
       fp_dict[fp_name] = fp_param['fp_func'](mol)
    # 通过四种指纹对输入的分子进行靶点预测
    # 遍历所有可能的靶点预测
    for idx, chembl_id in enumerate(target_list):
        # print idx
        target_result = dict()
        # 遍历四种指纹
        for fp_name, fp_param in FP_PARAM.iteritems():
            if fp_name == 'maccs':
                mol_fp = fp_dict[fp_name]
                target_mol_pkl = cPickle.load(open(os.path.join(TARGET_FOLDER_BASE, fp_name, chembl_id), 'r'))
                # 计算rs值是否是大于指纹设定的阙值
                rs = raw_score(target_mol_pkl, mol_fp, fp_param['tc'])
                if rs:
                    # 如果rs值大于指纹设定的阙值,则对某种靶点具有活性计算z_score pvalue
                    zscore = z_score(rs, len(target_mol_pkl), fp_param['mean'], fp_param['sd'], fp_param['sd_exp'])
                    pvalue = p_value(zscore)
                    target_result[fp_name] = pvalue

        if target_result:
            target_result['chembl_id'] = chembl_id
            result.append(target_result)
            # 返回预测的靶点result 字典形式返回
    return result


def dock_status(work_name, status):
    """
    分子对接任务开始时修改任务执行状态
    :param work_name: 输入任务id 号确认修改任务状态
    :param status: 输入需要修改成为的状态
    :return: None 只修改状态无返回值
    """
    work = Dock.objects.get(work_name=work_name)
    work.status = status
    work.save()


def dock_status2(work_name, status):
    """
    分子对接任务开始时修改任务执行状态
    :param work_name: 输入任务id 号确认修改任务状态
    :param status: 输入需要修改成为的状态
    :return: None 只修改状态无返回值
    """
    work = AutoDock2.objects.get(work_name=work_name)
    work.status = status
    work.save()


def dock_out(work_name, out_path):
    """
    保存分子对接任务结果
    :param work_name: 任务id
    :param out_path: 输出文件路径
    :return:
    """
    work = Dock.objects.get(work_name=work_name)
    work.out_path = out_path
    work.save()


def dock_out2(work_name, out_path):
    """
    保存分子对接任务结果
    :param work_name: 任务id
    :param out_path: 输出文件路径
    :return:
    """
    work = AutoDock2.objects.get(work_name=work_name)
    work.out_path = out_path
    work.save()


def dock_affinity(work_name, affinity):
    """
    保存分子对接任务的打分值
    :param work_name: 输入任务id 号确认任务
    :param affinity: 输入需要保存的分子对接打分值
    :return: None 只保存分子对接任务的打分值,无返回
    """
    work = Dock.objects.get(work_name=work_name)
    work.affinity = affinity
    work.save()


def dock_affinity2(work_name, affinity):
    """
    保存分子对接任务的打分值
    :param work_name: 输入任务id 号确认任务
    :param affinity: 输入需要保存的分子对接打分值
    :return: None 只保存分子对接任务的打分值,无返回
    """
    work = AutoDock2.objects.get(work_name=work_name)
    work.affinity = affinity
    work.save()


def screen_status(work_name, status):
    work = VirScreen.objects.get(work_name=work_name)
    work.status = status
    work.save()


def screen_status2(work_name, status):
    work = VirtualScreen2.objects.get(work_name=work_name)
    work.status = status
    work.save()
# def VirScreen_status(work_name, status):
#     work = VirScreen.objects.get(work_name=work_name)
#     work.status = status
#     work.save()


def gbsa_status(work_name, status):
    work = Gbsa.objects.get(work_name=work_name)
    work.status = status
    work.save()


def gen_user_db_qt_smiles(input_file, output_path):
    """
    处理用户上传的数据库
    :param input_file:
    :param output_path:
    :return:
    """
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    sdf = Chem.SDMolSupplier(input_file.encode('utf-8'))
    sdf = [n for n in sdf if n is not None]
    smiles = []
    for mol in sdf:
        name = mol.GetProp('_Name'.encode('utf-8'))
        name = name + '.pdb'.encode('utf-8')
        mol_path = os.path.join(output_path, name)
        try:
            Chem.MolToPDBFile(mol=mol, filename=mol_path.encode('utf-8'))
            smile = Chem.MolToSmiles(mol)
            smiles.append([name, smile])
        except IOError:
            pass
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


@shared_task
def perform_screen(work_name, user_target, center_x, center_y, center_z, size_x, size_y, size_z, mol_db, pdb_path, res_path):
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
    screen_status(work_name, status='computing')
    email = UserProfile.objects.filter(id=VirtualScreen.objects.filter(pdb_file=pdb_path.replace('/home/wz/pywork/try27/drug/media', '')).values()[0]['user_id']).values()[0]['email']
    screen_out = 'screen_out.csv'
    res = os.path.join(res_path, 'res')

    if not os.path.exists(res):
        os.mkdir(res)

    ligand_db = os.path.join(drugdb, mol_db)

    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r %s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(pdbqt):
        os.system("mv %s %s" % (pdbqt, res_path))
    if len(user_target) > 1:
        user_target = user_target.split(';')
        target_list = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))
        smile_file = os.path.join(ligand_db, 'smiles.csv')
        df = pd.read_csv(smile_file, header=None, encoding='utf-8')
        smile_data = df.values.tolist()
        curr_proc = mp.current_process()
        curr_proc.daemon = False
        p = mp.Pool(processes=mp.cpu_count())
        curr_proc.daemon = True
        pool_lst = []
        for ligand in smile_data:
            smiles = ligand[1]
            targets = p.apply_async(pred, args=(smiles, target_list))
            pool_lst.append([ligand[0], targets])
        p.close()
        p.join()
        pool_lst = [[n[0], n[1].get()] for n in pool_lst]
        for target in pool_lst:
            if target[1]:
                pred_target = []
                for pred_tar in target[1]:
                    pred_target.append(pred_tar['chembl_id'])
                same_target = [l for l in pred_target if l in user_target]
                if same_target:
                    ligand = target[0].split('.')[0] + '.pdbqt'
                    ligand_path = os.path.join(ligand_db, ligand)
                    os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                              " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                                      ligand_path, center_x, center_y,
                                                                                      center_z, size_x, size_y, size_z))
                    os.system("mv %s %s" % (ligand_path.split('.')[0]+'_out.pdbqt', res))
                    os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                            ligand_path.split('/')[-1].split('.')[0]+'_out.pdbqt')))
    else:
        ligand_lst = os.listdir(ligand_db)
        ligand_lst = [v for v in ligand_lst if v.endswith('pdbqt')]
        for ligand_ in ligand_lst:
            ligand_path = os.path.join(ligand_db, ligand_)
            os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                      " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                              ligand_path, center_x, center_y,
                                                                              center_z, size_x, size_y, size_z))
            os.system("mv %s %s" % (ligand_path.split('.')[0] + '_out.pdbqt', res))
            os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                                                                     ligand_path.split(
                                                                                                         '/')[-1].split(
                                                                                                         '.')[
                                                                                                         0] + '_out.pdbqt')))
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
            insert_lst.append(Screen(work_name=work_name, screen_cat='screen', affinity=min(med),
                                     path=os.path.join(res.split('media/')[1], out[:-2])))
    if screen_res:
        Screen.objects.bulk_create(insert_lst)
        arr = np.array(screen_res)
        df = pd.DataFrame(arr, columns=['id', 'Affinity (kcal/mol)'])
        df = df.sort_values("Affinity (kcal/mol)", ascending=False)
        df.to_csv(screen_out, index=False)
        os.system("mv %s %s" % (screen_out, res))
        email_status(email, res_path + '/res/screen_out.csv')
    else:
        email_status(email, '')
    screen_status(work_name=work_name, status='completed')


@shared_task
def perform_screen2(work_name, user_target, mol_db, pdb_path, resi_path, res_path):
    """
    用户指定几个残基进行筛选
    :param work_name:
    :param mol_db:
    :param pdb_path:
    :param lig_path:
    :param res_path:
    :return:
    """
    screen_status2(work_name, status='computing')
    email = UserProfile.objects.filter(id=VirtualScreen2.objects.filter(pdb_file=pdb_path.replace('/home/wz/pywork/try27/drug/media', '')).values()[0]['user_id']).values()[0]['email']
    with open(resi_path, 'r') as f:
        lines = f.readlines()
    lines = [n.rstrip() for n in lines if len(n) > 1]

    x, y, z = [], [], []
    for n in lines:
        x.append(float(n[30:38]))
        y.append(float(n[38:46]))
        z.append(float(n[46:54]))
    center_x = float('%.3f' % (sum(x)/len(x)))
    center_y = float('%.3f' % (sum(y)/len(y)))
    center_z = float('%.3f' % (sum(z)/len(z)))
    size_x = max(x) - min(x)
    size_y = max(y) - min(y)
    size_z = max(z) - min(z)

    screen_out = 'screen_out.csv'
    res = os.path.join(res_path, 'res')

    if not os.path.exists(res):
        os.mkdir(res)

    ligand_db = os.path.join(drugdb, mol_db)

    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r /%s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(pdbqt):
        os.system("mv %s %s" % (pdbqt, res_path))

    if len(user_target) > 1:
        user_target = user_target.split(';')
        target_list = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))
        smile_file = os.path.join(ligand_db, 'smiles.csv')
        df = pd.read_csv(smile_file, header=None, encoding='utf-8')
        smile_data = df.values.tolist()
        curr_proc = mp.current_process()
        curr_proc.daemon = False
        p = mp.Pool(processes=mp.cpu_count())
        curr_proc.daemon = True
        pool_lst = []
        for ligand in smile_data:
            smiles = ligand[1]
            targets = p.apply_async(pred, args=(smiles, target_list))
            pool_lst.append([ligand[0], targets])
        p.close()
        p.join()
        pool_lst = [[n[0], n[1].get()] for n in pool_lst]
        for target in pool_lst:
            if target[1]:
                pred_target = []
                for pred_tar in target[1]:
                    pred_target.append(pred_tar['chembl_id'])
                same_target = [l for l in pred_target if l in user_target]
                if same_target:
                    ligand = target[0].split('.')[0] + '.pdbqt'
                    ligand_path = ligand_db + '/' + ligand
                    os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                              " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                                      ligand_path, center_x, center_y, center_z,
                                                                                      size_x, size_y, size_z))
                    os.system("mv %s %s" % (ligand_path.split('.')[0]+'_out.pdbqt', res))
                    os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                            ligand_path.split('/')[-1].split('.')[0]+'_out.pdbqt')))
    else:
        ligand_lst = os.listdir(ligand_db)
        ligand_lst = [v for v in ligand_lst if v.endswith('pdbqt')]
        for ligand_ in ligand_lst:
            ligand_path = os.path.join(ligand_db, ligand_)
            os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                      " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                              ligand_path, center_x, center_y,
                                                                              center_z, size_x, size_y, size_z))
            os.system("mv %s %s" % (ligand_path.split('.')[0] + '_out.pdbqt', res))
            os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                                                                     ligand_path.split(
                                                                                                         '/')[-1].split(
                                                                                                         '.')[0] + '_out.pdbqt')))
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
            insert_lst.append(Screen(work_name=work_name, screen_cat='screen2', affinity=min(med),
                                     path=os.path.join(res.split('media/')[1], out[:-2])))
    if screen_res:
        Screen.objects.bulk_create(insert_lst)
        arr = np.array(screen_res)
        df = pd.DataFrame(arr, columns=['id', 'Affinity (kcal/mol)'])
        df = df.sort_values("Affinity (kcal/mol)", ascending=False)
        df.to_csv(screen_out, index=False)
        os.system("mv %s %s" % (screen_out, res))
        email_status(email, res_path + '/res/screen_out.csv')
    else:
        email_status(email, '')
    screen_status2(work_name=work_name, status='completed')


@shared_task
def perform_screen_user(work_name, user_target, center_x, center_y, center_z, size_x, size_y, size_z, user_db_name, pdb_path, res_path):
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
    screen_status(work_name, status='computing')
    email = UserProfile.objects.filter(id=VirtualScreen.objects.filter(pdb_file=pdb_path.replace('/home/wz/pywork/try27/drug/media/', '')).values()[0]['user_id']).values()[0]['email']
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
    if os.path.exists(pdbqt):
        os.system("mv %s %s" % (pdbqt, res_path))

    if len(user_target) > 1:
        user_target = user_target.split(';')
        target_list = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))
        smile_file = os.path.join(user_db, 'smiles.csv')
        df = pd.read_csv(smile_file, header=None, encoding='utf-8')
        smile_data = df.values.tolist()
        curr_proc = mp.current_process()
        curr_proc.daemon = False
        p = mp.Pool(processes=mp.cpu_count())
        curr_proc.daemon = True
        pool_lst = []
        for ligand in smile_data:
            smiles = ligand[1]
            targets = p.apply_async(pred, args=(smiles, target_list))
            pool_lst.append([ligand[0], targets])
        p.close()
        p.join()
        pool_lst = [[n[0], n[1].get()] for n in pool_lst]
        for target in pool_lst:
            if target[1]:
                pred_target = []
                for pred_tar in target[1]:
                    pred_target.append(pred_tar['chembl_id'])
                same_target = [l for l in pred_target if l in user_target]
                if same_target:
                    ligand = target[0].split('.')[0] + '.pdbqt'
                    ligand_path = os.path.join(user_db, ligand)
                    if os.path.exists(ligand_path):
                        os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                                  " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                                          ligand_path, center_x, center_y,
                                                                                          center_z, size_x, size_y, size_z))
                        os.system("mv %s %s" % (ligand_path.split('.')[0] + '_out.pdbqt', res))
                        os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                                ligand_path.split('/')[-1].split('.')[0] + '_out.pdbqt')))
    else:
        ligand_lst = os.listdir(user_db)
        ligand_lst = [v for v in ligand_lst if v.endswith('pdbqt')]
        for ligand_ in ligand_lst:
            ligand_path = os.path.join(user_db, ligand_)
            os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                      " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                              ligand_path, center_x, center_y,
                                                                              center_z, size_x, size_y, size_z))
            os.system("mv %s %s" % (ligand_path.split('.')[0] + '_out.pdbqt', res))
            os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                                                                     ligand_path.split(
                                                                                                         '/')[-1].split(
                                                                                                         '.')[0] + '_out.pdbqt')))
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
            insert_lst.append(Screen(work_name=work_name, screen_cat='screen', affinity=min(med),
                                     path=os.path.join(res.split('media/')[1],
                                                       out[:-2])))
    if screen_res:
        Screen.objects.bulk_create(insert_lst)
        arr = np.array(screen_res)
        df = pd.DataFrame(arr, columns=['id', 'Affinity (kcal/mol)'])
        df = df.sort_values("Affinity (kcal/mol)", ascending=False)
        df.to_csv(screen_out, index=False)
        os.system("mv %s %s" % (screen_out, res))
        email_status(email, res_path + '/res/screen_out.csv')
    else:
        email_status(email,'')
    screen_status(work_name, status='completed')
    # shutil.rmtree(res_path)


@shared_task
def perform_screen2_user(work_name, user_target, user_db_name, pdb_path, resi_path, res_path):
    """
    用户提供数据库以及残基进行筛选
    :param work_name:
    :param user_db_name:
    :param pdb_path:
    :param resi_path:
    :param res_path:
    :return:
    """
    screen_status2(work_name, status='computing')
    email = UserProfile.objects.filter(id=VirtualScreen2.objects.filter(pdb_file=pdb_path.replace('/home/wz/pywork/try27/drug/media/', '')).values()[0]['user_id']).values()[0]['email']

    with open(resi_path, 'r') as f:
        lines = f.readlines()
    lines = [n.rstrip() for n in lines if len(n) > 1]

    x, y, z = [], [], []
    for n in lines:
        x.append(float(n[30:38]))
        y.append(float(n[38:46]))
        z.append(float(n[46:54]))
    center_x = float('%.3f' % (sum(x)/len(x)))
    center_y = float('%.3f' % (sum(y)/len(y)))
    center_z = float('%.3f' % (sum(z)/len(z)))
    size_x = max(x) - min(x)
    size_y = max(y) - min(y)
    size_z = max(z) - min(z)

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
    if os.path.exists(pdbqt):
        os.system("mv %s %s" % (pdbqt, res_path))
    if len(user_target) > 1:
        user_target = user_target.split(';')
        target_list = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))
        smile_file = os.path.join(user_db, 'smiles.csv')
        df = pd.read_csv(smile_file, header=None, encoding='utf-8')
        smile_data = df.values.tolist()
        curr_proc = mp.current_process()
        curr_proc.daemon = False
        p = mp.Pool(processes=mp.cpu_count())
        curr_proc.daemon = True
        pool_lst = []
        for ligand in smile_data:
            smiles = ligand[1]
            targets = p.apply_async(pred, args=(smiles, target_list))
            pool_lst.append([ligand[0], targets])
        p.close()
        p.join()
        pool_lst = [[n[0], n[1].get()] for n in pool_lst]
        for target in pool_lst:
            if target[1]:
                pred_target = []
                for pred_tar in target[1]:
                    pred_target.append(pred_tar['chembl_id'])
                same_target = [l for l in pred_target if l in user_target]
                if same_target:
                    ligand = target[0].split('.')[0] + '.pdbqt'
                    ligand_path = os.path.join(user_db, ligand)
                    if os.path.exists(ligand_path):
                        os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                                  " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                                          ligand_path, center_x, center_y,
                                                                                          center_z, size_x, size_y, size_z))
                        os.system("mv %s %s" % (ligand_path.split('.')[0] + '_out.pdbqt', res))
                        os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                            ligand_path.split('/')[-1].split('.')[0]+'_out.pdbqt')))

    else:
        ligand_lst = os.listdir(user_db)
        ligand_lst = [v for v in ligand_lst if v.endswith('pdbqt')]
        for ligand_ in ligand_lst:
            ligand_path = os.path.join(user_db, ligand_)
            os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                      " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                              ligand_path, center_x, center_y,
                                                                              center_z, size_x, size_y, size_z))
            os.system("mv %s %s" % (ligand_path.split('.')[0] + '_out.pdbqt', res))
            os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                                                                     ligand_path.split(
                                                                                                         '/')[-1].split(
                                                                                                         '.')[0] + '_out.pdbqt')))
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
            insert_lst.append(Screen(work_name=work_name, screen_cat='screen2', affinity=min(med),
                                     path=os.path.join(res.split('media/')[1],
                                                       out[:-2])))
    if screen_res:
        Screen.objects.bulk_create(insert_lst)
        arr = np.array(screen_res)
        df = pd.DataFrame(arr, columns=['id', 'Affinity (kcal/mol)'])
        df = df.sort_values("Affinity (kcal/mol)", ascending=False)
        df.to_csv(screen_out, index=False)
        os.system("mv %s %s" % (screen_out, res))
        email_status(email, res_path + '/res/screen_out.csv')
    else:
        email_status(email, '')
    screen_status2(work_name, status='completed')
    # shutil.rmtree(res_path)


# 定义用户指定对接中心分子对接任务
@shared_task
def perform_dock(work_name, center_x, center_y, center_z, size_x, size_y, size_z, pdb_path, lig_path, res_path, email):
    """
    用户指定对接中心以及盒子大小
    :param work_name: 任务id
    :param center_x:  用户制定的对接中心以及盒子大小六个参数
    :param center_y:
    :param center_z:
    :param size_x:
    :param size_y:
    :param size_z:
    :param pdb_path: 用户上传的受体文件
    :param lig_path: 用户上传的配体文件
    :param res_path: 用户任务保存路径
    :param email: 用户的邮箱
    :return: None 只执行任务无返回
    """
    # 任务开始 修改任务状态为computing
    dock_status(work_name=work_name, status='computing')
    # email = UserProfile.objects.filter(id=AutoDock.objects.filter(pdb_file=pdb_path.replace('/home/wz/pywork/try27/drug/media/', '')).values()[0]['user_id']).values()[0]['email']
    # 创建任务保存路径
    res = os.path.join(res_path, 'res')
    if not os.path.exists(res):
        os.mkdir(res)
    # 处理受体文件
    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r %s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(pdbqt):
        os.system("mv %s %s" % (pdbqt, res_path))
    # 处理配体文件
    os.system("python %s/extra_apps/vina/prepare_ligand4.py -l %s -v" % (BASE_DIR, lig_path))
    ligqt = lig_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(ligqt):
        os.system("mv %s %s" % (ligqt, res_path))
    # 执行对接
    os.system("%s --receptor %s/%s --ligand %s/%s --center_x %s --center_y %s --center_z %s --size_x %s --size_y %s"
              " --size_z %s" % ('vina', res_path, pdbqt, res_path, ligqt, center_x, center_y, center_z,
                                size_x, size_y, size_z))
    outqt = ligqt.split('.')[0]+'_out.pdbqt'
    outqt_path = os.path.join(res_path, outqt)
    # 取出分子对接打分值
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
            affinity = str(min(med))
            dock_affinity(work_name=work_name, affinity=affinity)
        os.system("mv %s %s" % (outqt_path, res))
        os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res, outqt)))
        out_path = os.path.join(res.split('media/')[1], outqt[:-2])
        dock_out(work_name=work_name, out_path=out_path)
        email_status(email, '/home/wz/pywork/try27/drug/media/' + out_path.encode('raw_unicode_escape'))
    # 对接任务完成,把任务状态修改为completed
    dock_status(work_name=work_name, status='completed')
    shutil.rmtree(res_path)
    # 对接任务完成后发送email


# 定义用户指定关键残基分子对接任务
@shared_task
def perform_dock2(work_name, pdb_path, lig_path, resn, res_path, email):
    """
    :param work_name: 任务id
    :param pdb_path: 用户上传受体文件
    :param lig_path: 用户上传配体文件
    :param resn: 用户指定的关键残基
    :param res_path: 任务保存路径
    :param email: 用户的邮箱
    :return: None 只执行任务,无返回
    """
    # 任务开始 修改任务状态为computing
    # email = UserProfile.objects.filter(id=Dock.objects.filter(pdb_file=pdb_path.replace('/home/wz/pywork/try27/drug/media', '')).values()[0]['user_id']).values()[0]['email']
    dock_status(work_name=work_name, status='computing')
    res = '%s/res' % res_path
    # 创建任务保存路径
    if not os.path.exists(res):
        os.mkdir(res)
    # 根据用户指定的关键残基计算对接中心和盒子大小
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    lines = [n.rstrip() for n in lines if n.startswith("ATOM")]
    resn_lst = resn.split(";")
    x, y, z = [], [], []
    for n in lines:
        if n.split()[5] in resn_lst:
            x.append(float(n.split()[6]))
            y.append(float(n.split()[7]))
            z.append(float(n.split()[8]))
    center_x = float('%.3f' % (sum(x)/len(x)))
    center_y = float('%.3f' % (sum(y)/len(y)))
    center_z = float('%.3f' % (sum(z)/len(z)))
    size_x = max(x) - min(x)
    size_y = max(y) - min(y)
    size_z = max(z) - min(z)
    # 处理受体文件
    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r %s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(pdbqt):
        os.system("mv %s %s" % (pdbqt, res_path))
    # 处理配体文件
    os.system("python %s/extra_apps/vina/prepare_ligand4.py -l %s -v" % (BASE_DIR, lig_path))
    ligqt = lig_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(ligqt):
        os.system("mv %s %s" % (ligqt, res_path))
    # 执行分子对接任务
    os.system("%s --receptor %s/%s --ligand %s/%s --center_x %s --center_y %s --center_z %s --size_x %s --size_y %s"
              " --size_z %s" % ('vina', res_path, pdbqt, res_path, ligqt, center_x, center_y, center_z,
                                size_x, size_y, size_z))
    outqt = ligqt.split('.')[0]+'_out.pdbqt'
    outqt_path = os.path.join(res_path, outqt)
    # 取出分子对接打分值
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
            dock_affinity(work_name=work_name, affinity=affinity)
        os.system("mv %s %s" % (outqt_path, res))
        os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res, outqt)))
        out_path = os.path.join(res.split('media/')[1], outqt[:-2])
        dock_out(work_name=work_name, out_path=out_path)
        email_status(email, '/home/wz/pywork/try27/drug/media/' + out_path.encode('raw_unicode_escape'))
    # 分子对接任务完成,把任务状态修改为completed
    dock_status(work_name=work_name, status='completed')
    # 分子对接任务完成,给用户发送邮箱
    shutil.rmtree(res_path)


# 用户上传共晶结构下分子对接任务
@shared_task
def perform_dock3(work_name, pdb_path, lig_path, reference_path, res_path, email):
    """
    :param work_name: 任务id
    :param pdb_path: 用户上传的受体文件
    :param lig_path: 用户上传的配体文件
    :param reference_path: 用户上传的共晶结构文件
    :param res_path: 任务保存路径
    :param email: 用户的email
    :return:
    """
    # 任务开始 修改任务状态为computing
    dock_status(work_name=work_name, status='computing')
    # 创建任务保存路径
    res = '%s/res' % res_path
    if not os.path.exists(res):
        os.mkdir(res)
    # 根据用户上传的共晶结构计算对接中心和盒子大小
    with open(reference_path, 'r') as f:
        lines = f.readlines()
    lines = [n.rstrip() for n in lines if n.startswith("HETATM")]

    x, y, z = [], [], []
    for n in lines:
        if n.split()[3].startswith("MK"):
            x.append(float(n.split()[6]))
            y.append(float(n.split()[7]))
            z.append(float(n.split()[8]))
    center_x = float('%.3f' % (sum(x)/len(x)))
    center_y = float('%.3f' % (sum(y)/len(y)))
    center_z = float('%.3f' % (sum(z)/len(z)))
    size_x = max(x) - min(x)
    size_y = max(y) - min(y)
    size_z = max(z) - min(z)
    # 处理受体文件
    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r %s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(pdbqt):
        os.system("mv %s %s" % (pdbqt, res_path))
    # 处理配体文件
    os.system("python %s/extra_apps/vina/prepare_ligand4.py -l %s -v" % (BASE_DIR, lig_path))
    ligqt = lig_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(ligqt):
        os.system("mv %s %s" % (ligqt, res_path))
    # 执行分子对接任务
    os.system("%s --receptor %s/%s --ligand %s/%s --center_x %s --center_y %s --center_z %s --size_x %s --size_y %s"
              " --size_z %s" % ('vina', res_path, pdbqt, res_path, ligqt, center_x, center_y, center_z,
                                size_x, size_y, size_z))
    outqt = ligqt.split('.')[0]+'_out.pdbqt'
    outqt_path = os.path.join(res_path, outqt)
    # 取出分子对接任务的打分值
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
            dock_affinity(work_name=work_name, affinity=affinity)
        os.system("mv %s %s" % (outqt_path, res))
        os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res, outqt)))
        out_path = os.path.join(res.split('media/')[1], outqt[:-2])
        dock_out(work_name=work_name, out_path=out_path)
    # 分子对接任务完成,把任务状态修改为completed
    dock_status(work_name=work_name, status='completed')
    # 分子渡劫任务完成,给用户发送邮箱
    email_status(email, '/home/wz/pywork/try27/drug/media/' + out_path.encode('raw_unicode_escape'))


@shared_task
def perform_sea(work_name, smiles, email):
    """
    :param work_name:
    :param smiles:
    :return:
    """
    screen_status(work_name=work_name, status='computing')
    results = pred(smiles, target_list=TargetList)
    insert_lst = []
    if results:
        targets = []
        for n in results:
            targets.append(n['chembl_id'])
        tar_lst = ';'.join(targets)
        insert_lst.append(SeaTarget(work_name=work_name, smiles=smiles, target=tar_lst))
    else:
        insert_lst.append(SeaTarget(work_name=work_name, smiles=smiles, target='null'))
    SeaTarget.objects.bulk_create(insert_lst)
    screen_status(work_name=work_name, status='completed')
    email_status(email,'')


# 定义虚拟筛选任务函数,执行虚拟筛选任务
@shared_task
def perform_virscreen(work_name, user_target, mol_db, pdb_path, resn, res_path, email):
    """
    虚拟筛选任务
    :param work_name: 任务id
    :param user_target: 用户选择的靶点
    :param mol_db: 用户选择的数据库
    :param pdb_path: 用户上传的是受体文件
    :param resn: 用户选择的关键残疾确定活性位点
    :param res_path: 用户任务保存的路径
    :param email: 用户的email
    :return: None 只执行任务无返回
    """
    # 开始执行任务,把任务状态修改为computing
    screen_status(work_name=work_name, status='computing')
    # 创建保存结果的文件
    res = '%s/res' % res_path
    if not os.path.exists(res):
        os.mkdir(res)
    # 根据用户输入的关键残基确认活性位点中心和盒子大小
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    lines = [n.rstrip() for n in lines if n.startswith("ATOM")]
    resn_lst = resn.split(";")
    x, y, z = [], [], []
    for n in lines:
        if n.split()[5] in resn_lst:
            x.append(float(n.split()[6]))
            y.append(float(n.split()[7]))
            z.append(float(n.split()[8]))
    center_x = float('%.3f' % (sum(x)/len(x)))
    center_y = float('%.3f' % (sum(y)/len(y)))
    center_z = float('%.3f' % (sum(z)/len(z)))
    size_x = max(x) - min(x)
    size_y = max(y) - min(y)
    size_z = max(z) - min(z)

    screen_out = 'screen_out.csv'
    res = os.path.join(res_path, 'res')

    # 找到需要筛选的数据库
    ligand_db = os.path.join(drugdb, mol_db)

    # 处理受体文件生成晶格文件
    os.system("python %s/extra_apps/vina/prepare_receptor4.py -r /%s"
              " -A checkhydrogens " % (BASE_DIR, pdb_path))
    pdbqt = pdb_path.split('/')[-1].split('.')[0] + '.pdbqt'
    if os.path.exists(pdbqt):
        os.system("mv %s %s" % (pdbqt, res_path))
    # 获取所有靶点列表
    target_list = os.listdir(os.path.join(TARGET_FOLDER_BASE, 'maccs'))
    # 读取数据库中所有分子的smiles
    smile_file = os.path.join(ligand_db, 'smiles.csv')
    df = pd.read_csv(smile_file, header=None, encoding='utf-8')
    smile_data = df.values.tolist()
    # 使用多进程技术进行并行计算
    curr_proc = mp.current_process()
    curr_proc.daemon = False
    p = mp.Pool(processes=mp.cpu_count())
    curr_proc.daemon = True
    pool_lst = []
    for ligand in smile_data:
        smiles = ligand[1]
        # 把靶点预测任务添加到进程池中进行并行计算
        targets = p.apply_async(pred, args=(smiles, target_list))
        pool_lst.append([ligand[0], targets])
    p.close()
    p.join()
    # 并行计算完成之后取出并行计算任务的结果
    pool_lst = [[n[0], n[1].get()] for n in pool_lst]
    for target in pool_lst:
        if target[1]:
            pred_target = []
            for pred_tar in target[1]:
                pred_target.append(pred_tar['chembl_id'])
            same_target = [l for l in pred_target if l == user_target]
            # 判断预测出来的靶点与用户制定的靶点是否一致
            if same_target:
                    ligand = target[0].split('.')[0] + '.pdbqt'
                    ligand_path = ligand_db + '/' + ligand
                    # SEA算法预测完成后调用后分子对接软件进行分子对接
                    os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                              " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                                      ligand_path, center_x, center_y, center_z,
                                                                                      size_x, size_y, size_z))
                    os.system("mv %s %s" % (ligand_path.split('.')[0]+'_out.pdbqt', res))
                    os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                            ligand_path.split('/')[-1].split('.')[0]+'_out.pdbqt')))
    else:
        ligand_lst = os.listdir(ligand_db)
        ligand_lst = [v for v in ligand_lst if v.endswith('pdbqt')]
        for ligand_ in ligand_lst:
            ligand_path = os.path.join(ligand_db, ligand_)
            os.system("%s --receptor %s/%s --ligand %s --center_x %s --center_y %s"
                      " --center_z %s --size_x %s --size_y %s --size_z %s" % ('vina', res_path, pdbqt,
                                                                              ligand_path, center_x, center_y,
                                                                              center_z, size_x, size_y, size_z))
            os.system("mv %s %s" % (ligand_path.split('.')[0] + '_out.pdbqt', res))
            os.system("python %s/extra_apps/vina/pdbqt_to_pdb.py -f %s -v" % (BASE_DIR, os.path.join(res,
                                                                                                     ligand_path.split(
                                                                                                         '/')[-1].split(
                                                                                                         '.')[0] + '_out.pdbqt')))
    # 取出所有分子对接任务配体的多个构想打分值
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
            insert_lst.append(SeaVirScreen(work_name=work_name, affinity=min(med),
                                     path=os.path.join(res.split('media/')[1], out[:-2])))
    # 保存虚拟筛选任务的结果
    if screen_res:
        SeaVirScreen.objects.bulk_create(insert_lst)
        arr = np.array(screen_res)
        df = pd.DataFrame(arr, columns=['id', 'Affinity (kcal/mol)'])
        df = df.sort_values("Affinity (kcal/mol)", ascending=False)
        df.to_csv(screen_out, index=False)
        os.system("mv %s %s" % (screen_out, res))
    # 虚拟筛选任务完成之后把热舞状态修改为completed
    screen_status(work_name=work_name, status='completed')
    # 任务完成之后给用户发邮箱提示
    email_status(email, res_path + '/res/screen_out.csv')


@shared_task
def perform_gbsa(work_name, pdb_path, lig_path, complex_path, res_path, email):
    """
    :param work_name:
    :param pdb_path:
    :param lig_path:
    :param complex_path:
    :param res_path:
    :param email:
    :return:
    """
    gbsa_status(work_name=work_name, status='computing')
    res = os.path.join(res_path, 'res')
    if not os.path.exists(res):
        os.mkdir(res)
    # os.system()
    gbsa_status(work_name=work_name, status='completed')
    email_status(email)
