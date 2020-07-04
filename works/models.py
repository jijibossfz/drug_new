# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from datetime import datetime
from django.db import models
from django.contrib.auth import get_user_model

from drug import settings
from drug.settings import MEDIA_ROOT
from datetime import datetime
from django.db import models
from django.contrib.auth.models import AbstractUser
# User = get_user_model()


class UserProfile(AbstractUser):
    """
    用户
    """
    name = models.CharField(max_length=10, null=True, blank=True, verbose_name='username')
    mobile = models.CharField(max_length=11, verbose_name='电话')
    email = models.EmailField(max_length=100, null=True, blank=True, verbose_name='email')
    work_org = models.CharField(max_length=20, default='', verbose_name='work_org')
    research_dir = models.CharField(max_length=20, default='', verbose_name='research_dir')
    add_time = models.DateTimeField(default=datetime.now, verbose_name='add_time')

    class Meta:
        verbose_name = 'user'
        verbose_name_plural = 'user'

    def __unicode__(self):
        return self.username


class VerifyCode(models.Model):
    """
    验证码
    """
    code = models.CharField(max_length=6, verbose_name='code')
    email = models.CharField(max_length=100, verbose_name='email')
    add_time = models.DateTimeField(default=datetime.now, verbose_name='add_time')

    class Meta:
        verbose_name = 'verifycode'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.code


class Passwordreset(models.Model):
    email = models.EmailField(max_length=100, verbose_name='email')
    add_time = models.DateTimeField(default=datetime.now, verbose_name='add_time')


def upload_to(instance, filename):
    return '/'.join([MEDIA_ROOT, instance.user.username, instance.work_name, filename])


def dock_upload_to(instance, filename):
    return '/'.join(['dock', instance.user.username, instance.work_name, filename])


def dock2_upload_to(instance, filename):
    return '/'.join(['dock2', instance.user.username, instance.work_name, filename])


def screen_upload_to(instance, filename):
    return '/'.join(['screen', instance.user.username, instance.work_name, filename])


def screen2_upload_to(instance, filename):
    return '/'.join(['screen2', instance.user.username, instance.work_name, filename])


def gbsa_upload_to(instance, filename):
    return '/'.join(['gbsa', instance.user.username, instance.work_name, filename])


class Banner(models.Model):
    """
    首页轮播图
    """
    image = models.ImageField(upload_to='banner/', verbose_name="image")
    index = models.IntegerField(verbose_name="index")
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'banner'
        verbose_name_plural = verbose_name

    def __str__(self):
        return "banner"


class Product(models.Model):
    """
    服务内容
    """
    name = models.CharField(max_length=20, verbose_name='name')
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'product'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.name


class Target(models.Model):
    """
    靶点列表
    """
    target = models.CharField(max_length=20, verbose_name="target_name")

    class Meta:
        verbose_name = 'target_name'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return "靶点列表"


class AutoDock(models.Model):
    """
    分子对接  用户指点中心坐标以及盒子大小
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='user')
    work_name = models.CharField(max_length=20, verbose_name='work_name', unique=True)
    work_decs = models.CharField(max_length=100, default='', verbose_name='work_decs')
    size_x = models.FloatField(verbose_name='size_x')
    size_y = models.FloatField(verbose_name='size_y')
    size_z = models.FloatField(verbose_name='size_z')
    center_x = models.FloatField(verbose_name='center_x')
    center_y = models.FloatField(verbose_name='center_y')
    center_z = models.FloatField(verbose_name='center_z')
    pdb_file = models.FileField(upload_to=dock_upload_to, verbose_name='pdb_file')
    lig_file = models.FileField(upload_to=dock_upload_to, verbose_name='lig_file')
    price = models.IntegerField(default=10000, verbose_name="price")
    status = models.CharField(default='waiting', max_length=10, verbose_name='status')
    out_path = models.FileField(null=True, verbose_name='out_path')
    affinity = models.CharField(default='the position is unreasonable', max_length=100, verbose_name='affinity')
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'dock'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.user.username


class AutoDock2(models.Model):
    """
    分子对接  用户指点对接的残基
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='user')
    work_name = models.CharField(max_length=20, verbose_name='work_name', unique=True)
    work_decs = models.CharField(max_length=100, default='', verbose_name='work_decs')
    pdb_file = models.FileField(upload_to=dock2_upload_to, verbose_name='pdb_file')
    lig_file = models.FileField(upload_to=dock2_upload_to, verbose_name='lig_file')
    resi_file = models.FileField(upload_to=dock2_upload_to, verbose_name='resi_file')
    price = models.IntegerField(default=10000, verbose_name="price")
    status = models.CharField(default='waiting', max_length=10, verbose_name='status')
    out_path = models.FileField(null=True, verbose_name='out_path')
    affinity = models.CharField(default='the position is unreasonable', max_length=100,verbose_name='affinity')
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'dock2'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.user.username


class VirtualScreen(models.Model):
    """
    虚拟筛选
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='user')
    work_name = models.CharField(max_length=20, verbose_name='work_name', unique=True)
    work_decs = models.CharField(max_length=100, default='', verbose_name='work_decs')
    # mol_db = models.CharField(max_length=20, null=True, choices=(('zinc', 'zinc'), ('chembl', 'chembl'), ('wi', 'wi'),
    #                                                              ('user_db_file', 'user_db_file')),
    #                           verbose_name='mol_db')
    size_x = models.FloatField(verbose_name='size_x')
    size_y = models.FloatField(verbose_name='size_y')
    size_z = models.FloatField(verbose_name='size_z')
    center_x = models.FloatField(verbose_name='center_x')
    center_y = models.FloatField(verbose_name='center_y')
    center_z = models.FloatField(verbose_name='center_z')
    target = models.CharField(max_length=100000, verbose_name='target', default='')
    pdb_file = models.FileField(upload_to=screen_upload_to, verbose_name='pdb_file')
    user_db = models.FileField(upload_to=screen_upload_to, null=True, verbose_name='user_db')
    price = models.IntegerField(default=10000, verbose_name="price")
    status = models.CharField(default='waiting', max_length=10, verbose_name='status')
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'virtualscreen'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.user.username


class VirtualScreen2(models.Model):
    """
    虚拟筛选2
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='user')
    work_name = models.CharField(max_length=20, verbose_name='work_name', unique=True)
    work_decs = models.CharField(max_length=100, default='', verbose_name='work_decs')
    # mol_db = models.CharField(max_length=20, null=True, choices=(('zinc', 'zinc'), ('chembl', 'chembl'), ('wi', 'wi'),
    #                                                              ('user_db_file', 'user_db_file')),
    #                           default=0, verbose_name='mol_db')
    target = models.CharField(max_length=100000, verbose_name='target', default='')
    pdb_file = models.FileField(upload_to=screen2_upload_to, verbose_name='pdb_file')
    resi_file = models.FileField(upload_to=screen2_upload_to, verbose_name='resi_file')
    user_db = models.FileField(upload_to=screen2_upload_to, null=True, verbose_name='user_db')
    price = models.IntegerField(default=10000, verbose_name="price")
    status = models.CharField(default='waiting', max_length=10, verbose_name='status')
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'virtualscreen2'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.user.username


class Screen(models.Model):
    work_name = models.CharField(max_length=100, verbose_name='work_name')
    screen_cat = models.CharField(max_length=10, verbose_name='screen_cat')
    affinity = models.FloatField(verbose_name='affinity')
    path = models.FileField(max_length=100, verbose_name='out_file')

    class Meta:
        verbose_name = 'screen'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.work_name


# class VsBlast(models.Model):
#     """
#     VsleadBlast
#     """
#     user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='用户')
#     work_name = models.CharField(max_length=20, verbose_name='任务名称')
#     work_decs = models.CharField(max_length=100, default='', verbose_name='任务描述')
#     sequence = models.CharField(max_length=1000, verbose_name='蛋白序列')
#     protein_db = models.CharField(max_length=10, choices=((1, 'zinc'), (2, 'chembl')))
#     e_value = models.CharField(max_length=10, choices=((1, 0.00001), (2, 0.0001), (3, 0.001),
#                                                        (4, 0.01), (5, 0.01), (6, 1), (7, 10)),
#                                verbose_name='E值选择')
#     out_format = models.CharField(max_length=10, choices=((1, 'pariwise'), (2, 'XML blast output'),
#                                                           (3, 'tabular')), verbose_name='输出格式选择')
#     add_time = models.DateTimeField(default=datetime.now, verbose_name="添加时间")
#
#     class Meta:
#         verbose_name = 'VsleadBlast'
#         verbose_name_plural = verbose_name
#
#     def __unicode__(self):
#         return self.user.username


# class ReverseVirtualScreen(models.Model):
#     """
#     反向虚拟筛选
#     """
#     user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='用户')
#     work_name = models.CharField(max_length=20, verbose_name='任务名称')
#     work_decs = models.CharField(max_length=100, default='', verbose_name='任务描述')
#     target_db = models.CharField(max_length=10, choices=((1, 'zine'), (2, 'chembl')), verbose_name='靶点数据库')
#     mol_file = models.FileField(upload_to=upload_to, verbose_name='靶点文件')
#     top_n = models.IntegerField(verbose_name='返回结果数目')
#     add_time = models.DateTimeField(default=datetime.now, verbose_name="添加时间")
#
#     class Meta:
#         verbose_name = 'ReverseVirtualScreen'
#         verbose_name_plural = verbose_name
#
#     def __unicode__(self):
#         return self.user.username


# class Dynamic(models.Model):
#     """
#     动力学模拟
#     """
#     user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='用户')
#     work_name = models.CharField(max_length=20, verbose_name='任务名称')
#     work_decs = models.CharField(max_length=100, default='', verbose_name='任务描述')
#     conf_info = models.CharField(max_length=20, choices=((1, '通用计算'), (2, '计算结合能'), (3, '计算氢键'),
#                                                          (4, '聚类分析')), verbose_name='配置信息')
#     protein_file = models.FileField(upload_to=upload_to, verbose_name='蛋白文件')
#     mol_file = models.FileField(upload_to=upload_to, verbose_name='小分子文件')
#     conf_project = models.CharField(max_length=100, choices=((1, '有水'), (2, '无水'), (3, '自计算')),
#                                     verbose_name='配置项目')
#     s_file = models.FileField(upload_to=upload_to, default='', verbose_name='S信息文件')
#     lig_file = models.FileField(upload_to=upload_to, default='', verbose_name='lig文件')
#     frcmod_file = models.FileField(upload_to=upload_to, default='', verbose_name='frcmod文件')
#     res_num = models.IntegerField(default=0, verbose_name='氨基酸数目')
#     add_time = models.DateTimeField(default=datetime.now, verbose_name="添加时间")
#
#     class Meta:
#         verbose_name = 'Dynamic'
#         verbose_name_plural = verbose_name
#
#     def __unicode__(self):
#         return self.user.username


class Admet(models.Model):
    """
    ADMET预测
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='user')
    work_name = models.CharField(max_length=20, verbose_name='work_name')
    work_decs = models.CharField(max_length=100, default='', verbose_name='work_decs')
    mol_file = models.FileField(upload_to=upload_to, verbose_name='mol_file')
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'ADMET_prediction'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.user.username


class Dock(models.Model):
    """
    分子对接  用户指点中心坐标以及盒子大小
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='user')
    work_name = models.CharField(max_length=20, verbose_name='work_name', unique=True)
    algorithm = models.CharField(max_length=20, choices=(('AutoDock', 'AutoDock'), ('AutoDock-vina', 'AutoDock-vina'),
                                                         ('Glide', 'Glide')))
    size_x = models.FloatField(verbose_name='size_x', null=True)
    size_y = models.FloatField(verbose_name='size_y', null=True)
    size_z = models.FloatField(verbose_name='size_z', null=True)
    center_x = models.FloatField(verbose_name='center_x', null=True)
    center_y = models.FloatField(verbose_name='center_y', null=True)
    center_z = models.FloatField(verbose_name='center_z', null=True)
    resn = models.CharField(max_length=100, verbose_name='resn', null=True)
    pdb_file = models.FileField(upload_to=dock_upload_to, verbose_name='pdb_file')
    lig_file = models.FileField(upload_to=dock_upload_to, verbose_name='lig_file')
    reference_file = models.FileField(upload_to=dock_upload_to, verbose_name='reference_file', null=True)
    status = models.CharField(default='waiting', max_length=10, verbose_name='status')
    out_path = models.FileField(null=True, verbose_name='out_path')
    affinity = models.CharField(default='the position is unreasonable', max_length=100, verbose_name='affinity')
    email = models.CharField(max_length=100, verbose_name='email')
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'dock'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.user.username


class VirScreen(models.Model):
    """
    虚拟筛选
    """
    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='user')
    work_name = models.CharField(max_length=20, verbose_name='work_name', unique=True)
    mol_db = models.CharField(max_length=20, null=True, choices=(('zinc', 'zinc'), ('chembl', 'chembl'),
                                                                 ('drugbank', 'drugbank'),
                                                                 ('chinese-medicine', 'chinese-medicine'),
                                                                 ('taosu', 'taosu'), ('bailingwei', 'bailingwei'),
                                                                 ('jianqiao', 'jianqiao')), verbose_name='mol_db')
    target = models.CharField(max_length=100, verbose_name='target', null=True)
    smiles = models.CharField(max_length=200, verbose_name='smiles', null=True)
    pdb_file = models.FileField(upload_to=screen_upload_to, verbose_name='pdb_file', null=True)
    resn = models.CharField(max_length=100, verbose_name='resn', null=True)
    status = models.CharField(default='waiting', max_length=10, verbose_name='status')
    email = models.CharField(default="", max_length=100, verbose_name='email')
    add_time = models.DateTimeField(default=datetime.now, verbose_name="add_time")

    class Meta:
        verbose_name = 'VirScreen'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.user.username


class SeaTarget(models.Model):
    work_name = models.CharField(max_length=100, verbose_name='work_name')
    smiles = models.CharField(max_length=200, verbose_name='smiles')
    target = models.CharField(max_length=500, verbose_name='targets', default='null')

    class Meta:
        verbose_name = 'SeaTarget'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.work_name


class SeaVirScreen(models.Model):
    work_name = models.CharField(max_length=100, verbose_name='work_name')
    affinity = models.FloatField(verbose_name='affinity')
    path = models.FileField(max_length=100, verbose_name='out_file')

    class Meta:
        verbose_name = 'SeaVirScreen'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.work_name


class Gbsa(models.Model):
    user = models.ForeignKey(settings.AUTH_USER_MODEL, verbose_name='user')
    work_name = models.CharField(max_length=100, verbose_name='work_name')
    pdb_file = models.FileField(upload_to=gbsa_upload_to, verbose_name='pdb_file')
    lig_file = models.FileField(upload_to=gbsa_upload_to, verbose_name='lig_file')
    complex_file = models.FileField(upload_to=gbsa_upload_to, verbose_name='complex_file')
    status = models.CharField(max_length=20, default='waiting', verbose_name='status')
    email = models.CharField(max_length=100, verbose_name='email')
    add_time = models.DateTimeField(default=datetime.now, verbose_name='add_time')

    class Meta:
        verbose_name = 'gbsa'
        verbose_name_plural = verbose_name

    def __unicode__(self):
        return self.work_name
