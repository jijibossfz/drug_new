# coding=utf-8
import xadmin

from .models import Banner, Product, AutoDock, AutoDock2, VirtualScreen, VirtualScreen2, Screen


class BannerAdmin(object):
    list_display = ['image', 'index', 'add_time']
    search_fields = ['image', 'index', 'add_time']
    list_filter = ['image', 'index', 'add_time']


class ProductAdmin(object):
    list_display = ['name', 'add_time']
    search_fields = ['name', 'add_time']
    list_filter = ['name', 'add_time']


class AutoDoctAdmin(object):
    list_display = ['user', 'work_name', 'work_decs', 'size_x', 'size_y', 'size_z', 'center_x', 'center_y', 'center_z',
                    'pdb_file', 'lig_file', 'out_path', 'affinity', 'status', 'add_time']
    search_fields = ['user', 'work_name', 'work_decs', 'size_x', 'size_y', 'size_z', 'center_x', 'center_y', 'center_z',
                     'pdb_file', 'lig_file', 'out_path', 'affinity', 'status', 'add_time']
    list_filter = ['user', 'work_name', 'work_decs', 'size_x', 'size_y', 'size_z', 'center_x', 'center_y', 'center_z',
                   'pdb_file', 'lig_file', 'out_path', 'affinity', 'status', 'add_time']


class Autodock2Admin(object):
    list_display = ['user', 'work_name', 'work_decs', 'lig_file', 'pdb_file','out_path', 'affinity', 'status',
                    'add_time']
    search_fields = ['user', 'work_name', 'work_decs', 'lig_file', 'pdb_file', 'out_path', 'affinity', 'status',
                     'add_time']
    list_filter = ['user', 'work_name', 'work_decs', 'lig_file', 'pdb_file', 'out_path', 'affinity', 'status',
                   'add_time']


class VirtualScreenAdmin(object):
    list_display = ['user', 'work_name', 'work_decs', 'mol_db', 'user_db', 'target', 'size_x', 'size_y', 'size_z',
                    'center_x', 'center_y', 'center_z', 'pdb_file', 'status', 'add_time']
    search_fields = ['user', 'work_name', 'work_decs', 'mol_db', 'user_db', 'target', 'size_x', 'size_y', 'size_z',
                     'center_x', 'center_y', 'center_z', 'pdb_file', 'status', 'add_time']
    list_filter = ['user', 'work_name', 'work_decs', 'mol_db', 'user_db', 'target', 'size_x', 'size_y', 'size_z',
                   'center_x', 'center_y', 'center_z', 'pdb_file', 'status', 'add_time']


class VirtualScreen2Admin(object):
    list_display = ['user', 'work_name', 'work_decs', 'mol_db', 'user_db', 'target', 'pdb_file', 'status', 'add_time']
    search_fields = ['user', 'work_name', 'work_decs', 'mol_db', 'user_db', 'target', 'pdb_file', 'status', 'add_time']
    list_filter = ['user', 'work_name', 'work_decs', 'mol_db', 'user_db', 'target', 'pdb_file', 'status', 'add_time']


class ScreenAdmin(object):
    list_display = ['work_name', 'screen_cat', 'affinity', 'path']
    search_fields = ['work_name', 'screen_cat', 'affinity', 'path']
    list_filter = ['work_name', 'screen_cat', 'affinity', 'path']


xadmin.site.register(Banner, BannerAdmin)
xadmin.site.register(Product, ProductAdmin)
xadmin.site.register(AutoDock, AutoDoctAdmin)
xadmin.site.register(AutoDock2, Autodock2Admin)
xadmin.site.register(VirtualScreen, VirtualScreenAdmin)
xadmin.site.register(VirtualScreen2, VirtualScreen2Admin)
xadmin.site.register(Screen, ScreenAdmin)
