# coding=utf-8

from __future__ import unicode_literals

import rest_framework_jwt
from django.conf.urls import url, include
from drug.settings import MEDIA_ROOT
from django.views.static import serve
from rest_framework.documentation import include_docs_urls
from rest_framework.routers import DefaultRouter
from works.views import UserRegViewset, EmailCodeViewset, AutoDockOrderViewset, AutoDock2OrderViewset, \
    VirtualScreenOrderViewset, VirtualScreen2OrderViewset, UserinfoViewset, PasswordresetViewset
from works.views import BannerViewset, ProductViewset, AutoDockViewset, VirtualScreenViewset
from works.views import AdmetViewset
from works.views import AutoDock2Viewset, VirtualScreen2Viewset, ScreenViewset, TargetViewset, GbsaViewset, \
    GbsaOrderViewset, DockViewset, DockOrderViewset, VirScreenViewset, SeaTargetViewset, SeaVirScreenOrderViewset, \
    VirScreenOrderViewset
import xadmin
from rest_framework.authtoken import views
from rest_framework_jwt.views import obtain_jwt_token, refresh_jwt_token, verify_jwt_token

router = DefaultRouter()

router.register(r'docks', DockViewset, base_name='dock')
router.register(r'dockorder', DockOrderViewset, base_name='dockorder')
router.register(r'virscreens', VirScreenViewset, base_name='virscreen')
router.register(r'virscreenorder', VirScreenOrderViewset, base_name='virscreenorder')
router.register(r'seatargetorder', SeaTargetViewset, base_name='seatarget')
# router.register(r'seavirscreenorder', SeaVirScreenOrderViewset, base_name='seavirscreenorder')
# router.register(r'gbsa', GbsaViewset, base_name='gbsa')
# router.register(r'gbsaorder', GbsaOrderViewset, base_name='gbsaorder')
router.register(r'codes', EmailCodeViewset, base_name='验证码')
router.register(r'banners', BannerViewset, base_name='首页轮播图')
router.register(r'products', ProductViewset, base_name='服务内容')
# router.register(r'targets', TargetViewset, base_name='靶点列表')
# router.register(r'autodocks', AutoDockViewset, base_name='分子对接')
# router.register(r'autodock2s', AutoDock2Viewset, base_name='分子对接2')
# router.register(r'virtualscreens', VirtualScreenViewset, base_name='虚拟筛选')
# router.register(r'virtualscreen2s', VirtualScreen2Viewset, base_name='虚拟筛选2')
# router.register(r'vsleadblasts', VsBlastViewset, base_name='vsleadblast')
# router.register(r'reversevirtualscreens', ReverseVirtualScreenViewset, base_name='反向虚拟筛选')
# router.register(r'dynamics', DynamicViewset, base_name='动力学模拟')
# router.register(r'admets', AdmetViewset, base_name='admet预测')
router.register(r'registers', UserRegViewset, base_name='用户注册')
# router.register(r'autodockorders', AutoDockOrderViewset, base_name='autoduckorder')
# router.register(r'autodock2orders', AutoDock2OrderViewset, base_name='autoduck2order')
# router.register(r'virtualscreenorders', VirtualScreenOrderViewset, base_name='virtualscreenorders')
# router.register(r'virtualscreen2orders', VirtualScreen2OrderViewset, base_name='virtualscreen2orders')
router.register(r'users', UserinfoViewset, base_name="user")
router.register(r'passwordreset', PasswordresetViewset, base_name='passwordreset')
router.register(r'screens', ScreenViewset, base_name='screen')


urlpatterns = [
    url(r'', include(router.urls)),
    url(r'^xadmin/', xadmin.site.urls),
    url(r'^media/(?P<path>.*)$', serve, {"document_root": MEDIA_ROOT}),
    url(r'^docs/', include_docs_urls(title="高通量药物筛选平台")),
    url(r'^api-token-auth/', rest_framework_jwt.views.obtain_jwt_token),
    url(r'^api-token-refresh/', rest_framework_jwt.views.refresh_jwt_token),
    url(r'^login/', obtain_jwt_token),
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    # url(r'^api-token-auth/', obtain_jwt_token),
    # url(r'^api-token-refresh/', refresh_jwt_token),
    url(r'^api-token-verify/', verify_jwt_token),
    url(r'^rest-auth/', include('rest_auth.urls')),
    ]
