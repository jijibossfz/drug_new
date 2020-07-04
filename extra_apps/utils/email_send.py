# coding:utf-8
from __future__ import unicode_literals

import os
from email.header import make_header

from django.core.mail import send_mail, EmailMultiAlternatives

from drug import settings
from drug.settings import EMAIL_FROM


def email_send(email, code, send_type="register"):

    if send_type == "register":
        email_title = "邮箱注册验证码"
        email_body = "注册高通量药物筛选平台邮箱验证码为 %s" % code

        send_status = send_mail(subject=email_title, message=email_body,
                                from_email=EMAIL_FROM,
                                recipient_list=[email])
        return send_status
    elif send_type == "forget":
        email_title = "邮箱找回验证码"
        email_body = "重置高通量药物筛选平台邮箱验证码为 %s" % code

        send_status = send_mail(subject=email_title, message=email_body,
                                from_email=EMAIL_FROM, recipient_list=[email])
        return send_status


def email_status(email,file_path):
    email_title = '虚拟筛选服务器'
    email_body = '您的任务已经完成'
    from_email = settings.DEFAULT_FROM_EMAIL
    msg = EmailMultiAlternatives(email_title, email_body, from_email, [email])
    if len(file_path) > 0:
        text = open(file_path, 'rb').read()
        file_name = os.path.basename(file_path)
        b = make_header([(file_name, 'utf-8')]).encode('utf-8')
        msg.attach(b, text)
        msg.send()
    else:
        msg.send()
    # send_status = send_mail(subject=email_title, message=email_body, from_email=EMAIL_FROM, recipient_list=[email], fail_silently=False)
    # return send_status


def email_status_sea(email,insert_list):
    email_title = '虚拟筛选服务器'
    email_body = '您的任务已经完成:'+str(insert_list)
    from_email = settings.DEFAULT_FROM_EMAIL
    msg = EmailMultiAlternatives(email_title, email_body, from_email, [email])
    msg.send()