# coding=utf-8
from __future__ import unicode_literals
from django.core.mail import send_mail, EmailMultiAlternatives
from drug.settings import EMAIL_FROM

subject = 'test'
content = 'boy'
from_email = EMAIL_FROM
to_email = 'my_liulei@sina.com'


def send_emails(subject, content, from_email, to_email):
    msg = EmailMultiAlternatives(subject, content, from_email, [to_email])
    # msg.content_subtype = 'html'
    msg.attach_file('./email_send.py')
    msg.send()


send_emails(subject, content, from_email, to_email)
