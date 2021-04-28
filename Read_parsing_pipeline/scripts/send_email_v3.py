#!/usr/bin/python

# -------------------------------------------------
#
# Send email
#
# ------------------------------------------------

# To do:

import sys
import smtplib
import time
from datetime import datetime, timedelta
import imaplib
import email
#import html2text
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart


def send_email(email_from, email_to, msg, testMode):

    # Send the email via our own SMTP server.
    try:
        s = smtplib.SMTP('localhost')
        s.sendmail(email_from, email_to, msg.as_string())
        s.quit()

    except:
        print("Error: unable to send email")
        print(sys.exc_info()[0], "occured.")


def main(argv):

    toEmail = argv[0]
    subject = argv[1]
    message = argv[2]

    testMode = False

    # Create the container (outer) email message.
    me = 'prm88@cornell.edu'
    # you = ['6072802109@vtext.com', 'programmeratlarge@gmail.com']
    # msg['To'] = ', '.join(you)
    toEmailList = toEmail.split(',')
    subject = subject.strip('\"')
    body = message.strip('\"')
    attachment_type = 'plain'

    for you in toEmailList:
        if you.split('@')[0].isdigit():
            msg = MIMEText(body)
            msg['From'] = me
            msg['To'] = you
            msg['Subject'] = subject
        else:
            msg = MIMEMultipart()
            msg.preamble = subject
            msg['From'] = me
            msg['To'] = you
            msg['Subject'] = subject
            msg.attach(MIMEText(body, attachment_type))

        send_email(me, you, msg, testMode)


if __name__ == "__main__":
    main(sys.argv[1:])

