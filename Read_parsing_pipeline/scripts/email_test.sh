#!/usr/bin/bash
set -x #echo on

# Test send of email via Python
python send_email_v3.py prm88@cornell.edu "Test send email v03" "Message body"
python send_email_v3.py 6072802109@vtext.com "Test send email v03" "Text message body"

echo "Done!"
