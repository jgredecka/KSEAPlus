from datetime import timedelta
import os

# 6 MB file upload limit
MAX_CONTENT_LENGTH = 6 * 1024 * 1024
# session expiry
PERMANENT_SESSION_LIFETIME = timedelta(minutes=6)
# secret key
SECRET_KEY = '\xea\x99DkZN\xc1\x9f\x01n/\xa5\x91\x9f\xfc\x02-\x03\n\x8d\x0f\xd5\xcas'
DEBUG = True