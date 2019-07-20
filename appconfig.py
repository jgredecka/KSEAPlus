from datetime import timedelta
import os

# 6 MB file upload limit
MAX_CONTENT_LENGTH = 6 * 1024 * 1024
# session expiry
PERMANENT_SESSION_LIFETIME = timedelta(minutes=6)
# secret key
SECRET_KEY = os.urandom(24)
DEBUG = True
