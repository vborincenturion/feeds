#!/usr/bin/env python

import os
import urllib.request
import zipfile

# Create the "db" directory if it doesn't exist
if not os.path.exists('db'):
    os.makedirs('db')

# URL of the file to download
url = 'https://figshare.com/ndownloader/articles/22194535/versions/1'

# Name of the file to download
filename = '22194535.zip'

# Download the file from the URL and save it to the "db" directory
urllib.request.urlretrieve(url, os.path.join('db', filename))

# Unzip the downloaded file
with zipfile.ZipFile(os.path.join('db', filename), 'r') as zip_ref:
    zip_ref.extractall('db')

# Remove the ZIP file
os.remove(os.path.join('db', filename))
