#!/usr/bin/env python

import os
import urllib.request
import zipfile

# Create the "db" directory if it doesn't exist
if not os.path.exists('db'):
    os.makedirs('db')

# URL of the file to download
url = 'https://figshare.com/ndownloader/articles/22194535/versions/3'

# Name of the file to download
filename = '22194535.zip'

# Define the callback function to show the download progress
def reporthook(count, block_size, total_size):
    progress = count * block_size * 100 / total_size
    print(f"Downloaded {count * block_size} bytes out of {total_size} ({progress:.1f}%)", end="\r")

# Download the file from the URL and save it to the "db" directory
urllib.request.urlretrieve(url, os.path.join('db', filename), reporthook=reporthook)

# Unzip the downloaded file
with zipfile.ZipFile(os.path.join('db', filename), 'r') as zip_ref:
    zip_ref.extractall('db')

# Remove the ZIP file
os.remove(os.path.join('db', filename))
