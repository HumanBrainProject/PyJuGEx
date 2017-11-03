#!/usr/bin/env python
# -*- coding: utf-8 -*-

def upload(filename):
    client = get_hbp_service_client()
    collab_path = get_collab_storage_path()
    if not clients.storage.exists(collab_path + '/pmaps'):
        pmap_folder = client.storage.mkdir(collab_path + '/pmaps')
    else:
        pmap_folder = collab_path + '/pmaps'
    ufilename = collab_path + '/pmaps/' + filename
    if filename.split('.')[-1] == 'gz':
        client.storage.upload_file(filename, ufilename, 'image/gznii')
    elif filename.split('.')[-1] == 'nii':
        client.storage.upload_file(filename, ufilename, 'image/nii')
    else:
        print('Not a valid file format')
        exit()
