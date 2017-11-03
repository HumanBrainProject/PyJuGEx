#!/usr/bin/env python
# -*- coding: utf-8 -*-

def upload(filename):
    clients = get_hbp_service_client()
    collab_path = get_collab_storage_path()
    if not clients.storage.exists(collab_path + '/pmaps'):
        pmap_folder = clients.storage.mkdir(collab_path + '/pmaps')
    else:
        pmap_folder = collab_path + '/pmaps'
    print(pmap_folder)
    ufilename = collab_path + '/pmaps/' + filename
    if clients.storage.exists(ufilename):
        clients.storage.delete(ufilename)
    if filename.split('.')[-1] == 'gz':
        clients.storage.upload_file(filename, ufilename, 'image/gznii')
    elif filename.split('.')[-1] == 'nii':
        clients.storage.upload_file(filename, ufilename, 'image/nii')
    else:
        print('Not a valid file format')
        exit()
