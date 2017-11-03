#!/usr/bin/env python
# -*- coding: utf-8 -*-
import hbp_service_client.client
class collab_pyjugex:
    @classmethod
    def upload(cls, filename):
        client = get_hbp_service_client()
        collab_path = get_collab_storage_path()
        pmap_folder = collab_path + '/pmaps/'
        if not clients.storage.exists(pmap_folder):
            client.storage.mkdir(pmap_folder)
        ufilename = pmap_folder + filename
        if filename.split('.')[-1] == 'gz':
            client.storage.upload_file(filename, ufilename, 'image/gznii')
        elif filename.split('.')[-1] == 'nii':
            client.storage.upload_file(filename, ufilename, 'image/nii')
        else:
            print('Not a valid file format')
            exit()
