# -*- coding: utf-8 -*-

from pyArango.connection import *

class querydb:
    def __init__(self):
        self.conn = Connection(arangoURL='http://cudaws02.ime.kfa-juelich.de:8529', username='root', password='root_passwd')
