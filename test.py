#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
zscores = np.zeros((3, 3))
zscores = [[i*(j+1) for i in range(0, 3)] for j in range(0, 3)]
print(zscores)

dups = np.zeros((3, 3))
for i in range(0, 3):
    for j in range(0, 3):
        dups[j][i] = i*(j+1)

print(dups)

for i in range(0, 3):
    for j in range(0, 3):
        print(dups[i][j],' ',zscores[i][j])
