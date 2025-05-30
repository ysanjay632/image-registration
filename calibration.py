import sunpy.map
import glob
import numpy as np
import aiapy.data.sample as sample_data
from aiapy.calibrate import  register
from sunpy.io import write_file
import sunkit_image.coalignment as c
import astropy.units as u

import matplotlib.pyplot as plt
from sunpy.image import resample


def align_and_register(file_path):
    filelist = glob.glob(file_path + "/*.fits")
    map_list = [sunpy.map.Map(f) for f in filelist]
    maps = [[] for _ in range(len(map_list))]
    smap = [[] for _ in range(len(map_list))]

    for i in range(len(map_list)):
        smap[i] = register(map_list[i], missing=0, order=3, method='scipy')
        smap[i].meta['naxis1'] = 4096
        smap[i].meta['naxis2'] = 4096

        width_diff = 4096 - smap[i].data.shape[0]
        height_diff = 4096 - smap[i].data.shape[1]
        crpix1 = smap[i].meta['crpix1'] + width_diff / 2
        crpix2 = smap[i].meta['crpix2'] + height_diff / 2

        smap[i].meta['crpix1'] = crpix1
        smap[i].meta['crpix2'] = crpix2

        sm = resample.resample(smap[i].data, [4096, 4096])
        maps[i] = sunpy.map.Map(sm, smap[i].meta)
        print(i)

    map_sequence = sunpy.map.MapSequence(maps)
    data = [[] for _ in range(len(maps))]
  
    yshift = [[] for _ in range(len(maps))]
    xshift = [[] for _ in range(len(maps))]
    y = [[] for _ in range(len(maps))]
    x = [[] for _ in range(len(maps))]
    sy = [[] for _ in range(len(maps))]
    sx = [[] for _ in range(len(maps))]

    for i in range(len(maps)):
        data[i] = ((maps[i].data))
       

    temp = data[0][2200:3000, 2000:3000]  # can be chosen differently
    print(np.max(temp))

    for i in range(len(maps)):
        shift = c.calculate_shift(data[i], temp)
        print(shift[1], shift[0])
        yshift[i] = shift[0]
        xshift[i] = shift[1]

    for i in range(len(maps)):
        sx[i] = xshift[0] - xshift[i]
        sy[i] = yshift[0] - yshift[i]

    y = u.Quantity(sy, unit=u.pix)
    x = u.Quantity(sx, unit=u.pix)
    print(x, y)
    aligned = c.apply_shifts(map_sequence, y, x, clip=False)

    for i in range(len(aligned)):
        print(aligned[i].data.shape)
        print(maps[i].data.shape)

    for i in range(len(maps)):
        t = aligned[i].meta['DATE-OBS']
        t = t.replace(':', '_')
        file = f'AIA_1700_{t}.fits'

        write_file(file, aligned[i].data, aligned[i].meta, filetype= 'fits')
align_and_register(file_path)
