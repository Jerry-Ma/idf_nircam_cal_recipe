#!/usr/bin/env python


import sys
from astropy.table import Table
import numpy as np
from astropy.io import fits
from pathlib import Path
from astropy.wcs import WCS


scriptdir = Path(__file__).parent


if __name__ == "__main__":

    tbl = Table.read(scriptdir.joinpath("/tables/idf_zpcor_1130_pmap.ecsv"))

    print(tbl)

    for f in sys.argv[1:]:
        suffix = "_cal.fits"
        if not f.endswith(suffix):
            raise ValueError("invalid input file")
        hl = fits.open(f)
        hl.info()
        hdr = hl[0].header
        det = hdr['DETECTOR'].lower().strip()
        filt = hdr['FILTER'].lower().strip()
        mod = det[3:4]
        det = det[4:]
        print(mod, det, filt)
        if 'ep1' in Path(f).parent.name:
            ep = 1
        elif 'ep2' in Path(f).parent.name:
            ep = 2
        else:
            raise ValueError("cannot get ep")

        # get entry
        e = tbl[(tbl['det'] == det) & (tbl['filt'] == filt)][0]
        scale = 10 ** (e[f'dmag_{mod}'] / 2.5)

        print(f"apply scaling factor: {scale}")
        hl[1].data = hl[1].data * scale

        # apply slope correction if they > 0
        if e['slope'] != 0:
            k = e['slope']
            print(f"apply slope {k=} to {mod=} {det=} {filt=} {f}")
            wcsobj = WCS(hl[1].header)
            ps = wcsobj.proj_plane_pixel_scales()[0]
            print(ps)
            # import pdb
            # pdb.set_trace()
            ny, nx = hl[1].data.shape
            # the negative sign takes into account the filpped x axis
            if ep == 1:
                # kk = -abs(k) / 2
                kk = 0
            else:
                # dimmer at nx
                kk = -abs(k)
            gm = kk * np.arange(nx) * ps.to_value('arcmin')
            gscale = 10 ** (gm / 2.5)
            print(f"apply gradient scaling {gscale[0]} @0 {gscale[-1]} @nx")
            hl[1].data = hl[1].data * gscale.reshape((1, -1))
        fout = f.replace(suffix, '_cal_scaled.fits')
        hl.writeto(fout, overwrite=True)
