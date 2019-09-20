#!/bin/bash
for f in $* ; do
    imagetype=`gethead IMAGETYP $f`
    if echo $imagetype|grep -q OBJECT; then
	g=${f##*/}
	ln -s $f $g

        radec=`dct-radec.py $f|awk '{print $2, $3}'`
	ra=`echo $radec|awk '{print $1}'`
	dec=`echo $radec|awk '{print $2}'`
	binning=`gethead CCDSUM $f|awk '{print $1}'`
	scalelow=`calc "$binning * 0.12 - 0.01"`
	scalehigh=`calc "$binning * 0.12 + 0.01"`
	solve-field --no-verify --use-sextractor --skip-solved --scale-low=$scalelow --scale-high=$scalehigh --scale-units=arcsecperpix -t2 --radius=0.5 --ra=$ra --dec=$dec $g --depth=50
	if [ -e ${g%.fits}.solved ]; then
	    #cphead ${g%.fits}.wcs $f
	    python3 <<EOF
from astropy.io import fits
h = fits.getheader("${g%.fits}.wcs")
with fits.open("$f", mode='update') as hdu:
    hdu[0].header.update(h)
EOF
	fi
    fi
done
