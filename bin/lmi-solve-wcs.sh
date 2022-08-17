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
		scalelow=`python -c "$binning * 0.12 - 0.01" | tr -d '[:space:]'`
		scalehigh=`python -c "$binning * 0.12 + 0.01" | tr -d '[:space:]'`
		solve-field --no-verify --use-source-extractor --skip-solved --scale-low=$scalelow --scale-high=$scalehigh --scale-units=arcsecperpix -t0 --radius=0.5 --ra=$ra --dec=$dec $g --depth=50
		if [ -e ${g%.fits}.solved ]; then
			#cphead ${g%.fits}.wcs $f
			python3 <<EOF
from astropy.io import fits
h = fits.getheader("${g%.fits}.wcs")
with fits.open("$f", mode='update') as hdu:
	hdu[0].header.update(h)
	hdu[0].header["LMIWCS"] = (1, "WCS solution added by lmi-solve-wcs.sh")
EOF
		fi
    fi
done
