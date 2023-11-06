# dct-redux v0.8.1
Python scripts for reducing Lowell Discovery Telescope (formerly Discovery Channel Telescope) data, mostly focused on observations of comets.

## Status
Most, but not all, of these scripts are being routinely used with good results.  However, it is best to keep safe backups of your data to guard against errors.

## Requirements
* Python 3
* astropy
* ccdproc >=2.2
* sep
* DS9
* pyds9
* astroquery
* sbpy >= 0.4
* calviacat
* mskpy
* MontagePy (for lmi-stack-comets.py)

For file name normalization and FITS header fixing, just `astropy` is needed.

With the data in a `raw/` sub-directory, the scripts will create and use a directory tree similar to the following:

```
.
├── all.list
├── file-summary.txt
├── lmi-add-cat.log
├── lmi-catalog-calibrate.log
├── lmi-rx.log
├── msk-notes.txt
├── phot
│   ├── cal-standard-phot-colorcor.pdf
│   ├── cal-standard-phot-colorcor.png
│   ├── cal-standard-phot-extinction.pdf
│   ├── cal-standard-phot-extinction.png
│   ├── cal-standard-phot-oh.pdf
│   ├── cal-standard-phot-oh.png
│   ├── cal-standard-phot.txt
│   ├── catalog-cal-summary.txt
│   ├── catalog-extinction.txt
│   ├── cat.db
│   ├── lmi-catalog-calibrate.log
│   ├── lmi-standard-phot.log
│   ├── lmi-zeropoints.png
│   ├── standard-phot.txt
│   └── zeropoints
│       ├── lmi_20160212_0036_ppp.png
│       ├── ...
│       └── lmi_20160212_0275_ppp.png
├── ppp
│   ├── bias_20160212_x1to6144_y1to6160_2x2.fits
│   ├── lmi_20160212_0001_ppp.fits
│   ├── ...
│   ├── lmi_20160212_0275_ppp.fits
│   ├── skyflat_20160212_BC-2x2-D.fits
│   ├── skyflat_20160212_CN-2x2-D.fits
│   ├── skyflat_20160212_OH-2x2-D.fits
│   ├── skyflat_20160212_RC-2x2-D.fits
│   └── skyflat_20160212_SDSS-R-2x2-D.fits
├── raw
│   ├── fixes.txt
│   ├── lmi_20160212_0001_raw.fits
│   ├── ...
│   ├── lmi_20160212_0275_raw.fits
│   └── original-headers
│       └── lmi_20160212_0036_raw.header
├── sorted
│   ├── 15p_finlay
│   │   └── SDSS-R
│   │       ├── lmi_20160212_0173_ppp.fits -> ../../../ppp/lmi_20160212_0173_ppp.fits
│   │       ├── ...
│   │       └── lmi_20160212_0176_ppp.fits -> ../../../ppp/lmi_20160212_0176_ppp.fits
│   ├── 21p_giacobini-zinner
│   │   └── SDSS-R
│   │       └── ...
│   └── ...
├── target-summary.csv
└── target-summary.txt
```

## General steps

### LMI semi-automated data processing

The bash script `lmi-rx-to-calibration.sh` guides the user through most of the scripts with interactive prompts.  It requires a `raw/` directory with the raw data files to be processed.  It assumes all the data were taken from a single night.

### For any LDT data:
1. File name normalization.
   ```bash
   cd raw
   dct-filename-normalization.py lmi.*.fits
   ```
1. FITS header fixes.  Create a file of the fixes to be applied, and apply them:
   ```bash
   cat >header-fixes.txt <<EOF
   first frame, last frame, keyword, value
   11, 26,	object|objname, 96 Her
   27, 42, object|objname, HD 219188
   43, 56, object|objname, C/2014 Q2 (Lovejoy)
   57, 70, object|objname, 53P/Van Biesbroeck
   EOF
   dct-fits-header-fixes.py header-fixes.txt
   ```
### For LMI data

1. Create bias and flat field frames, apply to data to create partially processed products (PPP files):
   ```bash
   cd ..
   ls raw/*fits > all.list
   # edit bad files from all.list
   emacs all.list
   lmi-rx.py all.list
   ```

   Need to reuse a flat or bias from another night?  Make a symbolic link or copy the file to ppp/ and lmi-rx will find it.

   ```bash
   cd ppp
   ln -s ../../20210113/ppp/bias1_3x3.fits
   cd ..
   lmi-rx.py all.list
   ```

1. Add world coordinates: (**LMI's nominal WCS is fairly good now.  May not be necessary anymore?**)

	```bash
	mkdir wcs
	cd wcs
	lmi-solve-wcs.sh ../ppp/lmi*fits
	cd ..
	```
	   
1. Create a hierarchical path structure with the data sorted by target and filter.  This does not copy the data, but makes symbolic links to the `ppp/` directory.
   ```bash
   lmi-sort.py ppp
   ```
   
1. Add world coordinate system centered on moving targets ephemeris:

	```bash
	dct-add-moving-wcs.py sorted/TARGET/*/*
	```

1. Add world coordinate system centered on moving targets:

	```bash
	dct-add-moving-wcs.py sorted/TARGET/*/*
	```

1. Copy world coordinate system for frames without solutions with lmi-copy-wcs.  This tool is not routinely used at the moment.

1. Add photometry catalog:

	```bash
	lmi-add-cat.py ppp/lmi*fits
	```

1. Calibrate it:

   ```bash
   mkdir phot
   cd phot
   lmi-calibrate-catalog.py ../ppp/lmi*fits --plot
   ```

   Default is to calibrate g', r', i', z', VR, BC, RC to the RefCat2 photometric catalog (Tonry et al. 2018).  BC, RC, and VR are calibrated to g, r, and r, respectively.  To limit the filters:

   ```
   lmi-calibrate-catalog.py ../ppp/lmi*fits --plot --filter=SDSS-R --filter=VR
   ```

1. Analyze and plot the catalog calibration results.  This step helps find non-photometric data or periods during the night.

   ```bash
   # in the phot/ directory...
   lmi-fit-cat-zps.py ../ppp/
   ```

1. Find photometric standards and calibrate:

   ```bash
   # in the phot/ directory:
   lmi-standard-phot.py ../ppp/lmi*fits -o standard-phot.txt
   lmi-standard-calibration.py standard-phot.txt
   ```

   For more advanced calibration, consider copying and editing the standard-phot.txt file to suit your needs.


## Acknowledgements

Development of this code is supported by NASA Solar System Workings Program award ID 80NSSC21K0164, awarded to the University of Maryland.