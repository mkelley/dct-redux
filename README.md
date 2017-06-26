# dct-redux
Python scripts for reducing DCT data, mostly focused on comet data.

## Status
The scripts for early steps in the reduction, e.g., file name normalization and FITS header fixing, are stable, but the later steps are still under development.

## Requirements
* Python 3
* astropy
* ccdproc
* DS9
* pyds9
* callhorizons
* mskpy

For file name normalization and FITS header fixing, just `astropy` is needed.

## General steps, should apply to any DCT data:
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
## For LMI data:

1. Create bias and flat field frames, apply to data to create partially processed products (PPP files):
   ```bash
   cd ..
   ls raw/*fits > all.list
   # edit bad files from all.list
   emacs all.list
   lmi-rx.py all.list
   ```

1. Create a hierarchical path structure with the data sorted by target and filter.  This does not copy the data, but makes hard links to the `ppp/` directory.
   ```bash
   lmi-sort.py ppp
   ```

1. Center on targets or stars by hand.  Uses script from `mskpy`:
   ```bash
   center-target ppp/*fits
   ```

1. Update WCS, all files.  Targets coordinates are grabbed from SIMBAD or HORIZONS, as needed.
   ```bash
   lmi-update-wcs.py ppp/*.fits
   ```

1. Update WCS for those images that were not aligned on the target, or the target name is not found in SIMBAD:
   ```bash
   cd sorted
   lmi-update-wcs.py --celestial=197.230094,45.588819 c2013a1/*/*fits
   lmi-update-wcs.py --celestial=204.52143223,-8.85685162 giacobini-zinner/SDSS-R/*fits
   lmi-update-wcs.py --celestial=232.66475833,6.02027222 pg*/*/*fits
   ```

1. Add moving target WCS (WCS key "m"):
   ```bash
   for d in c2013a1 c2013x1 c2014s2 c2015v2 churyumov-gerasimenko clark finlay giacobini-zinner linear12 tempel1; do
     (cd $d; dct-add-moving-wcs.py */*fits)
   done
   ```

More later.
