#!/bin/bash -e

cat<<EOT
lmi-rx-to-calibration.sh - Reduce and calibrate LMI data.

This script is intended to be run in a directory with a single raw/
sub-directory that contains the data to be reduced.

EOT
while true;
do
  echo "Press"
  echo "  c to continue"
  echo "  e to exit script"
  read -rsn1 response
  case "$response" in
  "c") break;;
  "e") exit;;
  esac
done

cd raw
cat<<EOT

--------------------------
Raw filename normalization
--------------------------
EOT
while true;
do
  cat<<EOT
Press
  c to run and commit to changes
  s to skip this step
  e to exit script
  otherwise perform test run
EOT
  read -rsn1 response
  case "$response" in
  "c")
    dct-filename-normalization.py lmi*fits
    break
    ;;
  "s") break;;
  "e") exit;;
  *) dct-filename-normalization.py -n lmi*fits;;
  esac
done

cat<<EOT

-----------------
FITS header fixes
-----------------
EOT
while true;
do
  cat<<EOT
Edit raw/fixes.txt
  c to commit to changes
  t to create fixes.txt from template
  m to show basic file metadata
  s to skip this step
  e to exit script
  otherwise test changes in updated file
EOT
  read -rsn1 response
  echo
  case "$response" in
  "c")
    dct-fits-header-fixes.py fixes.txt
    break
    ;;
  "t") cat>fixes.txt<<EOF
first frame, last frame, keyword, value
# 11, 26,	object|objname, 96 Her
EOF
  echo
  echo "--> Created raw/fixes.txt"
  echo
  ;;
  "m")
    gethead IMAGETYP OBJECT RA DEC FILTERS lmi*.fits | less
    echo
    ;;
  "s") break;;
  "e") exit;;
  *) dct-fits-header-fixes.py -n fixes.txt;;
  esac
done

cat<<EOT

---------------------
Bias and flat correct
---------------------
EOT
cd ..
ls raw/lmi*fits > all.list
read -p "Edit bad files from all.list, then press enter to continue."
lmi-rx.py all.list

while true;
do
  # everything flat corrected?  then continue, else prompt user
  gethead obstype flatfile ppp/lmi*fits | grep OBJECT | grep -v flat || break
  echo <<EOT
$n files not corrected
  c to continue anyway
  e to exit script
  otherwise add missing calibration files to ppp and press any other key
EOT
  read -rsn1 response
  echo
  case "$response" in
  "c") break;;
  "e") exit;;
  *) lmi-rx.py all.list;;
  esac
done

dct-summary.py ppp/lmi*fits

mkdir wcs
cd wcs
lmi-solve-wcs.sh ../ppp/lmi*fits
cd ..

lmi-sort.py ppp
for t in sorted/{[cp][12]*,[1-9]*p_*}; do
  dct-add-moving-wcs.py $t/*/*
done


lmi-add-cat.py ppp/lmi*fits

mkdir phot
cd phot
lmi-calibrate-catalog.py ../ppp/lmi*fits --plot
lmi-standard-phot.py ../ppp/lmi*fits

## edit out bad photometry (something wrong with images... saturation?)
lmi-standard-calibration.py standard-phot.txt
lmi-plot-cat-zps.py
cd ..
