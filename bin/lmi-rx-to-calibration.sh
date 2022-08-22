#!/bin/bash -e


function filename_normalization()
{
  pushd .
  cd raw
  cat<<EOT


--------------------------
Raw filename normalization
--------------------------
EOT
  while true
  do
    cat<<EOT

Press
  t to perform test run
  c to run and commit to changes
  s to skip this step
  q to quit the script
EOT
    read -rsn1 response
    case "$response" in
    "t") dct-filename-normalization.py -n lmi*fits;;
    "c")
      dct-filename-normalization.py lmi*fits
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
  popd
}


function header_fixes() {
  pushd .
  cat<<EOT


-----------------
FITS header fixes
-----------------
EOT

  if [ ! -e fixes.txt ]
  then 
    cat>fixes.txt<<EOF
first frame, last frame, keyword, value
# 11, 26,	object|objname, 96 Her
EOF
    echo
    echo "--> Created raw/fixes.txt"
  fi

  while true
  do
    cat<<EOT

Press
  n to create a new raw/fixes.txt from template
  e to edit raw/fixes.txt with "$EDITOR"
  m to show basic file metadata
  t to test changes in updated file
  c to continue and commit to changes
  s to skip this step
  q to quit the script
EOT
    read -rsn1 response
    echo
    case "$response" in
    "n")
      cat>fixes.txt<<EOF
first frame, last frame, keyword, value
# 11, 26,	object|objname, 96 Her
EOF
      echo
      echo "--> Created raw/fixes.txt"
      echo
      ;;
    "e")
      $EDITOR fixes.txt
      ;;
    "m")
      gethead IMAGETYP OBJECT RA DEC FILTERS lmi*.fits | less
      echo
      ;;
    "t") dct-fits-header-fixes.py -n fixes.txt;;
    "c")
      dct-fits-header-fixes.py fixes.txt
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
  popd
}


function file_list() {
  cat<<EOT


---------------
Input file list
---------------
EOT

  if [ ! -e all.list ]
  then 
    ls raw/lmi*fits > all.list
    echo
    echo "--> Created all.list"
  fi

  while true
  do
    cat<<EOT

Press
  n to create a new all.list (an existing file will be backed up)
  e to edit bad files from all.list with "$EDITOR"
  c to continue and process all.list with lmi-rx.py
  s to skip this step
  q to quit the script
EOT
    read -rsn1 response
    echo
    case "$response" in
    "n")
      if [ -e all.list ]
      then
        cp -f --backup=numbered all.list all.list
        echo "--> Created backup of all.list"
      fi
      ls raw/lmi*fits > all.list
      echo "--> Created all.list"
      echo
      ;;
    "e")
      $EDITOR all.list
      ;;
    "c") 
      lmi-rx.py all.list
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
}


function bias_and_flat() {
  cat<<EOT


---------------------
Bias and flat correct
---------------------
EOT

  while true
  do
    # everything flat corrected?  then continue, else prompt user
    uncorrected=`gethead obstype flatfile ppp/lmi*fits | grep OBJECT | grep -v flat || echo "0"`
    if [ $uncorrected == "0" ]
    then
      echo "All files corrected"
      break
    fi

    cat <<EOT

Some files not corrected
  - If lmi-rx.py has already run, add missing calibration files to ppp before continuing

Press
  l to list files not yet corrected
  c to continue and process files with lmi-rx.py
  s to skip this step
  q to quit the script
EOT
    read -rsn1 response
    echo
    case "$response" in
    "l")
      gethead -a FLATFILE OBSTYPE OBJECT CCDSUM FILTERS FMDSTAT ppp/lmi*fits | grep OBJECT | awk '$2 == "___" { print; }'
      echo
      ;;
    "c") lmi-rx.py all.list;;
    "s") break;;
    "q") exit;;
    esac
  done
}


function summarize_files() {
  cat<<EOT


---------------
Summarize files
---------------
EOT

  while true
  do

    echo
    echo -n target-summary.txt
    [ -e target-summary.txt ] && echo " exists" || echo " does not exist"
    echo -n "file-summary.txt"
    [ -e file-summary.txt ] && echo " exists" || echo " does not exist"

    cat<<EOT

Press
  c to continue and (re)create the summaries
  s to skip this step
  q to quit the script
EOT
    read -rsn1 response
    echo
    case "$response" in
    "c")
      dct-summary.py ppp/lmi*fits
      echo
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
}


function wcs() {
  cat<<EOT


-----------------------
World coordinate system
-----------------------
EOT

  while true
  do

    unsolved=`gethead -a LMIWCS OBSTYPE ppp/lmi*fits | grep OBJECT | awk '$2 == "___" { print; }'`
    if [ -z "$unsolved" ]
    then
      echo "All files solved"
      break
    fi

    cat<<EOT

Some files unsolved

Press
  l to list files not yet solved
  c to continue and process files with lmi-solve-wcs.sh
  s to skip this step
  q to quit the script
EOT
    read -rsn1 response
    echo
    case "$response" in
    "l")
      gethead -a LMIWCS OBSTYPE OBJECT ppp/lmi*fits | grep OBJECT | awk '$2 == "___" { print; }'
      echo
      ;;
    "c")
      cd wcs
      lmi-solve-wcs.sh ../ppp/lmi*fits
      cd ..
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
}


function sort_files() {
  cat<<EOT


-------------
Sorting files
-------------
EOT
  lmi-sort.py ppp
}


function moving_wcs() {
  cat<<EOT


----------
Moving WCS
----------
EOT

  declare -A dirs
  for d in `(cd sorted && ls -dv1 *) | tail -n+2`; do
    # need || echo "" here else this command fails and the script terminates
    test=`expr match "$d" "^\(\([cp][12]\)\|\([1-9]\)\)" || echo ""`
    if [ -z "$test" ]
    then dirs[$d]="-"
    else dirs[$d]='+'
    fi
  done

  declare -A toggle=( ["+"]="-" ["-"]="+" )

  while true
  do
    for d in ${!dirs[@]}; do
      echo ${dirs[$d]} $d
    done | sort -n -k2

    cat<<EOT

Directories marked with "+" will be processed

Press
  t to toggle directory processing
  c to continue and add moving WCSs
  s to skip this step
  q to quit the script
EOT

    read -rsn1 response
    case "$response" in
    "t")
      echo -n "Enter directory name: "
      read -r d
      if [ ! -z "${dirs[$d]}" ]
      then
        dirs[$d]=${toggle[${dirs[$d]}]}
      else
        echo
        echo "Directory $d not found in sorted/"
      fi
      echo
      ;;
    "c")
      for d in ${!dirs[@]}
      do
        if [ "${dirs[$d]}" = "+" ]
        then
          dct-add-moving-wcs.py sorted/$d/*/*
        fi
      done
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
}

function photometric_catalogs() {
  cat<<EOT


--------------------
Photometric catalogs
--------------------
EOT

  while true
  do
    cat<<EOT

Press
  c to continue and add catalog FITS extensions, if missing
  r to continue, reprocessing images when catalogs are already defined
  s to skip this step
  q to quit the script
EOT

    read -rsn1 response
    case "$response" in
    "c")
      lmi-add-cat.py ppp/lmi*fits
      break
      ;;
    "r")
      lmi-add-cat.py --reprocess ppp/lmi*fits
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
}

function calibrate_catalogs() {
  cat<<EOT


------------------
Calibrate catalogs
------------------
EOT

  FETCH=default

  while true
  do
    cat<<EOT

Fetch: $FETCH

Press
  f to set catalog fetch option
  c to continue
  s to skip this step
  q to quit the script
EOT

    read -rsn1 response
    case "$response" in
    "f")
      echo "Catalog fetch option: as needed (default), for all images (all), or for none of them (none)"
      read -r FETCH
      ;;
    "c")
      cd phot
      [ -e catalog-extinction.txt ] && cp -f --backup=numbered catalog-extinction.txt catalog-extinction.txt
      lmi-calibrate-catalog.py ../ppp/lmi*fits --plot --fetch=$FETCH
      echo
      cat catalog-extinction.txt
      cd ..
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
}

function _standardphot() {
  [ -e standard-phot.txt ] && cp -f --backup=numbered standard-phot.txt standard-phot.txt && echo "Backed up standard-phot.txt"
  lmi-standard-phot.py --keep-all $*
}

function _zeropoints() {
  [ -e cal-standard-phot.txt ] && cp -f --backup=numbered cal-standard-phot.txt cal-standard-phot.txt
  lmi-standard-calibration.py standard-phot.txt
  lmi-plot-cat-zps.py
  cat cal-standard-phot.txt
}

function standard_stars() {
  cat<<EOT


-------------------------
Standard star calibration
-------------------------
EOT

  DMAX=2
  HBDMAX=20
  cd phot

  while true
  do
    cat<<EOT

First measure standard stars, then derive zero-points.

Offset tolerances:
  HB catalog: $HBDMAX"
  Others: $DMAX"

Press
  b to measure the brightest source near the star coordinates
  k to measure all sources within the offset tolerances
  t to set the offset tolerances
  e to manually edit standard-phot.txt (e.g., to remove bad data) with $EDITOR
  s to skip this step
  q to quit the script
EOT

    opts="-f --dmax=${DMAX}arcsec --hb-dmax=${HBDMAX}arcsec ../ppp/lmi*fits"

    read -rsn1 response
    case "$response" in
    "b")
      echo "Measure stars"
      _standardphot $opts
      echo "Derive zero-points"
      _zeropoints
      echo "See plots in phot/"
      echo
      echo "Edit bad photometry from standard-phot.txt as needed"
      ;;
    "k")
      echo "Measure stars"
      _standardphot --keep-all $opts
      echo "Derive zero-points"
      _zeropoints
      echo "See plots in phot/"
      echo
      echo "Edit bad photometry from standard-phot.txt as needed"
      ;;
    "t")
      echo -n "Enter HB catalog tolerance (arcsec): "
      read -r HBDMAX
      echo -n "Enter other catalogs tolerance (arcsec): "
      read -r DMAX
      echo
      ;;
    "e")
      [ -e standard-phot.txt ] && cp -f --backup=numbered standard-phot.txt standard-phot.txt && echo "Backed up standard-phot.txt"
      $EDITOR standard-phot.txt
      echo "Derive zero-points"
      _zeropoints
      echo "See plots in phot/"
      ;;
    "s") break;;
    "q") 
      cd ..
      exit
      ;;
    esac
  done
  cd ..
}

cat<<EOT
lmi-rx-to-calibration.sh - Reduce and calibrate LMI data.

This script is intended to be run in a directory with a single raw/
sub-directory that contains the data to be reduced.
EOT

if [ ! -n $EDITOR ]
then
  echo '\$EDITOR environment variable not set, using "nano"'
  EDITOR=nano
fi

while true
do
  cat<<EOT

Press
  0 to continue with "Filename normalization"
  1 to skip to "FITS header fixes"
  2 to skip to "Input file list"
  3 to skip to "Bias and flat correct"
  4 to skip to "Summarize files"
  5 to skip to "World coordinate system"
  6 to skip to "Sorting files"
  7 to skip to "Moving WCS"
  8 to skip to "Photometric catalogs"
  9 to skip to "Calibrate catalogs"
  a to skip to "Standard star calibration"
  c to continue
  q to quit the script
EOT
  read -rsn1 first_step
  case "$first_step" in
  "c")
    first_step="0"
    break;;
  [0-9a]) break;;
  "q") exit;;
  esac
done

mkdir -p wcs phot

# using hexadecimal comparison
[[ "0x$first_step" -le "0x0" ]] && filename_normalization
[[ "0x$first_step" -le "0x1" ]] && header_fixes
[[ "0x$first_step" -le "0x2" ]] && file_list
[[ "0x$first_step" -le "0x3" ]] && bias_and_flat
[[ "0x$first_step" -le "0x4" ]] && summarize_files
[[ "0x$first_step" -le "0x5" ]] && wcs
[[ "0x$first_step" -le "0x6" ]] && sort_files
[[ "0x$first_step" -le "0x7" ]] && moving_wcs
[[ "0x$first_step" -le "0x8" ]] && photometric_catalogs
[[ "0x$first_step" -le "0x9" ]] && calibrate_catalogs
[[ "0x$first_step" -le "0xA" ]] && standard_stars
