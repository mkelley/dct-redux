#!/bin/bash -e


function filename_normalization()
{
  pushd .
  cd raw
  cat<<EOT


-----------------------------
0. Raw filename normalization
-----------------------------
EOT
  while true
  do
    cat<<EOT

Press
  t to perform test run
  c to run and commit to changes
  s to skip to the next step
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
  cd raw
  cat<<EOT


--------------------
1. FITS header fixes
--------------------

EOT

  if [ -e fixes.txt ]
  then
    echo "fixes.txt already exists"
  else
    cat>fixes.txt<<EOF
first frame, last frame, keyword, value
# 11, 26,	object|objname, 96 Her
EOF
    echo "--> Created raw/fixes.txt"
  fi

  while true
  do
    cat<<EOT

Press
  n to create a new raw/fixes.txt from template
  e to edit fixes.txt with "$EDITOR"
  m to show basic file metadata
  t to test changes in updated file
  c to continue and commit to changes
  s to skip to the next step
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


------------------
2. Input file list
------------------

EOT

  if [ -e all.list ]
  then 
    echo "all.list already exists"
  else
    ls raw/lmi*fits > all.list
    echo "--> Created all.list"
  fi

  while true
  do
    cat<<EOT

Press
  n to create a new all.list (an existing file will be backed up)
  e to edit bad files from all.list with "$EDITOR"
  c to continue
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
    "c") break;;
    "q") exit;;
    esac
  done
}


function summarize_files() {
  cat<<EOT
------------------
3. Summarize files
------------------
EOT

  while true
  do

    echo
    echo -n "target-summary.txt"
    [ -e target-summary.txt ] && echo " exists" || echo " does not exist"
    echo -n "file-summary.txt"
    [ -e file-summary.txt ] && echo " exists" || echo " does not exist"

    cat<<EOT

Press
  c to continue and (re)create the summaries
  s to skip to the next step
  q to quit the script
EOT
    read -rsn1 response
    echo
    case "$response" in
    "c")
      dct-summary.py --list all.list
      echo
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
}

function bias_and_flat() {
  cat<<EOT


------------------------
4. Bias and flat correct
------------------------
EOT

  while true
  do
    echo
    # everything flat corrected?  then continue, else prompt user
    if [ -e ppp ]
    then
      uncorrected=`gethead obstype flatfile ppp/lmi*fits | grep OBJECT | grep -v flat || echo "0"`
      if [ "$uncorrected" == "0" ]
      then
        echo "All files corrected"
      else
        cat<<EOT
Some files not corrected
  - If lmi-rx.py has already run, add missing calibration files to ppp/ before continuing
EOT
      fi
    fi

    cat <<EOT
Press
  l to list files not yet corrected
  c to continue and process files in all.list with lmi-rx.py
  o to reprocess files with imagetyp = object
  a to reprocess all bias, flat, and object files
  s to skip to the next step
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
    "o") lmi-rx.py all.list --reprocess-data;;
    "a") lmi-rx.py all.list --reprocess-all;;
    "s") break;;
    "q") exit;;
    esac
  done
}


function wcs() {
  cat<<EOT


--------------------------
5. World coordinate system
--------------------------
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
  s to skip to the next step
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


----------------
6. Sorting files
----------------
EOT
  lmi-sort.py ppp
}


function moving_wcs() {
  cat<<EOT


-------------
7. Moving WCS
-------------
EOT

  declare -A dirs
  for d in `(cd sorted && ls -dv1 *) | tail -n+2`; do
    # need || echo "" here else this command fails and the script terminates
    test=`expr match "$d" "^\(\([acp][12]\)\|\([1-9]\)\)" || echo ""`
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
  s to skip to the next step
  q to quit the script
EOT

    read -rsn1 response
    case "$response" in
    "t")
      echo -n "Enter space-separated directory names: "
      read -r dlist
      if [ ! -z "$dlist" ]
      then
        for d in $dlist
        do
          if [ ! -z "${dirs[$d]}" ]
          then

            dirs[$d]=${toggle[${dirs[$d]}]}
          else
            echo
            echo "Directory $d not found in sorted/"
          fi
        done
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


-----------------------
8. Photometric catalogs
-----------------------
EOT

  while true
  do
    CATS=`fitsinfo ppp/lmi*.fits|awk '$2 == "CAT" { print; }'|wc -l`
    FILES=`gethead IMAGETYP ppp/lmi*fits|grep OBJECT|wc -l`
    echo
    if [ "$CATS" = "$FILES" ]
    then
      echo "All $FILES object files have catalogs."
    else
      echo "$CATS of $FILES object files have catalogs."
    fi

    cat<<EOT

Press
  c to continue and add catalog FITS extensions, if missing
  r to continue, reprocessing images when catalogs are already defined
  s to skip to the next step
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


---------------------
9. Calibrate catalogs
---------------------
EOT

  FETCH=default
  ZPOPTS=
  cd phot

  while true
  do
    cat<<EOT

Fetch option: $FETCH
Extinction fit options: $ZPOPTS

First calibrate the catalogs, then fit and plot the zeropoints.

Press
  f to set catalog fetch option
  c to calibrate the catalogs
  e to set extinction calculation options
  z to fit the zeropoints as a function of airmass
  s to skip to the next step
  q to quit the script
EOT

    read -rsn1 response
    case "$response" in
    "f")
      echo "Catalog fetch option: as needed (default), for all images (all), or for none of them (none)"
      read -r FETCH
      ;;
    "c")
      # backup and expire old files
      [ -e catalog-cal-summary.txt ] && cp -f --backup=numbered catalog-cal-summary.txt catalog-cal-summary.txt && rm catalog-cal-summary.txt
      [ -e catalog-extinction.txt ] && cp -f --backup=numbered catalog-extinction.txt catalog-extinction.txt && rm catalog-extinction.txt
      lmi-calibrate-catalog.py ../ppp/lmi*fits --plot --fetch=$FETCH
      ;;
    "e")
      echo "Options to pass to lmi-fit-cat-zps.py, e.g., --nrange or --arange:"
      read -r ZPOPTS
      ;;
    "z")
      [ -e catalog-extinction.txt ] && cp -f --backup=numbered catalog-extinction.txt catalog-extinction.txt && echo "catalog-extinction.txt backed up."
      lmi-fit-cat-zps.py ../ppp $ZPOPTS
      [ -e catalog-extinction.txt ] && cat catalog-extinction.txt || echo "Missing catalog extinction results!"
      ;;
    "s") break;;
    "q") exit;;
    esac
  done

  cd ..
}

function _standardphot() {
  [ -e standard-phot.txt ] && cp -f --backup=numbered standard-phot.txt standard-phot.txt && echo "Backed up standard-phot.txt"
  # backup and expire old standard extinction file, if it exists
  [ -e cal-standard-phot.txt ] && cp -f --backup=numbered cal-standard-phot.txt cal-standard-phot.txt && rm cal-standard-phot.txt
  lmi-standard-phot.py --keep-all $*
}

function standard_stars() {
  cat<<EOT


----------------------------
a. Standard star calibration
----------------------------
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
  z to fit the zerpoints
  s to skip to the next step
  q to quit the script
EOT

    opts="-f --dmax=${DMAX}arcsec --hb-dmax=${HBDMAX}arcsec ../ppp/lmi*fits"

    read -rsn1 response
    case "$response" in
    "b")
      echo "Measure stars"
      _standardphot $opts
      echo
      echo "Edit bad photometry from standard-phot.txt as needed"
      ;;
    "k")
      echo "Measure stars"
      _standardphot --keep-all $opts
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
      ;;
    "z")
      [ -e cal-standard-phot.txt ] && cp -f --backup=numbered cal-standard-phot.txt cal-standard-phot.txt
      lmi-standard-calibration.py standard-phot.txt
      [ -e cal-standard-phot.txt ] && cat cal-standard-phot.txt || echo "Missing standard calibration results!"
      echo
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

function comet_stacks() {
  cat<<EOT


---------------
7. Comet stacks
---------------
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
  c to continue and stack the data
  s to skip to the next step
  q to quit the script
EOT

    read -rsn1 response
    case "$response" in
    "t")
      echo -n "Enter space-separated directory names: "
      read -r dlist
      if [ ! -z "$dlist" ]
      then
        for d in $dlist
        do
          if [ ! -z "${dirs[$d]}" ]
          then

            dirs[$d]=${toggle[${dirs[$d]}]}
          else
            echo
            echo "Directory $d not found in sorted/"
          fi
        done
      fi
      echo
      ;;
    "c")
      for d in ${!dirs[@]}
      do
        if [ "${dirs[$d]}" = "+" ]
        then
          lmi-stack-comets.py sorted/$d/*/*
        fi
      done
      break
      ;;
    "s") break;;
    "q") exit;;
    esac
  done
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
  3 to skip to "Summarize files"
  4 to skip to "Bias and flat correct"
  5 to skip to "World coordinate system"
  6 to skip to "Sorting files"
  7 to skip to "Moving WCS"
  8 to skip to "Photometric catalogs"
  9 to skip to "Calibrate catalogs"
  a to skip to "Standard star calibration"
  b to skip to "Comet stacks"
  c to continue
  q to quit the script
EOT
  read -rsn1 first_step
  case "$first_step" in
  "c")
    first_step="0"
    break;;
  [0-9ab]) break;;
  "q") exit;;
  esac
done

mkdir -p wcs phot

# using hexadecimal comparison
[[ "0x$first_step" -le "0x0" ]] && filename_normalization
[[ "0x$first_step" -le "0x1" ]] && header_fixes
[[ "0x$first_step" -le "0x2" ]] && file_list
[[ "0x$first_step" -le "0x3" ]] && summarize_files
[[ "0x$first_step" -le "0x4" ]] && bias_and_flat
[[ "0x$first_step" -le "0x5" ]] && wcs
[[ "0x$first_step" -le "0x6" ]] && sort_files
[[ "0x$first_step" -le "0x7" ]] && moving_wcs
[[ "0x$first_step" -le "0x8" ]] && photometric_catalogs
[[ "0x$first_step" -le "0x9" ]] && calibrate_catalogs
[[ "0x$first_step" -le "0xA" ]] && standard_stars
[[ "0x$first_step" -le "0xB" ]] && comet_stacks
