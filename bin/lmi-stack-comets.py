#!/usr/bin/env python3
import os
import re
import warnings
import argparse
from glob import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.time import Time
from MontagePy.main import (
    mImgtbl,
    mProjExec,
    mGetHdr,
    mAdd,
    mOverlaps,
    mDiffFitExec,
    mBgModel,
    mBgExec,
)
from mskpy import gcentroid

parser = argparse.ArgumentParser(
    description="Stack comet images together.",
    epilog="Data are stacked by filter",
)
parser.add_argument("file", nargs="+", help="files to stack")
parser.add_argument(
    "--overwrite",
    "-f",
    action="store_true",
    help="overwrite previously saved files",
)
parser.add_argument(
    "--no-bg-match",
    dest="bg_match",
    action="store_false",
    help="disable background matching",
)
parser.add_argument(
    "--bg-subtract",
    action="store_true",
    help="subtract the 'BG' FITS extension, if available",
)
parser.add_argument(
    "--centroid",
    action="store_true",
    help="align on the centroid of the ephemeris position of the target",
)
parser.add_argument(
    "--box", default=21, type=int, help="use this box size for centroiding"
)
parser.add_argument(
    "--niter", default=3, type=int, help="centroiding iterations"
)
parser.add_argument(
    "--shrink",
    action="store_true",
    help="shrink the centroiding box with each iteration",
)
parser.add_argument(
    "--destination",
    "--dest",
    "--dir",
    default="stacks",
    help="files are written to this directory",
)

args = parser.parse_args()

if not os.path.exists(args.destination):
    os.mkdir(args.destination)

dest = os.path.abspath(args.destination)


class Config:
    file_template = "lmi_[0-9]{8}_[0-9]{4}_ppp.fits"
    comet_pat = (
        "([0-9]*[CPDX])(-[A-Z]+)?/(([12][0-9]{3} [A-Z]+[0-9]+) \((.*)\)|(.*))"
    )
    asteroid_pat = "(\(([0-9]+)\) (.*))|([12][0-9]{3} [A-Za-z]+[0-9]+)"


def centroid_offset(hdu):
    gx = hdu.header["CRPIX1M"]
    gy = hdu.header["CRPIX2M"]
    print(gx, gy)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        cy, cx = gcentroid(
            hdu.data,
            (gy, gx),
            box=args.box,
            niter=args.niter,
            shrink=args.shrink,
        )
    print(cx, cy)
    return cy - gy, cx - gx


def object2field(obj):
    if re.match(Config.comet_pat, obj) is not None:
        field = (
            obj.strip()
            .lower()
            .replace(" ", "")
            .replace("(", "_")
            .replace(")", "")
        )
        if field[0].isdigit():
            field = field.replace("/", "_")
        else:
            field = field.replace("/", "")
    elif re.match(Config.asteroid_pat, obj) is not None:
        m = re.findall(Config.asteroid_pat, obj)[0]
        if len(m[0]) == 0:
            field = m[3].replace(" ", "")
        else:
            field = m[2].lower().replace(" ", "")
    else:
        field = obj.lower().replace(" ", "").replace("/", "_")
        m = re.findall("\(\s*[^)]+\)", field)
        if len(m) > 0:
            for s in m:
                field = field.replace(s, "")
    return field


# collect image and filter names
files = {}
for f in sorted(args.file):
    h = fits.getheader(f)
    if h["FILTERS"] not in files:
        files[h["FILTERS"]] = []
    if h.get("FMDSTAT", "HOME") != "HOME":
        print(f"skipping {f}, FMDSTAT != HOME")
        continue
    files[h["FILTERS"]].append(f)

home = os.path.abspath(os.path.curdir)
for filt, ff in files.items():
    os.chdir(home)
    print(filt)
    print(ff)

    if os.path.exists("/tmp/lmi-stack-comets/input/"):
        os.system("rm -f /tmp/lmi-stack-comets/input/*")
    else:
        os.system("mkdir -p /tmp/lmi-stack-comets/input")

    if os.path.exists("/tmp/lmi-stack-comets/output/"):
        os.system("rm -f /tmp/lmi-stack-comets/output/projected/*")
        os.system("rmdir /tmp/lmi-stack-comets/output/projected")
        os.system("rm -f /tmp/lmi-stack-comets/output/corrected/*")
        os.system("rmdir /tmp/lmi-stack-comets/output/corrected")
        os.system("rm -f /tmp/lmi-stack-comets/output/*")
        os.system("rmdir /tmp/lmi-stack-comets/output")

    os.mkdir("/tmp/lmi-stack-comets/output")

    headers = []
    for f in ff:
        cf = "/tmp/lmi-stack-comets/input/" + f.split("/")[-1]
        with fits.open(f, mode="readonly") as hdu:
            mask = hdu["MASK"].data.astype(bool)
            hdu[0].data[mask] = np.nan
            if args.bg_subtract and "BG" in hdu:
                hdu[0].data = hdu[0].data - hdu["BG"].data

            hdu[0].data = hdu[0].data / hdu[0].header["EXPTIME"]
            hdu[0].header["BUNIT"] = "adu/s"

            hdu.writeto(cf, overwrite=True)
            headers.append(hdu[0].header)

    copied = sorted(glob("/tmp/lmi-stack-comets/input/*fits"))
    exptime = 0
    obsnums = []
    os.chdir("/tmp/lmi-stack-comets/output")
    t = []
    for f in copied:
        with fits.open(f, mode="update") as hdu:
            exptime += hdu[0].header["EXPTIME"]
            obsnums.append(hdu[0].header["OBSERNO"])

            if args.centroid:
                dy, dx = centroid_offset(hdu[0])
                print(f, "offset", dx, dy)
            else:
                dy, dx = 0, 0

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FITSFixedWarning)
                w = WCS(hdu[0].header, key="M")

            hdu[0].header.update(w.to_header(key=" ", relax=True))
            hdu[0].header["CRPIX1"] += dx
            hdu[0].header["CRPIX2"] += dy
            t.append(hdu[0].header["DATE"])

        if f == copied[0]:
            mGetHdr(f, "header")

    date = Time(Time(t).mjd.mean(), format="mjd").isot[:16]
    file_date = date.replace("-", "").replace(":", "")

    combines = []
    for combine in ["median", "mean"]:
        if combine == "mean" and len(ff) < 3:
            # median and mean are the same for 1 or 2 images
            continue

        outf = f"{dest}/{object2field(h['object'])}-{file_date}-{filt}-{combine}.fits"
        if os.path.exists(outf) and not args.overwrite:
            print("exists, skipping.")
            continue
        else:
            combines.append(combine)

    if len(combines) == 0:
        continue

    r = mImgtbl("/tmp/lmi-stack-comets/input", "input-images.tbl")
    print("mImgtbl", r)

    os.mkdir("/tmp/lmi-stack-comets/output/projected")
    r = mProjExec(
        "/tmp/lmi-stack-comets/input/",
        "input-images.tbl",
        "header",
        projdir="projected",
    )
    print("mProjExec", r)

    r = mImgtbl("projected", "projected-images.tbl")
    print("mImgtbl", r)

    if args.bg_match:
        r = mOverlaps("projected-images.tbl", "diffs.tbl")
        print(r)

        # Generate difference images and fit them.
        r = mDiffFitExec(
            "projected", "diffs.tbl", "header", "diffs", "diff-fits.tbl"
        )
        print(r)

        # Model the background corrections.
        r = mBgModel("projected-images.tbl", "diff-fits.tbl", "corrections.tbl")
        print(r)

        os.mkdir("/tmp/lmi-stack-comets/output/corrected")
        r = mBgExec(
            "projected",
            "projected-images.tbl",
            "corrections.tbl",
            "corrected",
        )
        print(r)

        r = mImgtbl("corrected", "corrected-images.tbl")
        print(r)

        path = "corrected"
        table = "corrected-images.tbl"
    else:
        path = "projected"
        table = "projected-images.tbl"

    for combine in combines:
        print("  -", combine, end=" ")
        outf = f"{object2field(h['object'])}-{file_date}-{filt}-{combine}.fits"
        coadd = {"median": 1, "mean": 0, "cov": 2}[combine]

        r = mAdd(path, table, "header", outf, coadd=coadd)
        print(r)

        with fits.open(outf, mode="update") as hdu:
            hdu[0].header["EXPTIME"] = exptime
            hdu[0].header["DATEMID"] = date
            hdu[0].header["MOSAIC"] = combine
            hdu[0].header["NIMAGES"] = len(ff), "number of images combined"
            hdu[0].header["OBSERNOS"] = (
                ",".join([str(h["OBSERNO"]) for h in headers]),
                "observation numbers combined",
            )
            hdu[0].header.add_history("Mosaicked: {}".format(obsnums))
            hdu[0].header.add_history(
                "EXPTIME, DATEMID updated for this mosaic."
            )

        os.system(f"cp -f /tmp/lmi-stack-comets/output/{outf} {dest}")
        os.system(
            f"cp -f /tmp/lmi-stack-comets/output/{outf[:-5]}_area.fits {dest}"
        )

        print(outf, "saved")
