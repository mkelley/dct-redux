#!/usr/bin/env python3
import os
import re
import sys
import logging
import argparse
from glob import glob
from datetime import datetime


class Config:
    file_template = "lmi_[0-9]{8}_[0-9]{4}_ppp.fits"
    comet_pat = (
        "([0-9]*[CPDX])(-[A-Z]+)?/(([12][0-9]{3} [A-Z]+[0-9]+) \((.*)\)|(.*))"
    )
    asteroid_pat = "(\(([0-9]+)\) (.*))|([12][0-9]{3} [A-Za-z]+[0-9]+)"


def file_filter(source):
    """Returns a function that tests the file name.

    Files must match Config.file_template.

    Parameters
    ----------
    source : string
      The source directory.

    """
    pat = os.sep.join([source, Config.file_template])

    def f(fn):
        if re.fullmatch(pat, fn) is None:
            return False
        return True

    return f


def summarize(files):
    from astropy.io import fits
    from astropy.table import Table
    from astropy.coordinates import Angle

    summary = Table(
        names=["file", "obstype", "object", "ra", "dec", "filter"],
        dtype=([(str, 128)] * 3 + [float, float, (str, 128)]),
    )
    for f in files:
        h = fits.getheader(f)
        ra = Angle(h["ra"], unit="hr")
        dec = Angle(h["dec"], unit="deg")
        summary.add_row(
            (
                f,
                h["obstype"],
                h["object"],
                ra.hourangle,
                dec.degree,
                h["filters"],
            )
        )

    return summary


def object2field(obj):
    if re.match(Config.comet_pat, row["object"]) is not None:
        field = (
            row["object"]
            .strip()
            .lower()
            .replace(" ", "")
            .replace("(", "_")
            .replace(")", "")
        )
        if field[0].isdigit():
            field = field.replace("/", "_")
        else:
            field = field.replace("/", "")
    elif re.match(Config.asteroid_pat, row["object"]) is not None:
        m = re.findall(Config.asteroid_pat, row["object"])[0]
        if len(m[0]) == 0:
            field = m[3].replace(" ", "")
        else:
            field = m[2].lower().replace(" ", "")
    else:
        field = row["object"].lower().replace(" ", "").replace("/", "_")
        m = re.findall("\(\s*[^)]+\)", field)
        if len(m) > 0:
            for s in m:
                field = field.replace(s, "")
    return field


######################################################################
parser = argparse.ArgumentParser(
    description="Sort files into field/filter directories."
)
parser.add_argument("source", help="name of the source directory")
parser.add_argument(
    "--target",
    default="sorted",
    help="name of the target directory [default: sorted]",
)
parser.add_argument("-v", action="store_true", help="verbose output")
args = parser.parse_args()

######################################################################
# setup logging
logger = logging.Logger("LMI sort")
logger.setLevel(logging.DEBUG if args.v else logging.INFO)

formatter = logging.Formatter("%(levelname)s: %(message)s")

console = logging.StreamHandler(sys.stdout)
console.setLevel(logging.DEBUG)
console.setFormatter(formatter)
logger.addHandler(console)

logger.info("#" * 70)
logger.info(datetime.now().isoformat())
logger.info("Command line: " + " ".join(sys.argv[1:]))

######################################################################
# parameter and file check
assert os.path.isdir(args.source)
logger.info("Finding FITS files.")
files = glob("{}/*.fits".format(args.source))
files = sorted(filter(file_filter(args.source), files))
summary = summarize(files)

if os.path.exists(args.target):
    assert os.path.isdir(args.target)
else:
    os.mkdir(args.target)
    logger.info("Created target directory: {}".format(args.target))

######################################################################
# sort
i = summary["obstype"] == "OBJECT"
linked = 0
exists = 0
for row in summary[i]:
    field = object2field(row["object"])
    path = os.sep.join([args.target, field])
    if os.path.exists(path):
        assert os.path.isdir(path)
    else:
        os.mkdir(path)
        logger.info("Created directory: {}".format(path))

    path = os.sep.join([path, row["filter"]])
    if os.path.exists(path):
        assert os.path.isdir(path)
    else:
        os.mkdir(path)
        logger.info("Created directory: {}".format(path))

    dest = os.sep.join([path, os.path.split(row["file"])[1]])
    src = row["file"]
    if os.path.exists(dest):
        if not os.path.samefile(dest, src):
            logger.warning(
                "Found file {}, but it is not a link of {}.  Skipping.".format(
                    dest, src
                )
            )
        exists += 1
    else:
        # all files are sorted/object/filter so add a few dots
        os.symlink(os.path.join("..", "..", "..", src), dest)
        logger.debug("Linked: {}".format(dest))
        linked += 1

logger.info(
    "Found {} files, linked {}, {} already existed.".format(
        linked + exists, linked, exists
    )
)
