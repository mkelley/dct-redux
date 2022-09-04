#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from sbpy.data import Names, Ephem, natural_sort_key

parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="+")
parser.add_argument(
    "--list",
    "-l",
    action="store_true",
    help="files are lists of files to summarize (lines starting with # are ignored)",
)

args = parser.parse_args()


def read_list(fn):
    # lists of files are relative to fn's location
    path = os.path.dirname(fn)
    with open(fn, "r") as inf:
        files = [
            os.path.join(path, line.strip()) for line in inf if line[0] != "#"
        ]
    print(files)
    return files


def get_ephemeris(target):
    objtype = Names.asteroid_or_comet(target)
    opts = dict(epochs=Time(h["DATE-OBS"]), id_type="designation")
    if objtype == "comet":
        parsed = Names.parse_comet(target)
        t = str(parsed.get("number", "")) + parsed["type"]
        if "desig" in parsed:
            t += "/" + parsed["desig"]

        opts.update(dict(no_fragments=True, closest_apparition=True))
    else:
        parsed = Names.parse_asteroid(target)
        t = parsed.get("number", parsed.get("desig"))

    if t is None:
        raise ValueError

    return Ephem.from_horizons(t, **opts)


# summarize files and save all headers
# OBSTYPE OBJECT CCDSUM FILTERS FMDSTAT AORGX_01 AENDX_01 AORGY_01 AENDY_01
if args.list:
    filenames = sum([read_list(fn) for fn in args.files], [])
else:
    filenames = args.files

if len(filenames) == 0:
    print("No files to summarize.", file=sys.stderr)
    exit(1)

headers = []
files = []
for f in sorted(filenames):
    h = fits.getheader(f)
    row = {"file": f}
    row.update(
        {
            k: h.get(k, None)
            for k in [
                "OBSTYPE",
                "OBJECT",
                "CCDSUM",
                "FILTERS",
                "FMDSTAT",
            ]
        }
    )
    row["x0"] = h["AORGX_01"]
    row["x1"] = h["AENDX_01"]
    row["y0"] = h["AORGY_01"]
    row["y1"] = h["AENDY_01"]
    files.append(row)
    headers.append(h)

# sort by object and date
headers = sorted(
    headers, key=lambda h: (h["object"].replace(" ", ""), h["date-obs"])
)
headers.append(None)  # triggers summary of last target

targets = {}
last_k = None
for h in headers:
    if h is None:
        k = None
    else:
        target = h["object"]
        k = target.replace(" ", "")

    if k in targets:
        filters.add(h["FILTERS"])
        airmass.append(h["airmass"])
        mjd.append(Time(h["DATE-OBS"]).mjd)
        continue
    elif last_k not in [k, None] or h is None:
        targets[last_k]["filters"] = " ".join(filters)
        date, time = Time(np.mean(mjd), format="mjd").iso.split()
        targets[last_k]["date"] = date
        targets[last_k]["time"] = time
        targets[last_k]["mjd"] = np.mean(mjd)
        targets[last_k]["airmass"] = np.round(np.mean(airmass), 1)

        if h is None:
            break

    date, time = h["DATE-OBS"].split("T")
    mjd = [Time(h["DATE-OBS"]).mjd]
    airmass = [h["airmass"]]
    filters = set([h["FILTERS"]])
    targets[k] = {"target": target, "ra": h["ra"], "dec": h["dec"]}
    last_k = k

    try:
        eph = get_ephemeris(target)
    except:
        continue

    rdot = eph["rdot"].value.filled(-999)[0]
    rh = eph["rh"].value.filled(-999)[0]
    delta = eph["delta"].value.filled(-999)[0]
    phase = eph["phase"].value.filled(-999)[0]

    targets[k].update(
        {
            "rh": np.round(np.sign(rdot) * rh, 2),
            "delta": np.round(delta, 2),
            "phase": np.round(phase, 1),
        }
    )

rows = sorted(
    list(targets.values()), key=lambda row: natural_sort_key(row["target"])
)
tab = Table(rows)[
    "target",
    "date",
    "time",
    "mjd",
    "ra",
    "dec",
    "airmass",
    "filters",
    "rh",
    "delta",
    "phase",
]

Table(files).write(
    "file-summary.txt", format="ascii.fixed_width_two_line", overwrite=True
)
tab.filled(-999).write(
    "target-summary.txt", format="ascii.fixed_width_two_line", overwrite=True
)
tab.filled(-999).write("target-summary.csv", overwrite=True)
os.system("cat target-summary.csv")
