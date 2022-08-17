#!/usr/bin/env python3
from astropy.modeling.models import Const2D, Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling import Fittable2DModel, Parameter
import sys
import logging
import argparse
import warnings
from datetime import datetime

import numpy as np
import scipy.ndimage as nd
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

import astropy
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.table import Table
from astropy.wcs import FITSFixedWarning
import photutils as pu
import sep


parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="*", help="files to process")
parser.add_argument("--reprocess", action="store_true")
parser.add_argument(
    "--verbose", "-v", action="store_true", help="verbose logging"
)
args = parser.parse_args()

######################################################################


class GaussianConst2D(Fittable2DModel):
    """A model for a 2D Gaussian plus a constant.

    Code from photutils (Copyright (c) 2011, Photutils developers).

    Parameters
    ----------
    constant : float
        Value of the constant.

    amplitude : float
        Amplitude of the Gaussian.

    x_mean : float
        Mean of the Gaussian in x.

    y_mean : float
        Mean of the Gaussian in y.

    x_stddev : float
        Standard deviation of the Gaussian in x. ``x_stddev`` and
        ``y_stddev`` must be specified unless a covariance matrix
        (``cov_matrix``) is input.

    y_stddev : float
        Standard deviation of the Gaussian in y. ``x_stddev`` and
        ``y_stddev`` must be specified unless a covariance matrix
        (``cov_matrix``) is input.

    theta : float, optional
        Rotation angle in radians. The rotation angle increases
        counterclockwise.

    """

    constant = Parameter(default=1)
    amplitude = Parameter(default=1)
    x_mean = Parameter(default=0)
    y_mean = Parameter(default=0)
    x_stddev = Parameter(default=1)
    y_stddev = Parameter(default=1)
    theta = Parameter(default=0)

    @staticmethod
    def evaluate(
        x, y, constant, amplitude, x_mean, y_mean, x_stddev, y_stddev, theta
    ):
        """Two dimensional Gaussian plus constant function."""

        model = Const2D(constant)(x, y) + Gaussian2D(
            amplitude, x_mean, y_mean, x_stddev, y_stddev, theta
        )(x, y)
        return model


def fit_2dgaussian(data):
    """Fit a 2D Gaussian plus a constant to a 2D image.

    Based on code from photutils (Copyright (c) 2011, Photutils developers).

    Parameters
    ----------
    data : array_like
        The 2D array of the image.

    Returns
    -------
    result : A `GaussianConst2D` model instance.
        The best-fitting Gaussian 2D model.

    """

    if np.ma.count(data) < 7:
        raise ValueError(
            "Input data must have a least 7 unmasked values to "
            "fit a 2D Gaussian plus a constant."
        )

    data.fill_value = 0.0
    data = data.filled()

    # Subtract the minimum of the data as a rough background estimate.
    # This will also make the data values positive, preventing issues with
    # the moment estimation in data_properties. Moments from negative data
    # values can yield undefined Gaussian parameters, e.g., x/y_stddev.
    data = data - np.min(data)
    guess_y, guess_x = np.array(data.shape) / 2

    init_amplitude = np.ptp(data)
    g_init = GaussianConst2D(
        constant=0,
        amplitude=init_amplitude,
        x_mean=guess_x,
        y_mean=guess_y,
        x_stddev=3,
        y_stddev=3,
        theta=0,
    )
    fitter = LevMarLSQFitter()
    y, x = np.indices(data.shape)
    with astropy.log.log_to_list() as log_list:
        gfit = fitter(g_init, x, y, data)

    return gfit


######################################################################
# setup logging
logger = logging.Logger("LMI Add Catalog")
logger.setLevel(logging.DEBUG)

# this allows logging to work when lmi-add-cat is run multiple times from
# ipython
if len(logger.handlers) == 0:
    formatter = logging.Formatter("%(levelname)s: %(message)s")
    level = logging.DEBUG if args.verbose else logging.INFO

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(level)
    console.setFormatter(formatter)
    logger.addHandler(console)

    logfile = logging.FileHandler("lmi-add-cat.log")
    logfile.setLevel(level)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

logger.info("#" * 70)
logger.info(datetime.now().isoformat())
logger.info("Command line: " + " ".join(sys.argv[1:]))

######################################################################
# suppress unnecessary warnings
warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", FITSFixedWarning)

######################################################################


def show_objects(im, objects):
    # plot background-subtracted image
    fig, ax = plt.subplots()
    m, s = np.mean(im), np.std(im)
    im = ax.imshow(
        im,
        interpolation="nearest",
        cmap="gray",
        vmin=m - s,
        vmax=m + s,
        origin="lower",
    )

    # plot an ellipse for each object
    for i in range(len(objects)):
        e = Ellipse(
            xy=(objects["x"][i], objects["y"][i]),
            width=6 * objects["a"][i],
            height=6 * objects["b"][i],
            angle=objects["theta"][i] * 180.0 / np.pi,
        )
        e.set_facecolor("none")
        e.set_edgecolor("red")
        ax.add_artist(e)


for f in args.files:
    if fits.getheader(f)["IMAGETYP"] != "OBJECT":
        continue

    logger.debug(f)

    with fits.open(f, mode="update") as hdu:
        im = hdu[0].data + 0
        h = hdu[0].header

        if h["IMAGETYP"].upper() != "OBJECT":
            logger.warning(
                f'Refusing to measure {f} with image type {h["imagetyp"]}.'
            )
            continue

        if "MASK" in hdu:
            mask = hdu["MASK"].data.astype(bool)
        else:
            mask = np.zeros_like(im, bool)

        if "cat" in hdu:
            if not args.reprocess:
                continue
            else:
                del hdu["cat"]

        det = np.zeros_like(mask)
        for iteration in range(3):
            bkg = sep.Background(im, mask=det | mask, bw=64, bh=64, fw=3, fh=3)

            # mask potential sources
            det = ((im - bkg) / bkg.globalrms) > 3
            #   remove isolated pixels
            det = nd.binary_closing(det)

        bkg = sep.Background(im, mask=det | mask, bw=64, bh=64, fw=3, fh=3)

        if "bg" in hdu:
            del hdu["bg"]
        hdu.append(fits.ImageHDU(bkg.back(), name="bg"))
        hdu["bg"].header["bg"] = bkg.globalback
        hdu["bg"].header["rms"] = bkg.globalrms

        data = im - bkg
        data[mask] = 0
        try:
            objects, labels = sep.extract(
                data, 3, err=bkg.globalrms, segmentation_map=True
            )
        except Exception as e:
            logger.error(f"{f}: Object detection failed - {str(e)}")
            continue

        hdu[0].header["ncat"] = len(objects), "number of objects in catalog"
        if len(objects) == 0:
            continue

        # show_objects(data, objects)

        # estimate seeing
        fwhms = []
        segmap = pu.SegmentationImage(labels)
        for i in np.random.choice(len(segmap.segments), 50):
            obj = segmap.segments[i].make_cutout(data, masked_array=True)
            try:
                g = fit_2dgaussian(obj)
            except:
                continue

            fwhm = np.mean((g.x_stddev.value, g.y_stddev.value)) * 2.35
            if fwhm < 1:
                continue
            fwhms.append(fwhm)

        fwhm = sigma_clipped_stats(fwhms)[1]
        rap = fwhm * 2 if np.isfinite(fwhm) else 10

        flux, fluxerr, flag = sep.sum_circle(
            data,
            objects["x"],
            objects["y"],
            rap,
            err=bkg.globalrms,
            gain=h["gain"],
        )

        kronrad, krflag = sep.kron_radius(
            data,
            objects["x"],
            objects["y"],
            objects["a"],
            objects["b"],
            objects["theta"],
            6.0,
        )
        krflux, krfluxerr, _flag = sep.sum_ellipse(
            data,
            objects["x"],
            objects["y"],
            objects["a"],
            objects["b"],
            np.minimum(objects["theta"], np.pi / 2.00001),
            2.5 * kronrad,
            subpix=1,
            err=bkg.globalrms,
            gain=h["gain"],
        )
        krflag |= _flag  # combine flags

        wcs = WCS(h)
        ra, dec = wcs.all_pix2world(objects["x"], objects["y"], 0)

        tab = Table(
            (
                objects["x"],
                objects["y"],
                ra,
                dec,
                flux,
                fluxerr,
                flag,
                objects["a"],
                objects["b"],
                objects["theta"],
                kronrad,
                krflux,
                krfluxerr,
                krflag,
            ),
            names=(
                "x",
                "y",
                "ra",
                "dec",
                "flux",
                "fluxerr",
                "flag",
                "a",
                "b",
                "theta",
                "kronrad",
                "krflux",
                "krfluxerr",
                "krflag",
            ),
        )
        if "cat" in hdu:
            del hdu["cat"]
        hdu.append(fits.BinTableHDU(tab, name="cat"))
        hdu["cat"].header["FWHM"] = fwhm, "estimated median FWHM"
        hdu["cat"].header["RADIUS"] = 2 * fwhm, "aperture photometry radius"

        logger.info(
            f"{f}: {len(tab)} objects, seeing = {fwhm:.1f}, background mean/rms = "
            f"{hdu['bg'].header['bg']:.1f}/{hdu['bg'].header['rms']:.1f}"
        )
