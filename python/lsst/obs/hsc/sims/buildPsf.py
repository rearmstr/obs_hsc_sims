#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 20082014 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import numpy

import lsst.daf.base
import lsst.pex.config
import lsst.pipe.base
import lsst.afw.table
import lsst.afw.image

import lsst.meas.algorithms

class BuildControlPsfConfig(lsst.pex.config.Config):
    warpKernel = lsst.pex.config.Field(
        dtype=str,
        doc="kernel to use when shifting the star image; see afw.math.makeWarpingKernel for choices",
        default="lanczos5"
        )
    size = lsst.pex.config.Field(
        dtype=int,
        doc="size of subimage of the postage PSF stamp",
        default=64
        )

class BuildControlPsfTask(lsst.pipe.base.Task):

    ConfigClass = BuildControlPsfConfig

    _DefaultName = "psf"

    def __init__(self, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)

    def run(self, dataRef, dataType):
        image = dataRef.get(dataType + "starfield_image", immediate=True)
        nx = self.config.size-1
        ny = self.config.size-1
        # extract only the lower-left image
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0,0), lsst.afw.geom.Extent2I(nx, ny))
        image = image[bbox].convertD()
        # shift by half a pixel in both directions, so it's centered on a pixel
        image = lsst.afw.math.offsetImage(image, -0.5, -0.5)
        assert image.getWidth() == nx
        assert image.getHeight() == ny
        kernel = lsst.afw.math.FixedKernel(image)
        return lsst.meas.algorithms.KernelPsf(kernel)

class BuildVariablePsfConfig(lsst.pex.config.Config):
    determiner = lsst.meas.algorithms.psfDeterminerRegistry.makeField(
        doc="PSF determination algorithm", default="pca"
        )

class BuildVariablePsfTask(lsst.pipe.base.Task):

    ConfigClass = BuildVariablePsfConfig

    _DefaultName = "psf"

    def __init__(self, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.makeSubtask("determiner")

    def run(self, dataRef):
        raise NotImplementedError()
