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
import lsst.pex.config as pexConfig
import lsst.pipe.base
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.deblender import SourceDeblendTask
import lsst.meas.extensions.shapeHSM
import lsst.meas.modelfit
import lsst.meas.extensions.photometryKron
from lsst.afw.geom import makeSkyWcs, makeCdMatrix


try:
    import scipy.spatial
    spatialAvailable = True
except ImportError:
    spatialAvailable = False

from .buildPsf import BuildControlPsfTask
bfdAvailable = True
try:
    import lsst.desc.old.bfd
except ImportError:
    bfdAvailable = False


class ProcessImageConfig(pexConfig.Config):
    psf = pexConfig.ConfigurableField(
        target=BuildControlPsfTask,
        doc="PSF determination"
    )
    measurement = pexConfig.ConfigurableField(target=SingleFrameMeasurementTask, doc="Source measurement")
    varianceBorderWidth = pexConfig.Field(
        dtype=int,
        default=1,
        doc=("Use a border of this many pixels around each postage stamp (combined over the full image)"
             " to compute the variance")
    )
    dataType = pexConfig.Field(
        dtype=str,
        default='',
        doc=("default data type")
    )
    maxObjects = pexConfig.Field(
        dtype=int,
        default=None,
        optional=True,
        doc=("If this is not None, then only process this many objects.")
    )
    numPerRow = pexConfig.Field(
        dtype=int,
        default=100,
        optional=True,
        doc=("How many objects per row")
    )
    doDetection = pexConfig.Field(
        dtype=bool,
        default=False,
        doc=("run source detection on the image")
    )
    detection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="detection algorithm",
    )
    deblend = pexConfig.ConfigurableField(
        target=SourceDeblendTask,
        doc="deblend algorithm",
    )
    addGridPoints = pexConfig.Field(
        dtype=bool,
        default=False,
        doc=("add footprints at the grid points if they are not in the detection list")
    )
    addGridDist = pexConfig.Field(
        dtype=float,
        default=6.,
        doc=("add objects grid points that do not have a source within this distance")
    )
    writeGridDist = pexConfig.Field(
        dtype=bool,
        default=False,
        doc=("Write out distance to closest grid point")
    )
    correlatedNoiseFile = pexConfig.Field(
        doc="File that contains list of k values and power spectrum to use with correlated noise.  Still"
        " needs to be scaled by noise variance.  For BFD only",
        dtype=str,
        default='',
    )
    doCorrelatedNoise = pexConfig.Field(
        doc="Run with correlated noise",
        dtype=bool,
        default=False,
    )
    recomputeVariance = pexConfig.Field(
        doc="recompute Variance after deblending",
        dtype=bool,
        default=False,
    )
    addTruthCols = pexConfig.ListField(
        dtype=str,
        default=[],
        optional=True,
        doc="List of columns to store in truth table"
    )
    truthPrefix = pexConfig.Field(
        dtype=str,
        default='truth',
        optional=True,
        doc="Prefix for truth values"
    )
    matchTruth = pexConfig.Field(
        doc="match the detected objects to truth",
        dtype=bool,
        default=False,
    )
    matchTol = pexConfig.Field(
        doc="match tolerance to truth if no size information available",
        dtype=float,
        default=3,
    )
    PIXEL_SCALE = pexConfig.Field(
        doc="pixel scale in arcsec/pixel",
        dtype=float,
        default=0.168,
    )

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)
        # self.measurement.slots.centroid = "centroid.sdss"
        # self.measurement.slots.instFlux = None
        # self.measurement.slots.modelFlux = "cmodel.flux"
        # self.measurement.slots.calibFlux = None
        # self.measurement.slots.apFlux = "flux.sinc"
        # self.measurement.plugins.names = ["shape.sdss", "flux.sinc", "flux.psf", "shape.hsm.regauss",
        #                                   "flux.kron", "cmodel", "classification.extendedness"]
        #                               ]
        # self.measurement.plugins.names |= lsst.meas.extensions.multiShapelet.algorithms
        # self.measurement.plugins['classification.extendedness'].fluxRatio = 0.985
        if bfdAvailable:
            self.measurement.plugins.names |= ["bfdKMoment"]
            self.measurement.plugins['bfdKMoment'].sigma = 2
            self.measurement.plugins['bfdKMoment'].useRecVariance = False
            self.measurement.plugins['bfdKMoment'].useTableVariance = True
            self.measurement.plugins['bfdKMoment'].shift = True
            self.measurement.plugins['bfdKMoment'].wIndex = 4
            self.measurement.plugins['bfdKMoment'].useNoisePs = False

        self.detection.thresholdValue = 5
        self.detection.reEstimateBackground = False


class ProcessImageTask(lsst.pipe.base.CmdLineTask):

    ConfigClass = ProcessImageConfig
    _DefaultName = "processImage"

    def __init__(self, **kwargs):
        super(ProcessImageTask, self).__init__(**kwargs)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = lsst.daf.base.PropertyList()
        self.makeSubtask("psf")
        self.makeSubtask("measurement", schema=self.schema, algMetadata=self.algMetadata)

        self.bfdAvailable = bfdAvailable

        if self.config.doDetection:
            self.makeSubtask("detection", schema=self.schema)
            self.makeSubtask("deblend", schema=self.schema)
            if self.config.writeGridDist:
                self.indexKey = self.schema.addField("index", type=int, doc="grid index of galaxy")
                self.gridDistKey = self.schema.addField("gridDist", type=float,
                                                        doc="distance to nearest grid point")

        self.truthKeys = {}
        for key in self.config.addTruthCols:
            self.truthKeys[key] = self.schema.addField(self.config.truthPrefix + "." + key, type=float,
                                                       doc=key+"from truth catalog")

        if self.config.matchTruth:
            self.multTruthMatchKey = self.schema.addField("multTruthMatch", type='Flag',
                                                          doc="multiple truth matches")
            self.noTruthMatchKey = self.schema.addField("noTruthMatch", type='Flag', doc="no truth matches")

        self.gridPositions = []
        self.gridBoxes = []

    def computeVariance(self, image):
        array = image.getArray()
        n = array.shape[0] // self.config.numPerRow
        assert n * self.config.numPerRow == array.shape[0]
        mask = numpy.zeros(array.shape, dtype=bool)
        for i in range(self.config.varianceBorderWidth):
            mask[i::n, :] = True
            mask[:, i::n] = True
            mask[n-i::n, :] = True
            mask[:, n-i::n] = True
        borderPixels = array[mask]
        return numpy.std(borderPixels, dtype=numpy.float64)**2

    def buildExposure(self, dataRef):

        image = dataRef.get(self.config.dataType + "image", immediate=True)
        exposure = afwImage.ExposureF(image.getBBox(afwImage.PARENT))
        exposure.getMaskedImage().getImage().getArray()[:, :] = image.getArray()

        variance = self.computeVariance(image)
        exposure.getMaskedImage().getVariance().set(variance)
        self.algMetadata.set('noise_variance', variance)

        exposure.setPsf(self.psf.run(dataRef, self.config.dataType))

        pixelScale = self.config.PIXEL_SCALE*afwGeom.arcseconds
        cdMatrix = makeCdMatrix(scale=pixelScale)
        crpix = afwGeom.Point2D(0, 0)
        crval = lsst.afw.coord.IcrsCoord(0.*afwGeom.degrees, 0.*afwGeom.degrees)
        wcs = makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cdMatrix)
        exposure.setWcs(wcs)

        calib = afwImage.Calib()
        calib.setFluxMag0(1e12)
        exposure.setCalib(calib)

        return exposure

    def buildSourceCatalog(self, imageBBox, dataRef):
        """Build an empty source catalog, using the provided sim catalog's position to generate
        square Footprints and its ID to set that of the source catalog.
        """
        sourceCat = afwTable.SourceCatalog(self.schema)
        sourceCat.getTable().setMetadata(self.algMetadata)
        simCat = dataRef.get(self.config.dataType + "epoch_catalog", immediate=True)
        xKey = simCat.schema.find('x').key
        yKey = simCat.schema.find('y').key
        idKey = simCat.schema.find('num').key
        self.inputTruthKeys = {}
        for key in self.config.addTruthCols:
            if key == 'star':
                continue
            self.inputTruthKeys[key] = simCat.schema.find(key).key

        n = imageBBox.getWidth() / self.config.numPerRow
        assert n * self.config.numPerRow == imageBBox.getWidth()
        dims = afwGeom.Extent2I(n, n)
        offset = afwGeom.Extent2I(int(numpy.min(simCat[xKey])), int(numpy.min(simCat[yKey])))

        if self.config.matchTruth:
            self.truePos = [(simRecord.get(xKey), simRecord.get(yKey)) for simRecord in simCat]
            self.truthTree = scipy.spatial.KDTree(self.truePos)
            self.simCat = simCat

        max_objects = len(simCat)
        if self.config.maxObjects is not None:
            max_objects = self.config.maxObjects

        for simRecord in simCat[:max_objects]:
            sourceRecord = sourceCat.addNew()
            sourceRecord.setId(simRecord.get(idKey))
            for key in self.config.addTruthCols:
                if key == 'star':
                    if simRecord.get('type') == 'star':
                        sourceRecord.set(self.truthKeys[key], 1)
                    else:
                        sourceRecord.set(self.truthKeys[key], 0)
                    continue
                sourceRecord.set(self.truthKeys[key], simRecord.get(self.inputTruthKeys[key]))
            position = afwGeom.Point2I(int(simRecord.get(xKey)), int(simRecord.get(yKey)))
            bbox = afwGeom.Box2I(position - offset, dims)
            mask = afwImage.Mask(bbox)
            mask.set(100)
            spans = afwGeom.SpanSet.fromMask(mask, 100)
            footprint = afwDet.Footprint(spans)
            # add dummy value of 100 for peak value
            footprint.addPeak(position.getX(), position.getY(), 100)
            self.gridPositions.append((position.getX(), position.getY()))
            self.gridBoxes.append(bbox)
            sourceRecord.setFootprint(footprint)

        if self.config.doCorrelatedNoise and bfdAvailable:
            data = numpy.genfromtxt(self.config.correlatedNoiseFile)
            sourceRecord.getTable().getMetadata().set('kData', data[:, 0])
            sourceRecord.getTable().getMetadata().set('psData', data[:, 1])

        return sourceCat

    def run(self, dataRef):

        exposure = self.buildExposure(dataRef)
        sourceCat = self.buildSourceCatalog(exposure.getBBox(afwImage.PARENT), dataRef)

        if self.config.doDetection:
            table = afwTable.SourceTable.make(self.schema)
            table.setMetadata(self.algMetadata)

            detections = self.detection.makeSourceCatalog(table, exposure)
            sourceCat = detections.sources
            if self.config.maxObjects is not None:
                sourceCat = sourceCat[:self.config.maxObjects]

            if self.config.addGridPoints:
                if spatialAvailable:
                    lastId = sourceCat[-1].getId() + 1

                    peakPosList = [(peak.getIx(), peak.getIy()) for src in sourceCat
                                   for peak in src.getFootprint().getPeaks()]

                    gridPos = numpy.array(self.gridPositions)
                    peakPos = numpy.array(peakPosList)
                    tree = scipy.spatial.cKDTree(peakPos)

                    for (x, y), box in zip(self.gridPositions, self.gridBoxes):
                        dist, index = tree.query([x, y])
                        gridDetection = False
                        if dist < self.config.addGridDist:
                            gridDetection = True

                        if gridDetection is False:
                            sourceRecord = sourceCat.addNew()
                            sourceRecord.setId(lastId+1)
                            lastId += 1
                            footprint = lsst.afw.detection.Footprint(box, exposure.getBBox())
                            footprint.addPeak(x, y, 100)
                            sourceRecord.setFootprint(footprint)
                else:
                    self.log.warn('Cannot add grid points requires scipy.spatial')

            if self.config.recomputeVariance:
                mask = exposure.getMaskedImage().getMask()
                mask &= (mask.getPlaneBitMask("DETECTED"))
                iso_mask = (mask.getArray() != (mask.getPlaneBitMask("DETECTED")))
                vals = exposure.getMaskedImage().getImage().getArray()[iso_mask]
                sigma_mad = 1.4826*numpy.median(numpy.abs(vals - numpy.median(vals)))
                good_mask = numpy.abs(vals < 10*sigma_mad)
                variance = numpy.var(vals[good_mask], dtype=numpy.float64)
                self.log.info('Computed variance after detection: %f' % variance)
                exposure.getMaskedImage().getVariance().set(variance)
                self.algMetadata.set('noise_variance', variance)

            self.deblend.run(exposure, sourceCat, exposure.getPsf())

        self.measurement.run(sourceCat, exposure)

        if self.config.doDetection:
            # Remove parents from catalog, only keeping the children
            mask = numpy.array([a.get('deblend_nChild') == 0 for a in sourceCat])
            sourceCat = sourceCat.subset(mask)

            if self.config.writeGridDist and spatialAvailable:
                gridPos = numpy.array(self.gridPositions)
                tree = scipy.spatial.cKDTree(gridPos)

                srcPosList = numpy.array([(src.getX(), src.getY()) for src in sourceCat])
                minDist, index = tree.query(srcPosList)

                # assume a square image
                stampSize = exposure.getMaskedImage().getImage().getWidth() / self.config.numPerRow
                for dist, src in zip(minDist, sourceCat):

                    xIndex = src.getX() // stampSize
                    yIndex = src.getY() // stampSize
                    index = xIndex*self.config.numPerRow + yIndex
                    src.set(self.gridDistKey, dist)
                    src.set(self.indexKey, int(index))
            else:
                if not spatialAvailable:
                    self.log.warn("Not computing distances or indexes because can't find scipy.spatial")

        if self.config.matchTruth:

            measPos = [(src.getX(), src.getY()) for src in sourceCat]
            measTree = scipy.spatial.KDTree(measPos)
            radii = []
            for ii, src in enumerate(sourceCat):
                x = src.getX()
                y = src.getY()

                radius = None
                if src.get('shape.sdss.flags') is False:
                    radius = src.get('shape.sdss').getTraceRadius()
                elif (src.get('bfd.flags') is False) or (radius > 10):
                    radius = 1./(numpy.sqrt(0.5*(src['bfd.moments'][3])/(src['bfd.moments'][0])))
                elif radius is None or radius > 10:
                    radius = src.get('shape.sdss.psf').getTraceRadius()
                else:
                    radius = src.get('shape.sdss.psf').getTraceRadius()

                radii.append(radius)
                # count how many truth objects overlap this source close to center
                indexes = self.truthTree.query_ball_point([x, y], self.config.matchTol*radius)
                if len(indexes) == 0:
                    src.set(self.noTruthMatchKey, True)
                elif len(indexes) == 1:
                    simRecord = self.simCat[indexes[0]]
                    for key in self.config.addTruthCols:
                        if key == 'star':
                            obj_type = simRecord.get('type')
                            val = 0
                            if obj_type == 'star':
                                val = 1
                            src.set(self.truthKeys[key], val)
                            continue
                        src.set(self.truthKeys[key], simRecord.get(self.inputTruthKeys[key]))
                else:

                    # For each potential truth see if this src is the closest one to it
                    use_indexes = []
                    for index in indexes:
                        xx, yy = (self.simCat['x'][index], self.simCat['y'][index])
                        minDist, mindex = measTree.query([xx, yy])
                        # This is the closest match and within a radius
                        if mindex == ii and minDist < radius*self.config.matchTol:
                            print('  closest', minDist, radius)
                            use_indexes.append(index)
                        else:
                            print('  match another')
                    if len(use_indexes) == 0:
                        src.set(self.noTruthMatchKey, True)
                        continue

                    if len(use_indexes) > 1:
                        src.set(self.multTruthMatchKey, True)

                    weights = self.simCat['flux'][use_indexes]
                    weights /= weights.sum()
                    for key in self.config.addTruthCols:
                        # for num just take the most luminous
                        if key == 'num':
                            max_index = use_indexes[numpy.argmax(self.simCat['flux'][use_indexes])]
                            src.set(self.truthKeys[key], self.simCat[max_index].get('num'))
                            continue
                        if key == 'star':
                            max_index = use_indexes[numpy.argmax(self.simCat['flux'][use_indexes])]
                            obj_type = self.simCat[max_index].get('type')
                            val = 0
                            if obj_type == 'star':
                                val = 1
                            src.set(self.truthKeys[key], val)
                            continue
                        if key == 'flux':
                            src.set(self.truthKeys[key], numpy.sum(self.simCat['flux'][use_indexes]))
                            continue
                        if key == 'g1':
                            src.set(self.truthKeys[key], numpy.sum(self.simCat['g1'][use_indexes]))
                            continue
                        if key == 'g2':
                            src.set(self.truthKeys[key], numpy.sum(self.simCat['g2'][use_indexes]))
                            continue
                        val = 0
                        for index, w in zip(use_indexes, weights):
                            simRecord = self.simCat[index]
                            val += w*simRecord.get(self.inputTruthKeys[key])

                        src.set(self.truthKeys[key], val)

        dataRef.put(sourceCat, self.config.dataType + "src")

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="image", level="image",
                               help="data ID, e.g. --id subfield=0")
        return parser

    def writeConfig(self, butler, clobber=False, doBackup=False):
        pass

    def writeSchemas(self, butler, clobber=False, doBackup=False):
        pass

    def writeMetadata(self, dataRef):
        pass

    def writeEupsVersions(self, butler, clobber=False, doBackup=False):
        pass
