#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
from lsst.daf.persistence import Butler
import desc_old_bfd as bfd

import glob
import os
import numpy
import pygmmis

from .processImage import ProcessImageConfig, ProcessImageTask

class Runner(pipeBase.TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        parentDir = parsedCmd.input
        while os.path.exists(os.path.join(parentDir, "_parent")):
            parentDir = os.path.realpath(os.path.join(parentDir, "_parent"))
        butler2 = Butler(root=os.path.join(parentDir, "rerun", parsedCmd.prior_rerun), calibRoot=parsedCmd.calib)
        idParser = parsedCmd.prior.__class__(parsedCmd.prior.level)
        idParser.idList = parsedCmd.prior.idList
        idParser.datasetType = parsedCmd.prior.datasetType
        butler = parsedCmd.butler
        parsedCmd.butler = butler2
        idParser.makeDataRefList(parsedCmd)
        parsedCmd.butler = butler
        return [ (parsedCmd.id.refList,  dict(priorRefList=idParser.refList, **kwargs))]

class ProcessBfdModelPriorConfig(ProcessImageConfig):
    invariantCovariance = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="will all the galaxies have the same covariance"
    )
    numProc = pexConfig.Field(
        dtype=int,
        default=1,
        optional=True,
        doc="number of processors to use"
    )
    sampleSeed = pexConfig.Field(
        dtype=int,
        default=0,
        optional=True,
        doc="random seed for selecting prior galaxies"
    )
    sampleFraction = pexConfig.Field(
        dtype=float,
        default=0.25,
        optional=True,
        doc="only add this fraction of galaxies from the prior files"
    )
    covFile = pexConfig.Field(
        dtype=str,
        default='./test.fits',
        optional=True,
        doc="file that contains the covariance matrices to know which to choose"
    )
    setStar = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="set if the object was a star"
    )
    rerunLabel = pexConfig.Field(
        dtype=str,
        default='b0',
        optional=False,
        doc="label for the prior"
    )
    generateFake = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="generate fake data with "
    )

class ProcessBfdModelPriorTask(ProcessImageTask):

    ConfigClass = ProcessBfdModelPriorConfig
    _DefaultName = "processBfdModelPrior"
    RunnerClass = Runner


    def __init__(self,  **kwargs):
        """Initialize the measurement task, including the modelfits catalog schema,
        the model, prior, and calib objects, and the fitter subtask.
        """

        ProcessImageTask.__init__(self, **kwargs)
        #self.schema =  afwTable.Schema()
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        # Should move these into C++?
        self.flagMomKey = self.schema.addField("bfd.flags.moment", doc="flag bad input on moments",
                                               type='Flag')
        self.notSelFluxKey = self.schema.addField("bfd.ns.flux", doc="not selected because of flux",
                                                  type='Flag')
        self.notSelVarKey = self.schema.addField("bfd.ns.var", doc="not selected because of variance",
                                                 type='Flag')
        self.varBinKey = self.schema.addField("bfd.var.bin", doc="variance bin",type=int)
        self.pqrKey = bfd.BfdPqrKey.addFields(self.schema, "bfd")
        self.momKey = self.schema.addField("moments",doc="moments",type="ArrayF",size=6)
        self.probKey = self.schema.addField("logL", doc="redshift used in pqr", type=float)

        self.log.info('Init')

    def getCovariances(self):
        cat = afwTable.BaseCatalog.readFits(self.config.covFile)
        covList = []
        for rec in cat:
            cov = rec.get('isoCov')
            label = rec.get('label')
            minVariance = rec.get('min')
            maxVariance = rec.get('max')
            cov = numpy.array(cov.reshape(6,6),dtype=numpy.float32)

            covList.append((cov,label,minVariance,maxVariance))

        return covList


    def loadPrior(self, priorRefList):
        self.prior = bfd.MomentPrior()
        self.log.info('Load prior')
        first = True
        for ii,priorRef in enumerate(priorRefList):
            self.log.info("Adding prior %s" % priorRef.dataId)
            try:
                cat = priorRef.get('prior', immediate=True)
                self.prior.addCatalog(cat, self.config.invariantCovariance,
                                      self.config.sampleFraction, self.config.sampleSeed)
                # Should be same for all prior catalogs
                if first:
                    self.cov = numpy.array(cat.getTable().getMetadata().getArrayDouble('COV')).reshape(6,6)
                    first=False
            except  Exception as e:
                print('Failed to read',e)
                continue
            if first is False:
                break

        self.fluxMin = self.prior.getFluxMin()
        self.fluxMax = self.prior.getFluxMax()
        self.varMin = self.prior.getVarMin()
        self.varMax = self.prior.getVarMax()

        bins = numpy.arange(0.05,1.25,0.05)
        self.varBin = numpy.digitize([(self.varMax + self.varMin)/2.], bins)[0]-1
        priorLabel = self.label#priorRefList[0].dataId['label']
        self.log.info("Creating mixture model")
        self.gmm = pygmmis.GMM()
        file_name = 'gmm_files/gmm_%s_%s.npz' % (self.config.rerunLabel, priorLabel)
        self.gmm.load(file_name)
        self.log.info("File loaded")

    def prepCatalog(self, inputs):
        """Prepare the prior and return the output catalog
        """
        outCat =  afwTable.SourceCatalog(self.schema)
        srcCat = inputs

        for srcRecord in srcCat:
            outRecord = outCat.addNew()
            #outRecord.setId(srcCat.get('id'))

        return outCat


    def run(self, dataRefList, priorRefList):
        """Main driver
        """
        self.log.info('running')
        self.label = priorRefList[0].dataId['label']
        self.log.info('ready to load')
        self.loadPrior(priorRefList)
        for dataRef in dataRefList:
            try:
                self.log.info("processing %s" % str(dataRef.dataId))
                src = dataRef.get('src', immediate=True)
                cov = src[src['bfd.flags']==False][0]['bfd.momentsCov'][0]

                if ( (cov > self.varMax or cov < self.varMin) and (self.varMax > 0 and self.varMin > 0) ):
                    print('Flux variance not within bounds, skipping %f,%f,%f'%(cov,self.varMin,self.varMax))
                    continue

                outCat = self.prepCatalog(src)
                self.runMeasure(src, outCat)
                self.log.info("Writing outputs")
                dataRef.dataId['label'] = self.label
                dataRef.put(outCat,'pqr')
            except Exception as e:
                 self.log.warn("Could not process %s" % e)
                 continue

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_argument("--prior_rerun", required=True, help="rerun for prior data")
        parser.add_id_argument(name="--id", datasetType="image", level="image",
                               help="data ID, e.g. --id subfield=0")
        parser.add_id_argument(name="--prior", datasetType="prior", level="image",
                               help="data ID, e.g. --id subfield=0 label=b1")
        return parser

    def runMeasure(self, sources, outCat):

        flags = sources.get('bfd.flags')
        flux = sources.get('bfd.moments')[:,0]
        noise = sources.get('bfd.momentsCov')[:,0]
        pqrKey = self.schema.find('bfd.pqr').key

        # Preslection cuts
        fail_pre_sel = flags == True
        [rec.set(self.flagMomKey,True) for rec in outCat[fail_pre_sel]]

        pre_sel = numpy.logical_not(fail_pre_sel)
        pre_frac =  1-1.*numpy.sum(fail_pre_sel)/len(sources)
        self.log.info('Fraction passing preselection %g' % pre_frac)

        # Flux selection
        flux_sel = numpy.logical_and.reduce((pre_sel,
                                             flux > self.fluxMin,
                                             flux < self.fluxMax
        ))

        not_flux_sel = numpy.logical_and.reduce((pre_sel,
                                                 numpy.logical_or(flux < self.fluxMin, flux > self.fluxMax)
        ))

        total = numpy.sum(flux_sel) + numpy.sum(not_flux_sel)

        if total == 0:
            sel_frac = 0
        else:
            sel_frac =  1.*numpy.sum(flux_sel)/(numpy.sum(flux_sel) + numpy.sum(not_flux_sel))


        [rec.set(self.varBinKey, self.varBin) for rec in outCat]
        [rec.set(self.momKey, src.get('bfd.moments')) for rec,src in zip(outCat,sources)]
        self.log.info('Remaining fraction passing flux selection %g / %g (total)' % (sel_frac,numpy.sum(flux_sel)/(1.*len(sources))))

        # Set flag key to false initially, a failure will be set in c++
        if self.config.generateFake is False:
            prob = self.gmm.logL(sources['bfd.moments'][flux_sel])
        else:
            prob = self.gmm.draw(int(numpy.sum(flux_sel)))
        [rec.set(self.probKey, float(p)) for rec,p in zip(outCat[flux_sel], prob)]



