#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import copy
import os
import weakref
import lsst.daf.persistence as dafPersist
from lsst.daf.persistence import Policy
from lsst.obs.base import ImageMapping, ExposureMapping, CalibrationMapping, DatasetMapping
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
from lsst.afw.fits import readMetadata
import lsst.log as lsstLog
import lsst.pex.policy as pexPolicy
import lsst.pex.exceptions as pexExcept
from lsst.utils import getPackageDir

__all__ = ["HscSimsMapper"]


class HscSimsMapper(dafPersist.Mapper):

    """CameraMapper is a base class for mappers that handle images from a
    camera and products derived from them.  This provides an abstraction layer
    between the data on disk and the code.
    Public methods: keys, queryMetadata, getDatasetTypes, map,
    canStandardize, standardize
    Mappers for specific data sources (e.g., CFHT Megacam, LSST
    simulations, etc.) should inherit this class.
    The CameraMapper manages datasets within a "root" directory. Note that
    writing to a dataset present in the input root will hide the existing
    dataset but not overwrite it.  See #2160 for design discussion.
    A camera is assumed to consist of one or more rafts, each composed of
    multiple CCDs.  Each CCD is in turn composed of one or more amplifiers
    (amps).  A camera is also assumed to have a camera geometry description
    (CameraGeom object) as a policy file, a filter description (Filter class
    static configuration) as another policy file, and an optional defects
    description directory.
    Information from the camera geometry and defects are inserted into all
    Exposure objects returned.
    The mapper uses one or two registries to retrieve metadata about the
    images.  The first is a registry of all raw exposures.  This must contain
    the time of the observation.  One or more tables (or the equivalent)
    within the registry are used to look up data identifier components that
    are not specified by the user (e.g. filter) and to return results for
    metadata queries.  The second is an optional registry of all calibration
    data.  This should contain validity start and end entries for each
    calibration dataset in the same timescale as the observation time.
    Subclasses will typically set MakeRawVisitInfoClass:
    MakeRawVisitInfoClass: a class variable that points to a subclass of
    MakeRawVisitInfo, a functor that creates an
    lsst.afw.image.VisitInfo from the FITS metadata of a raw image.
    Subclasses must provide the following methods:
    _extractDetectorName(self, dataId): returns the detector name for a CCD
    (e.g., "CFHT 21", "R:1,2 S:3,4") as used in the AFW CameraGeom class given
    a dataset identifier referring to that CCD or a subcomponent of it.
    _computeCcdExposureId(self, dataId): see below
    _computeCoaddExposureId(self, dataId, singleFilter): see below
    Subclasses may also need to override the following methods:
    _transformId(self, dataId): transformation of a data identifier
    from colloquial usage (e.g., "ccdname") to proper/actual usage (e.g., "ccd"),
    including making suitable for path expansion (e.g. removing commas).
    The default implementation does nothing.  Note that this
    method should not modify its input parameter.
    getShortCcdName(self, ccdName): a static method that returns a shortened name
    suitable for use as a filename. The default version converts spaces to underscores.
    _getCcdKeyVal(self, dataId): return a CCD key and value
    by which to look up defects in the defects registry.
    The default value returns ("ccd", detector name)
    _mapActualToPath(self, template, actualId): convert a template path to an
    actual path, using the actual dataset identifier.
    The mapper's behaviors are largely specified by the policy file.
    See the MapperDictionary.paf for descriptions of the available items.
    The 'exposures', 'calibrations', and 'datasets' subpolicies configure
    mappings (see Mappings class).
    Common default mappings for all subclasses can be specified in the
    "policy/{images,exposures,calibrations,datasets}.yaml" files. This provides
    a simple way to add a product to all camera mappers.
    Functions to map (provide a path to the data given a dataset
    identifier dictionary) and standardize (convert data into some standard
    format or type) may be provided in the subclass as "map_{dataset type}"
    and "std_{dataset type}", respectively.
    If non-Exposure datasets cannot be retrieved using standard
    daf_persistence methods alone, a "bypass_{dataset type}" function may be
    provided in the subclass to return the dataset instead of using the
    "datasets" subpolicy.
    Implementations of map_camera and bypass_camera that should typically be
    sufficient are provided in this base class.
    Notes
    -----
    TODO:
    - Handle defects the same was as all other calibration products, using the calibration registry
    - Instead of auto-loading the camera at construction time, load it from the calibration registry
    - Rewrite defects as AFW tables so we don't need pyfits to unpersist them; then remove all mention
      of pyfits from this package.
    """
    packageName = "obs_hsc_sims"

    # a class or subclass of MakeRawVisitInfo, a functor that makes an
    # lsst.afw.image.VisitInfo from the FITS metadata of a raw image
    # MakeRawVisitInfoClass = MakeRawVisitInfo

    # a class or subclass of PupilFactory
    # PupilFactoryClass = afwCameraGeom.PupilFactory

    def __init__(self, root=None, registry=None, calibRoot=None, calibRegistry=None,
                 provided=None, parentRegistry=None, repositoryCfg=None):
        """Initialize the CameraMapper.
        Parameters
        ----------
        policy : daf_persistence.Policy,
            Can also be pexPolicy.Policy, only for backward compatibility.
            Policy with per-camera defaults already merged.
        repositoryDir : string
            Policy repository for the subclassing module (obtained with
            getRepositoryPath() on the per-camera default dictionary).
        root : string, optional
            Path to the root directory for data.
        registry : string, optional
            Path to registry with data's metadata.
        calibRoot : string, optional
            Root directory for calibrations.
        calibRegistry : string, optional
            Path to registry with calibrations' metadata.
        provided : list of string, optional
            Keys provided by the mapper.
        parentRegistry : Registry subclass, optional
            Registry from a parent repository that may be used to look up
            data's metadata.
        repositoryCfg : daf_persistence.RepositoryCfg or None, optional
            The configuration information for the repository this mapper is
            being used with.
        """
        policyFile = Policy.defaultPolicyFile("obs_hsc_sims", "HscSimsMapper.yaml", "policy")
        policy = Policy(policyFile)

        dafPersist.Mapper.__init__(self)

        self.log = lsstLog.Log.getLogger("HscSimsMapper")

        if root:
            self.root = root
        elif repositoryCfg:
            self.root = repositoryCfg.root
        else:
            self.root = None
        if isinstance(policy, pexPolicy.Policy):
            policy = dafPersist.Policy(policy)

        repoPolicy = repositoryCfg.policy if repositoryCfg else None
        if repoPolicy is not None:
            policy.update(repoPolicy)

        # Don't load the default policy from obs_base
        # defaultPolicyFile = dafPersist.Policy.defaultPolicyFile("obs_base",
        #                                                        "MapperDictionary.paf",
        #                                                        "policy")
        # dictPolicy = dafPersist.Policy(defaultPolicyFile)
        # policy.merge(dictPolicy)

        # Levels
        self.levels = dict()
        if 'levels' in policy:
            levelsPolicy = policy['levels']
            for key in levelsPolicy.names(True):
                self.levels[key] = set(levelsPolicy.asArray(key))
        self.defaultLevel = policy['defaultLevel']
        self.defaultSubLevels = dict()
        if 'defaultSubLevels' in policy:
            self.defaultSubLevels = policy['defaultSubLevels']

        # Root directories
        if root is None:
            root = "."
        root = dafPersist.LogicalLocation(root).locString()

        self.rootStorage = dafPersist.Storage.makeFromURI(uri=root)

        # If the calibRoot is passed in, use that. If not and it's indicated in
        # the policy, use that. And otherwise, the calibs are in the regular
        # root.
        # If the location indicated by the calib root does not exist, do not
        # create it.
        calibStorage = None
        if calibRoot is not None:
            calibRoot = dafPersist.Storage.absolutePath(root, calibRoot)
            calibStorage = dafPersist.Storage.makeFromURI(uri=calibRoot,
                                                          create=False)
        else:
            calibRoot = policy.get('calibRoot', None)
            if calibRoot:
                calibStorage = dafPersist.Storage.makeFromURI(uri=calibRoot,
                                                              create=False)
        if calibStorage is None:
            calibStorage = self.rootStorage

        self.root = root

        # Registries
        self.registry = self._setupRegistry("registry", "exposure", registry, policy, "registryPath",
                                            self.rootStorage, searchParents=False,
                                            posixIfNoSql=True)
        if not self.registry:
            self.registry = parentRegistry
        needCalibRegistry = policy.get('needCalibRegistry', None)
        if needCalibRegistry:
            if calibStorage:
                self.calibRegistry = self._setupRegistry("calibRegistry", "calib", calibRegistry, policy,
                                                         "calibRegistryPath", calibStorage,
                                                         posixIfNoSql=False)  # NB never use posix for calibs
            else:
                raise RuntimeError(
                    "'needCalibRegistry' is true in Policy, but was unable to locate a repo at " +
                    "calibRoot ivar:%s or policy['calibRoot']:%s" %
                    (calibRoot, policy.get('calibRoot', None)))
        else:
            self.calibRegistry = None

        # Dict of valid keys and their value types
        self.keyDict = dict()

        self._initMappings(policy, self.rootStorage, calibStorage, provided=None)
        self._initWriteRecipes()

        # Camera geometry
        # #self.cameraDataLocation = None  # path to camera geometry config file
        # #self.camera = self._makeCamera(policy=policy, repositoryDir=repositoryDir)

        # Defect registry and root. Defects are stored with the camera and the registry is loaded from the
        # camera package, which is on the local filesystem.
        # #self.defectRegistry = None
        # #if 'defects' in policy:
        # #    self.defectPath = os.path.join(repositoryDir, policy['defects'])
        # #    defectRegistryLocation = os.path.join(self.defectPath, "defectRegistry.sqlite3")
        # #    self.defectRegistry = dafPersist.Registry.create(defectRegistryLocation)

        # Filter translation table
        self.filters = None

        # verify that the class variable packageName is set before attempting
        # to instantiate an instance
        # #if self.packageName is None:
        # #    raise ValueError('class variable packageName must not be None')

        # #self.makeRawVisitInfo = self.MakeRawVisitInfoClass(log=self.log)

    def _initMappings(self, policy, rootStorage=None, calibStorage=None, provided=None, use_default=True):
        """Initialize mappings
        For each of the dataset types that we want to be able to read, there are
        methods that can be created to support them:
        * map_<dataset> : determine the path for dataset
        * std_<dataset> : standardize the retrieved dataset
        * bypass_<dataset> : retrieve the dataset (bypassing the usual retrieval machinery)
        * query_<dataset> : query the registry
        Besides the dataset types explicitly listed in the policy, we create
        additional, derived datasets for additional conveniences, e.g., reading
        the header of an image, retrieving only the size of a catalog.
        Parameters
        ----------
        policy : `lsst.daf.persistence.Policy`
            Policy with per-camera defaults already merged
        rootStorage : `Storage subclass instance`
            Interface to persisted repository data.
        calibRoot : `Storage subclass instance`
            Interface to persisted calib repository data
        provided : `list` of `str`
            Keys provided by the mapper
        use_default : `bool`
            Load default camera mappings
        """
        # Sub-dictionaries (for exposure/calibration/dataset types)
        imgMappingPolicy = dafPersist.Policy(dafPersist.Policy.defaultPolicyFile(
            "obs_base", "ImageMappingDictionary.paf", "policy"))
        expMappingPolicy = dafPersist.Policy(dafPersist.Policy.defaultPolicyFile(
            "obs_base", "ExposureMappingDictionary.paf", "policy"))
        calMappingPolicy = dafPersist.Policy(dafPersist.Policy.defaultPolicyFile(
            "obs_base", "CalibrationMappingDictionary.paf", "policy"))
        dsMappingPolicy = dafPersist.Policy(dafPersist.Policy.defaultPolicyFile(
            "obs_base", "DatasetMappingDictionary.paf", "policy"))

        # Mappings
        mappingList = (
            ("images", imgMappingPolicy, ImageMapping),
            ("exposures", expMappingPolicy, ExposureMapping),
            ("calibrations", calMappingPolicy, CalibrationMapping),
            ("datasets", dsMappingPolicy, DatasetMapping)
        )
        self.mappings = dict()
        for name, defPolicy, cls in mappingList:
            if name in policy:
                datasets = policy[name]

                # Centrally-defined datasets
                defaultsPath = os.path.join(getPackageDir("obs_base"), "policy", name + ".yaml")
                if os.path.exists(defaultsPath) and use_default:
                    datasets.merge(dafPersist.Policy(defaultsPath))

                mappings = dict()
                setattr(self, name, mappings)
                for datasetType in datasets.names(True):
                    subPolicy = datasets[datasetType]
                    subPolicy.merge(defPolicy)

                    if not hasattr(self, "map_" + datasetType) and 'composite' in subPolicy:
                        def compositeClosure(dataId, write=False, mapper=None, mapping=None,
                                             subPolicy=subPolicy):
                            components = subPolicy.get('composite')
                            assembler = subPolicy['assembler'] if 'assembler' in subPolicy else None
                            disassembler = subPolicy['disassembler'] if 'disassembler' in subPolicy else None
                            python = subPolicy['python']
                            butlerComposite = dafPersist.ButlerComposite(assembler=assembler,
                                                                         disassembler=disassembler,
                                                                         python=python,
                                                                         dataId=dataId,
                                                                         mapper=self)
                            for name, component in components.items():
                                butlerComposite.add(id=name,
                                                    datasetType=component.get('datasetType'),
                                                    setter=component.get('setter', None),
                                                    getter=component.get('getter', None),
                                                    subset=component.get('subset', False),
                                                    inputOnly=component.get('inputOnly', False))
                            return butlerComposite
                        setattr(self, "map_" + datasetType, compositeClosure)
                        # for now at least, don't set up any other handling for this dataset type.
                        continue

                    if name == "calibrations":
                        mapping = cls(datasetType, subPolicy, self.registry, self.calibRegistry, calibStorage,
                                      provided=provided, dataRoot=rootStorage)
                    else:
                        mapping = cls(datasetType, subPolicy, self.registry, rootStorage, provided=provided)
                    self.keyDict.update(mapping.keys())
                    mappings[datasetType] = mapping
                    self.mappings[datasetType] = mapping
                    if not hasattr(self, "map_" + datasetType):
                        def mapClosure(dataId, write=False, mapper=weakref.proxy(self), mapping=mapping):
                            return mapping.map(mapper, dataId, write)
                        setattr(self, "map_" + datasetType, mapClosure)
                    if not hasattr(self, "query_" + datasetType):
                        def queryClosure(format, dataId, mapping=mapping):
                            return mapping.lookup(format, dataId)
                        setattr(self, "query_" + datasetType, queryClosure)
                    if hasattr(mapping, "standardize") and not hasattr(self, "std_" + datasetType):
                        def stdClosure(item, dataId, mapper=weakref.proxy(self), mapping=mapping):
                            return mapping.standardize(mapper, item, dataId)
                        setattr(self, "std_" + datasetType, stdClosure)

                    def setMethods(suffix, mapImpl=None, bypassImpl=None, queryImpl=None):
                        """Set convenience methods on CameraMapper"""
                        mapName = "map_" + datasetType + "_" + suffix
                        bypassName = "bypass_" + datasetType + "_" + suffix
                        queryName = "query_" + datasetType + "_" + suffix
                        if not hasattr(self, mapName):
                            setattr(self, mapName, mapImpl or getattr(self, "map_" + datasetType))
                        if not hasattr(self, bypassName):
                            if bypassImpl is None and hasattr(self, "bypass_" + datasetType):
                                bypassImpl = getattr(self, "bypass_" + datasetType)
                            if bypassImpl is not None:
                                setattr(self, bypassName, bypassImpl)
                        if not hasattr(self, queryName):
                            setattr(self, queryName, queryImpl or getattr(self, "query_" + datasetType))

                    # Filename of dataset
                    setMethods("filename", bypassImpl=lambda datasetType, pythonType, location, dataId:
                               [os.path.join(location.getStorage().root, p) for p in location.getLocations()])
                    # Metadata from FITS file
                    if subPolicy["storage"] == "FitsStorage":  # a FITS image
                        setMethods("md", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                   readMetadata(location.getLocationsWithRoot()[0]))

                        # Add support for configuring FITS compression
                        addName = "add_" + datasetType
                        if not hasattr(self, addName):
                            setattr(self, addName, self.getImageCompressionSettings)

                        if name == "exposures":
                            setMethods("wcs", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwGeom.makeSkyWcs(readMetadata(location.getLocationsWithRoot()[0])))
                            setMethods("calib", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwImage.Calib(readMetadata(location.getLocationsWithRoot()[0])))
                            setMethods("visitInfo",
                                       bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwImage.VisitInfo(readMetadata(location.getLocationsWithRoot()[0])))
                            setMethods("filter",
                                       bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwImage.Filter(readMetadata(location.getLocationsWithRoot()[0])))
                            setMethods("detector",
                                       mapImpl=lambda dataId, write=False:
                                           dafPersist.ButlerLocation(
                                               pythonType="lsst.afw.cameraGeom.CameraConfig",
                                               cppType="Config",
                                               storageName="Internal",
                                               locationList="ignored",
                                               dataId=dataId,
                                               mapper=self,
                                               storage=None,
                                           ),
                                       bypassImpl=lambda datasetType, pythonType, location, dataId:
                                           self.camera[self._extractDetectorName(dataId)]
                                       )
                            setMethods("bbox", bypassImpl=lambda dsType, pyType, location, dataId:
                                       afwImage.bboxFromMetadata(
                                           readMetadata(location.getLocationsWithRoot()[0], hdu=1)))

                        elif name == "images":
                            setMethods("bbox", bypassImpl=lambda dsType, pyType, location, dataId:
                                       afwImage.bboxFromMetadata(
                                           readMetadata(location.getLocationsWithRoot()[0])))

                    if subPolicy["storage"] == "FitsCatalogStorage":  # a FITS catalog
                        setMethods("md", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                   readMetadata(os.path.join(location.getStorage().root,
                                                             location.getLocations()[0]), hdu=1))

                    # Sub-images
                    if subPolicy["storage"] == "FitsStorage":
                        def mapSubClosure(dataId, write=False, mapper=weakref.proxy(self), mapping=mapping):
                            subId = dataId.copy()
                            del subId['bbox']
                            loc = mapping.map(mapper, subId, write)
                            bbox = dataId['bbox']
                            llcX = bbox.getMinX()
                            llcY = bbox.getMinY()
                            width = bbox.getWidth()
                            height = bbox.getHeight()
                            loc.additionalData.set('llcX', llcX)
                            loc.additionalData.set('llcY', llcY)
                            loc.additionalData.set('width', width)
                            loc.additionalData.set('height', height)
                            if 'imageOrigin' in dataId:
                                loc.additionalData.set('imageOrigin',
                                                       dataId['imageOrigin'])
                            return loc

                        def querySubClosure(key, format, dataId, mapping=mapping):
                            subId = dataId.copy()
                            del subId['bbox']
                            return mapping.lookup(format, subId)
                        setMethods("sub", mapImpl=mapSubClosure, queryImpl=querySubClosure)

                    if subPolicy["storage"] == "FitsCatalogStorage":
                        # Length of catalog
                        setMethods("len", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                   readMetadata(os.path.join(location.getStorage().root,
                                                             location.getLocations()[0]),
                                                hdu=1).get("NAXIS2"))

                        # Schema of catalog
                        if not datasetType.endswith("_schema") and datasetType + "_schema" not in datasets:
                            setMethods("schema", bypassImpl=lambda datasetType, pythonType, location, dataId:
                                       afwTable.Schema.readFits(os.path.join(location.getStorage().root,
                                                                             location.getLocations()[0])))

    def _search(self, path):
        """Search for path in the associated repository's storage.
        Parameters
        ----------
        path : string
            Path that describes an object in the repository associated with
            this mapper.
            Path may contain an HDU indicator, e.g. 'foo.fits[1]'. The
            indicator will be stripped when searching and so will match
            filenames without the HDU indicator, e.g. 'foo.fits'. The path
            returned WILL contain the indicator though, e.g. ['foo.fits[1]'].
        Returns
        -------
        string
            The path for this object in the repository. Will return None if the
            object can't be found. If the input argument path contained an HDU
            indicator, the returned path will also contain the HDU indicator.
        """
        return self.rootStorage.search(path)

    def backup(self, datasetType, dataId):
        """Rename any existing object with the given type and dataId.
        The CameraMapper implementation saves objects in a sequence of e.g.:
        - foo.fits
        - foo.fits~1
        - foo.fits~2
        All of the backups will be placed in the output repo, however, and will
        not be removed if they are found elsewhere in the _parent chain.  This
        means that the same file will be stored twice if the previous version was
        found in an input repo.
        """

        # Calling PosixStorage directly is not the long term solution in this
        # function, this is work-in-progress on epic DM-6225. The plan is for
        # parentSearch to be changed to 'search', and search only the storage
        # associated with this mapper. All searching of parents will be handled
        # by traversing the container of repositories in Butler.

        def firstElement(list):
            """Get the first element in the list, or None if that can't be done.
            """
            return list[0] if list is not None and len(list) else None

        n = 0
        newLocation = self.map(datasetType, dataId, write=True)
        newPath = newLocation.getLocations()[0]
        path = dafPersist.PosixStorage.search(self.root, newPath, searchParents=True)
        path = firstElement(path)
        oldPaths = []
        while path is not None:
            n += 1
            oldPaths.append((n, path))
            path = dafPersist.PosixStorage.search(self.root, "%s~%d" % (newPath, n), searchParents=True)
            path = firstElement(path)
        for n, oldPath in reversed(oldPaths):
            self.rootStorage.copyFile(oldPath, "%s~%d" % (newPath, n))

    def keys(self):
        """Return supported keys.
        Returns
        -------
        iterable
            List of keys usable in a dataset identifier
        """
        return iter(self.keyDict.keys())

    def getKeys(self, datasetType, level):
        """Return a dict of supported keys and their value types for a given dataset
        type at a given level of the key hierarchy.
        Parameters
        ----------
        datasetType :  `str`
            Dataset type or None for all dataset types.
        level :  `str` or None
            Level or None for all levels or '' for the default level for the
            camera.
        Returns
        -------
        `dict`
            Keys are strings usable in a dataset identifier, values are their
            value types.
        """

        # not sure if this is how we want to do this. what if None was intended?
        if level == '':
            level = self.getDefaultLevel()

        if datasetType is None:
            keyDict = copy.copy(self.keyDict)
        else:
            keyDict = self.mappings[datasetType].keys()
        if level is not None and level in self.levels:
            keyDict = copy.copy(keyDict)
            for l in self.levels[level]:
                if l in keyDict:
                    del keyDict[l]
        return keyDict

    def getDefaultLevel(self):
        return self.defaultLevel

    def getDefaultSubLevel(self, level):
        if level in self.defaultSubLevels:
            return self.defaultSubLevels[level]
        return None

    @classmethod
    def getPackageName(cls):
        """Return the name of the package containing this CameraMapper."""
        if cls.packageName is None:
            raise ValueError('class variable packageName must not be None')
        return cls.packageName

    @classmethod
    def getCameraName(cls):
        return 'HscSims'

    @classmethod
    def getPackageDir(cls):
        """Return the base directory of this package"""
        return getPackageDir(cls.getPackageName())

    def std_bfKernel(self, item, dataId):
        """Disable standardization for bfKernel
        bfKernel is a calibration product that is numpy array,
        unlike other calibration products that are all images;
        all calibration images are sent through _standardizeExposure
        due to CalibrationMapping, but we don't want that to happen to bfKernel
        """
        return item

###############################################################################
#
# Utility functions
#
###############################################################################

    def _setupRegistry(self, name, description, path, policy, policyKey, storage, searchParents=True,
                       posixIfNoSql=True):
        """Set up a registry (usually SQLite3), trying a number of possible
        paths.
        Parameters
        ----------
        name : string
            Name of registry.
        description: `str`
            Description of registry (for log messages)
        path : string
            Path for registry.
        policy : string
            Policy that contains the registry name, used if path is None.
        policyKey : string
            Key in policy for registry path.
        storage : Storage subclass
            Repository Storage to look in.
        searchParents : bool, optional
            True if the search for a registry should follow any Butler v1
            _parent symlinks.
        posixIfNoSql : bool, optional
            If an sqlite registry is not found, will create a posix registry if
            this is True.
        Returns
        -------
        lsst.daf.persistence.Registry
            Registry object
        """
        if path is None and policyKey in policy:
            path = dafPersist.LogicalLocation(policy[policyKey]).locString()
            if os.path.isabs(path):
                raise RuntimeError("Policy should not indicate an absolute path for registry.")
            if not storage.exists(path):
                newPath = storage.instanceSearch(path)

                newPath = newPath[0] if newPath is not None and len(newPath) else None
                if newPath is None:
                    self.log.warn("Unable to locate registry at policy path (also looked in root): %s",
                                  path)
                path = newPath
            else:
                self.log.warn("Unable to locate registry at policy path: %s", path)
                path = None

        # Old Butler API was to indicate the registry WITH the repo folder, New Butler expects the registry to
        # be in the repo folder. To support Old API, check to see if path starts with root, and if so, strip
        # root from path. Currently only works with PosixStorage
        try:
            root = storage.root
            if path and (path.startswith(root)):
                path = path[len(root + '/'):]
        except AttributeError:
            pass

        # determine if there is an sqlite registry and if not, try the posix registry.
        registry = None

        def search(filename, description):
            """Search for file in storage
            Parameters
            ----------
            filename : `str`
                Filename to search for
            description : `str`
                Description of file, for error message.
            Returns
            -------
            path : `str` or `None`
                Path to file, or None
            """
            result = storage.instanceSearch(filename)
            if result:
                return result[0]
            self.log.debug("Unable to locate %s: %s", description, filename)
            return None

        # Search for a suitable registry database
        if path is None:
            path = search("%s.pgsql" % name, "%s in root" % description)
        if path is None:
            path = search("%s.sqlite3" % name, "%s in root" % description)
        if path is None:
            path = search(os.path.join(".", "%s.sqlite3" % name), "%s in current dir" % description)

        if path is not None:
            if not storage.exists(path):
                newPath = storage.instanceSearch(path)
                newPath = newPath[0] if newPath is not None and len(newPath) else None
                if newPath is not None:
                    path = newPath
            localFileObj = storage.getLocalFile(path)
            self.log.info("Loading %s registry from %s", description, localFileObj.name)
            registry = dafPersist.Registry.create(localFileObj.name)
            localFileObj.close()
        elif not registry and posixIfNoSql:
            try:
                self.log.info("Loading Posix %s registry from %s", description, storage.root)
                registry = dafPersist.PosixRegistry(storage.root)
            except Exception:
                registry = None

        return registry

    def _transformId(self, dataId):
        """Generate a standard ID dict from a camera-specific ID dict.
        Canonical keys include:
        - amp: amplifier name
        - ccd: CCD name (in LSST this is a combination of raft and sensor)
        The default implementation returns a copy of its input.
        Parameters
        ----------
        dataId : `dict`
            Dataset identifier; this must not be modified
        Returns
        -------
        `dict`
            Transformed dataset identifier.
        """

        return dataId.copy()

    def _mapActualToPath(self, template, actualId):
        """Convert a template path to an actual path, using the actual data
        identifier.  This implementation is usually sufficient but can be
        overridden by the subclass.
        Parameters
        ----------
        template : `str`
            Template path
        actualId : `dict`
            Dataset identifier
        Returns
        -------
        `str`
            Pathname
        """

        try:
            transformedId = self._transformId(actualId)
            return template % transformedId
        except Exception as e:
            raise RuntimeError("Failed to format %r with data %r: %s" % (template, transformedId, e))

    def _setFilter(self, mapping, item, dataId):
        """Set the filter object in an Exposure.  If the Exposure had a FILTER
        keyword, this was already processed during load.  But if it didn't,
        use the filter from the registry.
        Parameters
        ----------
        mapping : `lsst.obs.base.Mapping`
            Where to get the filter from.
        item : `lsst.afw.image.Exposure`
            Exposure to set the filter in.
        dataId : `dict`
            Dataset identifier.
        """

        if not (isinstance(item, afwImage.ExposureU) or isinstance(item, afwImage.ExposureI) or
                isinstance(item, afwImage.ExposureF) or isinstance(item, afwImage.ExposureD)):
            return

        if item.getFilter().getId() != afwImage.Filter.UNKNOWN:
            return

        actualId = mapping.need(['filter'], dataId)
        filterName = actualId['filter']
        if self.filters is not None and filterName in self.filters:
            filterName = self.filters[filterName]
        item.setFilter(afwImage.Filter(filterName))

    # Default standardization function for exposures
    def _standardizeExposure(self, mapping, item, dataId, filter=True,
                             trimmed=True, setVisitInfo=True):
        """Default standardization function for images.
        This sets the Detector from the camera geometry
        and optionally set the Fiter. In both cases this saves
        having to persist some data in each exposure (or image).
        Parameters
        ----------
        mapping : `lsst.obs.base.Mapping`
            Where to get the values from.
        item : image-like object
            Can be any of lsst.afw.image.Exposure,
            lsst.afw.image.DecoratedImage, lsst.afw.image.Image
            or lsst.afw.image.MaskedImage
        dataId : `dict`
            Dataset identifier
        filter : `bool`
            Set filter? Ignored if item is already an exposure
        trimmed : `bool`
            Should detector be marked as trimmed?
        setVisitInfo : `bool`
            Should Exposure have its VisitInfo filled out from the metadata?
        Returns
        -------
        `lsst.afw.image.Exposure`
            The standardized Exposure.
        """
        try:
            item = exposureFromImage(item, dataId, mapper=self, logger=self.log, setVisitInfo=setVisitInfo)
        except Exception as e:
            self.log.error("Could not turn item=%r into an exposure: %s" % (repr(item), e))
            raise

        if mapping.level.lower() == "amp":
            self._setAmpDetector(item, dataId, trimmed)
        elif mapping.level.lower() == "ccd":
            self._setCcdDetector(item, dataId, trimmed)

        if filter:
            self._setFilter(mapping, item, dataId)

        return item

    def getRegistry(self):
        """Get the registry used by this mapper.
        Returns
        -------
        Registry or None
            The registry used by this mapper for this mapper's repository.
        """
        return self.registry

    def getImageCompressionSettings(self, datasetType, dataId):
        """Stuff image compression settings into a daf.base.PropertySet
        This goes into the ButlerLocation's "additionalData", which gets
        passed into the boost::persistence framework.
        Parameters
        ----------
        datasetType : `str`
            Type of dataset for which to get the image compression settings.
        dataId : `dict`
            Dataset identifier.
        Returns
        -------
        additionalData : `lsst.daf.base.PropertySet`
            Image compression settings.
        """
        mapping = self.mappings[datasetType]
        recipeName = mapping.recipe
        storageType = mapping.storage
        if storageType not in self._writeRecipes:
            return dafBase.PropertySet()
        if recipeName not in self._writeRecipes[storageType]:
            raise RuntimeError("Unrecognized write recipe for datasetType %s (storage type %s): %s" %
                               (datasetType, storageType, recipeName))
        recipe = self._writeRecipes[storageType][recipeName].deepCopy()
        seed = hash(tuple(dataId.items())) % 2**31
        for plane in ("image", "mask", "variance"):
            if recipe.exists(plane + ".scaling.seed") and recipe.get(plane + ".scaling.seed") == 0:
                recipe.set(plane + ".scaling.seed", seed)
        return recipe

    def _initWriteRecipes(self):
        """Read the recipes for writing files
        These recipes are currently used for configuring FITS compression,
        but they could have wider uses for configuring different flavors
        of the storage types. A recipe is referred to by a symbolic name,
        which has associated settings. These settings are stored as a
        `PropertySet` so they can easily be passed down to the
        boost::persistence framework as the "additionalData" parameter.
        The list of recipes is written in YAML. A default recipe and
        some other convenient recipes are in obs_base/policy/writeRecipes.yaml
        and these may be overridden or supplemented by the individual obs_*
        packages' own policy/writeRecipes.yaml files.
        Recipes are grouped by the storage type. Currently, only the
        ``FitsStorage`` storage type uses recipes, which uses it to
        configure FITS image compression.
        Each ``FitsStorage`` recipe for FITS compression should define
        "image", "mask" and "variance" entries, each of which may contain
        "compression" and "scaling" entries. Defaults will be provided for
        any missing elements under "compression" and "scaling".
        The allowed entries under "compression" are:
        * algorithm (string): compression algorithm to use
        * rows (int): number of rows per tile (0 = entire dimension)
        * columns (int): number of columns per tile (0 = entire dimension)
        * quantizeLevel (float): cfitsio quantization level
        The allowed entries under "scaling" are:
        * algorithm (string): scaling algorithm to use
        * bitpix (int): bits per pixel (0,8,16,32,64,-32,-64)
        * fuzz (bool): fuzz the values when quantising floating-point values?
        * seed (long): seed for random number generator when fuzzing
        * maskPlanes (list of string): mask planes to ignore when doing statistics
        * quantizeLevel: divisor of the standard deviation for STDEV_* scaling
        * quantizePad: number of stdev to allow on the low side (for STDEV_POSITIVE/NEGATIVE)
        * bscale: manually specified BSCALE (for MANUAL scaling)
        * bzero: manually specified BSCALE (for MANUAL scaling)
        A very simple example YAML recipe:
            FitsStorage:
              default:
                image: &default
                  compression:
                    algorithm: GZIP_SHUFFLE
                mask: *default
                variance: *default
        """
        recipesFile = os.path.join(getPackageDir("obs_base"), "policy", "writeRecipes.yaml")
        recipes = dafPersist.Policy(recipesFile)
        supplementsFile = os.path.join(self.getPackageDir(), "policy", "writeRecipes.yaml")
        validationMenu = {'FitsStorage': validateRecipeFitsStorage, }
        if os.path.exists(supplementsFile) and supplementsFile != recipesFile:
            supplements = dafPersist.Policy(supplementsFile)
            # Don't allow overrides, only supplements
            for entry in validationMenu:
                intersection = set(recipes[entry].names()).intersection(set(supplements.names()))
                if intersection:
                    raise RuntimeError("Recipes provided in %s section %s may not override those in %s: %s" %
                                       (supplementsFile, entry, recipesFile, intersection))
            recipes.update(supplements)

        self._writeRecipes = {}
        for storageType in recipes.names(True):
            if "default" not in recipes[storageType]:
                raise RuntimeError("No 'default' recipe defined for storage type %s in %s" %
                                   (storageType, recipesFile))
            self._writeRecipes[storageType] = validationMenu[storageType](recipes[storageType])


def exposureFromImage(image, dataId=None, mapper=None, logger=None, setVisitInfo=True):
    """Generate an Exposure from an image-like object
    If the image is a DecoratedImage then also set its WCS and metadata
    (Image and MaskedImage are missing the necessary metadata
    and Exposure already has those set)
    Parameters
    ----------
    image : Image-like object
        Can be one of lsst.afw.image.DecoratedImage, Image, MaskedImage or
        Exposure.
    Returns
    -------
    `lsst.afw.image.Exposure`
        Exposure containing input image.
    """
    metadata = None
    if isinstance(image, afwImage.MaskedImage):
        exposure = afwImage.makeExposure(image)
    elif isinstance(image, afwImage.DecoratedImage):
        exposure = afwImage.makeExposure(afwImage.makeMaskedImage(image.getImage()))
        metadata = image.getMetadata()
        try:
            wcs = afwGeom.makeSkyWcs(metadata, strip=True)
            exposure.setWcs(wcs)
        except pexExcept.TypeError as e:
            # raised on failure to create a wcs (and possibly others)
            if logger is None:
                logger = lsstLog.Log.getLogger("CameraMapper")
            logger.warn("wcs set to None; insufficient information found in metadata to create a valid wcs: "
                        "%s", e.args[0])

        exposure.setMetadata(metadata)
    elif isinstance(image, afwImage.Exposure):
        # Exposure
        exposure = image
        metadata = exposure.getMetadata()
    else:
        # Image
        exposure = afwImage.makeExposure(afwImage.makeMaskedImage(image))
    #
    # set VisitInfo if we can
    #
    if setVisitInfo and exposure.getInfo().getVisitInfo() is None:
        if metadata is not None:
            if mapper is None:
                if not logger:
                    logger = lsstLog.Log.getLogger("CameraMapper")
                logger.warn("I can only set the VisitInfo if you provide a mapper")
            else:
                exposureId = mapper._computeCcdExposureId(dataId)
                visitInfo = mapper.makeRawVisitInfo(md=metadata, exposureId=exposureId)

                exposure.getInfo().setVisitInfo(visitInfo)

    return exposure


def validateRecipeFitsStorage(recipes):
    """Validate recipes for FitsStorage
    The recipes are supplemented with default values where appropriate.
    TODO: replace this custom validation code with Cerberus (DM-11846)
    Parameters
    ----------
    recipes : `lsst.daf.persistence.Policy`
        FitsStorage recipes to validate.
    Returns
    -------
    validated : `lsst.daf.base.PropertySet`
        Validated FitsStorage recipe.
    Raises
    ------
    `RuntimeError`
        If validation fails.
    """
    # Schemas define what should be there, and the default values (and by the default
    # value, the expected type).
    compressionSchema = {
        "algorithm": "NONE",
        "rows": 1,
        "columns": 0,
        "quantizeLevel": 0.0,
    }
    scalingSchema = {
        "algorithm": "NONE",
        "bitpix": 0,
        "maskPlanes": ["NO_DATA"],
        "seed": 0,
        "quantizeLevel": 4.0,
        "quantizePad": 5.0,
        "fuzz": True,
        "bscale": 1.0,
        "bzero": 0.0,
    }

    def checkUnrecognized(entry, allowed, description):
        """Check to see if the entry contains unrecognised keywords"""
        unrecognized = set(entry.keys()) - set(allowed)
        if unrecognized:
            raise RuntimeError(
                "Unrecognized entries when parsing image compression recipe %s: %s" %
                (description, unrecognized))

    validated = {}
    for name in recipes.names(True):
        checkUnrecognized(recipes[name], ["image", "mask", "variance"], name)
        rr = dafBase.PropertySet()
        validated[name] = rr
        for plane in ("image", "mask", "variance"):
            checkUnrecognized(recipes[name][plane], ["compression", "scaling"],
                              name + "->" + plane)

            for settings, schema in (("compression", compressionSchema),
                                     ("scaling", scalingSchema)):
                prefix = plane + "." + settings
                if settings not in recipes[name][plane]:
                    for key in schema:
                        rr.set(prefix + "." + key, schema[key])
                    continue
                entry = recipes[name][plane][settings]
                checkUnrecognized(entry, schema.keys(), name + "->" + plane + "->" + settings)
                for key in schema:
                    value = type(schema[key])(entry[key]) if key in entry else schema[key]
                    rr.set(prefix + "." + key, value)
    return validated
