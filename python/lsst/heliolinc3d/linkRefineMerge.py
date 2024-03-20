import heliohypy
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes

class LinkRefineMergeConnections(lsst.pipe.base.PipelineTaskConnections,
                           dimensions=("instrument", "visitWindow"),
                           defaultTemplates={"timeSpan": "twoWeeks", "hypothesis": "mainBelt", "orbitType": "bound"}):
    visitTable = connectionTypes.Input(
        doc="visit stats plus observer coordinates and source indices",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "hlimageCatalog_augmented"
    )
    trackletSources = connectionTypes.Input(
        doc="sources that got included in tracklets",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "trackletSources"
     )
     summary = connectionTypes.Input(
        doc="one line summary of each linkage",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
         multiple=True,
        name = "{hypothesis}_hlclust_{orbitType}",
     )
     clust2det = connectionTypes.Input(
        doc="indices connecting linkages (clusters) to trackletSources",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
         multiple=True,
        name = "{hypothesis}_clust2det_{orbitType}",
     )
     summaryMerged = connectionTypes.Output(
        doc="one line summary of each linkage",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "{hypothesis}_hlclustMerged_{orbitType}",
     )
     clust2detMerged = connectionTypes.Output(
        doc="indices connecting linkages (clusters) to trackletSources",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "{hypothesis}_clust2detMerged_{orbitType}",
     )

class LinkRefineConfig(lsst.pex.config.Config):
  MJDref = lsst.pex.config.Field(
      dtype=float,
      default=0.0,
      doc="MJD or reference time. No sensible default is possible."
      )
  simptype = lsst.pex.config.Field(
      dtype=int,
      default=0,
      doc="Defines how simplex is constructed in the 2-D parameter space "
      "of geocentric distance at first detection (geodist1) and the last "
      "detection (geodist2). simptype=0 uses multiplicative scaling to "
      "create an approximately equilateral triangle. simptype=1 creates a "
      "simplex elongated along the direction defined by geodist1=geodist2. "
      "simptype=2 uses subtraction to create a precisely equilateral triangle."
      )
  ptpow = lsst.pex.config.Field(
      dtype=int,
      default=1,
      doc="Power to which we raise the number of unique detections, "
      "when calculating the cluster quality metric."
      )
  nightpow = lsst.pex.config.Field(
      dtype=int,
      default=1,
      doc="Power to which we raise the number of distinct observing nights, "
      "when calculating the cluster quality metric."
      )
  timepow = lsst.pex.config.Field(
      dtype=int,
      default=0,
      doc="Power to which we raise the total temporal span of the linkage, "
      "when calculating the cluster quality metric."
      )
  rmspow = lsst.pex.config.Field(
      dtype=int,
      default=2,
      doc="Power to which we raise the RMS astrometric residual "
      "when calculating the cluster quality metric."
      )
  maxrms = lsst.pex.config.Field(
      dtype=float,
      default=1.0e5
      doc="Maximum scaled RMS in km for a viable cluster, in km."
      )
  verbose = lsst.pex.config.Field(
      dtype=int,
      default=0,
      doc="Prints monitoring output."
      )

class LinkRefineMergeTask(lsst.pipe.base.PipelineTask):
    ConfigClass = LinkRefineConfig
    _DefaultName = "linkRefine"

    def run(self, visitTable, trackletSources, summary, clust2det):
        #Summary and clust2det will be lists of data references, rather than the data itself.
        """doc string 
           here
        """

        linkout = heliohypy.linkRefineHerget(self.config, visitTable, trackletSources, summary, clust2det)
        
        return lsst.pipe.base.Struct(summaryMerged=linkout[0],
                                     clust2detMerged=linkout[1]
                                     )
    
