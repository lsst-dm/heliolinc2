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
  clustrad = lsst.pex.config.Field(
      dtype=float,
      default=1.0e5
      doc="Clustering radius for the DBSCAN algorithm, in km."
      )
  dbscan_npt = lsst.pex.config.Field(
      dtype=int,
      default=3,
      doc="Number of points npt for the DBSCAN algorithm"
      )
  minobsnights = lsst.pex.config.Field(
      dtype=int,
      default=3,
      doc="Minimum number of distinct observing nights for a valid linkage"
      )
  mintimespan = lsst.pex.config.Field(
      dtype=float,
      default=1.0,
      doc="Minimum timespan for a valid linkage, in days"
      )
  mingeodist = lsst.pex.config.Field(
      dtype=float,
      default=0.10;
      doc="Geocentric distance (AU) at the center of the innermost distance bin"
      )
  maxgeodist = lsst.pex.config.Field(
      dtype=float,
      default=100.0
      doc="Minimum value in AU for the center of the outermost distance bin"
      )
  geologstep = lsst.pex.config.Field(
      dtype=float,
      default=1.5,
      doc="Factor by which distance increases from one bin to the next"
      )
  mingeoobs = lsst.pex.config.Field(
      dtype=float,
      default=0.0,
      doc="Minimum inferred geocentric distance for a valid tracklet"
      )
  minimpactpar = lsst.pex.config.Field(
      dtype=float,
      default=0.0,
      doc="Minimum inferred impact parameter (w.r.t Earth) for a valid tracklet"
      )
  use_univar = lsst.pex.config.Field(
      dtype=int,
      default=0,
      doc="Use the universal variable formulation of the Kepler equations, "
          "rather than the default fg function formulation"
      )
  max_v_inf = lsst.pex.config.Field(
      dtype=float,
      default=0.0,
      doc="Maximum value of v_infinity relative to the sun "
      "(must be greater than zero to probe interstellar orbits)"
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
    
