import heliohypy
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes

class MakeTrackletsConnections(lsst.pipe.base.PipelineTaskConnections, dimensions=("instrument", "visitWindow"), defaultTemplates={"timeSpan": "twoWeeks", "detType": "pointSource"}):
    sourceTable = connectionTypes.Input(
        doc="Table of unattributed sources",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "{detType}_hldetCatalog"
        )
    visitTable = connectionTypes.Input(
        doc="visit stats plus observer coordinates",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "{detType}_hlimageCatalog"
        )
    trackletSources = connectionTypes.Output(
        doc="sources that got included in tracklets",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "{detType}_trackletSources"
        )
     tracklets = connectionTypes.Output(
        doc="summary data for tracklets",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "{detType}_tracklets"
        )
     trk2det = connectionTypes.Output(
        doc="indices connecting tracklets to trackletSources",
        dimensions=("instrument", "visitWindow"),
        storageClass="SourceCatalog",
        name = "{detType}_trk2source"
        )

class MakeTrackletsConfig(lsst.pex.config.Config):
  mintrkpts = lsst.pex.config.Field(
      dtype=int,
      default=2,
      doc="minimum number of sources to qualify as a tracklet"
      )
  imagetimetol = lsst.pex.config.Field(
      dtype=float,
      default=getImageTimeTol()
      doc="Tolerance for matching image time, in days: e.g. 1 second"
      )
  maxvel = lsst.pex.config.Field(
      dtype=float,
      default=1.5,
      doc="Default max angular velocity in deg/day."
      )
  minvel = lsst.pex.config.Field(
      dtype=float,
      default=0,
      doc="Min angular velocity in deg/day"
      )
  minarc = lsst.pex.config.Field(
      dtype=float,
      default=0,
      doc="Min total angular arc in arcseconds."
      )
  maxtime = lsst.pex.config.Field(
      dtype=float,
      default=1.5/24,
      doc="Max inter-image time interval, in days."
      )
  mintime = lsst.pex.config.Field(
      dtype=float,
      default=1/86400,
      doc="Minimum inter-image time interval, in days."
      )
  imagerad = lsst.pex.config.Field(
      dtype=float,
      default=2.0,
      doc="radius from image center to most distant corner (deg)"
      )
  maxgcr = lsst.pex.config.Field(
      dtype=float,
      default=0.5,
      doc="Default maximum Great Circle Residual allowed for a valid tracklet (arcsec)"
      )
  forcerun = lsst.pex.config.Field(
      dtype=int,
      default=0,
      doc="Pushes through all but the immediately fatal errors."
      )
  verbose = lsst.pex.config.Field(
      dtype=int,
      default=0,
      doc="Prints monitoring output."
      )

class MakeTrackletsTask(lsst.pipe.base.PipelineTask):
    ConfigClass = MakeTrackletsConfig
    _DefaultName = "makeTracklets"

    def run(self, sourceTable, visitTable):
        """doc string 
           here
        """

        trackout = heliohypy.makeTracklets(self.config, sourceTable, visitTable)
        
        return lsst.pipe.base.Struct(trackletSources=trackout[0],
                                     tracklets=trackout[1],
                                     trk2det=trackout[2]
                                     )
    
