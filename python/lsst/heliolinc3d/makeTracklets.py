import heliohypy
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes
import solarsyst_dyn_geo as solardg
import pandas as pd
import numpy as np

class MakeTrackletsConnections(lsst.pipe.base.PipelineTaskConnections,
                               dimensions=["instrument"]):
    diaSourceTable = connectionTypes.Input(
        doc="Table of unattributed sources",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name = "sspDiaSourceInputs"
    )
    visitTable = connectionTypes.Input(
        doc="visit stats plus observer coordinates",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name = "sspVisitInputs"
    )
    trackletSources = connectionTypes.Output(
        doc="sources that got included in tracklets",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name = "sspTrackletSources"
    )
    tracklets = connectionTypes.Output(
        doc="summary data for tracklets",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name = "sspTracklets"
    )
    trk2det = connectionTypes.Output(
        doc="indices connecting tracklets to trackletSources",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name = "sspTrackletToSource"
    )

def getImageTimeTol():
  return 2./(24.*3600.)

class MakeTrackletsConfig(lsst.pipe.base.PipelineTaskConfig, pipelineConnections=MakeTrackletsConnections):
  mintrkpts = lsst.pex.config.Field(
      dtype=int,
      default=2,
      doc="minimum number of sources to qualify as a tracklet"
      )
  imagetimetol = lsst.pex.config.Field(
      dtype=float,
      default=getImageTimeTol(),
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
  timespan = lsst.pex.config.Field(
      dtype=float,
      default=14.0,
      doc="Default time to look back before most recent data (days)"
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


def df_to_numpy(df, dtype):
    sa = np.zeros(len(df), dtype=dtype)
    for col in df.columns:
        sa[col] = df[col]
    return sa

def det_to_numpy(df):
    return df_to_numpy(df, dtype=solardg.hldet)

def vis_to_numpy(df):
    return df_to_numpy(df, dtype=solardg.hlimage)

class MakeTrackletsTask(lsst.pipe.base.PipelineTask):
    ConfigClass = MakeTrackletsConfig
    _DefaultName = "makeTracklets"

    def run(self, diaSourceTable, visitTable):
        """doc string 
           here
        """

        # copy all config parameters from the Task's config object
        # to heliolinc's native config object.
        config = heliohypy.MakeTrackletsConfig()
        allvars = [ item for item in dir(heliohypy.MakeTrackletsConfig) if not item.startswith("__") ]
        for var in allvars:
            setattr(config, var, getattr(self.config, var))

        # convert dataframes to numpy array with dtypes that heliolinc expects
        diaSourceTable = det_to_numpy(diaSourceTable)
        visitTable     = vis_to_numpy(visitTable)

        (dets, tracklets, trac2det) = heliohypy.makeTracklets(config, diaSourceTable, visitTable)
        (dets, tracklets, trac2det) = map(pd.DataFrame, [dets, tracklets, trac2det])

        #Do something about trailed sources
        return lsst.pipe.base.Struct(trackletSources=dets,
                                     tracklets=tracklets,
                                     trk2det=trac2det
                                     )
    
