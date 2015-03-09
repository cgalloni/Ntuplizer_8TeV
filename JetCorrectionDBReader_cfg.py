import FWCore.ParameterSet.Config as cms
process = cms.Process("jectxt")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.StandardSequences.Services_cff')

# define your favorite global tag
process.GlobalTag.globaltag = 'FT_53_V21_AN4::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK5PF    = cms.EDAnalyzer('JetCorrectorDBReader',  
      # below is the communication to the database 
      payloadName    = cms.untracked.string('AK7PFchs'),
      # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
      # but it is recommended to use the GT name that you retrieved the files from.
      globalTag      = cms.untracked.string('FT_53_V21_AN4::All'),  
      printScreen    = cms.untracked.bool(False),
      createTextFile = cms.untracked.bool(True)
)
process.p = cms.Path(process.readAK5PF )
