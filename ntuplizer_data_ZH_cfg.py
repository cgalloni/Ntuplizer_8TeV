import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

runOnMC = False

process.load('CMGTools.Common.PAT.PATCMG_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuple')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(5)
)
##process.MessageLogger.cerr.FwkReport.reportEvery = 1000

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
options.maxEvents = -1

options.inputFiles = cms.untracked.vstring(
'dcap://t3se01.psi.ch/pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/pattuple/mc/aniello/JetEDBRTauTau_DATA/598e1284f39e7f866fbeb5af755bb49f/Jet__Run2012A-22Jan2013-v1__AOD_100_2_vFs.root',

)

options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound','SimpleJetCorrectionUncertainty')
                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('flatTuple.root')
                                   )

process.badEventFilter = cms.EDFilter("HLTHighLevel",
				      TriggerResultsTag =
				      cms.InputTag("TriggerResults","","PAT"),
				      HLTPaths =
				      cms.vstring('primaryVertexFilterPath',
						  'noscrapingFilterPath',
						  'hcalLaserEventFilterPath',
						  'HBHENoiseFilterPath',
						  'trackingFailureFilterPath',
						  'CSCTightHaloFilterPath',
						  'eeBadScFilterPath',
						  'EcalDeadCellTriggerPrimitiveFilterPath',
                                                  'tobtecfakesfilterPath' #bad reconstructed track events filter                                  
						  ),
				      eventSetupPathsKey = cms.string(''),
				      andOr = cms.bool(False), ## how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
				      throw = cms.bool(True)   ## throw exception on unknown path names
				      )
                                 

if runOnMC:
  jecLevelsCA8 = [
    'JEC/START53_V23::All_L1FastJet_AK7PFchs.txt',
    'JEC/START53_V23::All_L2Relative_AK7PFchs.txt',
    'JEC/START53_V23::All_L3Absolute_AK7PFchs.txt'
  ]
  jecLevelsAK5chs = [
    'JEC/START53_V23::All_L1FastJet_AK5PFchs.txt',
    'JEC/START53_V23::All_L2Relative_AK5PFchs.txt',
    'JEC/START53_V23::All_L3Absolute_AK5PFchs.txt'
  ]
  jecLevelsAK5 = [
    'JEC/START53_V23::All_L1FastJet_AK5PF.txt',
    'JEC/START53_V23::All_L2Relative_AK5PF.txt',
    'JEC/START53_V23::All_L3Absolute_AK5PF.txt'
  ]
else :
  jecLevelsCA8 = [
    'JEC/FT_53_V21_AN4::All_L1FastJet_AK7PFchs.txt',
    'JEC/FT_53_V21_AN4::All_L2Relative_AK7PFchs.txt',
    'JEC/FT_53_V21_AN4::All_L3Absolute_AK7PFchs.txt',
    'JEC/FT_53_V21_AN4::All_L2L3Residual_AK7PFchs.txt'
  ]
  jecLevelsAK5chs = [
    'JEC/FT_53_V21_AN4::All_L1FastJet_AK5PFchs.txt',
    'JEC/FT_53_V21_AN4::All_L2Relative_AK5PFchs.txt',
    'JEC/FT_53_V21_AN4::All_L3Absolute_AK5PFchs.txt',
    'JEC/FT_53_V21_AN4::All_L2L3Residual_AK5PFchs.txt'
  ]
  jecLevelsAK5 = [
    'JEC/FT_53_V21_AN4::All_L1FastJet_AK5PF.txt',
    'JEC/FT_53_V21_AN4::All_L2Relative_AK5PF.txt',
    'JEC/FT_53_V21_AN4::All_L3Absolute_AK5PF.txt',
    'JEC/FT_53_V21_AN4::All_L2L3Residual_AK5PF.txt'
  ]
      
if runOnMC is True:
  process.myPartons = cms.EDProducer("PartonSelector",
     src = cms.InputTag("genParticles"),
     withLeptons = cms.bool(False)
  )
  
  process.flavourByRef = cms.EDProducer("JetPartonMatcher",
     jets = cms.InputTag("patJetsCA8CHSprunedSubjets"),
     coneSizeToAssociate = cms.double(0.3),
     partons = cms.InputTag("myPartons")
  )

  process.flavourByVal = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("flavourByRef"),
    physicsDefinition = cms.bool(False)
  )

  process.ntuplizer = cms.EDAnalyzer('Ntuplizer',
                                   runOnMC = cms.bool(True),
                                   AK5jetColl = cms.InputTag("patJetsWithVarCHS"),
				   CA8jetColl = cms.InputTag("selectedPatJetsCA8CHSwithQJetsForBoostedTaus"),
				   CA8prunedjetColl = cms.InputTag("selectedPatJetsCA8CHSprunedForBoostedTaus"),
				   subjetsWithFlavor = cms.InputTag("flavourByVal"),
                                   eleColl = cms.InputTag("patElectronsWithTriggerBoosted"),
				   eleColl2 = cms.InputTag("patElectronsWithTrigger"),
				   muonsColl = cms.InputTag( "patMuonsWithTriggerBoosted" ),
				   muonsColl2 = cms.InputTag( "patMuonsWithTrigger" ),
				   tausColl = cms.InputTag("selectedPatTausBoosted"),
				   eleTausColl = cms.InputTag("selectedPatTausEleTau"),
				   muTausColl = cms.InputTag("selectedPatTausMuTau"),
				   METrawColl = cms.InputTag("patMETsRaw"),
				   METColl = cms.InputTag("patMetShiftCorrected"),
				   vtxColl = cms.InputTag("goodOfflinePrimaryVertices"),
				   vtxCollForMet = cms.InputTag("selectedVerticesForMEtCorr"),
				   AK5JetForMet = cms.InputTag("patJetsWithVar"),
				   corrMetPx = cms.string("+0.1166 + 0.0200*Nvtx"),
				   corrMetPy = cms.string("+0.2764 - 0.1280*Nvtx"),
				   rho = cms.InputTag("kt6PFJets", "rho"),
				   bs = cms.InputTag("offlineBeamSpot"),
				   conversions = cms.InputTag("allConversions"),
				   genpColl = cms.InputTag("genParticles"),
				   ak5GenJetsColl = cms.InputTag("ak5GenJets"),
				   ak5GenJetsNoNuColl = cms.InputTag("ak5GenJetsNoNu"),
				   ca8GenJetsColl = cms.InputTag("ca8GenJets"),
				   ca8GenJetsNoNuColl = cms.InputTag("ca8GenJetsNoNu"),
				   HLT = cms.InputTag("TriggerResults","","HLT"),
				   PUInfo = cms.InputTag("addPileupInfo"),
                                   jecCA8PayloadNames = cms.vstring( jecLevelsCA8 ),
                                   jecCA8UncName = cms.string('JEC/START53_V23::All_Uncertainty_AK7PFchs.txt'), 
                                   jecAK5chsPayloadNames = cms.vstring( jecLevelsAK5chs ),
                                   jecAK5chsUncName = cms.string('JEC/START53_V23::All_Uncertainty_AK5PFchs.txt'), 
				   jecAK5PayloadNames = cms.vstring( jecLevelsAK5 ),
				   jecpath = cms.string('/shome/cgalloni/TAU/CMSSW_5_3_13/src/EXOVVNtuplizer/Ntuplizer/'),
                                  )

  process.ca8GenJets = process.ak5GenJets.clone(doAreaFastjet = cms.bool(True),jetAlgorithm = cms.string('CambridgeAachen'),rParam = cms.double(0.8),jetPtMin = cms.double(25))
  process.p = cms.Path(process.badEventFilter*(process.myPartons+process.flavourByRef+process.flavourByVal+process.genParticlesForJets+process.ak5GenJets+process.ca8GenJets+process.ntuplizer))


if runOnMC is False:


  process.ntuplizer = cms.EDAnalyzer('Ntuplizer',
                   		   runOnMC = cms.bool(False),
                   		   AK5jetColl = cms.InputTag("patJetsWithVarCHS"),
		   		   CA8jetColl = cms.InputTag("selectedPatJetsCA8CHSwithQJetsForBoostedTaus"),
		   		   CA8prunedjetColl = cms.InputTag("selectedPatJetsCA8CHSprunedForBoostedTaus"),
                                   subjetsWithFlavor = cms.InputTag(""),
                                   eleColl = cms.InputTag("patElectronsWithTriggerBoosted"),
		   		   eleColl2 = cms.InputTag("patElectronsWithTrigger"),
		   		   muonsColl = cms.InputTag( "patMuonsWithTriggerBoosted" ),
		   		   muonsColl2 = cms.InputTag( "patMuonsWithTrigger" ),
		   		   tausColl = cms.InputTag("selectedPatTausBoosted"),
		   		   eleTausColl = cms.InputTag("selectedPatTausEleTau"),
		   		   muTausColl = cms.InputTag("selectedPatTausMuTau"),
		   		   METrawColl = cms.InputTag("patMETsRaw"),
		   		   METColl = cms.InputTag("patMetShiftCorrected"),
		   		   vtxColl = cms.InputTag("goodOfflinePrimaryVertices"),
				   vtxCollForMet = cms.InputTag("selectedVerticesForMEtCorr"),
				   AK5JetForMet = cms.InputTag("patJetsWithVar"),
				   corrMetPx = cms.string("+0.2661 + 0.3217*Nvtx"),
				   corrMetPy = cms.string("-0.2251 - 0.1747*Nvtx"),				   
		   		   rho = cms.InputTag("kt6PFJets", "rho"),
		   		   bs = cms.InputTag("offlineBeamSpot"),
		   		   conversions = cms.InputTag("allConversions"),
		   		   HLT = cms.InputTag("TriggerResults","","HLT"),
                                   jecCA8PayloadNames = cms.vstring( jecLevelsCA8 ),
                                   jecCA8UncName = cms.string('JEC/FT_53_V21_AN4::All_Uncertainty_AK7PFchs.txt'), 
                                   ##jecCA8UncName = cms.string('JEC/START53_V23::All_Uncertainty_AK7PFchs.txt'),
                                   jecAK5chsPayloadNames = cms.vstring( jecLevelsAK5chs ),
                                  ## jecAK5chsUncName = cms.string('JEC/START53_V23::All_Uncertainty_AK5PFchs.txt'),
                                   jecAK5chsUncName = cms.string('JEC/FT_53_V21_AN4::All_Uncertainty_AK5PFchs.txt'), 
				   jecAK5PayloadNames = cms.vstring( jecLevelsAK5 ),
				   jecpath = cms.string('/shome/cgalloni/TAU/CMSSW_5_3_13/src/EXOVVNtuplizer/Ntuplizer/'),
                                  )
  import FWCore.PythonUtilities.LumiList as LumiList
  import FWCore.ParameterSet.Types as CfgTypes
  process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
  JSONfile = '/shome/cgalloni/TAU/CMSSW_5_3_13/src/EXOVVNtuplizer/Ntuplizer/RunABCD_Cert.json'
  myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
  process.source.lumisToProcess.extend(myLumis)

  process.p = cms.Path(process.badEventFilter*process.ntuplizer)
