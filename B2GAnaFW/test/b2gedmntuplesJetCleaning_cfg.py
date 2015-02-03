### *****************************************************************************************
### Usage:
###
### cmsRun b2gedmntuples_cfg.py maxEvts=N 
###
### Default values for the options are set:
### maxEvts     = -1
### *****************************************************************************************
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')

options.register('maxEvts',
                 100,# default value: process all events
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.int,
                 'Number of events to process')

options.register('sample',
                 #'file:/afs/cern.ch/work/d/decosa/public/DMtt/miniAOD_Phys14.root',
                 #/TprimeJetToTH_allHdecays_M1200GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
                 '/store/mc/Phys14DR/TprimeJetToTH_allHdecays_M1200GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/94117DA2-009A-E411-9DFB-002590494CB2.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')

options.register('lheLabel',
                 'source',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'LHE module label')

options.register('outputLabel',
                 'B2GEDMNtuple.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Output label')

options.register('globalTag',
                 'PHYS14_25_V1',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Global Tag')

options.register('isData',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Is data?')

options.register('LHE',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Keep LHEProducts')

options.parseArguments()

if(options.isData):options.LHE = False

    
###inputTag labels
muLabel  = 'slimmedMuons'
elLabel  = 'slimmedElectrons'
jLabel = 'slimmedJets'
jLabelAK8 = 'slimmedJetsAK8'
ak8jetEILabel = 'selectedPatJetsAK8PFCHSEI' 
ak8jetLabel = 'selectedPatJetsAK8PFCHS' 
ak8subjetEILabel = 'selectedPatJetsAK8PFCHSEIPrunedSubjets' 
ak8subjetLabel = 'selectedPatJetsAK8PFCHSPrunedSubjets' 
pvLabel  = 'offlineSlimmedPrimaryVertices'
convLabel = 'reducedEgamma:reducedConversions'
particleFlowLabel = 'packedPFCandidates'    
metLabel = 'slimmedMETs'

triggerResultsLabel = "TriggerResults"
triggerSummaryLabel = "hltTriggerSummaryAOD"
hltMuonFilterLabel       = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoRhoFiltered0p15"
hltPathLabel             = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"
hltElectronFilterLabel  = "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8"

FileNames = [
    #/ZPrimeToTTJets_M3000GeV_W30GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM 
    #'/store/mc/Phys14DR/ZPrimeToTTJets_M3000GeV_W30GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/746BDDDC-3568-E411-BAA1-3417EBE2F0DF.root', 
    #'/store/mc/Phys14DR/ZPrimeToTTJets_M3000GeV_W30GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C691EB8B-0568-E411-AB94-00A0D1EEDDA8.root', 
    #'/store/mc/Phys14DR/ZPrimeToTTJets_M3000GeV_W30GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/14B4987B-2F68-E411-9FAD-00266CFAE740.root', 
    #'/store/mc/Phys14DR/ZPrimeToTTJets_M3000GeV_W30GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/24E57485-F467-E411-A9C6-F04DA275101A.root', 
    #'/store/mc/Phys14DR/ZPrimeToTTJets_M3000GeV_W30GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/4AB1C47F-DB67-E411-B18F-7845C4FC3A4C.root', 
    #'/store/mc/Phys14DR/ZPrimeToTTJets_M3000GeV_W30GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/D83DC464-EB67-E411-904F-3417EBE2F0DF.root', 
    #/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
    '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root'
    ]

process = cms.Process("b2gEDMNtuples")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('HLTrigReport')
### Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
### Number of maximum events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts) )
### Source file
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        #options.sample
        FileNames
        )
)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = options.globalTag 

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
#process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:startup_GRun')
#process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
#process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
#for pset in process.GlobalTag.toGet.value():
#    pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
##   Fix for multi-run processing:
#process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
#process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )
    
###
### Select leptons for jet cleaning
###
process.selectedMuons = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("slimmedMuons"), 
    cut = cms.string(
    #'''abs(eta)<2.5 && pt>10. &&
    #(pfIsolationR04().sumChargedHadronPt+
    #max(0.,pfIsolationR04().sumNeutralHadronEt+
    #pfIsolationR04().sumPhotonEt-
    #0.50*pfIsolationR04().sumPUPt))/pt < 0.20 && 
    #(isPFMuon && (isGlobalMuon || isTrackerMuon) )'''
    '''abs(eta)<2.5 && pt>10. &&
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''
    )
    )

process.selectedElectrons = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("slimmedElectrons"), 
    cut = cms.string(
    #'''abs(eta)<2.5 && pt>20. &&
    #gsfTrack.isAvailable() &&
    #gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2 &&
    #(pfIsolationVariables().sumChargedHadronPt+
    #max(0.,pfIsolationVariables().sumNeutralHadronEt+
    #pfIsolationVariables().sumPhotonEt-
    #0.5*pfIsolationVariables().sumPUPt))/pt < 0.15'''
    '''abs(eta)<2.5 && pt>20. &&
    gsfTrack.isAvailable() &&
    gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2''' 
    )
    )

### Do projections
process.pfCHS = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("packedPFCandidates"), 
    cut = cms.string("fromPV")
    )

# then remove the previously selected muons
process.pfNoMuonCHS =  cms.EDProducer("CandPtrProjector", 
    src = cms.InputTag("pfCHS"), 
    veto = cms.InputTag("selectedMuons")
    )
# then remove the previously selected electrons
process.pfNoElectronsCHS = cms.EDProducer("CandPtrProjector", 
    src = cms.InputTag("pfNoMuonCHS"), 
    veto =  cms.InputTag("selectedElectrons")
    )


## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("packedGenParticles"), 
    cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
    )

## Fat GenJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak8GenJetsNoNu = ak4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu")
    )

## Pruned fat GenJets (two jet collections are produced, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ak8GenJetsNoNuPruned = ak4GenJets.clone(
   SubJetParameters,
   rParam = cms.double(0.8),
   src = cms.InputTag("packedGenParticlesForJetsNoNu"),
   usePruning = cms.bool(True),
   writeCompound = cms.bool(True),
   jetCollInstanceName=cms.string("SubJets")
   )

## AK8 PF jets with EI 
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak8PFJetsCHSEI = ak4PFJets.clone(
   rParam = cms.double(0.8),
   src = cms.InputTag("pfNoElectronsCHS"),
   doAreaFastjet = cms.bool(True),
   jetPtMin = cms.double(50.)
   )

## Pruned AK8 PF jets with EI
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ak8PFJetsCHSEIPruned = ak5PFJetsPruned.clone(
   rParam = cms.double(0.8),
   src = cms.InputTag("pfNoElectronsCHS"),
   doAreaFastjet = cms.bool(True),
   writeCompound = cms.bool(True),
   jetCollInstanceName=cms.string("SubJets"),
   jetPtMin = cms.double(50.)
   )

## AK4 gen jets 
process.ak4GenJets = ak4GenJets.clone(
    src = 'packedGenParticlesForJetsNoNu'
    )

## AK4 PF jets with EI 
process.ak4PFJetsCHSEI = ak4PFJets.clone(
    src = 'pfNoElectronsCHS', 
    doAreaFastjet = True
    )

### End projections

## AK8 PF jets 
process.ak8PFJetsCHS = ak4PFJets.clone(
   rParam = cms.double(0.8),
   src = cms.InputTag("pfCHS"),
   doAreaFastjet = cms.bool(True),
   jetPtMin = cms.double(50.)
   )

## Pruned AK8 PF jets (two jet collections are produced, fat jets and subjets)
process.ak8PFJetsCHSPruned = ak5PFJetsPruned.clone(
   rParam = cms.double(0.8),
   src = cms.InputTag("pfCHS"),
   doAreaFastjet = cms.bool(True),
   writeCompound = cms.bool(True),
   jetCollInstanceName=cms.string("SubJets"),
   jetPtMin = cms.double(50.)
   )

from RecoJets.JetProducers.ak8PFJetsCHS_groomingValueMaps_cfi import *
process.ak8PFJetsCHSPrunedMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
    src = cms.InputTag("ak8PFJetsCHS"),
    matched = cms.InputTag("ak8PFJetsCHSPruned"),
    distMax = cms.double(0.8),
    value = cms.string('mass')
    )

process.ak8PFJetsCHSEIPrunedMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
    src = cms.InputTag("ak8PFJetsCHSEI"),
    matched = cms.InputTag("ak8PFJetsCHSEIPruned"),
    distMax = cms.double(0.8),
    value = cms.string('mass')
    )

#################################################
## Make PAT jets
#################################################

#for Inclusive Vertex Finder
process.load("RecoBTag/Configuration/RecoBTag_cff")
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
process.inclusiveVertexFinder.tracks = cms.InputTag("unpackedTracksAndVertices")
process.inclusiveVertexFinder.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
process.trackVertexArbitrator.tracks = cms.InputTag("unpackedTracksAndVertices")
process.trackVertexArbitrator.primaryVertices = cms.InputTag("unpackedTracksAndVertices")

#new input for impactParameterTagInfos, softleptons, IVF
process.impactParameterTagInfos.jetTracks = cms.InputTag("jetTracksAssociatorAtVertexSlimmedJetsAK8BTagged")
process.impactParameterTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.inclusiveVertexFinder.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
process.trackVertexArbitrator.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
process.softPFMuonsTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.softPFElectronsTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.softPFMuonsTagInfos.jets = cms.InputTag("patJetsSlimmedJetsAK8BTagged")
process.softPFElectronsTagInfos.jets = cms.InputTag("patJetsSlimmedJetsAK8BTagged") 
process.inclusiveSecondaryVertexFinderTagInfosV2 = process.inclusiveSecondaryVertexFinderTagInfos.clone()
process.inclusiveSecondaryVertexFinderTagInfosV2.trackSelection.qualityClass = cms.string('any')

## b-tag discriminators
bTagDiscriminators = [
    'trackCountingHighEffBJetTags',
    'trackCountingHighPurBJetTags',
    'jetProbabilityBJetTags',
    'jetBProbabilityBJetTags',
    'simpleSecondaryVertexHighEffBJetTags',
    'simpleSecondaryVertexHighPurBJetTags',
    'combinedSecondaryVertexBJetTags',
    'combinedInclusiveSecondaryVertexV2BJetTags'
    ]

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

addJetCollection(
    process,
    labelName = 'AK8PFCHSEI',
    jetSource = cms.InputTag('ak8PFJetsCHSEI'),
    algo = 'ak',  # needed for jet flavor clustering
    rParam = 0.8, # needed for jet flavor clustering
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    #svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNu')
    )

addJetCollection(
    process,
    labelName = 'AK8PFCHS',
    jetSource = cms.InputTag('ak8PFJetsCHS'),
    algo = 'ak',  # needed for jet flavor clustering
    rParam = 0.8, # needed for jet flavor clustering
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    #svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNu')
    )

getattr(process,'patJetPartons').particles = cms.InputTag('prunedGenParticles')
getattr(process,'patJetPartonMatchAK8PFCHSEI').matched = cms.InputTag('prunedGenParticles')
getattr(process,'patJetPartonMatchAK8PFCHS').matched = cms.InputTag('prunedGenParticles')

#if hasattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHS'):
#  getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHS').extSVCollection = cms.InputTag('slimmedSecondaryVertices')

getattr(process,'patJetsAK8PFCHSEI').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
getattr(process,'patJetsAK8PFCHSEI').addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD
getattr(process,'patJetsAK8PFCHS').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
getattr(process,'patJetsAK8PFCHS').addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD

process.jetTracksAssociatorAtVertexAK8PFCHSEI.tracks = cms.InputTag("unpackedTracksAndVertices")
process.jetTracksAssociatorAtVertexAK8PFCHS.tracks = cms.InputTag("unpackedTracksAndVertices")

### PATify pruned fat jets with EI
addJetCollection(
    process,
    labelName = 'AK8PFCHSEIPruned',
    jetSource = cms.InputTag('ak8PFJetsCHSEIPruned'),
    btagDiscriminators = ['None'],
    jetCorrections = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
    getJetMCFlavour = False # jet flavor disabled
    )
getattr(process,'patJetPartonMatchAK8PFCHSEIPruned').matched = cms.InputTag('prunedGenParticles')

### PATify pruned fat jets
addJetCollection(
    process,
    labelName = 'AK8PFCHSPruned',
    jetSource = cms.InputTag('ak8PFJetsCHSPruned'),
    btagDiscriminators = ['None'],
    jetCorrections = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
    getJetMCFlavour = False # jet flavor disabled
    )
getattr(process,'patJetPartonMatchAK8PFCHSPruned').matched = cms.InputTag('prunedGenParticles')

### PATify pruned subjets with EI
addJetCollection(
    process,
    labelName = 'AK8PFCHSEIPrunedSubjets',
    jetSource = cms.InputTag('ak8PFJetsCHSEIPruned','SubJets'),
    algo = 'ak',  # needed for subjet flavor clustering
    rParam = 0.8, # needed for subjet flavor clustering
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    #svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNuPruned','SubJets'),
    #Apparently only in CMSSW_7_3_X explicitJTA = True,  # needed for subjet b tagging
    #Apparently only in CMSSW_7_3_X svClustering = True, # needed for subjet b tagging
    #Apparently only in CMSSW_7_3_X fatJets=cms.InputTag('ak8PFJetsCHS'),             # needed for subjet flavor clustering
    #Apparently only in CMSSW_7_3_X groomedFatJets=cms.InputTag('ak8PFJetsCHSPruned') # needed for subjet flavor clustering
    )

if hasattr( process, 'jetTracksAssociatorAtVertex' + 'AK8PFCHSEIPrunedSubjets' ):
  process.jetTracksAssociatorAtVertexAK8PFCHSEIPrunedSubjets.tracks = cms.InputTag("unpackedTracksAndVertices")

getattr(process,'patJetPartonMatchAK8PFCHSEIPrunedSubjets').matched = cms.InputTag('prunedGenParticles')

getattr(process,'patJetsAK8PFCHSEIPrunedSubjets').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
getattr(process,'patJetsAK8PFCHSEIPrunedSubjets').addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD

process.selectedPatJetsAK8PFCHSEI.cut = cms.string("pt > 100 && abs(eta)<4.") 

### PATify pruned subjets
addJetCollection(
    process,
    labelName = 'AK8PFCHSPrunedSubjets',
    jetSource = cms.InputTag('ak8PFJetsCHSPruned','SubJets'),
    algo = 'ak',  # needed for subjet flavor clustering
    rParam = 0.8, # needed for subjet flavor clustering
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    #svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNuPruned','SubJets'),
    #Apparently only in CMSSW_7_3_X explicitJTA = True,  # needed for subjet b tagging
    #Apparently only in CMSSW_7_3_X svClustering = True, # needed for subjet b tagging
    #Apparently only in CMSSW_7_3_X fatJets=cms.InputTag('ak8PFJetsCHS'),             # needed for subjet flavor clustering
    #Apparently only in CMSSW_7_3_X groomedFatJets=cms.InputTag('ak8PFJetsCHSPruned') # needed for subjet flavor clustering
    )

if hasattr( process, 'jetTracksAssociatorAtVertex' + 'AK8PFCHSPrunedSubjets' ):
  process.jetTracksAssociatorAtVertexAK8PFCHSPrunedSubjets.tracks = cms.InputTag("unpackedTracksAndVertices")
#  from RecoJets.JetAssociationProducers.ak4aTA_cff import ak4JetTracksAssociatorExplicit
#  m = 'jetTracksAssociatorAtVertex' + 'AK8PFCHSPrunedSubjets'
#  print 'Switching ' + m + ' to explicit jet-track association'
#  setattr( process, m, ak4JetTracksAssociatorExplicit.clone(
#    jets = getattr(getattr(process,m),'jets'),
#    tracks = cms.InputTag("unpackedTracksAndVertices")
#    ) 
#    )

getattr(process,'patJetPartonMatchAK8PFCHSPrunedSubjets').matched = cms.InputTag('prunedGenParticles')
#if hasattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHSPrunedSubjets'):
#  getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHSPrunedSubjets').extSVCollection = cms.InputTag('slimmedSecondaryVertices')
getattr(process,'patJetsAK8PFCHSPrunedSubjets').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
getattr(process,'patJetsAK8PFCHSPrunedSubjets').addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD

process.selectedPatJetsAK8PFCHS.cut = cms.string("pt > 100 && abs(eta)<4.") 

## Establish references between PATified fat jets and subjets with EI using the BoostedJetMerger
process.selectedPatJetsAK8PFCHSEIPrunedPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsAK8PFCHSEIPruned"),
    subjetSrc=cms.InputTag("selectedPatJetsAK8PFCHSEIPrunedSubjets")
    )

## Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsAK8PFCHSPrunedPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsAK8PFCHSPruned"),
    subjetSrc=cms.InputTag("selectedPatJetsAK8PFCHSPrunedSubjets")
    )

from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

addJetCollection(
    process,
    postfix   = "",
    labelName = 'SlimmedJetsAK8BTagged',
    jetSource = cms.InputTag('slimmedJetsAK8'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'), 
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    jetCorrections = ('AK8PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags', 'combinedInclusiveSecondaryVertexV2BJetTags'],
    algo= 'AK', 
    rParam = 0.8
    )

process.patJetFlavourAssociation.jets = "slimmedJetsAK8"
process.patJetFlavourAssociation.rParam = 0.8 

#adjust MC matching
process.patJetGenJetMatchSlimmedJetsAK8BTagged.matched = "slimmedGenJets"
process.patJetPartonMatchSlimmedJetsAK8BTagged.matched = "prunedGenParticles"
process.patJetPartons.particles = "prunedGenParticles"

#adjust PV used for Jet Corrections
process.patJetCorrFactorsSlimmedJetsAK8BTagged.primaryVertices = "unpackedTracksAndVertices"

#adjust JTA cone size 
process.jetTracksAssociatorAtVertexSlimmedJetsAK8BTagged.coneSize = 0.8 

process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.combinedSecondaryVertex.trackMultiplicityMin = 1 #silly sv, uses un filtered tracks.. i.e. any pt

process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV()>0"))

### Add Nsubjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
process.Njettiness = Njettiness.clone(
    src = cms.InputTag("ak8PFJetsCHS"),
    cone = cms.double(0.8)
    )

process.patJetsAK8PFCHS.userData.userFloats.src += ['ak8PFJetsCHSPrunedMass','Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']

process.NjettinessEI = Njettiness.clone(
    src = cms.InputTag("ak8PFJetsCHSEI"),
    cone = cms.double(0.8)
    )

process.patJetsAK8PFCHSEI.userData.userFloats.src += ['ak8PFJetsCHSEIPrunedMass','NjettinessEI:tau1','NjettinessEI:tau2','NjettinessEI:tau3']


#$#$#$#$#$#$#$#$#$#
#   TOP TAG JETS  #
# patJetsCMSTopTagCHS
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
process.cmsTopTagCHS = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters.clone( src = cms.InputTag('chs'),
                           doAreaFastjet = cms.bool(True),
                           doRhoFastjet = cms.bool(False),
                           jetPtMin = cms.double(200.0)
                           ),
    AnomalousCellParameters,
    CATopJetParameters.clone( jetCollInstanceName = cms.string("SubJets"),
                              verbose = cms.bool(False),
                              algorithm = cms.int32(1), # 0 = KT, 1 = CA, 2 = anti-KT
                              tagAlgo = cms.int32(0), #0=legacy top
                              useAdjacency = cms.int32(2), # modified adjacency
                              centralEtaCut = cms.double(2.5), # eta for defining "central" jets
                              sumEtBins = cms.vdouble(0,1600,2600), # sumEt bins over which cuts vary. vector={bin 0 lower bound, bin 1 lower bound, ...}
                              rBins = cms.vdouble(0.8,0.8,0.8), # Jet distance paramter R. R values depend on sumEt bins.
                              ptFracBins = cms.vdouble(0.05,0.05,0.05), # minimum fraction of central jet pt for subjets (deltap)
                              deltarBins = cms.vdouble(0.19,0.19,0.19), # Applicable only if useAdjacency=1. deltar adjacency values for each sumEtBin
                              nCellBins = cms.vdouble(1.9,1.9,1.9),
                            ),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True)
    )

process.CATopTagInfos = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("cmsTopTagCHS"),
                                    TopMass = cms.double(171),
                                    TopMassMin = cms.double(0.),
                                    TopMassMax = cms.double(250.),
                                    WMass = cms.double(80.4),
                                    WMassMin = cms.double(0.0),
                                    WMassMax = cms.double(200.0),
                                    MinMassMin = cms.double(0.0),
                                    MinMassMax = cms.double(200.0),
                                    verbose = cms.bool(False)
                                    )
addJetCollection(
    process,
    labelName = 'CMSTopTagCHS',
    jetSource = cms.InputTag('cmsTopTagCHS'),
    jetCorrections = ('AK8PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag("unpackedTracksAndVertices"),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False
    )

process.patJetPartonMatchCMSTopTagCHS.matched='prunedGenParticles'
process.patJetCorrFactorsCMSTopTagCHS.primaryVertices = "unpackedTracksAndVertices"
process.patJetGenJetMatchCMSTopTagCHS.matched = 'slimmedGenJets'
process.patJetPartonMatchCMSTopTagCHS.matched = 'prunedGenParticles'
#process.jetTracksAssociatorAtVertexCMSTopTagCHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('cmsTopTagCHS'), coneSize = 0.8)
process.secondaryVertexTagInfosCMSTopTagCHS.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCMSTopTagCHS.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCMSTopTagCHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCMSTopTagCHS.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCMSTopTagCHS.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCMSTopTagCHS.jetTagComputer = cms.string('combinedSecondaryVertexCMSTopTagCHS')
process.patJetsCMSTopTagCHS.addTagInfos = True
process.patJetsCMSTopTagCHS.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopTagInfos')
    )

addJetCollection(
    process,
    labelName = 'CMSTopTagCHSSubjets',
    jetSource = cms.InputTag('cmsTopTagCHS','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag("unpackedTracksAndVertices"),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )
process.patJetPartonMatchCMSTopTagCHSSubjets.matched='prunedGenParticles'
process.patJetCorrFactorsCMSTopTagCHSSubjets.primaryVertices = "unpackedTracksAndVertices"
process.patJetGenJetMatchCMSTopTagCHSSubjets.matched = 'slimmedGenJets'
process.patJetPartonMatchCMSTopTagCHSSubjets.matched = 'prunedGenParticles'

process.patJetsCMSTopTagCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsCMSTopTagCHS" ),
    subjetSrc=cms.InputTag("patJetsCMSTopTagCHSSubjets")
    )

### Selected leptons and jets
process.skimmedPatMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag(muLabel),
    cut = cms.string("pt > 10 && abs(eta) < 2.5 && (isPFMuon && (isGlobalMuon || isTrackerMuon) )")
    )

process.skimmedPatElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag(elLabel),
    cut = cms.string("pt > 20 && abs(eta) < 2.5 && gsfTrack.isAvailable() && gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2")
    )

process.skimmedPatMET = cms.EDFilter(
    "PATMETSelector",
    src = cms.InputTag(metLabel),
    cut = cms.string("")
    )


process.skimmedPatJets = cms.EDFilter(
    "PATJetSelector",
    src = cms.InputTag(jLabel),
    cut = cms.string(" pt > 25 && abs(eta) < 4.")
    )

process.skimmedPatJetsAK8EI = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag(ak8jetEILabel),
    cut = cms.string("pt > 100 && abs(eta) < 4.")    
    )


process.skimmedPatJetsAK8 = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag(ak8jetLabel),
    cut = cms.string("pt > 100 && abs(eta) < 4.")    
    )

process.skimmedPatSubJetsAK8EI = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag(ak8subjetEILabel),
    cut = cms.string("")
    )

process.skimmedPatSubJetsAK8 = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag(ak8subjetLabel),
    cut = cms.string("")
    )

process.skimmedCMSTOPTAGSubJets = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag("selectedPatJetsCMSTopTagCHSSubjets"),
    cut = cms.string("pt > 1")
)


### Asking for at least 2 jets satisfying the selection above 
process.jetFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("skimmedPatJets"),
    minNumber = cms.uint32(2),
    filter = cms.bool(True)
    )

process.muonUserData = cms.EDProducer(
    'MuonUserData',
    muonLabel = cms.InputTag("skimmedPatMuons"),
    pv        = cms.InputTag(pvLabel),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltMuonFilter  = cms.InputTag(hltMuonFilterLabel),
    hltPath            = cms.string("HLT_IsoMu40_eta2p1_v11"),
    hlt2reco_deltaRmax = cms.double(0.1),
    # mainROOTFILEdir    = cms.string("../data/")
    )

process.jetUserData = cms.EDProducer(
    'JetUserData',
    pv        = cms.InputTag(pvLabel),
    jetLabel  = cms.InputTag("skimmedPatJets"),
    packedjetLabel  = cms.InputTag(""),
    subjetLabel  = cms.InputTag(""),
    doSubjets = cms.bool(False),
    elLabel   = cms.InputTag("skimmedPatElectrons"), 
    muLabel   = cms.InputTag("skimmedPatMuons"), 
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2),
    )

process.ak8jetEIUserData = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag("selectedPatJetsAK8PFCHSEI"),
    pv        = cms.InputTag(pvLabel),
    elLabel   = cms.InputTag("skimmedPatElectrons"), 
    muLabel   = cms.InputTag("skimmedPatMuons"), 
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2),
    doSubjets = cms.bool(True),
    packedjetLabel  = cms.InputTag("selectedPatJetsAK8PFCHSEIPrunedPacked"),
    subjetLabel  = cms.InputTag("selectedPatJetsAK8PFCHSEIPrunedSubjets"),
)


process.ak8jetUserData = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag("selectedPatJetsAK8PFCHS"),
#    jetLabel  = cms.InputTag("selectedPatJetsAK8PFCHSPrunedPacked"),
    pv        = cms.InputTag(pvLabel),
    elLabel   = cms.InputTag("skimmedPatElectrons"), 
    muLabel   = cms.InputTag("skimmedPatMuons"), 
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2),
    doSubjets = cms.bool(True),
       packedjetLabel  = cms.InputTag("selectedPatJetsAK8PFCHSPrunedPacked"),
       subjetLabel  = cms.InputTag("selectedPatJetsAK8PFCHSPrunedSubjets"),
)

process.ak8subjetsEIUserData = cms.EDProducer(
    'JetUserData',
    pv        = cms.InputTag(pvLabel),
    jetLabel  = cms.InputTag("skimmedPatSubJetsAK8EI"),
    packedjetLabel  = cms.InputTag(""),
    subjetLabel  = cms.InputTag(""),
    doSubjets = cms.bool(False),
    elLabel   = cms.InputTag("skimmedPatElectrons"), 
    muLabel   = cms.InputTag("skimmedPatMuons"), 
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2)
    )

process.subjetUserDataAK8 = cms.EDProducer(
    'JetUserData',
    pv        = cms.InputTag(pvLabel),
    jetLabel  = cms.InputTag("skimmedPatSubJetsAK8EI"),
    packedjetLabel  = cms.InputTag(""),
    subjetLabel  = cms.InputTag(""),
    doSubjets = cms.bool(False),
    elLabel   = cms.InputTag("skimmedPatElectrons"), 
    muLabel   = cms.InputTag("skimmedPatMuons"), 
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2)
    )


process.cmstoptagjetUserData = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag("patJetsCMSTopTagCHS"),
    packedjetLabel  = cms.InputTag("patJetsCMSTopTagCHSPacked"),
    subjetLabel  = cms.InputTag("patJetsCMSTopTagCHSSubjets"),
    pv        = cms.InputTag(pvLabel),
    elLabel   = cms.InputTag("skimmedPatElectrons"), 
    muLabel   = cms.InputTag("skimmedPatMuons"), 
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2),
    doSubjets          = cms.bool(True)
)

process.electronUserData = cms.EDProducer(
    'ElectronUserData',
    eleLabel = cms.InputTag("skimmedPatElectrons"),
    pv        = cms.InputTag(pvLabel),
    conversion        = cms.InputTag(convLabel),
    triggerResults = cms.InputTag(triggerResultsLabel),
    triggerSummary = cms.InputTag(triggerSummaryLabel),
    hltElectronFilter  = cms.InputTag(hltElectronFilterLabel),  ##trigger matching code to be fixed!
    hltPath             = cms.string("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"),
    #electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-veto"),
    #electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-tight"),
    )



from PhysicsTools.CandAlgos.EventShapeVars_cff import *
process.eventShapePFVars = pfEventShapeVars.clone()
process.eventShapePFVars.src = cms.InputTag(particleFlowLabel)

process.eventShapePFJetVars = pfEventShapeVars.clone()
process.eventShapePFJetVars.src = cms.InputTag("skimmedPatJets")

process.centrality = cms.EDProducer("CentralityUserData",
    src = cms.InputTag("skimmedPatJets")
    )                                    

### Including ntuplizer 
process.load("B2GAnaFW.B2GAnaFW.b2gedmntuples_cff")

process.options.allowUnscheduled = cms.untracked.bool(True)


### definition of Analysis sequence
process.analysisPath = cms.Path(
    process.selectedPatJetsAK8PFCHSEI + 
    process.selectedPatJetsAK8PFCHSEIPrunedSubjets + 
    process.selectedPatJetsAK8PFCHSEIPrunedPacked + 
    process.selectedPatJetsAK8PFCHS + 
    process.selectedPatJetsAK8PFCHSPrunedSubjets + 
    process.selectedPatJetsAK8PFCHSPrunedPacked + 
    process.skimmedPatElectrons +
    process.skimmedPatMuons +
    process.skimmedPatJets +
    process.skimmedPatJetsAK8EI +
    process.skimmedPatSubJetsAK8EI +
    process.skimmedPatJetsAK8 +
    process.skimmedPatSubJetsAK8 +
    process.skimmedPatMET +
    process.eventShapePFVars +
    process.eventShapePFJetVars +
    process.centrality
    )

#process.analysisPath+=process.jetFilter

#process.analysisPath+=process.egmGsfElectronIDSequence
process.analysisPath+=process.muonUserData
process.analysisPath+=process.jetUserData
#process.analysisPath+=process.jetUserDataAK8
#process.analysisPath+=process.subjetUserDataAK8
process.analysisPath+=process.electronUserData
process.analysisPath+=process.genPart
process.analysisPath+=process.muons
process.analysisPath+=process.electrons
process.analysisPath+=process.jetsAK4
process.analysisPath+=process.jetsAK8EI
process.analysisPath+=process.subjetsAK8EI
process.analysisPath+=process.jetsAK8
process.analysisPath+=process.subjetsAK8
process.analysisPath+=process.met

### Creating the filter path to use in order to select events
process.filterPath = cms.Path(
    process.jetFilter
    )

### keep info from LHEProducts if they are stored in PatTuples
if(options.LHE):
  process.LHEUserData = cms.EDProducer("LHEUserData",
  lheLabel = cms.InputTag(options.lheLabel)
  )
  process.analysisPath+=process.LHEUserData
  process.edmNtuplesOut.outputCommands+=('keep *_*LHE*_*_*',)
  process.edmNtuplesOut.outputCommands+=('keep LHEEventProduct_*_*_*',)
### end LHE products     

process.edmNtuplesOut.outputCommands+=('keep *_generator_*_*',)
process.edmNtuplesOut.fileName=options.outputLabel

process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('filterPath')
    )


process.fullPath = cms.Schedule(
    process.analysisPath,
    process.filterPath
    )

process.endPath = cms.EndPath(process.edmNtuplesOut)

#open('B2GEntupleFileDump.py','w').write(process.dumpPython())
