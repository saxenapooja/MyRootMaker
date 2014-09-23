import FWCore.ParameterSet.Config as cms
process = cms.Process("ROOTMAKER")

isData = False
skim = 128
year = 2012
period = 'Summer12'
useMarkov = False
if isData:
	globalTag = 'FT_53_V6C_AN4::All'
else:
	globalTag = 'START53_V23::All'

# Disable SVfit
doSVFit = False

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1 #0000
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')

process.GlobalTag.globaltag = cms.string(globalTag)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
#                'file:/storage/9/burgmeier/test/_SingleMu_Run2012A-13Jul2012-v1_AOD/FE749D99-45D2-E111-983D-1CC1DE1CF196.root'
#                'file:/portal/ekpcms5/home/burgmeier/CMSSW_4_2_8_TrackIso/src/pickevents5.root'
#                'file:/storage/9/burgmeier/test/_TauPlusX_Run2012A-PromptReco-v1_AOD/FA4A682D-E28C-E111-9347-001D09F2B30B.root'
#		 'file:/storage/9/burgmeier/test/_WH_ZH_TTH_HToTauTau_M-125_8TeV-pythia6-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/DCB0A3A6-73EF-E111-913D-00215E21D57C.root'
#		 'file:///storage/5/burgmeier/test/0292097E-6379-E211-8C11-E0CB4E29C4F6.root'
#	"/store/data/Run2012D/BJetPlusX/AOD/PromptReco-v1/000/203/773/C8EEC2FC-5E0A-E211-8058-001D09F26509.root"
#	"/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v1/000/161/008/8468EAFF-C855-E011-BB51-003048F117B4.root"
#/pnfs/desy.de/cms/tier2/store/mc/Summer12_DR53X/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/52ACA0C7-53E9-E111-A786-002618943919.root
	"/store/mc/Summer12_DR53X/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/52ACA0C7-53E9-E111-A786-002618943919.root"
#     'file:pickevents.root'
	)
#     noEventSort = cms.untracked.bool(True),
#     duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

# Load various tools before PAT
import os
cmssw_version = os.environ["CMSSW_VERSION"].replace("CMSSW_","")
if cmssw_version.startswith('4'): #earlier value was 4
        process.load("JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_42X_cff")
else:
        process.load("JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_cff")
	process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
	process.load("RecoJets.JetProducers.PileupJetID_cfi")
	
	process.load('CommonTools.ParticleFlow.pfMuons_cff')
	from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFMuonIso
	process.muIsoSequence = setupPFMuonIso(process, 'pfSelectedMuons')

	process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
	from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
	process.elIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
	
# Produce rho for PU correction
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True

# Load PAT stuff
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# We PAT taus for SVfit. However, to create PAT taus we also need all the
# other objects. So we cannot remove them from the PAT sequence any longer.
#process.patDefaultSequence.remove(getattr(process, 'selectedPatCandidates'))
#process.patDefaultSequence.remove(getattr(process, 'cleanPatCandidates'))
process.patDefaultSequence.remove(getattr(process, 'countPatCandidates'))

#process.patCandidates.remove(getattr(process, 'makePatElectrons'))
#process.patCandidates.remove(getattr(process, 'makePatMuons'))
#process.patCandidates.remove(getattr(process, 'makePatTaus'))
#process.patCandidates.remove(getattr(process, 'makePatPhotons'))
#process.patCandidates.remove(getattr(process, 'makePatMETs'))
from PhysicsTools.PatAlgos.tools.coreTools import *
if isData:
	removeMCMatching(process, ['All'], "", None, [])
	#removeMCMatching(process, ['All'], "", False)
	jetCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']#, 'L5Flavor', 'L7Parton']
else:
	jetCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']#, 'L5Flavor', 'L7Parton']

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,
	cms.InputTag('ak5PFJets'),
	outputModules = [],
	doJTA = True,
	doBTagging = True,
	jetCorrLabel = ('AK5PF', jetCorrLevels),
	doType1MET = False, # is only used for calo jets
	genJetCollection = cms.InputTag('ak5GenJets'),
	doJetID = True,
	jetIdLabel = 'ak5')
process.patJets.addTagInfos = cms.bool(True)
#process.patJetCorrFactorsAK5PF.rho = cms.InputTag("kt6PFJets", "rho")

# Use HPS for PAT taus.
from PhysicsTools.PatAlgos.tools.tauTools import *
# Only required for pre-52X:
switchToPFTauHPS(process)

# SVfit
process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.patElectrons.electronIDSources          = cms.PSet(mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"))
#process.patMuons.pfMuonSource = cms.InputTag('') # just to make sure this does not used only isolated muons
process.cleanPatElectronsForSVfit               = process.cleanPatElectrons.clone()
process.cleanPatElectronsForSVfit.preselection  = cms.string('pt > 24 & abs(eta) < 2.1 & electronID("mvaNonTrigV0") > 0.9')
process.cleanPatElectronsForSVfit.finalCut      = cms.string('')
process.cleanPatMuonsForSVfit                   = process.cleanPatMuons.clone()
process.cleanPatMuonsForSVfit.preselection      = cms.string('pt > 24 & abs(eta) < 2.1 & isGlobalMuon()')
process.cleanPatMuonsForSVfit.finalCut          = cms.string('')
process.cleanPatTausForSVfit                    = process.cleanPatTaus.clone()
process.cleanPatTausForSVfit.finalCut           = cms.string('')
process.cleanPatTausForSVfit.preselection       = cms.string('pt > 20 & abs(eta) < 2.3 & tauID("decayModeFinding") > 0.5 & (tauID("againstMuonLoose") > 0.5 || tauID("againstMuonLoose2")) & (tauID("againstElectronLoose") > 0.5)')

# Type-1 MET
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.patMETCorrections)

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.patMETsPFType1     = process.patMETs.clone(
	metSource          = cms.InputTag('pfType1CorrectedMet'),
	addMuonCorrections = cms.bool(False),
	genMETSource       = cms.InputTag('genMetTrue'))

# Enable MVA MET
process.isolatedPatElectrons              = process.cleanPatElectrons.clone()
process.isolatedPatElectrons.preselection = cms.string('pt > 20 & abs(eta) < 2.1 & electronID("mvaNonTrigV0") > 0.905 &  (chargedHadronIso() + max(0.0, neutralHadronIso() + photonIso() - 0.5 * puChargedHadronIso()))/pt < 0.1')
process.isolatedPatElectrons.finalCut     = cms.string('')
process.isolatedPatMuons                  = process.cleanPatMuons.clone()
process.isolatedPatMuons.preselection     = cms.string('pt > 20 & abs(eta) < 2.1 & isGlobalMuon() &  (chargedHadronIso() + max(0.0, neutralHadronIso() + photonIso() - 0.5 * puChargedHadronIso()))/pt < 0.1')
process.isolatedPatMuons.finalCut         = cms.string('')
process.isolatedPatTaus                   = process.cleanPatTaus.clone()
# Don't use the new discriminators here so that we have the same definition as the others
#process.isolatedPatTaus.preselection = cms.string('pt > 20 & abs(eta) < 2.3 & tauID("decayModeFinding") > 0.5 & (tauID("againstMuonLoose") > 0.5 || tauID("againstMuonLoose2")) & (tauID("againstElectronLoose") > 0.5 || tauID("againstElectronMVA") > 0.5 || tauID("againstElectronLooseMVA3") > 0.5) & (tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 || tauID("byLooseIsolationMVA2") > 0.5 || tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5)')
process.isolatedPatTaus.preselection      = cms.string('pt > 20 && abs(eta) < 2.3 && tauID("decayModeFinding") > 0.5 && tauID("againstMuonLoose") > 0.5 && tauID("againstElectronLoose") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5')
process.isolatedPatTaus.finalCut          = cms.string('')

if isData:
	process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
else:
	process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3")
process.pfMEtMVA.srcLeptons = cms.VInputTag('isolatedPatElectrons', 'isolatedPatMuons', 'isolatedPatTaus')
process.pfMEtMVA.verbosity = cms.int32(0)

process.patPFMetByMVA = process.patMETs.clone(
	metSource          = cms.InputTag('pfMEtMVA'),
	addMuonCorrections = cms.bool(False),
	genMETSource       = cms.InputTag('genMetTrue'))

process.pfMEtMVACov = cms.EDProducer("PFMEtSignCovMatrixUnembedder",
				     src = cms.InputTag("patPFMetByMVA")
				     )

# SVfit
process.load("TauAnalysis.CandidateTools.diTauPairProduction_cff")
process.diTau = process.selectedDiTauPairs2ndTauElectronVetoCumulative.clone()
process.diTau.srcLeg1 = cms.InputTag("cleanPatTausForSVfit")
process.diTau.srcLeg2 = cms.InputTag("cleanPatTausForSVfit")

process.diTau.srcMET      = cms.InputTag("patPFMetByMVA")
process.diTau.doPFMEtSign = cms.bool(False)
process.diTau.srcPrimaryVertex = cms.InputTag("offlinePrimaryVertices")
process.diTau.dRmin12     = cms.double(0.5)
process.diTau.doSVreco    = cms.bool(False)  #False, to run Standalone SVFit
'''
process.diTau.nSVfit.psKine_MEt_int.config.event.srcMEt = cms.InputTag("patPFMetByMVA")
process.diTau.nSVfit.psKine_MEt_int.config.event.likelihoodFunctions[0].srcMEtCovMatrix = cms.InputTag("pfMEtMVACov")
process.diTau.nSVfit.psKine_MEt_int.config.event.srcPrimaryVertex = cms.InputTag("offlinePrimaryVertices")
if useMarkov:
	process.diTau.nSVfit.psKine_MEt_int.algorithm = cms.PSet(
		pluginName = cms.string("nSVfitAlgorithmByIntegration2"),
		pluginType = cms.string("NSVfitAlgorithmByIntegration2"),
		markovChainOptions = cms.PSet(
		mode = cms.string("Metropolis"),
		initMode = cms.string("Gaus"),
		numIterBurnin = cms.uint32(10000),
		numIterSampling = cms.uint32(100000),
		numIterSimAnnealingPhase1 = cms.uint32(2000),
		numIterSimAnnealingPhase2 = cms.uint32(6000),
		T0 = cms.double(15.),
		alpha = cms.double(0.999),
		numChains = cms.uint32(1),
		numBatches = cms.uint32(1),
		L = cms.uint32(1),
		epsilon0 = cms.double(1.e-2),
		nu = cms.double(0.71)
		),
		max_or_median = cms.string("max"),
		verbosity = cms.int32(0)
		)
delattr(process.diTau.nSVfit, "psKine_MEt_logM_fit")
'''

process.selectedDiTau = cms.EDFilter("PATDiTauPairSelector",
				     src = cms.InputTag("diTau"),
				     cut = cms.string("dR12 > 0.5")
				     )
process.selectedDiTauFilter = cms.EDFilter("CandViewCountFilter",
					   src = cms.InputTag("selectedDiTau"),
					   minNumber = cms.uint32(1),
					   maxNumber = cms.uint32(999),
					   )

#import copy
#process.MuTauTauPairsWHInt   = copy.deepcopy(process.nSVfitProducerByIntegrationWH)

##################################################
# Muon isolation in 42X
from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedHadrons_cfi import pfAllChargedHadrons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi import pfAllNeutralHadrons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import pfAllPhotons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllMuons_cfi import pfAllMuons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllElectrons_cfi import pfAllElectrons
from CommonTools.ParticleFlow.Isolation.tools_cfi import *

pfAllChargedParticlesPdgIds = []
pfAllChargedParticlesPdgIds.extend(pfAllChargedHadrons.pdgId.value())
pfAllChargedParticlesPdgIds.extend(pfAllMuons.pdgId.value())
pfAllChargedParticlesPdgIds.extend(pfAllElectrons.pdgId.value())

process.pfChargedHadrons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfNoPileUpIso'), pdgId = cms.vint32(pfAllChargedHadrons.pdgId.value()))
process.pfChargedParticles = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfNoPileUpIso'), pdgId = cms.vint32(pfAllChargedParticlesPdgIds))
process.pfNeutralHadrons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('particleFlow'), pdgId = cms.vint32(pfAllNeutralHadrons.pdgId.value()))
process.pfPhotons = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('particleFlow'), pdgId = cms.vint32(pfAllPhotons.pdgId.value()))
process.pfPU = cms.EDFilter("PdgIdPFCandidateSelector", src = cms.InputTag('pfPileUpIso'), pdgId = cms.vint32(pfAllChargedParticlesPdgIds))

process.pfMuIsoDepositChargedHadrons = isoDepositReplace("muons", "pfChargedHadrons")
process.pfMuIsoDepositChargedParticles = isoDepositReplace("muons", "pfChargedParticles")
process.pfMuIsoDepositNeutralHadrons = isoDepositReplace("muons", "pfNeutralHadrons")
process.pfMuIsoDepositPhotons = isoDepositReplace("muons", "pfPhotons")
process.pfMuIsoDepositPU = isoDepositReplace("muons", "pfPU")

process.pfElIsoDepositChargedHadrons = isoDepositReplace("gsfElectrons", "pfChargedHadrons")
process.pfElIsoDepositChargedParticles = isoDepositReplace("gsfElectrons", "pfChargedParticles")
process.pfElIsoDepositNeutralHadrons = isoDepositReplace("gsfElectrons", "pfNeutralHadrons")
process.pfElIsoDepositPhotons = isoDepositReplace("gsfElectrons", "pfPhotons")
process.pfElIsoDepositPU = isoDepositReplace("gsfElectrons", "pfPU")

process.isoSequence = cms.Sequence(process.pfChargedHadrons*process.pfChargedParticles*process.pfNeutralHadrons*process.pfPhotons*process.pfPU*process.pfMuIsoDepositChargedHadrons*process.pfMuIsoDepositChargedParticles*process.pfMuIsoDepositNeutralHadrons*process.pfMuIsoDepositPhotons*process.pfMuIsoDepositPU*process.pfElIsoDepositChargedHadrons*process.pfElIsoDepositChargedParticles*process.pfElIsoDepositNeutralHadrons*process.pfElIsoDepositPhotons*process.pfElIsoDepositPU)

#####################################
# PAT Muon isolation
process.muPFIsoValueCharged04PFIso.deposits[0].vetos = cms.vstring('0.0001', 'Threshold(0.0)')
process.muPFIsoValueChargedAll04PFIso.deposits[0].vetos = cms.vstring('0.0001', 'Threshold(0.0)')
process.muPFIsoValueGamma04PFIso.deposits[0].vetos = cms.vstring('0.01', 'Threshold(0.5)')
process.muPFIsoValueNeutral04PFIso.deposits[0].vetos = cms.vstring('0.01', 'Threshold(0.5)')
process.muPFIsoValuePU04PFIso.deposits[0].vetos = cms.vstring('0.01', 'Threshold(0.5)')
process.patMuons.useParticleFlow = cms.bool(True)
process.muonMatch.src = cms.InputTag('pfSelectedMuons')
process.patMuons.pfMuonSource = cms.InputTag('pfSelectedMuons')
process.patMuons.muonSource = cms.InputTag('') # reco muons are not used -- reset input tag so we don't even try
process.patMuons.embedTrack = cms.bool(True) # Embed tracks for muon ID cuts to be done offline
process.patMuons.isolationValues = cms.PSet(
      pfChargedHadrons = cms.InputTag('muPFIsoValueCharged04PFIso'),
      pfChargedAll = cms.InputTag('muPFIsoValueChargedAll04PFIso'),
      pfPUChargedHadrons = cms.InputTag('muPFIsoValuePU04PFIso'),
      pfNeutralHadrons = cms.InputTag('muPFIsoValueNeutral04PFIso'),
      pfPhotons = cms.InputTag('muPFIsoValueGamma04PFIso'))
#process.selectedPatMuons.cut = cms.string('isGlobalMuon & pt > 10 & abs(eta) < 2.5')

#####################################
# PAT Electron isolation
process.elPFIsoValueCharged04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)', 'EcalBarrel:ConeVeto(0.010)', 'Threshold(0.0)')
process.elPFIsoValueChargedAll04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)', 'EcalBarrel:ConeVeto(0.010)', 'Threshold(0.0)')
process.elPFIsoValueGamma04NoPFIdPFIso.deposits[0].vetos = cms.vstring('ConeVeto(0.08)', 'Threshold(0.0)')
process.elPFIsoValueNeutral04NoPFIdPFIso.deposits[0].vetos = cms.vstring('Threshold(0.0)')
process.elPFIsoValuePU04NoPFIdPFIso.deposits[0].vetos = cms.vstring('Threshold(0.0)')
process.patElectrons.electronIDSources = cms.PSet(
      mvaTrigV0 = cms.InputTag("mvaTrigV0"),
      mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"))
process.patElectrons.isolationValuesNoPFId = cms.PSet(
      pfChargedHadrons = cms.InputTag('elPFIsoValueCharged04NoPFIdPFIso'),
      pfChargedAll = cms.InputTag('elPFIsoValueChargedAll04NoPFIdPFIso'),
      pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU04NoPFIdPFIso'),
      pfNeutralHadrons = cms.InputTag('elPFIsoValueNeutral04NoPFIdPFIso'),
      pfPhotons = cms.InputTag('elPFIsoValueGamma04NoPFIdPFIso'))
process.elPFIsoValueCharged04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)', 'EcalBarrel:ConeVeto(0.010)', 'Threshold(0.0)')
process.elPFIsoValueChargedAll04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)', 'EcalBarrel:ConeVeto(0.010)', 'Threshold(0.0)')
process.elPFIsoValueGamma04PFIdPFIso.deposits[0].vetos = cms.vstring('ConeVeto(0.08)', 'Threshold(0.0)')
process.elPFIsoValueNeutral04PFIdPFIso.deposits[0].vetos = cms.vstring('Threshold(0.0)')
process.elPFIsoValuePU04PFIdPFIso.deposits[0].vetos = cms.vstring('Threshold(0.0)')
process.patElectrons.isolationValues = cms.PSet(
      pfChargedHadrons = cms.InputTag('elPFIsoValueCharged04PFIdPFIso'),
      pfChargedAll = cms.InputTag('elPFIsoValueChargedAll04PFIdPFIso'),
      pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU04PFIdPFIso'),
      pfNeutralHadrons = cms.InputTag('elPFIsoValueNeutral04PFIdPFIso'),
      pfPhotons = cms.InputTag('elPFIsoValueGamma04PFIdPFIso'))

#process.selectedPatElectrons.cut = cms.string('pt > 10 & abs(eta) < 2.5')

##################################################
# Main
process.makeroottree = cms.EDAnalyzer("RootMaker",
IsData = cms.untracked.bool(isData),
GenSomeParticles = cms.untracked.bool(not isData and (skim & 128) == 0),
GenAllParticles = cms.untracked.bool( not isData and (skim & 128) != 0),
GenMet = cms.untracked.bool(not isData),
Year = cms.untracked.uint32(year),
Period = cms.untracked.string(period),
Skim = cms.untracked.uint32(skim),
TrackCollectionTag = cms.InputTag("generalTracks"),
PVTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
Trigger = cms.untracked.bool(True),
HLTriggerSelection = cms.untracked.vstring(),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),

RecCaloMet = cms.untracked.bool(False),
RecTCMet = cms.untracked.bool(False),
RecPFMet = cms.untracked.bool(True),

RecTrack = cms.untracked.bool(False),
RecTrackPtMin = cms.untracked.double(20.),
RecTrackEtaMax = cms.untracked.double(2.5),
RecTrackNum = cms.untracked.int32(100000),
RecTrackFilterPtMin = cms.untracked.double(20.),

RecSuperCluster = cms.untracked.bool(True),
RecSuperClusterBasicCluster = cms.untracked.bool(True),
RecSuperClusterHit = cms.untracked.bool(True),

RecMuon = cms.untracked.bool(True),
RecMuonPtMin = cms.untracked.double(20.),
RecMuonTrackIso = cms.untracked.double(1000000),
RecMuonEtaMax = cms.untracked.double(2.5),
RecMuonNum = cms.untracked.int32(0),
RecMuonHLTriggerMatching = cms.untracked.vstring(
'HLT_IsoMu24_v.*:hltSingleMuIsoL3IsoFiltered24,hltSingleMuL2QualIsoL3IsoFiltered24,hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15',
'HLT_IsoMu24_eta2p1_v.*:hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f24L3IsoFiltered,hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10,hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15',
'HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v.*:hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoFiltered10,hltL3crIsoL1sMu14erORMu16erL1f0L2f14QL3f17QL3crIsoRhoFiltered0p15,hltOverlapFilterIsoMu18LooseIsoPFTau20',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltL3crIsoL1sMu14erORMu16erL1f0L2f14QL3f17QL3crIsoRhoFiltered0p15,hltOverlapFilterIsoMu17LooseIsoPFTau20',
'HLT_IsoMu15_eta2p1_L1ETM20_v.*:hltL3crIsoL1sMu12Eta2p1L1f0L2f12QL3f15QL3crIsoRhoFiltered0p15,hltL3crIsoL1sMu12Eta2p1L1f0L2f12QL3f15QL3crIsoFiltered10',
'HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v.*:hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15,hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4',
'HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk1_eta2p1_v.*:hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15,hltOverlapFilterIsoMu18PFTau25TrackPt1Prong4'
),

RecPhoton = cms.untracked.bool(True),
RecPhotonHLTriggerMatching = cms.untracked.vstring(),
RecPhotonPtMin = cms.untracked.double(0.),
RecPhotonEtaMax = cms.untracked.double(0.),
RecPhotonNum = cms.untracked.int32(100000),

RecElectron = cms.untracked.bool(True),
RecElectronHLTriggerMatching = cms.untracked.vstring(
'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.*:hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter,hltEle27CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter', # 5e32
'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.*:hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter,hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter', # 1e33
'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.*:hltEle42CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter,hltEle42CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter', # 1e33
'HLT_Ele52_CaloIdVT_TrkIdT_v.*:hltEle52CaloIdVTTrkIdTDphiFilter', # 1.4e33
'HLT_Ele65_CaloIdVT_TrkIdT_v.*:hltEle65CaloIdVTTrkIdTDphiFilter', # 2e33, 3e33
'HLT_Ele80_CaloIdVT_TrkIdT_v.*:hltEle80CaloIdVTTrkIdTDphiFilter', # 5e33
'HLT_Ele80_CaloIdVT_GsfTrkIdT_v.*:hltEle80CaloIdVTGsfTrkIdTDPhiFilter', # 7e33, 8e33
'HLT_Ele90_CaloIdVT_GsfTrkIdT_v.*:hltEle90CaloIdVTGsfTrkIdTDPhiFilter', # 7e33, 8e33
'HLT_Ele27_WP80_v.*:hltEle27WP80TrackIsoFilter', # unprescaled since 2012
'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v.*:hltEle15CaloIdVTTrkIdTDphiFilter,hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter,hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter',
'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v.*:hltEle15CaloIdVTTrkIdTDphiFilter,hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter,hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter',
'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v.*:hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter', # only for very short period
'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v.*:hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter',
'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v.*:hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter', # MC Trigger
'HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v.*:hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20',
'HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v.*:hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1IsoEG18OrEG20',
'HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v.*:hltEle22WP90RhoTrackIsoFilter'
),
RecElectronPtMin = cms.untracked.double(0.),
RecElectronTrackIso = cms.untracked.double(1000000.),
RecElectronEta = cms.untracked.double(0.),
RecElectronNum = cms.untracked.int32(100000),

RecTau = cms.untracked.bool(True),
RecTauHLTriggerMatching = cms.untracked.vstring(
'HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu18LooseIsoPFTau20',
'HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v.*:hltOverlapFilterIsoMu17LooseIsoPFTau20',
'HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v.*:hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4',
'HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk1_eta2p1_v.*:hltOverlapFilterIsoMu18PFTau25TrackPt1Prong4',
'HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle20LooseIsoPFTau20',
'HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v.*:hltOverlapFilterIsoEle20WP90LooseIsoPFTau20',
'HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_v.*:hltDoublePFTau35TrackPt1MediumIsolationProng4Dz02',
'HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_v.*:hltDoublePFTau35TrackPt5MediumIsolationProng4Dz02',
'HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v.*:hltDoublePFTau30TrackPt1MediumIsolationProng4Dz02',
'HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v.*:hltDoublePFTau25TrackPt5MediumIsolationProng4Dz02',
'HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v.*:hltDoublePFTau30TrackPt5MediumIsolationProng4Dz02',
),

RecTauDiscriminators = cms.untracked.vstring(
'decayModeFinding',
'decayModeFindingNewDMs',
'decayModeFindingOldDMs',
#'againstElectronVLooseMVA2',
'againstElectronLoose',
#'againstElectronLooseMVA2',
#'againstElectronLooseMVA3',
#'againstElectronMVA',
'againstElectronMedium',
#'againstElectronMediumMVA2',
#'againstElectronMediumMVA3',
'againstElectronTight',
#'againstElectronTightMVA2',
#'againstElectronTightMVA3',
#'againstElectronVTightMVA3',
'againstElectronDeadECAL',
'againstElectronLooseMVA5',
#'againstElectronMVA5category',
#'againstElectronMVA5raw',
'againstElectronMediumMVA5',
'againstElectronTightMVA5',
'againstElectronVLooseMVA5',
'againstElectronVTightMVA5',
'againstMuonLoose',
'againstMuonMedium',
'againstMuonTight',
'againstMuonLoose2',
'againstMuonMedium2',
'againstMuonTight2',
'againstMuonLooseMVA',
#'againstMuonMVAraw',
'againstMuonMediumMVA',
'againstMuonTight3',
'againstMuonTightMVA',
#'byIsolationMVA3newDMwLTraw',
#'byIsolationMVA3newDMwoLTraw',
#'byIsolationMVA3oldDMwLTraw',
#'byIsolationMVA3oldDMwoLTraw',
'byLooseIsolationMVA3newDMwLT',
'byLooseIsolationMVA3newDMwoLT',
'byLooseIsolationMVA3oldDMwLT',
'byLooseIsolationMVA3oldDMwoLT',
'byMediumIsolationMVA3newDMwLT',
'byMediumIsolationMVA3newDMwoLT',
'byMediumIsolationMVA3oldDMwLT',
'byMediumIsolationMVA3oldDMwoLT',
'byTightIsolationMVA3newDMwLT',
'byTightIsolationMVA3newDMwoLT',
'byTightIsolationMVA3oldDMwLT',
'byTightIsolationMVA3oldDMwoLT',
'byVLooseIsolationMVA3newDMwLT',
'byVLooseIsolationMVA3newDMwoLT',
'byVLooseIsolationMVA3oldDMwLT',
'byVLooseIsolationMVA3oldDMwoLT',
'byVTightIsolationMVA3newDMwLT',
'byVTightIsolationMVA3newDMwoLT',
'byVTightIsolationMVA3oldDMwLT',
'byVTightIsolationMVA3oldDMwoLT',
'byVVTightIsolationMVA3newDMwLT',
'byVVTightIsolationMVA3newDMwoLT',
'byVVTightIsolationMVA3oldDMwLT',
'byVVTightIsolationMVA3oldDMwoLT',
'byLooseCombinedIsolationDeltaBetaCorr',
'byLooseCombinedIsolationDeltaBetaCorr3Hits',
'byLooseIsolation',
#'byLooseIsolationDeltaBetaCorr',
#'byLooseIsolationMVA',
#'byLooseIsolationMVA2',
'byMediumCombinedIsolationDeltaBetaCorr',
'byMediumCombinedIsolationDeltaBetaCorr3Hits',
#'byMediumIsolation',
#'byMediumIsolationDeltaBetaCorr',
#'byMediumIsolationMVA',
#'byMediumIsolationMVA2',
'byTightCombinedIsolationDeltaBetaCorr',
'byTightCombinedIsolationDeltaBetaCorr3Hits',
#'byTightIsolation',
#'byTightIsolationDeltaBetaCorr',
#'byTightIsolationMVA',
#'byTightIsolationMVA2',
'byVLooseCombinedIsolationDeltaBetaCorr',
#'byVLooseIsolation',
#'byVLooseIsolationDeltaBetaCorr'
),

RecTauTauPairs = cms.untracked.bool(False),
RecMuTauTauPairs = cms.untracked.bool(False),
RecElTauTauPairs = cms.untracked.bool(False),

RecAK5CaloJet = cms.untracked.bool(False),
RecAK5CaloPtMin = cms.untracked.double(20.),
RecAK5CaloEtaMax = cms.untracked.double(2.4),
RecAK5CaloNum = cms.untracked.int32(100000),
RecAK5CaloFilterPtMin = cms.untracked.double(20.),

RecAK5JPTJet = cms.untracked.bool(False),
RecAK5JPTPtMin = cms.untracked.double(20.),
RecAK5JPTEtaMax = cms.untracked.double(2.4),
RecAK5JPTNum = cms.untracked.int32(100000),
RecAK5JPTFilterPtMin = cms.untracked.double(20.),

RecAK5PFJet = cms.untracked.bool(True),
RecJetHLTriggerMatching = cms.untracked.vstring(
'HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v.*:hltTripleL2Jets30eta3',
'HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v.*:hltTripleL2Jets30eta3',
'HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v.*:hltTripleL2Jets30eta3',
),

RecAK5PFPtMin = cms.untracked.double(20.),
RecAK5PFEtaMax = cms.untracked.double(2.4),
RecAK5PFNum = cms.untracked.int32(100000),
RecAK5PFFilterPtMin = cms.untracked.double(10.),
## only for MuOnia
RecSecVertices = cms.untracked.bool(False),
RecVertexTRKChi2 = cms.untracked.double(5),
RecVertexTRKHitsMin = cms.untracked.int32(6),
RecVertexChi2 = cms.untracked.double(3),
RecVertexSig2D = cms.untracked.double(15),
RecKaonMasswin = cms.untracked.double(0.05),
RecLambdaMasswin = cms.untracked.double(0.02),

RecMuSecVertices = cms.untracked.bool(False),
RecAllConversion = cms.untracked.bool(False),
)
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('GluGluToHToTauTau_M-125_full.root')
)

# Re-run tau ID algorithms so that discrimintors required for PAT taus are
# created which might not exist in AOD (yet).
process.tau_step = cms.Path(process.PFTau)

# PU Jet ID
process.pileupJetIdProducer.residualsTxt = 'MyRootMaker/MyRootMaker/test/RootTree.py'
process.pileupJetIdProducerChs.residualsTxt = 'MyRootMaker/MyRootMaker/test/RootTree.py'
process.pileupJetIdProducer.jets = cms.InputTag("ak5PFJets")
process.pileupJetIdProducer.applyJec = cms.bool(True)
process.pileupJetIdProducer.inputIsCorrected = cms.bool(False)

# Rho
#process.kt6PFJetsAK5PF.Rho_EtaMax = cms.double(4.4)
process.kt6PFJets.Rho_EtaMax = cms.double(4.4)

if doSVFit:
	process.svfitSequence = cms.Sequence(process.cleanPatTausForSVfit*process.diTau*process.selectedDiTau*process.selectedDiTauFilter)
else:
	process.svfitSequence = cms.Sequence(process.cleanPatTausForSVfit)
	
process.p = cms.Path(
	# Jet algo for PU correction
	process.kt6PFJets*
	process.ak5PFJets*
	# MVA electron ID
	process.mvaTrigV0*
	process.mvaNonTrigV0*
	# Lepton isolation:
	process.pfNoPileUpSequence *
	process.pfParticleSelectionSequence *
	process.pfAllMuons *
	process.pfMuonsFromVertex *
	process.pfSelectedMuons *
	process.isoSequence *
	# PAT objects:
	process.muIsoSequence *
	process.elIsoSequence *
	process.patDefaultSequence*
	# MVA MET:
	process.isolatedPatElectrons*
	process.isolatedPatMuons*
	process.isolatedPatTaus*
	process.calibratedAK5PFJetsForPFMEtMVA*
	process.pfMEtMVA*
	process.patPFMetByMVA*
	process.pfMEtMVACov*
	# Pileup Jet ID
	process.pileupJetIdProducer*
	# Type-1 MET:
	process.producePFMETCorrections*
	process.patMETsPFType1*
	# SVfit:
	process.svfitSequence*
	# The analyzer
	process.makeroottree
)
