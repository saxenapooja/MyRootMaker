#include "MyRootMaker/MyRootMaker/interface/RootMaker.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
//#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisFwd.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisBaseFwd.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisByIntegration.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesisBase.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesis.h"
#include <DataFormats/RecoCandidate/interface/IsoDepositVetos.h>
//#include "AnalysisDataFormats/TauAnalysis/interface/PFMEtSignCovMatrix.h"
//#include "RecoJets/JetProducers/interface/PileupJetIdentifier.h"
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

using namespace reco;

typedef std::vector<NSVfitEventHypothesisByIntegration> NSVfitEventHypothesisByIntegrationCollection;

// http://root.cern.ch/root/html/ROOT__Math__SMatrix_double_3_3_-p1MatRepSym_double_3___.html#ROOT__Math__SMatrix_double_3_3_-p1MatRepSym_double_3___:_SMatrix_double_3_3_ROOT__Math__MatRepSym_double_3___
// http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/SVectorDoc.html
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;  // Standard Matrix representation for a general 3x3 matrix of type double.
typedef ROOT::Math::SVector<double, 3> SVector3; //// SVector: vector of size 3 

static const unsigned int SKIM_MUTAUTAU   = (1 << 0);     //1  : mu+tau+tau
static const unsigned int SKIM_ETAUTAU    = (1 << 1);     //2  : e+tau+tau  
static const unsigned int SKIM_MUMU       = (1 << 2);     //4  : mu+mu
static const unsigned int SKIM_EE         = (1 << 3);     //8  : e+e
static const unsigned int SKIM_ETAU       = (1 << 4);     //16 : e+tau
static const unsigned int SKIM_ALL        = (1 << 5);     //32 : all
static const unsigned int SKIM_EMU        = (1 << 6);     //64 : e+mu
static const unsigned int SKIM_TAUTAU     = (1 << 7);     //128: tau+tau


//to set the values from parameter
RootMaker::RootMaker(const edm::ParameterSet& iConfig) :  
  cdata(iConfig.getUntrackedParameter<bool>("IsData", false)),
  cgen(iConfig.getUntrackedParameter<bool>("GenSomeParticles", false)),
  cgenallparticles(iConfig.getUntrackedParameter<bool>("GenAllParticles", false)),
  ctrigger(iConfig.getUntrackedParameter<bool>("Trigger", false)),
  cbeamspot(iConfig.getUntrackedParameter<bool>("BeamSpot", false)),
  crectrack(iConfig.getUntrackedParameter<bool>("RecTrack", false)),
  crecprimvertex(iConfig.getUntrackedParameter<bool>("RecPrimVertex", false)),
  crecsupercluster(iConfig.getUntrackedParameter<bool>("RecSuperCluster", false)),
  crecsuperclustermember(iConfig.getUntrackedParameter<bool>("RecSuperClusterBasicCluster", false)),
  crecsuperclusterhit(iConfig.getUntrackedParameter<bool>("RecSuperClusterHit", false)),
  crecmuon(iConfig.getUntrackedParameter<bool>("RecMuon", false)),
  crectau(iConfig.getUntrackedParameter<bool>("RecTau", false)),
  crecmutautaupairs(iConfig.getUntrackedParameter<bool>("RecMuTauTauPairs", false)),
  creceltautaupairs(iConfig.getUntrackedParameter<bool>("RecElTauTauPairs", false)),
  crecelectron(iConfig.getUntrackedParameter<bool>("RecElectron", false)),
  crecphoton(iConfig.getUntrackedParameter<bool>("RecPhoton", false)),
  crecallconversion(iConfig.getUntrackedParameter<bool>("RecAllConversion", false)),
  crecak5calojet(iConfig.getUntrackedParameter<bool>("RecAK5CaloJet", false)),
  crecak5jptjet(iConfig.getUntrackedParameter<bool>("RecAK5JPTJet", false)),
  crecak5pfjet(iConfig.getUntrackedParameter<bool>("RecAK5PFJet", false)),
  creccalomet(iConfig.getUntrackedParameter<bool>("RecCaloMet", false)),
  crectcmet(iConfig.getUntrackedParameter<bool>("RecTCMet", false)),
  crecpfmet(iConfig.getUntrackedParameter<bool>("RecPFMet", false)),
  crecsecvertices(iConfig.getUntrackedParameter<bool>("RecSecVertices", false)),
  crecmusecvertices(iConfig.getUntrackedParameter<bool>("RecMuSecVertices", false)),
  cHLTriggerNamesSelection(iConfig.getUntrackedParameter<vector<string> >("HLTriggerSelection")),
  cTriggerProcess(iConfig.getUntrackedParameter<string>("TriggerProcess", "HLT")),
  cMuPtMin(iConfig.getUntrackedParameter<double>("RecMuonPtMin", 0.)),
  cMuTrackIso(iConfig.getUntrackedParameter<double>("RecMuonTrackIso", 100000.)),
  cMuEtaMax(iConfig.getUntrackedParameter<double>("RecMuonEtaMax", 1000000.)),
  cMuHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecMuonHLTriggerMatching")),
  cMuNum(iConfig.getUntrackedParameter<int>("RecMuonNum", 0)),
  cElPtMin(iConfig.getUntrackedParameter<double>("RecElectronPtMin", 0.)),
  cElTrackIso(iConfig.getUntrackedParameter<double>("RecElectronTrackIso", 1000000.)),
  cElEtaMax(iConfig.getUntrackedParameter<double>("RecElectronEtaMax", 1000000.)),
  cElHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecElectronHLTriggerMatching")),
  cElNum(iConfig.getUntrackedParameter<int>("RecElectronNum", 0)),
  cTauHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecTauHLTriggerMatching")),
  cTauDiscriminators(iConfig.getUntrackedParameter<vector<string> >("RecTauDiscriminators")),
  cTrackFilterPtMin(iConfig.getUntrackedParameter<double>("RecTrackFilterPtMin", 0.)),
  cTrackPtMin(iConfig.getUntrackedParameter<double>("RecTrackPtMin", 0.)),
  cTrackEtaMax(iConfig.getUntrackedParameter<double>("RecTrackEtaMax", 1000000.)),
  cTrackNum(iConfig.getUntrackedParameter<int>("RecTrackNum", 0)),
  cPhotonPtMin(iConfig.getUntrackedParameter<double>("RecPhotonPtMin", 0.)),
  cPhotonEtaMax(iConfig.getUntrackedParameter<double>("RecPhotonEtaMax", 1000000.)),
  cPhotonHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecPhotonHLTriggerMatching")),
  cPhotonNum(iConfig.getUntrackedParameter<int>("RecPhotonNum", 0)),
  cAK5CaloFilterPtMin(iConfig.getUntrackedParameter<double>("RecAK5CaloFilterPtMin", 0.)),
  cAK5CaloPtMin(iConfig.getUntrackedParameter<double>("RecAK5CaloPtMin", 0.)),
  cAK5CaloEtaMax(iConfig.getUntrackedParameter<double>("RecAK5CaloEtaMax", 1000000.)),
  cAK5CaloNum(iConfig.getUntrackedParameter<int>("RecAK5CaloNum", 0)),
  cAK5JPTFilterPtMin(iConfig.getUntrackedParameter<double>("RecAK5JPTFilterPtMin", 0.)),
  cAK5JPTPtMin(iConfig.getUntrackedParameter<double>("RecAK5JPTPtMin", 0.)),
  cAK5JPTEtaMax(iConfig.getUntrackedParameter<double>("RecAK5JPTEtaMax", 1000000.)),
  cAK5JPTNum(iConfig.getUntrackedParameter<int>("RecAK5JPTNum", 0)),
  cAK5PFFilterPtMin(iConfig.getUntrackedParameter<double>("RecAK5PFFilterPtMin", 0.)),
  cAK5PFPtMin(iConfig.getUntrackedParameter<double>("RecAK5PFPtMin", 0.)),
  cAK5PFEtaMax(iConfig.getUntrackedParameter<double>("RecAK5PFEtaMax", 1000000.)),
  cJetHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecJetHLTriggerMatching")),
  cAK5PFNum(iConfig.getUntrackedParameter<int>("RecAK5PFNum", 0)),
  cVertexTRKChi2(iConfig.getUntrackedParameter<double>("RecVertexTRKChi2", 10.)),
  cVertexTRKHitsMin(iConfig.getUntrackedParameter<int>("RecVertexTRKHitsMin", 6)),
  cVertexChi2(iConfig.getUntrackedParameter<double>("RecVertexChi2", 5.)),
  cVertexSig2D(iConfig.getUntrackedParameter<double>("RecVertexSig2D", 15.)),
  cKaonMassWindow(iConfig.getUntrackedParameter<double>("RecVertexKaonMassWin", 0.05)),
  cLambdaMassWindow(iConfig.getUntrackedParameter<double>("RecVertexLambdaMassWin", 0.02)),
  cYear(iConfig.getUntrackedParameter<unsigned int>("Year")),
  cPeriod(iConfig.getUntrackedParameter<std::string>("Period")),
  cSkim(iConfig.getUntrackedParameter<unsigned int>("Skim")),
  TrackCollectionTag_(iConfig.getParameter<edm::InputTag>("TrackCollectionTag")),
  PVTag_(iConfig.getParameter<edm::InputTag>("PVTag")),
  propagatorWithMaterial(0)
{
  if(cYear != 2011 && cYear != 2012)
    throw cms::Exception("RootMaker") << "Invalid Year, only 2011 and 2012 are allowed!";
  if(cPeriod != "Summer11" && cPeriod != "Fall11" && cPeriod != "Summer12")
    throw cms::Exception("RootMaker") << "Invalid period, only Summer11, Fall11 and Summer12 are allowed!";
  
  //vector<string>
  bdisclabel.push_back("trackCountingHighPurBJetTags");          //trackCount_highPure_BJets_tags
  bdisclabel.push_back("trackCountingHighEffBJetTags");          //trackCount_highEff_BJets_tags
  bdisclabel.push_back("combinedSecondaryVertexBJetTags");       //combined_secVtx_BJets_tags
  bdisclabel.push_back("combinedSecondaryVertexMVABJetTags");    //combined_secVtx_MVA_BJets_tags
  bdisclabel.push_back("simpleSecondaryVertexHighPurBJetTags");  //simple_secVtx_highPure_BJets_tags
  bdisclabel.push_back("simpleSecondaryVertexHighEffBJetTags");  //simple_secVtx_highEff_BJets_tags 
  
  double barrelRadius = 129.;  //p81, p50, ECAL TDR
  double endcapZ      = 320.5; // fig 3.26, p81, ECAL TDR
  Surface::RotationType rot;
  ecalBarrel         = Cylinder::build(Surface::PositionType(0, 0, 0), rot, barrelRadius);
  ecalNegativeEtaEndcap = Plane::build(Surface::PositionType(0, 0, -endcapZ), rot);
  ecalPositiveEtaEndcap = Plane::build(Surface::PositionType(0, 0, endcapZ), rot);
  
  const char* cmssw_base = getenv("CMSSW_BASE");
  if(!cmssw_base) throw cms::Exception("No CMSSW runtime environment found");
  std::string prefix = std::string(cmssw_base) + "/src/";

  //Muon
  std::vector<std::string> muonid_weightfiles;
  std::vector<std::string> muoniso_weightfiles;
  muonid_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Tracker_V0_BDTG.weights.xml");
  muonid_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Global_V0_BDTG.weights.xml");

  muoniso_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
  muoniso_weightfiles.push_back(prefix + "Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");

  fMuonIDMVA.initialize("MuonID_BDTG"           , MuonMVAEstimator::kID      , true, muonid_weightfiles );
  fMuonIsoMVA.initialize("MuonIso_BDTG_IsoRings", MuonMVAEstimator::kIsoRings, true, muoniso_weightfiles);
  //fMuonIsoMVA.SetPrintMVADebug(kTRUE);

  //electron
  vector<string> eleiso_weightfiles;
  //eleiso_weightfiles.push_back("prefix + UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIsoMVA_V0_BarrelPt5To10_V0_BDTG.weights.xml");
  //eleiso_weightfiles.push_back("prefix + UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIsoMVA_V0_EndcapPt5To10_V0_BDTG.weights.xml");
  //eleiso_weightfiles.push_back("prefix + UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIsoMVA_V0_BarrelPt10ToInf_V0_BDTG.weights.xml");
  //eleiso_weightfiles.push_back("prefix + UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIsoMVA_V0_EndcapPt10ToInf_V0_BDTG.weights.xml");
  eleiso_weightfiles.push_back(prefix + "UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_BarrelPt5To10.weights.xml");
  eleiso_weightfiles.push_back(prefix + "UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_EndcapPt5To10.weights.xml");
  eleiso_weightfiles.push_back(prefix + "UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_BarrelPt10ToInf.weights.xml");
  eleiso_weightfiles.push_back(prefix + "UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_EndcapPt10ToInf.weights.xml");
  
  /*fElectronIsoMVA.initialize("EleIso_BDTG_IsoRings",
    EGammaMvaEleEstimator::kIsoRings,
    kTRUE,
    eleiso_weightfiles);
    fElectronIsoMVA.SetPrintMVADebug(kFALSE);*/
}//RootMaker::RootMaker(const edm::ParameterSet& iConfig)


//destructor
RootMaker::~RootMaker(){
  if(propagatorWithMaterial != 0){ delete propagatorWithMaterial;}
  
}


void RootMaker::beginJob(){
  edm::Service<TFileService> FS;
  tree = FS->make<TTree>("AC1B", "AC1B", 1);
  nEvents = FS->make<TH1D>("nEvents", "nEvents", 2, -0.5, +1.5);
  
  tree->Branch("errors", &errors, "errors/i");
  tree->Branch("event_nr", &event_nr, "event_nr/i");
  tree->Branch("event_run", &event_run, "event_run/i");
  tree->Branch("event_timeunix", &event_timeunix, "event_timeunix/i");
  tree->Branch("event_timemicrosec", &event_timemicrosec, "event_timemicrosec/i");
  tree->Branch("event_luminosityblock", &event_luminosityblock, "event_luminosityblock/i");
  tree->Branch("trigger_level1bits", &trigger_level1bits, "trigger_level1bits[8]/b");
  tree->Branch("trigger_level1", &trigger_level1, "trigger_level1[128]/b");
  tree->Branch("trigger_HLT", &trigger_HLT, "trigger_HLT[128]/b");
  
  tree->Branch("beamspot_x", &beamspot_x, "beamspot_x/F");
  tree->Branch("beamspot_y", &beamspot_y, "beamspot_y/F");
  tree->Branch("beamspot_z", &beamspot_z, "beamspot_z/F");
  tree->Branch("beamspot_xwidth", &beamspot_xwidth, "beamspot_xwidth/F");
  tree->Branch("beamspot_ywidth", &beamspot_ywidth, "beamspot_ywidth/F");
  tree->Branch("beamspot_zsigma", &beamspot_zsigma, "beamspot_zsigma/F");
  tree->Branch("beamspot_cov", &beamspot_cov, "beamspot_cov[6]/F");
  
  tree->Branch("track_count", &track_count, "track_count/i"); 
  tree->Branch("track_px", track_px, "track_px[track_count]/F");
  tree->Branch("track_py", track_py, "track_py[track_count]/F");
  tree->Branch("track_pz", track_pz, "track_pz[track_count]/F");
  tree->Branch("track_outerx", track_outerx, "track_outerx[track_count]/F");
  tree->Branch("track_outery", track_outery, "track_outery[track_count]/F");
  tree->Branch("track_outerz", track_outerz, "track_outerz[track_count]/F");
  tree->Branch("track_closestpointx", track_closestpointx, "track_closestpointx[track_count]/F");
  tree->Branch("track_closestpointy", track_closestpointy, "track_closestpointy[track_count]/F");
  tree->Branch("track_closestpointz", track_closestpointz, "track_closestpointz[track_count]/F");
  tree->Branch("track_chi2", track_chi2, "track_chi2[track_count]/F");
  tree->Branch("track_ndof", track_ndof, "track_ndof[track_count]/F");
  tree->Branch("track_dxy", track_dxy, "track_dxy[track_count]/F");
  tree->Branch("track_dxyerr", track_dxyerr, "track_dxyerr[track_count]/F");
  tree->Branch("track_dz", track_dz, "track_dz[track_count]/F");
  tree->Branch("track_dzerr", track_dzerr, "track_dzerr[track_count]/F");
  tree->Branch("track_dedxharmonic2", track_dedxharmonic2, "track_dedxharmonic2[track_count]/F");
  tree->Branch("track_charge", track_charge, "track_charge[track_count]/I");
  tree->Branch("track_nhits", track_nhits, "track_nhits[track_count]/b");
  tree->Branch("track_npixelhits", track_npixelhits, "track_npixelhits[track_count]/b");
  tree->Branch("track_nmissinghits", track_nmissinghits, "track_nmissinghits[track_count]/b");
  tree->Branch("track_npixellayers", track_npixellayers, "track_npixellayers[track_count]/b");
  tree->Branch("track_nstriplayers", track_nstriplayers, "track_nstriplayers[track_count]/b");
  
  tree->Branch("primvertex_count", &primvertex_count, "primvertex_count/i"); 
  tree->Branch("primvertex_x", &primvertex_x, "primvertex_x/F");
  tree->Branch("primvertex_y", &primvertex_y, "primvertex_y/F");
  tree->Branch("primvertex_z", &primvertex_z, "primvertex_z/F");
  tree->Branch("primvertex_chi2", &primvertex_chi2, "primvertex_chi2/F");
  tree->Branch("primvertex_ndof", &primvertex_ndof, "primvertex_ndof/F");
  tree->Branch("primvertex_ptq", &primvertex_ptq, "primvertex_pdf/F");
  tree->Branch("primvertex_ntracks", &primvertex_ntracks, "primvertex_ntracks/I");
  tree->Branch("primvertex_cov", primvertex_cov, "primvertex_cov[6]/F");
  
  tree->Branch("supercluster_count", &supercluster_count, "supercluster_count/i");
  tree->Branch("supercluster_e", supercluster_e, "supercluster_e[supercluster_count]/F");
  tree->Branch("supercluster_x", supercluster_x, "supercluster_x[supercluster_count]/F");
  tree->Branch("supercluster_y", supercluster_y, "supercluster_y[supercluster_count]/F");
  tree->Branch("supercluster_z", supercluster_z, "supercluster_z[supercluster_count]/F");
  tree->Branch("supercluster_rawe", supercluster_rawe, "supercluster_rawe[supercluster_count]/F");
  tree->Branch("supercluster_phiwidth", supercluster_phiwidth, "supercluster_phiwidth[supercluster_count]/F");
  tree->Branch("supercluster_etawidth", supercluster_etawidth, "supercluster_etawidth[supercluster_count]/F");
  tree->Branch("supercluster_nbasiccluster", supercluster_nbasiccluster, "supercluster_nbasiccluster[supercluster_count]/I");
  tree->Branch("supercluster_basicclusterbegin", supercluster_basicclusterbegin, "supercluster_basicclusterbegin[supercluster_count]/I");
  tree->Branch("supercluster_esclusterbegin", supercluster_esclusterbegin, "supercluster_esclusterbegin[supercluster_count]/I");
  
  tree->Branch("supercluster_basiccluster_count", &supercluster_basiccluster_count, "supercluster_basiccluster_count/i");
  tree->Branch("supercluster_basiccluster_e", supercluster_basiccluster_e, "supercluster_basiccluster_e[supercluster_basiccluster_count]/F");
  tree->Branch("supercluster_basiccluster_x", supercluster_basiccluster_x, "supercluster_basiccluster_x[supercluster_basiccluster_count]/F");
  tree->Branch("supercluster_basiccluster_y", supercluster_basiccluster_y, "supercluster_basiccluster_y[supercluster_basiccluster_count]/F");
  tree->Branch("supercluster_basiccluster_z", supercluster_basiccluster_z, "supercluster_basiccluster_z[supercluster_basiccluster_count]/F");
  tree->Branch("supercluster_basiccluster_nhit", supercluster_basiccluster_nhit, "supercluster_basiccluster_nhit[supercluster_basiccluster_count]/I");
  tree->Branch("supercluster_basiccluster_hitbegin", supercluster_basiccluster_hitbegin, "supercluster_basiccluster_hitbegin[supercluster_basiccluster_count]/I");
  
  tree->Branch("supercluster_basiccluster_hit_count", &supercluster_basiccluster_hit_count, "supercluster_basiccluster_hit_count/i");
  tree->Branch("supercluster_basiccluster_hit_e", supercluster_basiccluster_hit_e, "supercluster_basiccluster_hit_e[supercluster_basiccluster_hit_count]/F");
  tree->Branch("supercluster_basiccluster_hit_x", supercluster_basiccluster_hit_x, "supercluster_basiccluster_hit_x[supercluster_basiccluster_hit_count]/F");
  tree->Branch("supercluster_basiccluster_hit_y", supercluster_basiccluster_hit_y, "supercluster_basiccluster_hit_y[supercluster_basiccluster_hit_count]/F");
  tree->Branch("supercluster_basiccluster_hit_z", supercluster_basiccluster_hit_z, "supercluster_basiccluster_hit_z[supercluster_basiccluster_hit_count]/F");
  
  tree->Branch("supercluster_escluster_count", &supercluster_escluster_count, "supercluster_escluster_count/i");
  tree->Branch("supercluster_escluster_e", supercluster_escluster_e, "supercluster_escluster_e[supercluster_escluster_count]/F");
  tree->Branch("supercluster_escluster_x", supercluster_escluster_x, "supercluster_escluster_x[supercluster_escluster_count]/F");
  tree->Branch("supercluster_escluster_y", supercluster_escluster_y, "supercluster_escluster_y[supercluster_escluster_count]/F");
  tree->Branch("supercluster_escluster_z", supercluster_escluster_z, "supercluster_escluster_z[supercluster_escluster_count]/F");
  tree->Branch("supercluster_escluster_nhit", supercluster_escluster_nhit, "supercluster_escluster_nhit[supercluster_escluster_count]/I");
  tree->Branch("supercluster_escluster_hitbegin", supercluster_escluster_hitbegin, "supercluster_escluster_hitbegin[supercluster_escluster_count]/I");
  
  tree->Branch("supercluster_escluster_hit_count", &supercluster_escluster_hit_count, "supercluster_escluster_hit_count/i");
  tree->Branch("supercluster_escluster_hit_e", supercluster_escluster_hit_e, "supercluster_escluster_hit_e[supercluster_escluster_hit_count]/F");
  tree->Branch("supercluster_escluster_hit_x", supercluster_escluster_hit_x, "supercluster_escluster_hit_x[supercluster_escluster_hit_count]/F");
  tree->Branch("supercluster_escluster_hit_y", supercluster_escluster_hit_y, "supercluster_escluster_hit_y[supercluster_escluster_hit_count]/F");
  tree->Branch("supercluster_escluster_hit_z", supercluster_escluster_hit_z, "supercluster_escluster_hit_z[supercluster_escluster_hit_count]/F");
  
  tree->Branch("muon_count", &muon_count, "muon_count/i");
  tree->Branch("muon_px", muon_px, "muon_px[muon_count]/F");
  tree->Branch("muon_py", muon_py, "muon_py[muon_count]/F");
  tree->Branch("muon_pz", muon_pz, "muon_pz[muon_count]/F");
  tree->Branch("muon_pterror", muon_pterror, "muon_pterror[muon_count]/F");
  tree->Branch("muon_chi2", muon_chi2, "muon_chi2[muon_count]/F");
  tree->Branch("muon_ndof", muon_ndof, "muon_ndof[muon_count]/F");
  tree->Branch("muon_nchambers", muon_nchambers, "muon_nchambers[muon_count]/i");
  tree->Branch("muon_innertrack_px", muon_innertrack_px, "muon_innertrack_px[muon_count]/F");
  tree->Branch("muon_innertrack_py", muon_innertrack_py, "muon_innertrack_py[muon_count]/F");
  tree->Branch("muon_innertrack_pz", muon_innertrack_pz, "muon_innertrack_pz[muon_count]/F");
  tree->Branch("muon_innertrack_outerx", muon_innertrack_outerx, "muon_innertrack_outerx[muon_count]/F");
  tree->Branch("muon_innertrack_outery", muon_innertrack_outery, "muon_innertrack_outery[muon_count]/F");
  tree->Branch("muon_innertrack_outerz", muon_innertrack_outerz, "muon_innertrack_outerz[muon_count]/F");
  tree->Branch("muon_innertrack_closestpointx", muon_innertrack_closestpointx, "muon_innertrack_closestpointx[muon_count]/F");
  tree->Branch("muon_innertrack_closestpointy", muon_innertrack_closestpointy, "muon_innertrack_closestpointy[muon_count]/F");
  tree->Branch("muon_innertrack_closestpointz", muon_innertrack_closestpointz, "muon_innertrack_closestpointz[muon_count]/F");
  tree->Branch("muon_innertrack_chi2", muon_innertrack_chi2, "muon_innertrack_chi2[muon_count]/F");
  tree->Branch("muon_innertrack_ndof", muon_innertrack_ndof, "muon_innertrack_ndof[muon_count]/F");
  tree->Branch("muon_innertrack_dxy", muon_innertrack_dxy, "muon_innertrack_dxy[muon_count]/F");
  tree->Branch("muon_innertrack_dxyerr", muon_innertrack_dxyerr, "muon_innertrack_dxyerr[muon_count]/F");
  tree->Branch("muon_innertrack_dz", muon_innertrack_dz, "muon_innertrack_dz[muon_count]/F");
  tree->Branch("muon_innertrack_dzerr", muon_innertrack_dzerr, "muon_innertrack_dzerr[muon_count]/F");
  tree->Branch("muon_innertrack_dedxharmonic2", muon_innertrack_dedxharmonic2, "muon_innertrack_dedxharmonic2[muon_count]/F");
  tree->Branch("muon_innertrack_charge", muon_innertrack_charge, "muon_innertrack_charge[muon_count]/I");
  tree->Branch("muon_innertrack_nhits", muon_innertrack_nhits, "muon_innertrack_nhits[muon_count]/b");
  tree->Branch("muon_innertrack_npixelhits", muon_innertrack_npixelhits, "muon_innertrack_npixelhits[muon_count]/b");
  tree->Branch("muon_innertrack_nmissinghits", muon_innertrack_nmissinghits, "muon_innertrack_nmissinghits[muon_count]/b");
  tree->Branch("muon_innertrack_npixellayers", muon_innertrack_npixellayers, "muon_innertrack_npixellayers[muon_count]/b");
  tree->Branch("muon_innertrack_nstriplayers", muon_innertrack_nstriplayers, "muon_innertrack_nstriplayers[muon_count]/b");
  tree->Branch("muon_outertrack_px", muon_outertrack_px, "muon_outertrack_px[muon_count]/F");
  tree->Branch("muon_outertrack_py", muon_outertrack_py, "muon_outertrack_py[muon_count]/F");
  tree->Branch("muon_outertrack_pz", muon_outertrack_pz, "muon_outertrack_pz[muon_count]/F");
  tree->Branch("muon_outertrack_hits", muon_outertrack_hits, "muon_outertrack_hits[muon_count]/b");
  tree->Branch("muon_outertrack_missinghits", muon_outertrack_missinghits, "muon_outertrack_missinghits[muon_count]/b");
  tree->Branch("muon_outertrack_chi2", muon_outertrack_chi2, "muon_outertrack_chi2[muon_count]/F");
  tree->Branch("muon_outertrack_ndof", muon_outertrack_ndof, "muon_outertrack_ndof[muon_count]/F");
  tree->Branch("muon_isolationr3track", muon_isolationr3track, "muon_isolationr3track[muon_count]/F");
  tree->Branch("muon_isolationr3ntrack", muon_isolationr3ntrack, "muon_isolationr3ntrack[muon_count]/I");
  tree->Branch("muon_isolationr3ecal", muon_isolationr3ecal, "muon_isolationr3ecal[muon_count]/F");
  tree->Branch("muon_isolationr3hcal", muon_isolationr3hcal, "muon_isolationr3hcal[muon_count]/F");
  tree->Branch("muon_pfisolationr3_sumchargedhadronpt", muon_pfisolationr3_sumchargedhadronpt, "muon_pfisolationr3_sumchargedhadronpt[muon_count]/F");
  tree->Branch("muon_pfisolationr3_sumchargedparticlept", muon_pfisolationr3_sumchargedparticlept, "muon_pfisolationr3_sumchargedparticlept[muon_count]/F");
  tree->Branch("muon_pfisolationr3_sumneutralhadronet", muon_pfisolationr3_sumneutralhadronet, "muon_pfisolationr3_sumneutralhadronet[muon_count]/F");
  tree->Branch("muon_pfisolationr3_sumphotonet", muon_pfisolationr3_sumphotonet, "muon_pfisolationr3_sumphotonet[muon_count]/F");
  tree->Branch("muon_pfisolationr3_sumPUpt", muon_pfisolationr3_sumPUpt, "muon_pfisolationr3_sumPUpt[muon_count]/F");
  tree->Branch("muon_pfisolationr4_sumchargedhadronpt", muon_pfisolationr4_sumchargedhadronpt, "muon_pfisolationr4_sumchargedhadronpt[muon_count]/F");
  tree->Branch("muon_pfisolationr4_sumchargedparticlept", muon_pfisolationr4_sumchargedparticlept, "muon_pfisolationr4_sumchargedparticlept[muon_count]/F");
  tree->Branch("muon_pfisolationr4_sumneutralhadronet", muon_pfisolationr4_sumneutralhadronet, "muon_pfisolationr4_sumneutralhadronet[muon_count]/F");
  tree->Branch("muon_pfisolationr4_sumphotonet", muon_pfisolationr4_sumphotonet, "muon_pfisolationr4_sumphotonet[muon_count]/F");
  tree->Branch("muon_pfisolationr4_sumPUpt", muon_pfisolationr4_sumPUpt, "muon_pfisolationr4_sumPUpt[muon_count]/F");
  tree->Branch("muon_ecalenergy", muon_ecalenergy, "muon_ecalenergy[muon_count]/F");
  tree->Branch("muon_hcalenergy", muon_hcalenergy, "muon_hcalenergy[muon_count]/F");
  tree->Branch("muon_charge", muon_charge, "muon_charge[muon_count]/I");
  tree->Branch("muon_numchambers", muon_numchambers, "muon_numchambers[muon_count]/I");
  tree->Branch("muon_numchamberswithsegments", muon_numchamberswithsegments, "muon_numchamberswithsegments[muon_count]/I");
  tree->Branch("muon_type", muon_type, "muon_type[muon_count]/i");
  tree->Branch("muon_trigger", muon_trigger, "muon_trigger[muon_count]/i");
  tree->Branch("muon_trackermuonquality", muon_trackermuonquality, "muon_trackermuonquality[muon_count]/i");
  tree->Branch("muon_mva_id", muon_mva_id, "muon_mva_id[muon_count]/F");
  tree->Branch("muon_mva_iso", muon_mva_iso, "muon_mva_iso[muon_count]/F");
  
  tree->Branch("ak5calojet_count", &ak5calojet_count, "ak5calojet_count/i");
  tree->Branch("ak5calojet_e", ak5calojet_e, "ak5calojet_e[ak5calojet_count]/F");
  tree->Branch("ak5calojet_px", ak5calojet_px, "ak5calojet_px[ak5calojet_count]/F");
  tree->Branch("ak5calojet_py", ak5calojet_py, "ak5calojet_py[ak5calojet_count]/F");
  tree->Branch("ak5calojet_pz", ak5calojet_pz, "ak5calojet_pz[ak5calojet_count]/F");
  tree->Branch("ak5calojet_hadronicenergy", ak5calojet_hadronicenergy, "ak5calojet_hadronicenergy[ak5calojet_count]/F");
  tree->Branch("ak5calojet_emenergy", ak5calojet_emenergy, "ak5calojet_emenergy[ak5calojet_count]/F");
  tree->Branch("ak5calojet_energycorr", ak5calojet_energycorr, "ak5calojet_energycorr[ak5calojet_count]/F");
  tree->Branch("ak5calojet_fhpd", ak5calojet_fhpd, "ak5calojet_fhpd[ak5calojet_count]/F");
  tree->Branch("ak5calojet_restrictedemf", ak5calojet_restrictedemf, "ak5calojet_restrictedemf[ak5calojet_count]/F");
  tree->Branch("ak5calojet_btag", ak5calojet_btag, "ak5calojet_btag[ak5calojet_count][6]/F");
  tree->Branch("ak5calojet_n90", ak5calojet_n90, "ak5calojet_n90[ak5calojet_count]/i");
  tree->Branch("ak5calojet_n60", ak5calojet_n60, "ak5calojet_n60[ak5calojet_count]/i");
  
  tree->Branch("ak5jptjet_count", &ak5jptjet_count, "ak5jptjet_count/i");
  tree->Branch("ak5jptjet_e", ak5jptjet_e, "ak5jptjet_e[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_px", ak5jptjet_px, "ak5jptjet_px[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_py", ak5jptjet_py, "ak5jptjet_py[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_pz", ak5jptjet_pz, "ak5jptjet_pz[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_hadronicenergy", ak5jptjet_hadronicenergy, "ak5jptjet_hadronicenergy[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_chargedhadronicenergy", ak5jptjet_chargedhadronicenergy, "ak5jptjet_chargedhadronicenergy[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_emenergy", ak5jptjet_emenergy, "ak5jptjet_emenergy[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_chargedemenergy", ak5jptjet_chargedemenergy, "ak5jptjet_chargedemenergy[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_chargedmulti", ak5jptjet_chargedmulti, "ak5jptjet_chargedmulti[ak5jptjet_count]/i");	
  tree->Branch("ak5jptjet_energycorr", ak5jptjet_energycorr, "ak5jptjet_energycorr[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_fhpd", ak5jptjet_fhpd, "ak5jptjet_fhpd[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_restrictedemf", ak5jptjet_restrictedemf, "ak5jptjet_restrictedemf[ak5jptjet_count]/F");
  tree->Branch("ak5jptjet_btag", ak5jptjet_btag, "ak5jptjet_btag[ak5jptjet_count][6]/F");
  tree->Branch("ak5jptjet_n90", ak5jptjet_n90, "ak5jptjet_n90[ak5jptjet_count]/i");
  
  tree->Branch("ak5pfjet_count", &ak5pfjet_count, "ak5pfjet_count/i");
  tree->Branch("ak5pfjet_e", ak5pfjet_e, "ak5pfjet_e[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_px", ak5pfjet_px, "ak5pfjet_px[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_py", ak5pfjet_py, "ak5pfjet_py[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_pz", ak5pfjet_pz, "ak5pfjet_pz[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_hadronicenergy", ak5pfjet_hadronicenergy, "ak5pfjet_hadronicenergy[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_chargedhadronicenergy", ak5pfjet_chargedhadronicenergy, "ak5pfjet_chargedhadronicenergy[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_emenergy", ak5pfjet_emenergy, "ak5pfjet_emenergy[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_chargedemenergy", ak5pfjet_chargedemenergy, "ak5pfjet_chargedemenergy[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_chargedmulti", ak5pfjet_chargedmulti, "ak5pfjet_chargedmulti[ak5pfjet_count]/i");	
  tree->Branch("ak5pfjet_neutralmulti", ak5pfjet_neutralmulti, "ak5pfjet_neutralmulti[ak5pfjet_count]/i");	
  tree->Branch("ak5pfjet_energycorr", ak5pfjet_energycorr, "ak5pfjet_energycorr[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_btag", ak5pfjet_btag, "ak5pfjet_btag[ak5pfjet_count][6]/F");
  tree->Branch("ak5pfjet_trigger", ak5pfjet_trigger, "ak5pfjet_trigger[ak5pfjet_count]/i");
  tree->Branch("ak5pfjet_sv_count", ak5pfjet_sv_count, "ak5pfjet_sv_count[ak5pfjet_count]/i");
  tree->Branch("ak5pfjet_sv_x", ak5pfjet_sv_x, "ak5pfjet_sv_x[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_sv_y", ak5pfjet_sv_y, "ak5pfjet_sv_y[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_sv_z", ak5pfjet_sv_z, "ak5pfjet_sv_z[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_sv_dx", ak5pfjet_sv_dx, "ak5pfjet_sv_dx[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_sv_dy", ak5pfjet_sv_dy, "ak5pfjet_sv_dy[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_sv_dz", ak5pfjet_sv_dz, "ak5pfjet_sv_dz[ak5pfjet_count]/F");
  //tree->Branch("ak5pfjet_pu_jet_cut_loose", ak5pfjet_pu_jet_cut_loose, "ak5pfjet_pu_jet_cut_loose[ak5pfjet_count]/O");
  //tree->Branch("ak5pfjet_pu_jet_cut_medium", ak5pfjet_pu_jet_cut_medium, "ak5pfjet_pu_jet_cut_medium[ak5pfjet_count]/O");
  //tree->Branch("ak5pfjet_pu_jet_cut_tight", ak5pfjet_pu_jet_cut_tight, "ak5pfjet_pu_jet_cut_tight[ak5pfjet_count]/O");
  //tree->Branch("ak5pfjet_pu_jet_cut_mva", ak5pfjet_pu_jet_cut_mva, "ak5pfjet_pu_jet_cut_mva[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_pu_jet_simple_loose", ak5pfjet_pu_jet_simple_loose, "ak5pfjet_pu_jet_simple_loose[ak5pfjet_count]/O");
  tree->Branch("ak5pfjet_pu_jet_simple_medium", ak5pfjet_pu_jet_simple_medium, "ak5pfjet_pu_jet_simple_medium[ak5pfjet_count]/O");
  tree->Branch("ak5pfjet_pu_jet_simple_tight", ak5pfjet_pu_jet_simple_tight, "ak5pfjet_pu_jet_simple_tight[ak5pfjet_count]/O");
  tree->Branch("ak5pfjet_pu_jet_simple_mva", ak5pfjet_pu_jet_simple_mva, "ak5pfjet_pu_jet_simple_mva[ak5pfjet_count]/F");
  tree->Branch("ak5pfjet_pu_jet_full_loose", ak5pfjet_pu_jet_full_loose, "ak5pfjet_pu_jet_full_loose[ak5pfjet_count]/O");
  tree->Branch("ak5pfjet_pu_jet_full_medium", ak5pfjet_pu_jet_full_medium, "ak5pfjet_pu_jet_full_medium[ak5pfjet_count]/O");
  tree->Branch("ak5pfjet_pu_jet_full_tight", ak5pfjet_pu_jet_full_tight, "ak5pfjet_pu_jet_full_tight[ak5pfjet_count]/O");
  tree->Branch("ak5pfjet_pu_jet_full_mva", ak5pfjet_pu_jet_full_mva, "ak5pfjet_pu_jet_full_mva[ak5pfjet_count]/F");

  tree->Branch("electron_count", &electron_count, "electron_count/i");
  tree->Branch("electron_px", electron_px, "electron_px[electron_count]/F");
  tree->Branch("electron_py", electron_py, "electron_py[electron_count]/F");
  tree->Branch("electron_pz", electron_pz, "electron_pz[electron_count]/F");
  tree->Branch("electron_trackchi2", electron_trackchi2, "electron_trackchi2[electron_count]/F");
  tree->Branch("electron_trackndof", electron_trackndof, "electron_trackndof[electron_count]/F");
  tree->Branch("electron_outerx", electron_outerx, "electron_outerx[electron_count]/F");
  tree->Branch("electron_outery", electron_outery, "electron_outery[electron_count]/F");
  tree->Branch("electron_outerz", electron_outerz, "electron_outerz[electron_count]/F");
  tree->Branch("electron_closestpointx", electron_closestpointx, "electron_closestpointx[electron_count]/F");
  tree->Branch("electron_closestpointy", electron_closestpointy, "electron_closestpointy[electron_count]/F");
  tree->Branch("electron_closestpointz", electron_closestpointz, "electron_closestpointz[electron_count]/F");
  tree->Branch("electron_esuperclusterovertrack", electron_esuperclusterovertrack, "electron_esuperclusterovertrack[electron_count]/F");
  tree->Branch("electron_eseedclusterovertrack", electron_eseedclusterovertrack, "electron_eseedclusterovertrack[electron_count]/F");
  tree->Branch("electron_deltaetasuperclustertrack", electron_deltaetasuperclustertrack, "electron_deltaetasuperclustertrack[electron_count]/F");
  tree->Branch("electron_deltaphisuperclustertrack", electron_deltaphisuperclustertrack, "electron_deltaphisuperclustertrack[electron_count]/F");
  tree->Branch("electron_e1x5", electron_e1x5, "electron_e1x5[electron_count]/F");
  tree->Branch("electron_e2x5", electron_e2x5, "electron_e2x5[electron_count]/F");
  tree->Branch("electron_e5x5", electron_e5x5, "electron_e5x5[electron_count]/F");
  tree->Branch("electron_sigmaetaeta", electron_sigmaetaeta, "electron_sigmaetaeta[electron_count]/F");
  tree->Branch("electron_sigmaietaieta", electron_sigmaietaieta, "electron_sigmaietaieta[electron_count]/F");
  tree->Branch("electron_ehcaloverecal", electron_ehcaloverecal, "electron_ehcaloverecal[electron_count]/F");
  tree->Branch("electron_ehcaloverecaldepth1", electron_ehcaloverecaldepth1, "electron_ehcaloverecaldepth1[electron_count]/F");
  tree->Branch("electron_ehcaloverecaldepth2", electron_ehcaloverecaldepth2, "electron_ehcaloverecaldepth2[electron_count]/F");
  tree->Branch("electron_isolationr3track", electron_isolationr3track, "electron_isolationr3track[electron_count]/F");
  tree->Branch("electron_isolationr3ecal", electron_isolationr3ecal, "electron_isolationr3ecal[electron_count]/F");
  tree->Branch("electron_isolationr3hcal", electron_isolationr3hcal, "electron_isolationr3hcal[electron_count]/F");
  tree->Branch("electron_isolationr4track", electron_isolationr4track, "electron_isolationr4track[electron_count]/F");
  tree->Branch("electron_isolationr4ecal", electron_isolationr4ecal, "electron_isolationr4ecal[electron_count]/F");
  tree->Branch("electron_isolationr4hcal", electron_isolationr4hcal, "electron_isolationr4hcal[electron_count]/F");
  tree->Branch("electron_pfisolationr3_sumchargedhadronpt", electron_pfisolationr3_sumchargedhadronpt, "electron_pfisolationr3_sumchargedhadronpt[electron_count]/F");
  tree->Branch("electron_pfisolationr3_sumchargedparticlept", electron_pfisolationr3_sumchargedparticlept, "electron_pfisolationr3_sumchargedparticlept[electron_count]/F");
  tree->Branch("electron_pfisolationr3_sumneutralhadronet", electron_pfisolationr3_sumneutralhadronet, "electron_pfisolationr3_sumneutralhadronet[electron_count]/F");
  tree->Branch("electron_pfisolationr3_sumphotonet", electron_pfisolationr3_sumphotonet, "electron_pfisolationr3_sumphotonet[electron_count]/F");
  tree->Branch("electron_pfisolationr3_sumPUpt", electron_pfisolationr3_sumPUpt, "electron_pfisolationr3_sumPUpt[electron_count]/F");
  tree->Branch("electron_pfisolationr4_sumchargedhadronpt", electron_pfisolationr4_sumchargedhadronpt, "electron_pfisolationr4_sumchargedhadronpt[electron_count]/F");
  tree->Branch("electron_pfisolationr4_sumchargedparticlept", electron_pfisolationr4_sumchargedparticlept, "electron_pfisolationr4_sumchargedparticlept[electron_count]/F");
  tree->Branch("electron_pfisolationr4_sumneutralhadronet", electron_pfisolationr4_sumneutralhadronet, "electron_pfisolationr4_sumneutralhadronet[electron_count]/F");
  tree->Branch("electron_pfisolationr4_sumphotonet", electron_pfisolationr4_sumphotonet, "electron_pfisolationr4_sumphotonet[electron_count]/F");
  tree->Branch("electron_pfisolationr4_sumPUpt", electron_pfisolationr4_sumPUpt, "electron_pfisolationr4_sumPUpt[electron_count]/F");
  tree->Branch("electron_nhits", electron_nhits, "electron_nhits[electron_count]/b");
  tree->Branch("electron_npixelhits", electron_npixelhits, "electron_npixelhits[electron_count]/b"); 
  tree->Branch("electron_nmissinghits", electron_nmissinghits, "electron_nmissinghits[electron_count]/b");
  tree->Branch("electron_npixellayers", electron_npixellayers, "electron_npixellayers[electron_count]/b");
  tree->Branch("electron_nstriplayers", electron_nstriplayers, "electron_nstriplayers[electron_count]/b");
  tree->Branch("electron_dxy", electron_dxy, "electron_dxy[electron_count]/F");
  tree->Branch("electron_dxyerr", electron_dxyerr, "electron_dxyerr[electron_count]/F");
  tree->Branch("electron_dz", electron_dz, "electron_dz[electron_count]/F");
  tree->Branch("electron_dzerr", electron_dzerr, "electron_dzerr[electron_count]/F"); 
  tree->Branch("electron_convdist", electron_convdist, "electron_convdist[electron_count]/F");
  tree->Branch("electron_convdcot", electron_convdcot, "electron_convdcot[electron_count]/F");
  tree->Branch("electron_convradius", electron_convradius, "electron_convradius[electron_count]/F");
  tree->Branch("electron_gapinfo", electron_gapinfo, "electron_gapinfo[electron_count]/i");
  tree->Branch("electron_chargeinfo", electron_chargeinfo, "electron_chargeinfo[electron_count]/i");
  tree->Branch("electron_fbrems", electron_fbrems, "electron_fbrems[electron_count]/F");
  tree->Branch("electron_numbrems", electron_numbrems, "electron_numbrems[electron_count]/I");
  tree->Branch("electron_charge", electron_charge, "electron_charge[electron_count]/I");
  tree->Branch("electron_superclusterindex", electron_superclusterindex, "electron_superclusterindex[electron_count]/I");
  tree->Branch("electron_info", electron_info, "electron_info[electron_count]/b");
  tree->Branch("electron_trigger", electron_trigger, "electron_trigger[electron_count]/i");
  tree->Branch("electron_mva_id_trig", electron_mva_id_trig, "electron_mva_id_trig[electron_count]/F");
  tree->Branch("electron_mva_id_nontrig", electron_mva_id_nontrig, "electron_mva_id_nontrig[electron_count]/F");
  tree->Branch("electron_mva_iso", electron_mva_iso, "electron_mva_iso[electron_count]/F");
  tree->Branch("electron_has_conversion", electron_has_conversion, "electron_has_conversion[electron_count]/O");
  
  tree->Branch("photon_count", &photon_count, "photon_count/i");
  tree->Branch("photon_px", photon_px, "photon_px[photon_count]/F");
  tree->Branch("photon_py", photon_py, "photon_py[photon_count]/F");
  tree->Branch("photon_pz", photon_pz, "photon_pz[photon_count]/F");
  tree->Branch("photon_e1x5", photon_e1x5, "photon_e1x5[photon_count]/F");
  tree->Branch("photon_e2x5", photon_e2x5, "photon_e2x5[photon_count]/F");
  tree->Branch("photon_e3x3", photon_e3x3, "photon_e3x3[photon_count]/F");
  tree->Branch("photon_e5x5", photon_e5x5, "photon_e5x5[photon_count]/F");
  tree->Branch("photon_sigmaetaeta", photon_sigmaetaeta, "photon_sigmaetaeta[photon_count]/F");
  tree->Branch("photon_sigmaietaieta", photon_sigmaietaieta, "photon_sigmaietaieta[photon_count]/F");
  tree->Branch("photon_ehcaloverecal", photon_ehcaloverecal, "photon_ehcaloverecal[photon_count]/F");
  tree->Branch("photon_ehcaloverecaldepth1", photon_ehcaloverecaldepth1, "photon_ehcaloverecaldepth1[photon_count]/F");
  tree->Branch("photon_ehcaloverecaldepth2", photon_ehcaloverecaldepth2, "photon_ehcaloverecaldepth2[photon_count]/F");
  tree->Branch("photon_maxenergyxtal", photon_maxenergyxtal, "photon_maxenergyxtal[photon_count]/F");
  tree->Branch("photon_isolationr3track", photon_isolationr3track, "photon_isolationr3track[photon_count]/F");
  tree->Branch("photon_isolationr3trackhollow", photon_isolationr3trackhollow, "photon_isolationr3trackhollow[photon_count]/F");
  tree->Branch("photon_isolationr3ntrack", photon_isolationr3ntrack, "photon_isolationr3ntrack[photon_count]/i");
  tree->Branch("photon_isolationr3ntrackhollow", photon_isolationr3ntrackhollow, "photon_isolationr3ntrackhollow[photon_count]/i");
  tree->Branch("photon_isolationr3ecal", photon_isolationr3ecal, "photon_isolationr3ecal[photon_count]/F");
  tree->Branch("photon_isolationr3hcal", photon_isolationr3hcal, "photon_isolationr3hcal[photon_count]/F");
  tree->Branch("photon_isolationr4track", photon_isolationr4track, "photon_isolationr4track[photon_count]/F");
  tree->Branch("photon_isolationr4trackhollow", photon_isolationr4trackhollow, "photon_isolationr4trackhollow[photon_count]/F");
  tree->Branch("photon_isolationr4ntrack", photon_isolationr4ntrack, "photon_isolationr4ntrack[photon_count]/i");
  tree->Branch("photon_isolationr4ntrackhollow", photon_isolationr4ntrackhollow, "photon_isolationr4ntrackhollow[photon_count]/i");
  tree->Branch("photon_isolationr4ecal", photon_isolationr4ecal, "photon_isolationr4ecal[photon_count]/F");
  tree->Branch("photon_isolationr4hcal", photon_isolationr4hcal, "photon_isolationr4hcal[photon_count]/F");
  tree->Branch("photon_superclusterindex", photon_superclusterindex, "photon_superclusterindex[photon_count]/I");
  tree->Branch("photon_info", photon_info, "photon_info[photon_count]/b");
  tree->Branch("photon_gapinfo", photon_gapinfo, "photon_gapinfo[photon_count]/i");
  tree->Branch("photon_trigger", photon_trigger, "photon_trigger[photon_count]/i");
  tree->Branch("photon_conversionbegin", photon_conversionbegin, "photon_conversionbegin[photon_count]/i");
  
  tree->Branch("conversion_count", &conversion_count, "conversion_count/i");
  tree->Branch("conversion_info", conversion_info, "conversion_info[conversion_count]/b");
  tree->Branch("conversion_vx", conversion_vx, "conversion_vx[conversion_count]/F");
  tree->Branch("conversion_vy", conversion_vy, "conversion_vy[conversion_count]/F");
  tree->Branch("conversion_vz", conversion_vz, "conversion_vz[conversion_count]/F");
  tree->Branch("conversion_chi2", conversion_chi2, "conversion_chi2[conversion_count]/F");
  tree->Branch("conversion_ndof", conversion_ndof, "conversion_ndof[conversion_count]/F");
  tree->Branch("conversion_cov", conversion_cov, "conversion_cov[conversion_count][6]/F");
  tree->Branch("conversion_mvaout", conversion_mvaout, "conversion_mvaout[conversion_count]/F");
  tree->Branch("conversion_trackecalpointx", conversion_trackecalpointx, "conversion_trackecalpointx[conversion_count][2]/F");
  tree->Branch("conversion_trackecalpointy", conversion_trackecalpointy, "conversion_trackecalpointy[conversion_count][2]/F");
  tree->Branch("conversion_trackecalpointz", conversion_trackecalpointz, "conversion_trackecalpointz[conversion_count][2]/F");
  tree->Branch("conversion_trackpx", conversion_trackpx, "conversion_trackpx[conversion_count][2]/F");
  tree->Branch("conversion_trackpy", conversion_trackpy, "conversion_trackpy[conversion_count][2]/F");
  tree->Branch("conversion_trackpz", conversion_trackpz, "conversion_trackpz[conversion_count][2]/F");
  tree->Branch("conversion_trackclosestpointx", conversion_trackclosestpointx, "conversion_trackclosestpointx[conversion_count][2]/F");
  tree->Branch("conversion_trackclosestpointy", conversion_trackclosestpointy, "conversion_trackclosestpointy[conversion_count][2]/F");
  tree->Branch("conversion_trackclosestpointz", conversion_trackclosestpointz, "conversion_trackclosestpointz[conversion_count][2]/F");
  tree->Branch("conversion_trackchi2", conversion_trackchi2, "conversion_trackchi2[conversion_count][2]/F");
  tree->Branch("conversion_trackndof", conversion_trackndof, "conversion_trackndof[conversion_count][2]/F");
  tree->Branch("conversion_trackdxy", conversion_trackdxy, "conversion_trackdxy[conversion_count][2]/F");
  tree->Branch("conversion_trackdxyerr", conversion_trackdxyerr, "conversion_trackdxyerr[conversion_count][2]/F");
  tree->Branch("conversion_trackdz", conversion_trackdz, "conversion_trackdz[conversion_count][2]/F");
  tree->Branch("conversion_trackdzerr", conversion_trackdzerr, "conversion_trackdzerr[conversion_count][2]/F");
  tree->Branch("conversion_trackcharge", conversion_trackcharge, "conversion_trackcharge[conversion_count][2]/I");
  tree->Branch("conversion_tracknhits", conversion_tracknhits, "conversion_tracknhits[conversion_count][2]/b");
  tree->Branch("conversion_tracknmissinghits", conversion_tracknmissinghits, "conversion_tracknmissinghits[conversion_count][2]/b");
  tree->Branch("conversion_tracknpixelhits", conversion_tracknpixelhits, "conversion_tracknpixelhits[conversion_count][2]/b");
  tree->Branch("conversion_tracknpixellayers", conversion_tracknpixellayers, "conversion_tracknpixellayers[conversion_count][2]/b");
  tree->Branch("conversion_tracknstriplayers", conversion_tracknstriplayers, "conversion_tracknstriplayers[conversion_count][2]/b");
  
  tree->Branch("allconversion_count", &allconversion_count, "allconversion_count/i");
  tree->Branch("allconversion_info", allconversion_info, "allconversion_info[allconversion_count]/b");
  tree->Branch("allconversion_vx", allconversion_vx, "allconversion_vx[allconversion_count]/F");
  tree->Branch("allconversion_vy", allconversion_vy, "allconversion_vy[allconversion_count]/F");
  tree->Branch("allconversion_vz", allconversion_vz, "allconversion_vz[allconversion_count]/F");
  tree->Branch("allconversion_chi2", allconversion_chi2, "allconversion_chi2[allconversion_count]/F");
  tree->Branch("allconversion_ndof", allconversion_ndof, "allconversion_ndof[allconversion_count]/F");
  tree->Branch("allconversion_cov", allconversion_cov, "allconversion_cov[allconversion_count][6]/F");
  tree->Branch("allconversion_mvaout", allconversion_mvaout, "allconversion_mvaout[allconversion_count]/F");
  tree->Branch("allconversion_trackecalpointx", allconversion_trackecalpointx, "allconversion_trackecalpointx[allconversion_count][2]/F");
  tree->Branch("allconversion_trackecalpointy", allconversion_trackecalpointy, "allconversion_trackecalpointy[allconversion_count][2]/F");
  tree->Branch("allconversion_trackecalpointz", allconversion_trackecalpointz, "allconversion_trackecalpointz[allconversion_count][2]/F");
  tree->Branch("allconversion_trackpx", allconversion_trackpx, "allconversion_trackpx[allconversion_count][2]/F");
  tree->Branch("allconversion_trackpy", allconversion_trackpy, "allconversion_trackpy[allconversion_count][2]/F");
  tree->Branch("allconversion_trackpz", allconversion_trackpz, "allconversion_trackpz[allconversion_count][2]/F");
  tree->Branch("allconversion_trackclosestpointx", allconversion_trackclosestpointx, "allconversion_trackclosestpointx[allconversion_count][2]/F");
  tree->Branch("allconversion_trackclosestpointy", allconversion_trackclosestpointy, "allconversion_trackclosestpointy[allconversion_count][2]/F");
  tree->Branch("allconversion_trackclosestpointz", allconversion_trackclosestpointz, "allconversion_trackclosestpointz[allconversion_count][2]/F");
  tree->Branch("allconversion_trackchi2", allconversion_trackchi2, "allconversion_trackchi2[allconversion_count][2]/F");
  tree->Branch("allconversion_trackndof", allconversion_trackndof, "allconversion_trackndof[allconversion_count][2]/F");
  tree->Branch("allconversion_trackdxy", allconversion_trackdxy, "allconversion_trackdxy[allconversion_count][2]/F");
  tree->Branch("allconversion_trackdxyerr", allconversion_trackdxyerr, "allconversion_trackdxyerr[allconversion_count][2]/F");
  tree->Branch("allconversion_trackdz", allconversion_trackdz, "allconversion_trackdz[allconversion_count][2]/F");
  tree->Branch("allconversion_trackdzerr", allconversion_trackdzerr, "allconversion_trackdzerr[allconversion_count][2]/F");
  tree->Branch("allconversion_trackcharge", allconversion_trackcharge, "allconversion_trackcharge[allconversion_count][2]/I");
  tree->Branch("allconversion_tracknhits", allconversion_tracknhits, "allconversion_tracknhits[allconversion_count][2]/b");
  tree->Branch("allconversion_tracknmissinghits", allconversion_tracknmissinghits, "allconversion_tracknmissinghits[allconversion_count][2]/b");
  tree->Branch("allconversion_tracknpixelhits", allconversion_tracknpixelhits, "allconversion_tracknpixelhits[allconversion_count][2]/b");
  tree->Branch("allconversion_tracknpixellayers", allconversion_tracknpixellayers, "allconversion_tracknpixellayers[allconversion_count][2]/b");
  tree->Branch("allconversion_tracknstriplayers", allconversion_tracknstriplayers, "allconversion_tracknstriplayers[allconversion_count][2]/b");

  tree->Branch("tau_count", &tau_count, "tau_count/i");
  tree->Branch("tau_e", tau_e, "tau_e[tau_count]/F");
  tree->Branch("tau_px", tau_px, "tau_px[tau_count]/F");
  tree->Branch("tau_py", tau_py, "tau_py[tau_count]/F");
  tree->Branch("tau_pz", tau_pz, "tau_pz[tau_count]/F");
  tree->Branch("tau_isolationneutralspt", tau_isolationneutralspt, "tau_isolationneutralspt[tau_count]/F");
  tree->Branch("tau_isolationneutralsnum", tau_isolationneutralsnum, "tau_isolationneutralsnum[tau_count]/i");
  tree->Branch("tau_isolationchargedpt", tau_isolationchargedpt, "tau_isolationchargedpt[tau_count]/F");
  tree->Branch("tau_isolationchargednum", tau_isolationchargednum, "tau_isolationchargednum[tau_count]/i");
  tree->Branch("tau_isolationgammapt", tau_isolationgammapt, "tau_isolationgammapt[tau_count]/F");
  tree->Branch("tau_isolationgammanum", tau_isolationgammanum, "tau_isolationgammanum[tau_count]/i");
  tree->Branch("tau_leadpfchargedhadrcandecalenergy", tau_leadpfchargedhadrcandecalenergy, "tau_leadpfchargedhadrcandecalenergy[tau_count]/F");
  tree->Branch("tau_leadpfchargedhadrcandhcalenergy", tau_leadpfchargedhadrcandhcalenergy, "tau_leadpfchargedhadrcandhcalenergy[tau_count]/F");
  tree->Branch("tau_leadpfchargedhadrcandp", tau_leadpfchargedhadrcandp, "tau_leadpfchargedhadrcandp[tau_count]/F");
  tree->Branch("tau_dxy", tau_dxy, "tau_dxy[tau_count]/F");
  tree->Branch("tau_dz", tau_dz, "tau_dz[tau_count]/F");
  tree->Branch("tau_ip3d", tau_ip3d, "tau_ip3d[tau_count]/F");
  tree->Branch("tau_ip3dSig", tau_ip3dSig, "tau_ip3dSig[tau_count]/F");
  tree->Branch("tau_nprongs", tau_nprongs, "tau_nprongs[tau_count]/i");
  tree->Branch("tau_charge", tau_charge, "tau_charge[tau_count]/I");
  tree->Branch("tau_emfraction", tau_emfraction, "tau_emfraction[tau_count]/F");
  tree->Branch("tau_newemfraction", tau_newemfraction, "tau_newemfraction[tau_count]/F");
  tree->Branch("tau_hcaltotoverplead", tau_hcaltotoverplead, "tau_hcaltotoverplead[tau_count]/F");
  tree->Branch("tau_hcal3x3overplead", tau_hcal3x3overplead, "tau_hcal3x3overplead[tau_count]/F");
  tree->Branch("tau_ecalstripsumeoverplead", tau_ecalstripsumeoverplead, "tau_ecalstripsumeoverplead[tau_count]/F");
  tree->Branch("tau_bremsrecoveryeoverplead", tau_bremsrecoveryeoverplead, "tau_bremsrecoveryeoverplead[tau_count]/F");
  tree->Branch("tau_calocomp", tau_calocomp, "tau_calocomp[tau_count]/F");
  tree->Branch("tau_segcomp", tau_segcomp, "tau_segcomp[tau_count]/F");
  tree->Branch("tau_dishps", tau_dishps, "tau_dishps[tau_count]/l");

  tree->Branch("tau_againstelectronmva5raw", tau_againstelectronmva5raw, "tau_againstelectronmva5raw[tau_count]/F");
  tree->Branch("tau_byIsolationmva3newDMwoLTraw", tau_byIsolationmva3newDMwoLTraw, "tau_byIsolationmva3newDMwoLTraw[tau_count]/F");
  tree->Branch("tau_byIsolationmva3newDMwLTraw", tau_byIsolationmva3newDMwLTraw, "tau_byIsolationmva3newDMwLTraw[tau_count]/F");
  tree->Branch("tau_againstelectronmva5raw", tau_againstelectronmva5raw, "tau_againstelectronmva5raw[tau_count]/F");
  tree->Branch("tau_againstelectronVLoosemva5", tau_againstelectronVLoosemva5, "tau_againstelectronVLoosemva5[tau_count]/F");
  tree->Branch("tau_againstelectronLoosemva5", tau_againstelectronLoosemva5, "tau_againstelectronLoosemva5[tau_count]/F");
  tree->Branch("tau_againstelectronMediummva5", tau_againstelectronMediummva5, "tau_againstelectronMediummva5[tau_count]/F");
  tree->Branch("tau_againstelectronTightmva5", tau_againstelectronTightmva5, "tau_againstelectronTightmva5[tau_count]/F");
  tree->Branch("tau_againstelectronDeadECAL", tau_againstelectronDeadECAL, "tau_againstelectronDeadECAL[tau_count]/F");

  tree->Branch("tau_againstelectronmva5category", tau_againstelectronmva5category, "tau_againstelectronmva5category[tau_count]/F");
  tree->Branch("tau_bycombinedisolationdeltabetacorrraw3hits", tau_bycombinedisolationdeltabetacorrraw3hits, "tau_bycombinedisolationdeltabetacorrraw3hits[tau_count]/F");

  tree->Branch("tau_trigger", tau_trigger, "tau_trigger[tau_count]/i");
  tree->Branch("tau_L1trigger_match", tau_L1trigger_match, "tau_L1trigger_match[tau_count]/O");
  tree->Branch("tau_signalPFChargedHadrCands_size", tau_signalPFChargedHadrCands_size, "tau_signalPFChargedHadrCands_size[tau_count]/i");
  tree->Branch("tau_signalPFGammaCands_size", tau_signalPFGammaCands_size, "tau_signalPFGammaCands_size[tau_count]/i");
  tree->Branch("tau_genTaudecayMode", (void*)tau_genTaudecayMode, "tau_genTaudecayMode[tau_count]/C");

  // ditau study
  tree->Branch("ditau_Index",&ditau_Index,"ditau_Index/i");
  tree->Branch("ditau_leg1_index", ditau_leg1_index,"ditau_leg1_index[ditau_Index]/i");
  tree->Branch("ditau_leg2_index", ditau_leg2_index,"ditau_leg2_index[ditau_Index]/i");

  tree->Branch("ditau_VtxX", &ditau_VtxX,"ditau_VtxX/F");
  tree->Branch("ditau_VtxY", &ditau_VtxY,"ditau_VtxY/F");
  tree->Branch("ditau_VtxZ", &ditau_VtxZ,"ditau_VtxZ/F");
  tree->Branch("ditau_VtxXErr", &ditau_VtxXErr,"ditau_VtxXErr/F");
  tree->Branch("ditau_VtxYErr", &ditau_VtxYErr,"ditau_VtxYErr/F");
  tree->Branch("ditau_VtxZErr", &ditau_VtxZErr,"ditau_VtxZErr/F");

  tree->Branch("ditau_reFitVtxX", ditau_reFitVtxX,"ditau_reFitVtxX[ditau_Index]/F");
  tree->Branch("ditau_reFitVtxY", ditau_reFitVtxY,"ditau_reFitVtxY[ditau_Index]/F");
  tree->Branch("ditau_reFitVtxZ", ditau_reFitVtxZ,"ditau_reFitVtxZ[ditau_Index]/F");
  tree->Branch("ditau_reFitVtxXErr", ditau_reFitVtxXErr,"ditau_reFitVtxXErr[ditau_Index]/F");
  tree->Branch("ditau_reFitVtxYErr", ditau_reFitVtxYErr,"ditau_reFitVtxYErr[ditau_Index]/F");
  tree->Branch("ditau_reFitVtxZErr", ditau_reFitVtxZErr,"ditau_reFitVtxZErr[ditau_Index]/F");
  tree->Branch("ditau_reFitVtxRho", ditau_reFitVtxRho,"ditau_reFitVtxRho[ditau_Index]/i");

  tree->Branch("ditau_leg1_dxy", ditau_leg1_dxy,"ditau_leg1_dxy[ditau_Index]/F");
  tree->Branch("ditau_leg1_dz", ditau_leg1_dz,"ditau_leg1_dz[ditau_Index]/F");
  tree->Branch("ditau_leg1_dxyErr", ditau_leg1_dxyErr,"ditau_leg1_dxyErr[ditau_Index]/F");
  tree->Branch("ditau_leg1_dzErr", ditau_leg1_dzErr,"ditau_leg1_dzErr[ditau_Index]/F");
  tree->Branch("ditau_leg1_dxy_OPV", ditau_leg1_dxy_OPV,"ditau_leg1_dxy_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg1_dz_OPV", ditau_leg1_dz_OPV,"ditau_leg1_dz_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg1_dxyErr_OPV", ditau_leg1_dxyErr_OPV,"ditau_leg1_dxyErr_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg1_dzErr_OPV", ditau_leg1_dzErr_OPV,"ditau_leg1_dzErr_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg1_X", ditau_leg1_X,"ditau_leg1_X[ditau_Index]/F");
  tree->Branch("ditau_leg1_Y", ditau_leg1_Y,"ditau_leg1_Y[ditau_Index]/F");
  tree->Branch("ditau_leg1_Z", ditau_leg1_Z,"ditau_leg1_Z[ditau_Index]/F");
  tree->Branch("ditau_leg1_X_OPV", ditau_leg1_X_OPV,"ditau_leg1_X_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg1_Y_OPV", ditau_leg1_Y_OPV,"ditau_leg1_Y_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg1_Z_OPV", ditau_leg1_Z_OPV,"ditau_leg1_Z_OPV[ditau_Index]/F");

  tree->Branch("ditau_leg2_dxy", ditau_leg2_dxy,"ditau_leg2_dxy[ditau_Index]/F");
  tree->Branch("ditau_leg2_dz", ditau_leg2_dz,"ditau_leg2_dz[ditau_Index]/F");
  tree->Branch("ditau_leg2_dxyErr", ditau_leg2_dxyErr,"ditau_leg2_dxyErr[ditau_Index]/F");
  tree->Branch("ditau_leg2_dzErr", ditau_leg2_dzErr,"ditau_leg2_dzErr[ditau_Index]/F");
  tree->Branch("ditau_leg2_dxy_OPV", ditau_leg2_dxy_OPV,"ditau_leg2_dxy_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg2_dz_OPV", ditau_leg2_dz_OPV,"ditau_leg2_dz_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg2_dxyErr_OPV", ditau_leg2_dxyErr_OPV,"ditau_leg2_dxyErr_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg2_dzErr_OPV", ditau_leg2_dzErr_OPV,"ditau_leg2_dzErr_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg2_X", ditau_leg2_X,"ditau_leg2_X[ditau_Index]/F");
  tree->Branch("ditau_leg2_Y",ditau_leg2_Y,"ditau_leg2_Y[ditau_Index]/F");
  tree->Branch("ditau_leg2_Z", ditau_leg2_Z,"ditau_leg2_Z[ditau_Index]/F");
  tree->Branch("ditau_leg2_X_OPV", ditau_leg2_X_OPV,"ditau_leg2_X_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg2_Y_OPV", ditau_leg2_Y_OPV,"ditau_leg2_Y_OPV[ditau_Index]/F");
  tree->Branch("ditau_leg2_Z_OPV", ditau_leg2_Z_OPV,"ditau_leg2_Z_OPV[ditau_Index]/F");

  tree->Branch("mutautaupair_count", &mutautaupair_count, "mutautaupair_count/i");
  tree->Branch("mutautaupair_leg1_px", mutautaupair_leg1_px, "mutautaupair_leg1_px[mutautaupair_count]/F");
  tree->Branch("mutautaupair_leg1_py", mutautaupair_leg1_py, "mutautaupair_leg1_py[mutautaupair_count]/F");
  tree->Branch("mutautaupair_leg1_pz", mutautaupair_leg1_pz, "mutautaupair_leg1_pz[mutautaupair_count]/F");
  tree->Branch("mutautaupair_leg1_energy", mutautaupair_leg1_energy, "mutautaupair_leg1_energy[mutautaupair_count]/F");
  tree->Branch("mutautaupair_leg2_px", mutautaupair_leg2_px, "mutautaupair_leg2_px[mutautaupair_count]/F");
  tree->Branch("mutautaupair_leg2_py", mutautaupair_leg2_py, "mutautaupair_leg2_py[mutautaupair_count]/F");
  tree->Branch("mutautaupair_leg2_pz", mutautaupair_leg2_pz, "mutautaupair_leg2_pz[mutautaupair_count]/F");
  tree->Branch("mutautaupair_leg2_energy", mutautaupair_leg2_energy, "mutautaupair_leg2_energy[mutautaupair_count]/F");
  tree->Branch("mutautaupair_mu_px", mutautaupair_mu_px, "mutautaupair_mu_px[mutautaupair_count]/F");
  tree->Branch("mutautaupair_mu_py", mutautaupair_mu_py, "mutautaupair_mu_py[mutautaupair_count]/F");
  tree->Branch("mutautaupair_mu_pz", mutautaupair_mu_pz, "mutautaupair_mu_pz[mutautaupair_count]/F");
  tree->Branch("mutautaupair_mu_energy", mutautaupair_mu_energy, "mutautaupair_mu_energy[mutautaupair_count]/F");
  tree->Branch("mutautaupair_svfit_int_valid", mutautaupair_svfit_int_valid, "mutautaupair_svfit_int_valid[mutautaupair_count]/O");
  tree->Branch("mutautaupair_svfit_mass_int", mutautaupair_svfit_mass_int, "mutautaupair_svfit_mass_int[mutautaupair_count]/F");
  tree->Branch("mutautaupair_svfit_mass_int_err_up", mutautaupair_svfit_mass_int_err_up, "mutautaupair_svfit_mass_int_err_up[mutautaupair_count]/F");
  tree->Branch("mutautaupair_svfit_mass_int_err_down", mutautaupair_svfit_mass_int_err_down, "mutautaupair_svfit_mass_int_err_down[mutautaupair_count]/F");
  tree->Branch("mutautaupair_dca2d", mutautaupair_dca2d, "mutautaupair_dca2d[mutautaupair_count]/F");
  tree->Branch("mutautaupair_dca2d_err", mutautaupair_dca2d_err, "mutautaupair_dca2d_err[mutautaupair_count]/F");
  tree->Branch("mutautaupair_dca3d", mutautaupair_dca3d, "mutautaupair_dca3d[mutautaupair_count]/F");
  tree->Branch("mutautaupair_dca3d_err", mutautaupair_dca3d_err, "mutautaupair_dca3d_err[mutautaupair_count]/F");
  
  tree->Branch("eltautaupair_count", &eltautaupair_count, "eltautaupair_count/i");
  tree->Branch("eltautaupair_leg1_px", eltautaupair_leg1_px, "eltautaupair_leg1_px[eltautaupair_count]/F");
  tree->Branch("eltautaupair_leg1_py", eltautaupair_leg1_py, "eltautaupair_leg1_py[eltautaupair_count]/F");
  tree->Branch("eltautaupair_leg1_pz", eltautaupair_leg1_pz, "eltautaupair_leg1_pz[eltautaupair_count]/F");
  tree->Branch("eltautaupair_leg1_energy", eltautaupair_leg1_energy, "eltautaupair_leg1_energy[eltautaupair_count]/F");
  tree->Branch("eltautaupair_leg2_px", eltautaupair_leg2_px, "eltautaupair_leg2_px[eltautaupair_count]/F");
  tree->Branch("eltautaupair_leg2_py", eltautaupair_leg2_py, "eltautaupair_leg2_py[eltautaupair_count]/F");
  tree->Branch("eltautaupair_leg2_pz", eltautaupair_leg2_pz, "eltautaupair_leg2_pz[eltautaupair_count]/F");
  tree->Branch("eltautaupair_leg2_energy", eltautaupair_leg2_energy, "eltautaupair_leg2_energy[eltautaupair_count]/F");
  tree->Branch("eltautaupair_el_px", eltautaupair_el_px, "eltautaupair_el_px[eltautaupair_count]/F");
  tree->Branch("eltautaupair_el_py", eltautaupair_el_py, "eltautaupair_el_py[eltautaupair_count]/F");
  tree->Branch("eltautaupair_el_pz", eltautaupair_el_pz, "eltautaupair_el_pz[eltautaupair_count]/F");
  tree->Branch("eltautaupair_el_energy", eltautaupair_el_energy, "eltautaupair_el_energy[eltautaupair_count]/F");
  tree->Branch("eltautaupair_svfit_int_valid", eltautaupair_svfit_int_valid, "eltautaupair_svfit_int_valid[eltautaupair_count]/O");
  tree->Branch("eltautaupair_svfit_mass_int", eltautaupair_svfit_mass_int, "eltautaupair_svfit_mass_int[eltautaupair_count]/F");
  tree->Branch("eltautaupair_svfit_mass_int_err_up", eltautaupair_svfit_mass_int_err_up, "eltautaupair_svfit_mass_int_err_up[eltautaupair_count]/F");
  tree->Branch("eltautaupair_svfit_mass_int_err_down", eltautaupair_svfit_mass_int_err_down, "eltautaupair_svfit_mass_int_err_down[eltautaupair_count]/F");
  tree->Branch("eltautaupair_dca2d", eltautaupair_dca2d, "eltautaupair_dca2d[eltautaupair_count]/F");
  tree->Branch("eltautaupair_dca2d_err", eltautaupair_dca2d_err, "eltautaupair_dca2d_err[eltautaupair_count]/F");
  tree->Branch("eltautaupair_dca3d", eltautaupair_dca3d, "eltautaupair_dca3d[eltautaupair_count]/F");
  tree->Branch("eltautaupair_dca3d_err", eltautaupair_dca3d_err, "eltautaupair_dca3d_err[eltautaupair_count]/F");
  
  tree->Branch("tau_ak5pfjet_e", tau_ak5pfjet_e, "tau_ak5pfjet_e[tau_count]/F");
  tree->Branch("tau_ak5pfjet_px", tau_ak5pfjet_px, "tau_ak5pfjet_px[tau_count]/F");
  tree->Branch("tau_ak5pfjet_py", tau_ak5pfjet_py, "tau_ak5pfjet_py[tau_count]/F");
  tree->Branch("tau_ak5pfjet_pz", tau_ak5pfjet_pz, "tau_ak5pfjet_pz[tau_count]/F");
  tree->Branch("tau_chargedbegin", tau_chargedbegin, "tau_chargedbegin[tau_count]/i");
  tree->Branch("tau_charged_count", &tau_charged_count, "tau_charged_count/i");
  tree->Branch("tau_charged_px", tau_charged_px, "tau_charged_px[tau_charged_count]/F");
  tree->Branch("tau_charged_py", tau_charged_py, "tau_charged_py[tau_charged_count]/F");
  tree->Branch("tau_charged_pz", tau_charged_pz, "tau_charged_pz[tau_charged_count]/F");
  tree->Branch("tau_charged_outerx", tau_charged_outerx, "tau_charged_outerx[tau_charged_count]/F");
  tree->Branch("tau_charged_outery", tau_charged_outery, "tau_charged_outery[tau_charged_count]/F");
  tree->Branch("tau_charged_outerz", tau_charged_outerz, "tau_charged_outerz[tau_charged_count]/F");
  tree->Branch("tau_charged_closestpointx", tau_charged_closestpointx, "tau_charged_closestpointx[tau_charged_count]/F");
  tree->Branch("tau_charged_closestpointy", tau_charged_closestpointy, "tau_charged_closestpointy[tau_charged_count]/F");
  tree->Branch("tau_charged_closestpointz", tau_charged_closestpointz, "tau_charged_closestpointz[tau_charged_count]/F");
  tree->Branch("tau_charged_chi2", tau_charged_chi2, "tau_charged_chi2[tau_charged_count]/F");
  tree->Branch("tau_charged_ndof", tau_charged_ndof, "tau_charged_ndof[tau_charged_count]/F");
  tree->Branch("tau_charged_dxy", tau_charged_dxy, "tau_charged_dxy[tau_charged_count]/F");
  tree->Branch("tau_charged_dxyerr", tau_charged_dxyerr, "tau_charged_dxyerr[tau_charged_count]/F");
  tree->Branch("tau_charged_dz", tau_charged_dz, "tau_charged_dz[tau_charged_count]/F");
  tree->Branch("tau_charged_dzerr", tau_charged_dzerr, "tau_charged_dzerr[tau_charged_count]/F");
  tree->Branch("tau_charged_dedxharmonic2", tau_charged_dedxharmonic2, "tau_charged_dedxharmonic2[tau_charged_count]/F");
  tree->Branch("tau_charged_charge", tau_charged_charge, "tau_charged_charge[tau_charged_count]/I");
  tree->Branch("tau_charged_nhits", tau_charged_nhits, "tau_charged_nhits[tau_charged_count]/b");
  tree->Branch("tau_charged_nmissinghits", tau_charged_nmissinghits, "tau_charged_nmissinghits[tau_charged_count]/b");
  tree->Branch("tau_charged_npixelhits", tau_charged_npixelhits, "tau_charged_npixelhits[tau_charged_count]/b");
  tree->Branch("tau_charged_npixellayers", tau_charged_npixellayers, "tau_charged_npixellayers[tau_charged_count]/b");
  tree->Branch("tau_charged_nstriplayers", tau_charged_nstriplayers, "tau_charged_nstriplayers[tau_charged_count]/b");
  
  tree->Branch("ak5pfjet_rho", &ak5pfjet_rho, "ak5pfjet_rho/F");
  tree->Branch("ak5pfjet_sigma", &ak5pfjet_sigma, "ak5pfjet_sigma/F");
  
  tree->Branch("calomet_ex", &calomet_ex, "calomet_ex/F");
  tree->Branch("calomet_ey", &calomet_ey, "calomet_ey/F");
  tree->Branch("calomet_exmuons", &calomet_exmuons, "calomet_exmuons/F");
  tree->Branch("calomet_eymuons", &calomet_eymuons, "calomet_eymuons/F");
  
  tree->Branch("tcmet_ex", &tcmet_ex, "tcmet_ex/F");
  tree->Branch("tcmet_ey", &tcmet_ey, "tcmet_ey/F");
  
  tree->Branch("pfmet_ex", &pfmet_ex, "pfmet_ex/F");
  tree->Branch("pfmet_ey", &pfmet_ey, "pfmet_ey/F");
  tree->Branch("pfmettype1_ex", &pfmettype1_ex, "pfmettype1_ex/F");
  tree->Branch("pfmettype1_ey", &pfmettype1_ey, "pfmettype1_ey/F");
  tree->Branch("pfmetmva_ex", &pfmetmva_ex, "pfmetmva_ex/F");
  tree->Branch("pfmetmva_ey", &pfmetmva_ey, "pfmetmva_ey/F");
  tree->Branch("pfmetsigxx", &pfmetsigxx, "pfmetsigxx/F");
  tree->Branch("pfmetsigxy", &pfmetsigxy, "pfmetsigxy/F");
  tree->Branch("pfmetsigyx", &pfmetsigyx, "pfmetsigyx/F");
  tree->Branch("pfmetsigyy", &pfmetsigyy, "pfmetsigyy/F");
  tree->Branch("secvertices_count", &secvertices_count, "secvertices_count/i");
  tree->Branch("secvertices_vx", secvertices_vx, "secvertices_vx[secvertices_count]/F");
  tree->Branch("secvertices_vy", secvertices_vy, "secvertices_vy[secvertices_count]/F");
  tree->Branch("secvertices_vz", secvertices_vz, "secvertices_vz[secvertices_count]/F");
  tree->Branch("secvertices_chi2", secvertices_chi2, "secvertices_chi2[secvertices_count]/F");
  tree->Branch("secvertices_ndof", secvertices_ndof, "secvertices_ndof[secvertices_count]/F");
  tree->Branch("secvertices_cov", secvertices_cov, "secvertices_cov[secvertices_count][6]/F");
  tree->Branch("secvertices_track_px", secvertices_track_px, "secvertices_track_px[secvertices_count][2]/F");
  tree->Branch("secvertices_track_py", secvertices_track_py, "secvertices_track_py[secvertices_count][2]/F");
  tree->Branch("secvertices_track_pz", secvertices_track_pz, "secvertices_track_pz[secvertices_count][2]/F");
  tree->Branch("secvertices_track_closestpointx", secvertices_track_closestpointx, "secvertices_track_closestpointx[secvertices_count][2]/F");
  tree->Branch("secvertices_track_closestpointy", secvertices_track_closestpointy, "secvertices_track_closestpointy[secvertices_count][2]/F");
  tree->Branch("secvertices_track_closestpointz", secvertices_track_closestpointz, "secvertices_track_closestpointz[secvertices_count][2]/F");
  tree->Branch("secvertices_track_chi2", secvertices_track_chi2, "secvertices_track_chi2[secvertices_count][2]/F");
  tree->Branch("secvertices_track_ndof", secvertices_track_ndof, "secvertices_track_ndof[secvertices_count][2]/F");
  tree->Branch("secvertices_track_dxy", secvertices_track_dxy, "secvertices_track_dxy[secvertices_count][2]/F");
  tree->Branch("secvertices_track_dxyerr", secvertices_track_dxyerr, "secvertices_track_dxyerr[secvertices_count][2]/F");
  tree->Branch("secvertices_track_dz", secvertices_track_dz, "secvertices_track_dz[secvertices_count][2]/F");
  tree->Branch("secvertices_track_dzerr", secvertices_track_dzerr, "secvertices_track_dzerr[secvertices_count][2]/F");
  tree->Branch("secvertices_track_dedxharmonic2", secvertices_track_dedxharmonic2, "secvertices_track_dedxharmonic2[track_count][2]/F");
  tree->Branch("secvertices_track_charge", secvertices_track_charge, "secvertices_track_charge[secvertices_count][2]/I");
  tree->Branch("secvertices_track_nhits", secvertices_track_nhits, "secvertices_track_nhits[secvertices_count][2]/b");
  tree->Branch("secvertices_track_nmissinghits", secvertices_track_nmissinghits, "secvertices_track_nmissinghits[secvertices_count][2]/b");
  tree->Branch("secvertices_track_npixelhits", secvertices_track_npixelhits, "secvertices_track_npixelhits[secvertices_count][2]/b");
  tree->Branch("secvertices_track_npixellayers", secvertices_track_npixellayers, "secvertices_track_npixellayers[secvertices_count][2]/b");
  tree->Branch("secvertices_track_nstriplayers", secvertices_track_nstriplayers, "secvertices_track_nstriplayers[secvertices_count][2]/b");
  
  tree->Branch("musecvertices_count", &musecvertices_count, "musecvertices_count/i");
  tree->Branch("musecvertices_vx", musecvertices_vx, "musecvertices_vx[musecvertices_count]/F");
  tree->Branch("musecvertices_vy", musecvertices_vy, "musecvertices_vy[musecvertices_count]/F");
  tree->Branch("musecvertices_vz", musecvertices_vz, "musecvertices_vz[musecvertices_count]/F");
  tree->Branch("musecvertices_chi2", musecvertices_chi2, "musecvertices_chi2[musecvertices_count]/F");
  tree->Branch("musecvertices_ndof", musecvertices_ndof, "musecvertices_ndof[musecvertices_count]/F");
  tree->Branch("musecvertices_cov", musecvertices_cov, "musecvertices_cov[musecvertices_count][6]/F");
  tree->Branch("musecvertices_track_px", musecvertices_track_px, "musecvertices_track_px[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_py", musecvertices_track_py, "musecvertices_track_py[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_pz", musecvertices_track_pz, "musecvertices_track_pz[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_closestpointx", musecvertices_track_closestpointx, "musecvertices_track_closestpointx[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_closestpointy", musecvertices_track_closestpointy, "musecvertices_track_closestpointy[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_closestpointz", musecvertices_track_closestpointz, "musecvertices_track_closestpointz[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_chi2", musecvertices_track_chi2, "musecvertices_track_chi2[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_ndof", musecvertices_track_ndof, "musecvertices_track_ndof[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_dxy", musecvertices_track_dxy, "musecvertices_track_dxy[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_dxyerr", musecvertices_track_dxyerr, "musecvertices_track_dxyerr[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_dz", musecvertices_track_dz, "musecvertices_track_dz[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_dzerr", musecvertices_track_dzerr, "musecvertices_track_dzerr[musecvertices_count][2]/F");
  tree->Branch("musecvertices_track_dedxharmonic2", musecvertices_track_dedxharmonic2, "musecvertices_track_dedxharmonic2[track_count][2]/F");
  tree->Branch("musecvertices_track_charge", musecvertices_track_charge, "musecvertices_track_charge[musecvertices_count][2]/I");
  tree->Branch("musecvertices_track_nhits", musecvertices_track_nhits, "musecvertices_track_nhits[musecvertices_count][2]/b");
  tree->Branch("musecvertices_track_nmissinghits", musecvertices_track_nmissinghits, "musecvertices_track_nmissinghits[musecvertices_count][2]/b");
  tree->Branch("musecvertices_track_npixelhits", musecvertices_track_npixelhits, "musecvertices_track_npixelhits[musecvertices_count][2]/b");
  tree->Branch("musecvertices_track_npixellayers", musecvertices_track_npixellayers, "musecvertices_track_npixellayers[musecvertices_count][2]/b");
  tree->Branch("musecvertices_track_nstriplayers", musecvertices_track_nstriplayers, "musecvertices_track_nstriplayers[musecvertices_count][2]/b");
  
  tree->Branch("genweight", &genweight, "genweight/F");
  tree->Branch("genid1", &genid1, "genid1/F");
  tree->Branch("genx1", &genx1, "genx1/F");
  tree->Branch("genid2", &genid2, "genid2/F");
  tree->Branch("genx2", &genx2, "genx2/F");
  tree->Branch("genScale", &genScale, "genScale/F");
  
  tree->Branch("numpileupinteractionsminus", &numpileupinteractionsminus, "numpileupinteractionsminus/I");
  tree->Branch("numpileupinteractions", &numpileupinteractions, "numpileupinteractions/I");
  tree->Branch("numpileupinteractionsplus", &numpileupinteractionsplus, "numpileupinteractionsplus/I");
  tree->Branch("numtruepileupinteractions", &numtruepileupinteractions, "numtruepileupinteractions/F");
  
  tree->Branch("genparticles_count", &genparticles_count, "genparticles_count/i");
  tree->Branch("genparticles_e", genparticles_e, "genparticles_e[genparticles_count]/F");
  tree->Branch("genparticles_px", genparticles_px, "genparticles_px[genparticles_count]/F");
  tree->Branch("genparticles_py", genparticles_py, "genparticles_py[genparticles_count]/F");
  tree->Branch("genparticles_pz", genparticles_pz, "genparticles_pz[genparticles_count]/F");
  tree->Branch("genparticles_vx", genparticles_vx, "genparticles_vx[genparticles_count]/F");
  tree->Branch("genparticles_vy", genparticles_vy, "genparticles_vy[genparticles_count]/F");
  tree->Branch("genparticles_vz", genparticles_vz, "genparticles_vz[genparticles_count]/F");
  tree->Branch("genparticles_pdgid", genparticles_pdgid, "genparticles_pdgid[genparticles_count]/I");
  tree->Branch("genparticles_status", genparticles_status, "genparticles_status[genparticles_count]/I");
  tree->Branch("genparticles_info", genparticles_info, "genparticles_info[genparticles_count]/i");

  tree->Branch("genallparticles_count", &genallparticles_count, "genallparticles_count/i");
  tree->Branch("genallparticles_e", genallparticles_e, "genallparticles_e[genallparticles_count]/F");
  tree->Branch("genallparticles_px", genallparticles_px, "genallparticles_px[genallparticles_count]/F");
  tree->Branch("genallparticles_py", genallparticles_py, "genallparticles_py[genallparticles_count]/F");
  tree->Branch("genallparticles_pz", genallparticles_pz, "genallparticles_pz[genallparticles_count]/F");
  tree->Branch("genallparticles_vx", genallparticles_vx, "genallparticles_vx[genallparticles_count]/F");
  tree->Branch("genallparticles_vy", genallparticles_vy, "genallparticles_vy[genallparticles_count]/F");
  tree->Branch("genallparticles_vz", genallparticles_vz, "genallparticles_vz[genallparticles_count]/F");
  tree->Branch("genallparticles_pdgid", genallparticles_pdgid, "genallparticles_pdgid[genallparticles_count]/I");
  tree->Branch("genallparticles_charge", genallparticles_charge, "genallparticles_charge[genallparticles_count]/I");
  tree->Branch("genallparticles_status", genallparticles_status, "genallparticles_status[genallparticles_count]/I");
  tree->Branch("genallparticles_motherbeg", genallparticles_motherbeg, "genallparticles_motherbeg[genallparticles_count]/i");
  tree->Branch("genallparticles_daughterbeg", genallparticles_daughterbeg, "genallparticles_daughterbeg[genallparticles_count]/i");
  
  tree->Branch("genallparticlesmother_count", &genallparticlesmother_count, "genallparticlesmother_count/i");
  tree->Branch("genallparticles_mothers", genallparticles_mothers, "genallparticles_mothers[genallparticlesmother_count]/i");
  
  tree->Branch("genallparticlesdaughter_count", &genallparticlesdaughter_count, "genallparticlesdaughter_count/i");
  tree->Branch("genallparticles_daughters", genallparticles_daughters, "genallparticles_daughters[genallparticlesdaughter_count]/i");

  tree->Branch("genmetcalo_ex", &genmetcalo_ex, "genmetcalo_ex/F");
  tree->Branch("genmetcalo_ey", &genmetcalo_ey, "genmetcalo_ey/F");
  tree->Branch("genmettrue_ex", &genmettrue_ex, "genmettrue_ex/F");
  tree->Branch("genmettrue_ey", &genmettrue_ey, "genmettrue_ey/F");
  
  lumitree = FS->make<TTree>("AC1Blumi", "AC1Blumi", 1);
  
  lumitree->Branch("lumi_run", &lumi_run, "lumi_run/i");
  lumitree->Branch("lumi_block", &lumi_block, "lumi_block/i");
  
  if(cdata)
    {
      lumitree->Branch("lumi_value", &lumi_value, "lumi_value/F");
      lumitree->Branch("lumi_valueerr", &lumi_valueerr, "lumi_valueerr/F");
      lumitree->Branch("lumi_livefrac", &lumi_livefrac, "lumi_livefrac/F");
      lumitree->Branch("lumi_deadfrac", &lumi_deadfrac, "lumi_deadfrac/F");
      lumitree->Branch("lumi_quality", &lumi_quality, "lumi_quality/i");
      lumitree->Branch("lumi_eventsprocessed", &lumi_eventsprocessed, "lumi_eventsprocessed/i");
      lumitree->Branch("lumi_eventsfiltered", &lumi_eventsfiltered, "lumi_eventsfiltered/i");
      lumitree->Branch("lumi_hltprescaletable", &lumi_hltprescaletable, "lumi_hltprescaletable/i");
      lumitree->Branch("lumi_l1algoprescaletable", &lumi_l1algoprescaletable, "lumi_l1algoprescaletable/i");
      lumitree->Branch("lumi_l1techprescaletable", &lumi_l1techprescaletable, "lumi_l1techprescaletable/i");
    }
  
  runtree = FS->make<TTree>("AC1Brun", "AC1Brun", 1);
  runtree->Branch("run_number", &run_number, "run_number/i");
  runtree->Branch("run_hltcount", &run_hltcount, "run_hltcount/i");
  runtree->Branch("run_hltnames", run_hltnames, "run_hltnames/C");
  runtree->Branch("run_hltmunames", run_hltmunames, "run_hltmunames/C");
  runtree->Branch("run_hltelnames", run_hltelnames, "run_hltelnames/C");
  runtree->Branch("run_hlttaunames", run_hlttaunames, "run_hlttaunames/C");
  runtree->Branch("run_hltphotonnames", run_hltphotonnames, "run_hltphotonnames/C");
  runtree->Branch("run_hltjetnames", run_hltjetnames, "run_hltjetnames/C");
  runtree->Branch("run_taudiscriminators", run_taudiscriminators, "run_taudiscriminators/C");
  runtree->Branch("run_hltprescaletablescount", &run_hltprescaletablescount, "run_hltprescaletablescount/i");
  runtree->Branch("run_hltprescaletables", run_hltprescaletables, "run_hltprescaletables[run_hltprescaletablescount]/i");
  runtree->Branch("run_l1algocount", &run_l1algocount, "run_l1algocount/i");
  runtree->Branch("run_l1algoprescaletablescount", &run_l1algoprescaletablescount, "run_l1algoprescaletablescount/i");
  runtree->Branch("run_l1algoprescaletables", run_l1algoprescaletables, "run_l1algoprescaletables[run_l1algoprescaletablescount]/i");
  runtree->Branch("run_l1techcount", &run_l1techcount, "run_l1techcount/i");
  runtree->Branch("run_l1techprescaletablescount", &run_l1techprescaletablescount, "run_l1techprescaletablescount/i");
  runtree->Branch("run_l1techprescaletables", run_l1techprescaletables, "run_l1techprescaletables[run_l1techprescaletablescount]/i");
}//void RootMaker::beginJob()



//http://www.boost.org/doc/libs/1_55_0/libs/regex/doc/html/boost_regex/introduction_and_overview.html
void AddTriggerList(const edm::Run& run, 
		    const HLTConfigProvider& hltConfig,
                    const std::vector<std::pair<boost::regex, std::string> >& regexes,
                    std::vector<std::pair<std::string, std::string> >& foundTriggers,
                    std::string& triggerNames   )
{
  for(std::size_t i = 0; i < hltConfig.size(); ++i) {
    
    // In some early 2011 runs saveTagsModules() does not give the
    // modules which actually save tags, but an empty list. So check all
    // the modules instead.
    
    //const std::vector<std::string> saveTagsModules = hltConfig.saveTagsModules(i);
    const std::vector<std::string> saveTagsModules = hltConfig.moduleLabels(i);
    
    for(std::size_t j = 0; j < regexes.size(); ++j) {
      boost::cmatch what;
      
      if(boost::regex_match(hltConfig.triggerName(i), regexes[j].first))
	//if(boost::regex_match(HLTConfiguration.triggerName(i).c_str(), what, muonregexes[j].first) && muontriggers.size() < 32)
	{
	  // Check for filter
	  if(regexes[j].second.size() != 0)
	      {
		// First, check for leg
		std::string::size_type legpos = regexes[j].second.find('|');
		std::string leg, filters;
		if(legpos != std::string::npos)
		  {
		    leg     = regexes[j].second.substr(0, legpos); 
		    filters = regexes[j].second.substr(legpos + 1); 
		  }
		else
		  {
		    filters = regexes[j].second;
		  }
		
		std::vector<std::string> strs;
		boost::split(strs, filters, boost::is_any_of(","));
		bool foundFilter = false;
		std::string filter;
		
		for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)
		  {
		    for(std::size_t l = 0; l < strs.size() && !foundFilter; ++l)
		      {
			if(saveTagsModules[k] == strs[l])
			  {
			    filter = saveTagsModules[k];
			    foundFilter = true;
			  }
		      }
		  }//for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)
		
		if(!foundFilter)
		  {
		    for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)
		      std::cout << saveTagsModules[k] << std::endl;
		    throw cms::Exception("RootMaker") << "Did not find filter for trigger " << hltConfig.triggerName(i) << " in run " << run.run() << std::endl;
		  }
		
		if(leg.empty())
		  triggerNames += hltConfig.triggerName(i) + string(":") + filter + string(" ");
		else
		  triggerNames += hltConfig.triggerName(i) + string("|") + leg + string(":") + filter + string(" ");
		
		foundTriggers.push_back(std::make_pair(hltConfig.triggerName(i), filter));
	      }//if(regexes[j].second.size() != 0)
	  
	  else
	    {
	      triggerNames += hltConfig.triggerName(i) + string(" ");
	      foundTriggers.push_back(std::make_pair(hltConfig.triggerName(i), saveTagsModules.back()));
	    }
	}
    }//for(std::size_t j = 0; j < regexes.size(); ++j)
  }//for(std::size_t i = 0; i < hltConfig.size(); ++i)
}// std::string& triggerNames


void RootMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  if(propagatorWithMaterial != 0){ delete propagatorWithMaterial;}
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  propagatorWithMaterial = new PropagatorWithMaterial(alongMomentum, 0.10566, &(*magneticField));
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", TTrackBuilder);
  
  run_number = iRun.run();
  //L1 prescales
  edm::ESHandle<L1GtPrescaleFactors> l1GtPfAlgo;
  iSetup.get<L1GtPrescaleFactorsAlgoTrigRcd>().get(l1GtPfAlgo);
  
  unsigned numl1algo = (l1GtPfAlgo.product()->gtPrescaleFactors())[0].size();
  unsigned numl1algotables = (l1GtPfAlgo.product()->gtPrescaleFactors()).size();
  
  run_l1algoprescaletablescount = numl1algo*numl1algotables;
  run_l1algocount = numl1algo;
  if(l1GtPfAlgo.isValid())
    {
      for(unsigned i = 0 ; i < numl1algotables ; i++)
	{
	  for(unsigned j = 0 ; j < numl1algo ; j++)
	    {
	      run_l1algoprescaletables[j+numl1algo*i] = (l1GtPfAlgo.product()->gtPrescaleFactors())[i][j];
	    }
	}
    }
  
  edm::ESHandle<L1GtPrescaleFactors> l1GtPfTech;
  iSetup.get<L1GtPrescaleFactorsTechTrigRcd>().get(l1GtPfTech);
  
  unsigned numl1tech = (l1GtPfTech.product()->gtPrescaleFactors())[0].size();
  unsigned numl1techtables = (l1GtPfTech.product()->gtPrescaleFactors()).size();
  
  run_l1techprescaletablescount = numl1tech*numl1techtables;
  run_l1techcount = numl1tech;
  if(l1GtPfTech.isValid())
    {
      for(unsigned i = 0 ; i < numl1techtables ; i++)
	{
	  for(unsigned j = 0 ; j < numl1tech ; j++)
	    {
	      run_l1techprescaletables[j+numl1tech*i] = (l1GtPfTech.product()->gtPrescaleFactors())[i][j];
	    }
	}
    }
  
  //HLT names and prescales
  muontriggers.clear();
  electrontriggers.clear();
  tautriggers.clear();
  photontriggers.clear();
  jettriggers.clear();
  HLTriggerIndexSelection.clear();
  
  bool changed = true;
  if(ctrigger)    HLTConfiguration.init(iRun, iSetup, cTriggerProcess, changed);
  
  vector<boost::regex> trigregexes;
  for(unsigned i = 0 ; i < cHLTriggerNamesSelection.size() ; i++)
    {
      trigregexes.push_back(boost::regex(cHLTriggerNamesSelection[i].c_str()));
    }
  vector<pair<boost::regex, string> > muonregexes;
  for(unsigned i = 0 ; i < cMuHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cMuHLTriggerMatching[i], boost::is_any_of(":"));
      if(strs.size() == 1) strs.push_back(string(""));
      muonregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }
  vector<pair<boost::regex, string> > electronregexes;
  for(unsigned i = 0 ; i < cElHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cElHLTriggerMatching[i], boost::is_any_of(":"));
      if(strs.size() == 1) strs.push_back(string(""));
      electronregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }
  vector<pair<boost::regex, string> > tauregexes;
  for(unsigned i = 0 ; i < cTauHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cTauHLTriggerMatching[i], boost::is_any_of(":"));	
      if(strs.size() == 1) strs.push_back(string(""));
      tauregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }
	vector<pair<boost::regex, string> > photonregexes;
	for(unsigned i = 0 ; i < cPhotonHLTriggerMatching.size() ; i++)
	  {
	    vector<string> strs;
	    boost::split(strs, cPhotonHLTriggerMatching[i], boost::is_any_of(":"));	
	    if(strs.size() == 1) strs.push_back(string(""));
	    photonregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
	  }
	vector<pair<boost::regex, string> > jetregexes;
	for(unsigned i = 0 ; i < cJetHLTriggerMatching.size() ; i++)
	  {
	    vector<string> strs;
	    boost::split(strs, cJetHLTriggerMatching[i], boost::is_any_of(":"));	
	    if(strs.size() == 1) strs.push_back(string(""));
	    jetregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
	  }
	
	run_hltcount = HLTConfiguration.size();
	string allnames;
	string allmuonnames;
	string allelectronnames;
	string alltaunames;
	string allphotonnames;
	string alljetnames;
	
	AddTriggerList(iRun, HLTConfiguration, muonregexes,     muontriggers,     allmuonnames);
	AddTriggerList(iRun, HLTConfiguration, electronregexes, electrontriggers, allelectronnames);
	AddTriggerList(iRun, HLTConfiguration, tauregexes,      tautriggers,      alltaunames);
	AddTriggerList(iRun, HLTConfiguration, photonregexes,   photontriggers,   allphotonnames);
	AddTriggerList(iRun, HLTConfiguration, jetregexes,      jettriggers,      alljetnames);

	
	for(unsigned i = 0; i < HLTConfiguration.size(); ++i)
	  {
	    for(unsigned j = 0 ; j < trigregexes.size() ; j++)
	      {
		if(boost::regex_match(HLTConfiguration.triggerName(i), trigregexes[j]))
		  {
		    HLTriggerIndexSelection.push_back(i);
		  }
	      }
	    allnames += HLTConfiguration.triggerName(i) + string(" ");
	  }
	if(allnames.size() > 20000)
	  {
	    throw cms::Exception("Too long trigger names");
	    //allnames = string("TOOLONGTRIGGERNAMES::ERROR ");
	  }
	if(muontriggers.size() > 32) throw cms::Exception("RootMaker") << "Too many muon triggers!" << std::endl;
	if(electrontriggers.size() > 32) throw cms::Exception("RootMaker") << "Too many electron triggers!" << std::endl;
	if(tautriggers.size() > 32) throw cms::Exception("RootMaker") << "Too many tau triggers!" << std::endl;
	if(photontriggers.size() > 32) throw cms::Exception("RootMaker") << "Too many photon triggers!" << std::endl;
	if(jettriggers.size() > 32) throw cms::Exception("RootMaker") << "Too many jet triggers!" << std::endl;
	
	strcpy(run_hltnames, allnames.c_str());
	strcpy(run_hltmunames, allmuonnames.c_str());
	strcpy(run_hltelnames, allelectronnames.c_str());
	strcpy(run_hlttaunames, alltaunames.c_str());
	strcpy(run_hltphotonnames, allphotonnames.c_str());
	strcpy(run_hltjetnames, alljetnames.c_str());
	string alltaudiscriminators;
	if(cTauDiscriminators.size() > sizeof(tau_dishps[0])*8)
	  throw cms::Exception("Too many tau discrimintors selected");
	for(unsigned i = 0 ; i < cTauDiscriminators.size() ; i++)
	  {
	    alltaudiscriminators += cTauDiscriminators[i] + string(" ");
	  }
	strcpy(run_taudiscriminators, alltaudiscriminators.c_str());
	
	run_hltprescaletablescount = HLTConfiguration.prescaleSize()*HLTConfiguration.size();
	for(unsigned j = 0 ; j < HLTConfiguration.prescaleSize() ; j++)
	  {
	    for(unsigned i = 0 ; i < HLTConfiguration.size() ; i++)
	      {
		run_hltprescaletables[i+HLTConfiguration.size()*j] = HLTConfiguration.prescaleValue(j, HLTConfiguration.triggerName(i));
	      }
	  }	
	runtree->Fill();
}

void RootMaker::beginLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup)
{
  lumi_run = iLumiBlock.run();
  lumi_block = iLumiBlock.luminosityBlock();
  
  if(cdata)
    {
      //edm::Handle<LumiSummary> lumiSummary;
      //iLumiBlock.getByLabel(edm::InputTag("lumiProducer"), lumiSummary);
      lumi_value = 0;//lumiSummary->avgInsDelLumi();
      lumi_valueerr = 0;//lumiSummary->avgInsDelLumiErr();
      lumi_livefrac = 0;//lumiSummary->lumiSecQual();
      lumi_deadfrac = 0;//lumiSummary->deadFrac();
      lumi_quality = 0;//lumiSummary->liveFrac();
      lumi_eventsprocessed = 0;
      lumi_eventsfiltered = 0;
    }
}

void RootMaker::endLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup)
{
  lumitree->Fill();
}


void RootMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(doDebug)  cout<<"inside the analyze function"<< endl;
  
  track_count = 0;
  primvertex_count = 0;
  supercluster_count = 0;
  supercluster_basiccluster_count = 0;
  supercluster_basiccluster_hit_count = 0;
  supercluster_escluster_count = 0;
  supercluster_escluster_hit_count = 0;
  muon_count = 0;
  tau_count = 0;
  ditau_Index = 0;
  tau_charged_count = 0;
  ak5calojet_count = 0;
  ak5jptjet_count = 0;
  ak5pfjet_count = 0;
  electron_count = 0;
  photon_count = 0;
  conversion_count = 0;
  secvertices_count = 0;
  musecvertices_count = 0;
  allconversion_count = 0;
  genallparticles_count = 0;
  genparticles_count = 0;
  genallparticlesmother_count = 0;
  genallparticlesdaughter_count = 0;
  errors = 0;
  bool takeevent = false;
  diTauCounter = -1;

  nEvents->Fill(0);
  pv_position = math::XYZPoint(0.,0.,0.);
  
  lumi_eventsprocessed++;
  
  event_nr      = iEvent.id().event();
  event_run      = iEvent.id().run();
  event_timeunix = iEvent.time().unixTime();
  event_timemicrosec = iEvent.time().microsecondOffset();
  event_luminosityblock = iEvent.getLuminosityBlock().luminosityBlock();
  
  
  // L1TriggerBits
  // https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h?v=CMSSW_6_2_0_SLHC2#04
  edm::Handle<L1GlobalTriggerReadoutRecord> L1trigger;
  iEvent.getByLabel(edm::InputTag("gtDigis"), L1trigger);
  
  const TechnicalTriggerWord& L1triggerbits = L1trigger->technicalTriggerWord();
  for(int i  = 0  ; i < 8 ; i++) trigger_level1bits[i] = 0;
  
  for(unsigned i = 0 ; i < min(unsigned(L1triggerbits.size()), unsigned(64)) ; i++)
    trigger_level1bits[i/8] |= (Byte_t)L1triggerbits[i] << (i % 8);  // bitwise OR -> | 
  
  //trigger_level1bits[i/8] = trigger_level1bits[i/8]  | (Byte_t) L1triggerbits[i] << (i % 8);  // bitwise OR -> | 
  
  
  //L1TriggerAlgos
  const DecisionWord& L1triggeralgos = L1trigger->decisionWord();
  for(int i = 0  ; i < 128 ; i++){trigger_level1[i] = 0;}
  for(unsigned i = 0 ; i < min(unsigned(L1triggeralgos.size()), unsigned(1024)) ; i++)
    {
      trigger_level1[i/8] |= (Byte_t)L1triggeralgos[i] << (i%8);
    }
  lumi_l1techprescaletable = (L1trigger->gtFdlWord()).gtPrescaleFactorIndexTech();
  lumi_l1algoprescaletable = (L1trigger->gtFdlWord()).gtPrescaleFactorIndexAlgo();	
  lumi_hltprescaletable = -1;
  
  //HLTriggerResults
  iEvent.getByLabel(edm::InputTag("TriggerResults", "", cTriggerProcess), HLTrigger);
  for(int i = 0  ; i < 128 ; i++){trigger_HLT[i] = 0;}
  
  assert(HLTrigger->size() <= 1024);
  for(unsigned i = 0 ; i < min(unsigned(HLTrigger->size()), unsigned(1024)) ; i++)
    {
      trigger_HLT[i/8] |= (Byte_t)HLTrigger->accept(i) << (i%8);
      if(HLTrigger->accept(i) && find(HLTriggerIndexSelection.begin(), HLTriggerIndexSelection.end(), i) != HLTriggerIndexSelection.end())
	{ takeevent = true;}
    }
  
  //TriggerEvent for matching
  iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", cTriggerProcess), HLTriggerEvent);
  //cout<<"PASSED!"<<endl;
  if(HLTriggerIndexSelection.size() == 0 || !ctrigger)
    {
      takeevent = true;
    }
  
  if(!takeevent) return;
  
  if(cbeamspot)
    {
      edm::Handle<BeamSpot> TheBeamSpot;
      iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), TheBeamSpot);
      if(TheBeamSpot.isValid())
	{
	  beamspot_x = TheBeamSpot->x0();
	  beamspot_y = TheBeamSpot->y0();
	  beamspot_z = TheBeamSpot->z0();
	  beamspot_xwidth = TheBeamSpot->BeamWidthX();
	  beamspot_ywidth = TheBeamSpot->BeamWidthY();
	  beamspot_zsigma = TheBeamSpot->sigmaZ();
	  beamspot_cov[0] = TheBeamSpot->covariance(0,0); 
	  beamspot_cov[1] = TheBeamSpot->covariance(0,1); 
	  beamspot_cov[2] = TheBeamSpot->covariance(0,2); 
	  beamspot_cov[3] = TheBeamSpot->covariance(1,1); 
	  beamspot_cov[4] = TheBeamSpot->covariance(1,2); 
	  beamspot_cov[5] = TheBeamSpot->covariance(2,2); 
	  pv_position = math::XYZPoint(TheBeamSpot->x0(), TheBeamSpot->y0(), TheBeamSpot->z0());
	}
      else
	{
	  beamspot_x = 0.;
	  beamspot_y = 0.;
	  beamspot_z = 0.;
	  beamspot_xwidth = 0.;
	  beamspot_ywidth = 0.;
	  beamspot_zsigma = 0.;
	  beamspot_cov[0] = 0.;
	  beamspot_cov[1] = 0.;
	  beamspot_cov[2] = 0.;
	  beamspot_cov[3] = 0.;
	  beamspot_cov[4] = 0.;
	  beamspot_cov[5] = 0.;
	}
    }
  
  if(crecprimvertex)
    {
      edm::Handle<VertexCollection> Vertex;
      iEvent.getByLabel(edm::InputTag("offlinePrimaryVertices"), Vertex);
      if(Vertex.isValid())
	{
	  for(unsigned i = 0 ; i < Vertex->size(); i++)
	    {
	      if((*Vertex)[i].isValid() && !(*Vertex)[i].isFake())
		{
		  if(primvertex_count == 0)
		    {
		      primvertex_x = (*Vertex)[i].x();
		      primvertex_y = (*Vertex)[i].y();
		      primvertex_z = (*Vertex)[i].z();
		      primvertex_chi2 = (*Vertex)[i].chi2();
		      primvertex_ndof = (*Vertex)[i].ndof();
		      primvertex_ntracks = (*Vertex)[i].tracksSize();
		      primvertex_cov[0] = (*Vertex)[i].covariance(0,0); // xError()
		      primvertex_cov[1] = (*Vertex)[i].covariance(0,1); 
		      primvertex_cov[2] = (*Vertex)[i].covariance(0,2);
		      primvertex_cov[3] = (*Vertex)[i].covariance(1,1); // yError()
		      primvertex_cov[4] = (*Vertex)[i].covariance(1,2);
		      primvertex_cov[5] = (*Vertex)[i].covariance(2,2); // zError()
		      Float_t ptq = 0.;
		      for(Vertex::trackRef_iterator it = (*Vertex)[i].tracks_begin() ; it != (*Vertex)[i].tracks_end() ; ++it)
			{
			  ptq += (*it)->pt() * (*it)->pt();
			}
		      primvertex_ptq = ptq;
		      
		      pv_position = (*Vertex)[0].position();
		      primvertex = (*Vertex)[0];
		    }
		  
		  primvertex_count++;
		}
	    }
	}
    }
  
  if(crecsupercluster)
    {
      bool error = false;
      edm::Handle<SuperClusterCollection> SCbarrel;
      iEvent.getByLabel(edm::InputTag("correctedHybridSuperClusters"), SCbarrel);
      edm::Handle<SuperClusterCollection> SCendcap;
      iEvent.getByLabel(edm::InputTag("correctedMulti5x5SuperClustersWithPreshower"), SCendcap);
      
      edm::ESHandle<CaloGeometry> caloGeo;
      edm::Handle<EcalRecHitCollection> barrelHits;
      edm::Handle<EcalRecHitCollection> endcapHits;
      edm::Handle<EcalRecHitCollection> esHits;
      if(crecsuperclusterhit)
	{
	  iSetup.get<CaloGeometryRecord>().get(caloGeo);
	  iEvent.getByLabel(edm::InputTag("ecalRecHit", "EcalRecHitsEB"), barrelHits);
	  iEvent.getByLabel(edm::InputTag("ecalRecHit", "EcalRecHitsEE"), endcapHits);
	  iEvent.getByLabel(edm::InputTag("ecalPreshowerRecHit", "EcalRecHitsES"), esHits);
	}
      
      if(SCbarrel.isValid())
	{
	  for(unsigned i = 0 ; i < SCbarrel->size() ; i++)
	    {
	      supercluster_e[supercluster_count] = (*SCbarrel)[i].energy();
	      supercluster_x[supercluster_count] = (*SCbarrel)[i].x();
	      supercluster_y[supercluster_count] = (*SCbarrel)[i].y();
	      supercluster_z[supercluster_count] = (*SCbarrel)[i].z();
	      supercluster_rawe[supercluster_count] = (*SCbarrel)[i].rawEnergy();
	      supercluster_phiwidth[supercluster_count] = (*SCbarrel)[i].phiWidth();
	      supercluster_etawidth[supercluster_count] = (*SCbarrel)[i].etaWidth();
	      supercluster_nbasiccluster[supercluster_count] = (*SCbarrel)[i].clustersSize();
	      if(crecsuperclustermember)
		{
		  supercluster_basicclusterbegin[supercluster_count] = supercluster_basiccluster_count;
		  for(CaloCluster_iterator CC = (*SCbarrel)[i].clustersBegin() ; CC != (*SCbarrel)[i].clustersEnd() ; CC++)
		    {
		      supercluster_basiccluster_e[supercluster_basiccluster_count] = (*CC)->energy();
		      supercluster_basiccluster_x[supercluster_basiccluster_count] = (*CC)->x();
		      supercluster_basiccluster_y[supercluster_basiccluster_count] = (*CC)->y();
		      supercluster_basiccluster_z[supercluster_basiccluster_count] = (*CC)->z();
		      supercluster_basiccluster_nhit[supercluster_basiccluster_count] = (*CC)->size();
		      if(crecsuperclusterhit && barrelHits.isValid())
			{
			  supercluster_basiccluster_hitbegin[supercluster_basiccluster_count] = supercluster_basiccluster_hit_count;
			  const std::vector< std::pair<DetId, float> >& hits = (*CC)->hitsAndFractions();
			  for(unsigned m = 0 ; m < hits.size() ; m++)
			    { 
			      float energy = -1;
			      for(unsigned u = 0 ; u < barrelHits->size() ; u++)
				{
				  if((*barrelHits)[u].id() == hits[m].first) {energy = (*barrelHits)[u].energy(); break;}
				} 
			      supercluster_basiccluster_hit_e[supercluster_basiccluster_hit_count] = energy*hits[m].second;
			      supercluster_basiccluster_hit_x[supercluster_basiccluster_hit_count] = caloGeo->getPosition(hits[m].first).x();
			      supercluster_basiccluster_hit_y[supercluster_basiccluster_hit_count] = caloGeo->getPosition(hits[m].first).y();
			      supercluster_basiccluster_hit_z[supercluster_basiccluster_hit_count] = caloGeo->getPosition(hits[m].first).z();
			      supercluster_basiccluster_hit_count++;
			      if(supercluster_basiccluster_hit_count == M_superclusterhitmaxcount) {error = true; break;}
			    }
			}
		      supercluster_basiccluster_count++;
		      if(error || supercluster_basiccluster_count == M_superclustermembermaxcount){error = true; break;}
		    }
		}
	      supercluster_count++;
	      if(error || supercluster_count == M_superclustermaxcount){cerr << "Error filling SuperClusters. They are missing." << endl; errors |= 1<<12; break;}
	    }
	}
      
      //EndCap Clusters
      if(SCendcap.isValid() && !error)
	{
	  for(unsigned i = 0 ; i < SCendcap->size() ; i++)
	    {
	      supercluster_e[supercluster_count] = (*SCendcap)[i].energy();
	      supercluster_x[supercluster_count] = (*SCendcap)[i].x();
	      supercluster_y[supercluster_count] = (*SCendcap)[i].y();
	      supercluster_z[supercluster_count] = (*SCendcap)[i].z();
	      supercluster_rawe[supercluster_count] = (*SCendcap)[i].rawEnergy();
	      supercluster_phiwidth[supercluster_count] = (*SCendcap)[i].phiWidth();
	      supercluster_etawidth[supercluster_count] = (*SCendcap)[i].etaWidth();
	      supercluster_nbasiccluster[supercluster_count] = (*SCendcap)[i].clustersSize();
	      if(crecsuperclustermember)
		{
		  supercluster_basicclusterbegin[supercluster_count] = supercluster_basiccluster_count;
		  for(CaloCluster_iterator CC = (*SCendcap)[i].clustersBegin() ; CC != (*SCendcap)[i].clustersEnd() ; CC++)
		    {
		      supercluster_basiccluster_e[supercluster_basiccluster_count] = (*CC)->energy();
		      supercluster_basiccluster_x[supercluster_basiccluster_count] = (*CC)->x();
		      supercluster_basiccluster_y[supercluster_basiccluster_count] = (*CC)->y();
		      supercluster_basiccluster_z[supercluster_basiccluster_count] = (*CC)->z();
		      supercluster_basiccluster_nhit[supercluster_basiccluster_count] = (*CC)->size();
		      if(crecsuperclusterhit && endcapHits.isValid())
			{
			  supercluster_basiccluster_hitbegin[supercluster_basiccluster_count] = supercluster_basiccluster_hit_count;
			  const std::vector< std::pair<DetId, float> >& hits = (*CC)->hitsAndFractions();
			  for(unsigned m = 0 ; m < hits.size() ; m++)
			    { 
			      float energy = -1;
			      for(unsigned u = 0 ; u < barrelHits->size() ; u++)
				{
				  if((*endcapHits)[u].id() == hits[m].first) {energy = (*endcapHits)[u].energy(); break;}
				} 
			      supercluster_basiccluster_hit_e[supercluster_basiccluster_hit_count] = energy*hits[m].second;
			      supercluster_basiccluster_hit_x[supercluster_basiccluster_hit_count] = caloGeo->getPosition(hits[m].first).x();
			      supercluster_basiccluster_hit_y[supercluster_basiccluster_hit_count] = caloGeo->getPosition(hits[m].first).y();
			      supercluster_basiccluster_hit_z[supercluster_basiccluster_hit_count] = caloGeo->getPosition(hits[m].first).z();
			      supercluster_basiccluster_hit_count++;
			      if(supercluster_basiccluster_hit_count == M_superclusterhitmaxcount) {error = true; break;}
			    }
			}
		      supercluster_basiccluster_count++;
		      if(error || supercluster_basiccluster_count == M_superclustermembermaxcount){error = true; break;}
		    }
		  
		  //filling Preshower
		  supercluster_esclusterbegin[supercluster_count] = supercluster_escluster_count;
		  for(CaloCluster_iterator CC = (*SCendcap)[i].preshowerClustersBegin() ; CC != (*SCendcap)[i].preshowerClustersEnd() ; CC++)
		    {
		      supercluster_escluster_e[supercluster_escluster_count] = (*CC)->energy();
		      supercluster_escluster_x[supercluster_escluster_count] = (*CC)->x();
		      supercluster_escluster_y[supercluster_escluster_count] = (*CC)->y();
		      supercluster_escluster_z[supercluster_escluster_count] = (*CC)->z();
		      supercluster_escluster_nhit[supercluster_escluster_count] = (*CC)->size();
		      if(crecsuperclusterhit && esHits.isValid())
			{
			  supercluster_escluster_hitbegin[supercluster_escluster_count] = supercluster_escluster_hit_count;
			  const std::vector< std::pair<DetId, float> >& hits = (*CC)->hitsAndFractions();
			  for(unsigned m = 0 ; m < hits.size() ; m++)
			    { 
			      float energy = -1;
			      for(unsigned u = 0 ; u < barrelHits->size() ; u++)
				{
				  if((*barrelHits)[u].id() == hits[m].first) {energy = (*esHits)[u].energy(); break;}
				} 
			      supercluster_escluster_hit_e[supercluster_escluster_hit_count] = energy*hits[m].second;
			      supercluster_escluster_hit_x[supercluster_escluster_hit_count] = caloGeo->getPosition(hits[m].first).x();
			      supercluster_escluster_hit_y[supercluster_escluster_hit_count] = caloGeo->getPosition(hits[m].first).y();
			      supercluster_escluster_hit_z[supercluster_escluster_hit_count] = caloGeo->getPosition(hits[m].first).z();
			      supercluster_escluster_hit_count++;
			      if(supercluster_escluster_hit_count == M_superclusterhitmaxcount) {error = true; break;}
			    }
			}
		      supercluster_escluster_count++;
		      if(error || supercluster_escluster_count == M_superclustermembermaxcount){error = true; break;}
		    }
		  
		}
	      supercluster_count++;
	      if(error || supercluster_count == M_superclustermaxcount){cerr << "Error filling SuperClusters. They are missing." << endl; errors |= 1<<12; 
		break;
	      }
	    }
	}
      
    }
  
  if(doDebug)	cout<<" passed the supercluster loop"<< endl;	
  takeevent = false;
  
  int haveElectrons = 0;
  int haveMuons = 0;
  int haveTaus = 0;
  if(crectrack)
    {
      takeevent = AddTracks(iEvent) || takeevent;
    }
  if(crecmuon)
    {
      haveMuons = AddMuons(iEvent);
      takeevent = haveMuons || takeevent;
    }
  if(crecelectron)
    {
      haveElectrons = AddElectrons(iEvent, iSetup);
      takeevent = haveElectrons || takeevent;
    }
  if(crecphoton)
    {
      takeevent = AddPhotons(iEvent) || takeevent;
    }
  if(crectau)
    {
      haveTaus = AddTaus(iEvent, iSetup);
      takeevent = haveTaus || takeevent;
    }
  if(crecmutautaupairs)
    {
      takeevent = AddMuTauTauPairs(iEvent) || takeevent;
    }
  if(creceltautaupairs)
    {
      takeevent = AddElTauTauPairs(iEvent) || takeevent;
    }
  if(crecak5calojet)
    {
      takeevent = AddAK5CaloJets(iEvent, iSetup) || takeevent;
    }
  if(crecak5jptjet)
    {
      takeevent = AddAK5JPTJets(iEvent, iSetup) || takeevent;
    }
  if(crecak5pfjet)
    {
      takeevent = AddAK5PFJets(iEvent, iSetup) || takeevent;
    }
  if(crecmusecvertices)
    {
      AddMuVertices(iEvent);
    }
  if(crecallconversion)
    {
      AddAllConversions(iEvent);
    }
  
  bool doSkim = false;
  if( (cSkim & SKIM_MUTAUTAU) != 0 && haveTaus >= 2 && haveMuons >= 1)
    doSkim = true;
  if( (cSkim & SKIM_ETAUTAU) != 0 && haveTaus >= 2 && haveElectrons >= 1)
    doSkim = true;
  if( (cSkim & SKIM_MUMU) != 0 && haveMuons >= 2)
    doSkim = true;
  if( (cSkim & SKIM_EE) != 0 && haveElectrons >= 2)
    doSkim = true;
  if( (cSkim & SKIM_ETAU) != 0 && haveTaus >= 1 && haveElectrons >= 1)
    doSkim = true;
  if( (cSkim & SKIM_ALL) != 0)
    doSkim = true;
  if( (cSkim & SKIM_EMU) != 0 && haveMuons >= 1 && haveElectrons >= 1)
    doSkim = true;
  if( (cSkim & SKIM_TAUTAU) != 0 && haveTaus >= 2)
    doSkim = true;
  
  if(!doSkim) return;
  takeevent = true;
  
  edm::Handle<double> rho;
  iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho", "ROOTMAKER"), rho);
  ak5pfjet_rho = *rho;
  
  edm::Handle<double> sigma;
  iEvent.getByLabel(edm::InputTag("kt6PFJets", "sigma", "ROOTMAKER"), sigma);
  ak5pfjet_sigma = *sigma;
  
#if 0
  std::cout << "Rho: " << *rho << std::endl;
  
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByLabel(edm::InputTag("isolatedPatTaus"), taus);
  std::cout << "isolatedPatTaus:" << std::endl;
  for(unsigned int i = 0; i < taus->size(); ++i)
    std::cout << "  Tau " << i << ": pt=" << (*taus)[i].pt() << ", eta=" << (*taus)[i].eta() << ", phi=" << (*taus)[i].phi() << std::endl;
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByLabel(edm::InputTag("isolatedPatElectrons"), electrons);
  std::cout << "isolatedPatElectrons:" << std::endl;
  for(unsigned int i = 0; i < electrons->size(); ++i)
    std::cout << "  Electron " << i << ": pt=" << (*electrons)[i].pt() << ", eta=" << (*electrons)[i].eta() << ", phi=" << (*electrons)[i].phi() << std::endl;
#endif
  
  if(creccalomet)
    {
      edm::Handle<reco::CaloMETCollection> caloMet;
      iEvent.getByLabel(edm::InputTag("met"), caloMet);
      if(caloMet.isValid() && caloMet->size() > 0)
	{
	  calomet_ex = (*caloMet)[0].px();
	  calomet_ey = (*caloMet)[0].py();
	}
      else
	{
	  errors |= 1<<22;
	}
      edm::Handle<reco::CaloMETCollection> caloMetMu;
      iEvent.getByLabel(edm::InputTag("corMetGlobalMuons"), caloMetMu);
      if(caloMet.isValid() && caloMet->size() > 0)
	{
	  calomet_exmuons = (*caloMetMu)[0].px();
	  calomet_eymuons = (*caloMetMu)[0].py();
	}
      else
	{
	  errors |= 1<<23;
	}
    }
  
  if(crectcmet)
    {
      edm::Handle<reco::METCollection> tcMet;
      iEvent.getByLabel(edm::InputTag("tcMet"), tcMet);
      if(tcMet.isValid() && tcMet->size() > 0)
	{
	  tcmet_ex = (*tcMet)[0].px();
	  tcmet_ey = (*tcMet)[0].py();
	}
      else
	{
	  errors |= 1<<20;
	}
    }
  
  if(crecpfmet)
    {
      edm::Handle<reco::PFMETCollection> pfMet;
      iEvent.getByLabel(edm::InputTag("pfMet"), pfMet);
      if(pfMet.isValid() && pfMet->size() > 0)
	{
	  pfmet_ex = (*pfMet)[0].px();
	  pfmet_ey = (*pfMet)[0].py();
	}
      else
	{
	  errors |= 1<<21;
	}
      
      edm::Handle<pat::METCollection> pfType1Met;
      iEvent.getByLabel(edm::InputTag("patMETsPFType1"), pfType1Met);
      
      assert(pfType1Met->size() > 0);
      pfmettype1_ex = (*pfType1Met)[0].px();
      pfmettype1_ey = (*pfType1Met)[0].py();
      
      edm::Handle<pat::METCollection> pfMetMVA;
      iEvent.getByLabel(edm::InputTag("patPFMetByMVA"), pfMetMVA);
      
      //std::cout << "MVAMET: " << (*pfMetMVA)[0].pt() << ", phi=" << (*pfMetMVA)[0].phi() << std::endl;
      
      assert(pfMetMVA->size() > 0);
      pfmetmva_ex = (*pfMetMVA)[0].px();
      pfmetmva_ey = (*pfMetMVA)[0].py();
      
      /*		edm::Handle<PFMEtSignCovMatrix> MetSignMatrix;
	iEvent.getByLabel(edm::InputTag("pfMEtSignCovMatrix"), MetSignMatrix);*/
      pfmetsigxx = 0.0f;//(*MetSignMatrix)(0,0);
      pfmetsigxy = 0.0f;//(*MetSignMatrix)(0,1);
      pfmetsigyx = 0.0f;//(*MetSignMatrix)(1,0);
      pfmetsigyy = 0.0f;//(*MetSignMatrix)(1,1);
    }
  
  genweight = 1.;
  numpileupinteractionsminus = -1;
  numpileupinteractions = -1;
  numpileupinteractionsplus = -1;
  numtruepileupinteractions = -1.0f;
  
  if(cgen || cgenallparticles)
    {
      edm::Handle<GenEventInfoProduct> HEPMC;
      iEvent.getByLabel(edm::InputTag("generator"), HEPMC);
      if(HEPMC.isValid())
	{
	  genweight = HEPMC->weight();
	  genid1 = HEPMC->pdf()->id.first;
	  genx1 = HEPMC->pdf()->x.second;
	  genid2 = HEPMC->pdf()->id.first;
	  genx2 = HEPMC->pdf()->x.second;
	  genScale = HEPMC->qScale();
	}
      
      edm::Handle<vector<PileupSummaryInfo> > PUInfo;
      iEvent.getByLabel(edm::InputTag("addPileupInfo"), PUInfo);
      
      if(PUInfo.isValid())
	{
	  for(vector<PileupSummaryInfo>::const_iterator PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI)
	    {
	      int BX = PVI->getBunchCrossing();
	      if(BX == -1)
		{ 
		  numpileupinteractionsminus = PVI->getPU_NumInteractions();
		}
	      else if(BX == 0)
		{ 
		  numpileupinteractions = PVI->getPU_NumInteractions();
		}
	      else if(BX == 1)
		{ 
		  numpileupinteractionsplus = PVI->getPU_NumInteractions();
		}
	      
	      numtruepileupinteractions = PVI->getTrueNumInteractions();
	    }
	}
      
      edm::Handle<GenMETCollection> GenMetCalo;
      iEvent.getByLabel(edm::InputTag("genMetCalo"), GenMetCalo);
      edm::Handle<GenMETCollection> GenMetTrue;
      iEvent.getByLabel(edm::InputTag("genMetTrue"), GenMetTrue);
      if(GenMetCalo.isValid() && GenMetCalo->size() > 0)
	{
	  genmetcalo_ex = (*GenMetCalo)[0].px();
	  genmetcalo_ey = (*GenMetCalo)[0].py();
	}
      else
	{
	  genmetcalo_ex = 0.;
	  genmetcalo_ey = 0.;
	  errors |= 1<<18;
	}	
      if(GenMetTrue.isValid() && GenMetTrue->size() > 0)
	{
	  genmettrue_ex = (*GenMetTrue)[0].px();
	  genmettrue_ey = (*GenMetTrue)[0].py();
	}
      else
	{
	  genmettrue_ex = 0.;
	  genmettrue_ey = 0.;
	  errors |= 1<<19;
	}
      
    }	
  
  // Generated particles information
  if(cgenallparticles)
    {
      edm::Handle<GenParticleCollection> GenParticles;
      iEvent.getByLabel(edm::InputTag("genPlusSimParticles"), GenParticles);
      if(!GenParticles.isValid())
	{
	  iEvent.getByLabel(edm::InputTag("genParticles"), GenParticles);
	}
      
      if(GenParticles.isValid())
	{
	  for(unsigned i = 0 ; i < GenParticles->size() ; i++)
	    {
	      genallparticles_e[genallparticles_count] = (*GenParticles)[i].energy();
	      genallparticles_px[genallparticles_count] = (*GenParticles)[i].px();
	      genallparticles_py[genallparticles_count] = (*GenParticles)[i].py();
	      genallparticles_pz[genallparticles_count] = (*GenParticles)[i].pz();
	      genallparticles_vx[genallparticles_count] = (*GenParticles)[i].vx();
	      genallparticles_vy[genallparticles_count] = (*GenParticles)[i].vy();
	      genallparticles_vz[genallparticles_count] = (*GenParticles)[i].vz();
	      genallparticles_pdgid[genallparticles_count] = (*GenParticles)[i].pdgId();
	      genallparticles_status[genallparticles_count] = (*GenParticles)[i].status();
	      genallparticles_charge[genallparticles_count] = (*GenParticles)[i].charge();
	      
	      genallparticles_count++;
	      if(genallparticles_count == M_genallparticlesmaxcount){
		cerr << "Number of genallparticles > M_genallparticlesmaxcount. They are missing." << endl; 
		errors |= 1<<15; 
		break;
	      }
	      
	    }
	  
	  for(unsigned i = 0 ; i < GenParticles->size() ; i++)
	    {
	      genallparticles_motherbeg[i] = genallparticlesmother_count;
	      genallparticles_daughterbeg[i] = genallparticlesdaughter_count;
	      
	      for(unsigned j = 0 ; j < (*GenParticles)[i].numberOfMothers() ; j++)
		{
		  genallparticles_mothers[genallparticlesmother_count] = FindGenParticle((*GenParticles)[i].mother(j));
		  genallparticlesmother_count++;
		  if(genallparticlesmother_count == M_genmotherdaughtermaxcount){break;}
		}
	      
	      for(unsigned j = 0 ; j < (*GenParticles)[i].numberOfDaughters() ; j++)
		{
		  genallparticles_daughters[genallparticlesdaughter_count] = FindGenParticle((*GenParticles)[i].daughter(j));
		  genallparticlesdaughter_count++;
		  if(genallparticlesdaughter_count == M_genmotherdaughtermaxcount){break;}
		}
	      
	      if(genallparticlesmother_count >= M_genmotherdaughtermaxcount){cerr << "Too many mothers" << endl; errors |= 1<<16; break;}
	      if(genallparticlesdaughter_count >= M_genmotherdaughtermaxcount){cerr << "Too many daughters" << endl; errors |= 1<<17; break;}
	    }
	}
    }
  
  // loop 2
  if(cgen)
    {
      edm::Handle<GenParticleCollection> GenParticles;
      iEvent.getByLabel(edm::InputTag("genPlusSimParticles"), GenParticles);
      if(!GenParticles.isValid())
	{
	  iEvent.getByLabel(edm::InputTag("genParticles"), GenParticles);
	}
      
      if(GenParticles.isValid())
	{
	  for(unsigned i = 0 ; i < GenParticles->size() ; i++)
	    {
	      bool fill = false;
	      UInt_t info = 0;
	      if(abs((*GenParticles)[i].pdgId()) == 13)
		{
		  fill = true;
		  if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0;}
		  if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;}
		  if(HasAnyMother(&(*GenParticles)[i], 6) > 0) {info |= 1<<2;}
		}
	      else if(abs((*GenParticles)[i].pdgId()) == 11)
		{
		  fill = true;
		  if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0;}
		  if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;}
		  if(HasAnyMother(&(*GenParticles)[i], 6) > 0) {info |= 1<<2;}
		}
	      else if(abs((*GenParticles)[i].pdgId()) == 15)
		{
		  fill = true;
		  if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0;}
		  if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;}
		  if(HasAnyMother(&(*GenParticles)[i], 6) > 0) {info |= 1<<2;}
		}
	      else if(abs((*GenParticles)[i].pdgId()) <= 6 && (*GenParticles)[i].pt() > 20.)
		{
		  fill = true;
		}
	      else if(abs((*GenParticles)[i].pdgId()) == 21 && (*GenParticles)[i].pt() > 20.)
		{
		  fill = true;
		}
	      else if(abs((*GenParticles)[i].pdgId()) == 22 && (*GenParticles)[i].pt() > 8. && HasAnyMother(&(*GenParticles)[i], 111) == 0)
		{
		  fill = true;
		  if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0;}
		  if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;}
		}
	      
	      // Save all W/Z/H bosons
	      else if(abs((*GenParticles)[i].pdgId()) == 23 || abs((*GenParticles)[i].pdgId()) == 24 || abs((*GenParticles)[i].pdgId()) == 25)
		{
		  fill = true;
		}
	      else if(abs((*GenParticles)[i].pdgId()) == 111 && (*GenParticles)[i].pt() > 8. && HasAnyMother(&(*GenParticles)[i], 111) == 0)
		{
		  fill = true;
		}
	      if(fill)
		{
		  genparticles_e[genparticles_count] = (*GenParticles)[i].energy();
		  genparticles_px[genparticles_count] = (*GenParticles)[i].px();
		  genparticles_py[genparticles_count] = (*GenParticles)[i].py();
		  genparticles_pz[genparticles_count] = (*GenParticles)[i].pz();
		  genparticles_vx[genparticles_count] = (*GenParticles)[i].vx();
		  genparticles_vy[genparticles_count] = (*GenParticles)[i].vy();
		  genparticles_vz[genparticles_count] = (*GenParticles)[i].vz();
		  genparticles_pdgid[genparticles_count] = (*GenParticles)[i].pdgId();
		  genparticles_status[genparticles_count] = (*GenParticles)[i].status();
		  genparticles_info[genparticles_count] = info;
		  genparticles_count++;
		}
	    } // for(unsigned i = 0 ; i < GenParticles->size() ; i++)
	} // if(GenParticles.isValid())
    } // if(cgen)
  
  if(doDebug)  cout<<"takeevent boolean is :"<<takeevent << endl;
  if(takeevent)
    {
      lumi_eventsfiltered++;
      nEvents->Fill(1);
      tree->Fill();
    }
} //void RootMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)



Int_t RootMaker::HasAnyMother(const GenParticle* particle, int id)
{
  vector<unsigned> bknummother;
  vector<const GenParticle*> bkparticle;
  bknummother.reserve(10);
  bkparticle.reserve(10);
  int level = 0;
  bkparticle.push_back(particle);
  bknummother.push_back(0);
  
  unsigned j = 0;
  while(true)
    {
      if(j == bkparticle[level]->numberOfMothers())
	{
	  level--;
	  if(level == -1){return(0);}
	  j = bknummother[level];
	  bkparticle.resize(level+1);
	  bknummother.resize(level+1);
	  continue;
	}
      
      if(bkparticle[level]->mother(j)->pdgId() == id) return(2);
      if(abs(bkparticle[level]->mother(j)->pdgId()) == abs(id)) return(1);
      
      if(bkparticle[level]->mother(j)->numberOfMothers() > 0)
	{
	  bknummother[level] = j+1;
	  bkparticle.push_back(dynamic_cast<const GenParticle*>(bkparticle[level]->mother(j)));
	  bknummother.push_back(0);
	  j = 0;
	  level++;
	  continue;
	}
      j++;
    }
  return(0);
}

void RootMaker::endJob()
{
}

UInt_t RootMaker::FindGenParticle(const Candidate* particle)
{
  for(unsigned i = 0 ; i < genallparticles_count ; i++)
    {
      if(particle->pdgId() == genallparticles_pdgid[i]     && 
	 particle->status() == genallparticles_status[i]   && 
	 float(particle->energy()) == genallparticles_e[i] && 
	 float(particle->px()) == genallparticles_px[i]    && 
	 float(particle->py()) == genallparticles_py[i]    && 
	 float(particle->pz()) == genallparticles_pz[i])
	{
	  return(i);
	}
    }
	return(genallparticles_count);
}

math::XYZPoint RootMaker::PositionOnECalSurface(TransientTrack& trTrack)
{
	math::XYZPoint ecalPosition(0.,0.,0.);
	const FreeTrajectoryState myTSOS = trTrack.initialFreeState();
	TrajectoryStateOnSurface stateAtECal = propagatorWithMaterial->propagate(myTSOS, *ecalBarrel);

	if(stateAtECal.isValid() && stateAtECal.globalPosition().eta() > 1.479) 
	{
		stateAtECal= propagatorWithMaterial->propagate(myTSOS, *ecalPositiveEtaEndcap);
	}    

	if(stateAtECal.isValid() && stateAtECal.globalPosition().eta() < -1.479) 
	{
		stateAtECal= propagatorWithMaterial->propagate(myTSOS, *ecalNegativeEtaEndcap);
	}    

	if(stateAtECal.isValid())
	{
		ecalPosition = stateAtECal.globalPosition();
	}
	return(ecalPosition);
}

static bool isPFMuon(const reco::Muon& muon, const std::vector<reco::PFCandidate>& pfCandidates)
{
	for(std::vector<reco::PFCandidate>::const_iterator iter = pfCandidates.begin(); iter != pfCandidates.end(); ++iter)
		if(abs(iter->pdgId()) == 13 && iter->muonRef().get() == &muon)
			return true;
	return false;
}

int RootMaker::AddMuons(const edm::Event& iEvent)
{
	edm::Handle<double> hRho;
	edm::InputTag tag("kt6PFJets","rho");
	iEvent.getByLabel(tag,hRho);

	edm::Handle<std::vector<reco::Vertex> > hVertices;
	iEvent.getByLabel(edm::InputTag("offlinePrimaryVertices"), hVertices);
	const reco::Vertex& primaryVtx = (*hVertices)[0];

	edm::Handle<std::vector<reco::PFCandidate> > hPfCandidates;
	iEvent.getByLabel(edm::InputTag("particleFlow"), hPfCandidates);
	const std::vector<reco::PFCandidate>& pfCandidates = *hPfCandidates;

	std::vector<reco::Muon> IdentifiedMuons;
	std::vector<reco::GsfElectron> IdentifiedElectrons;

	edm::Handle<std::vector<reco::GsfElectron> > inElectrons;
	if (!iEvent.getByLabel(edm::InputTag("gsfElectrons"), inElectrons)) throw cms::Exception("gsfElectrons not found");

	//edm::Handle<pat::MuonCollection> Muons;
	//iEvent.getByLabel(edm::InputTag("patMuons"), Muons);
	edm::Handle<MuonCollection> Muons;
	iEvent.getByLabel(edm::InputTag("muons"), Muons);

	for (std::vector<reco::GsfElectron>::const_iterator iE = inElectrons->begin(); iE != inElectrons->end(); ++iE)
	{
	        double electronTrackZ = 0;
	        if (iE->gsfTrack().isNonnull())
	        {
	                electronTrackZ = iE->gsfTrack()->dz(primaryVtx.position());
	        }
	        else if (iE->closestCtfTrackRef().isNonnull())
	        {
	                electronTrackZ = iE->closestCtfTrackRef()->dz(primaryVtx.position());
	        }
	        if(fabs(electronTrackZ) > 0.2)
	                continue;
	        if(fabs(iE->superCluster()->eta())<1.479)
	        {
	                if(iE->pt() > 20)
	                {
	                        if(iE->sigmaIetaIeta() > 0.01) continue;
	                        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
	                        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
	                        if(iE->hadronicOverEm() > 0.15) continue;
	                }
	                else
	                {
	                        if(iE->sigmaIetaIeta() > 0.012) continue;
	                        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
	                        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
	                        if(iE->hadronicOverEm() > 0.15) continue;
	                }
	        }
	        else
	        {
	                if(iE->pt() > 20)
	                {
	                        if(iE->sigmaIetaIeta() > 0.03) continue;
	                        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
	                        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8) continue;
	                }
	                else
	                {
	                        if(iE->sigmaIetaIeta() > 0.032) continue;
	                        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
	                        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8) continue;
	                }
	        }
	        IdentifiedElectrons.push_back(*iE);
	}

	for (MuonCollection::const_iterator iM = Muons->begin(); iM != Muons->end(); ++iM)
	{
	  if(!(iM->innerTrack().isNonnull()))
	    continue;
	  if(!(iM->isGlobalMuon() || iM->isTrackerMuon()))
	    continue;
	  if(iM->innerTrack()->numberOfValidHits() < 11 )
	    continue;
	  IdentifiedMuons.push_back(*iM);
	}

	MuonEffectiveArea::MuonEffectiveAreaTarget currentAreaTarget = MuonEffectiveArea::kMuEANoCorr;
	if (cYear == 2011)
	  {
	    if (iEvent.isRealData())
	      currentAreaTarget = MuonEffectiveArea::kMuEAData2011;
	    else
	      {
		if (cPeriod == "Summer11")
		  currentAreaTarget = MuonEffectiveArea::kMuEASummer11MC;
		if (cPeriod == "Fall11")
		  currentAreaTarget = MuonEffectiveArea::kMuEAFall11MC;
	        }
	  }
	if (cYear == 2012)
	  {
	    if (iEvent.isRealData())
	      currentAreaTarget = MuonEffectiveArea::kMuEAData2012;
	  }
	
	edm::Handle<edm::ValueMap<reco::IsoDeposit> > isoDepsChargedHadrons, isoDepsChargedParticles, isoDepsNeutralHadrons, isoDepsPhotons, isoDepsPU;
	iEvent.getByLabel(edm::InputTag("pfMuIsoDepositChargedHadrons"), isoDepsChargedHadrons);
	iEvent.getByLabel(edm::InputTag("pfMuIsoDepositChargedParticles"), isoDepsChargedParticles);
	iEvent.getByLabel(edm::InputTag("pfMuIsoDepositNeutralHadrons"), isoDepsNeutralHadrons);
	iEvent.getByLabel(edm::InputTag("pfMuIsoDepositPhotons"), isoDepsPhotons);
	iEvent.getByLabel(edm::InputTag("pfMuIsoDepositPU"), isoDepsPU);

	int nGoodMuons = 0;
	if(Muons.isValid())
	{
		for(unsigned i = 0 ; i < Muons->size() ; i++){
			muon_px[i] = (*Muons)[i].px();
			muon_py[i] = (*Muons)[i].py();
			muon_pz[i] = (*Muons)[i].pz();
			if((*Muons)[i].globalTrack().isNonnull())
			{
				muon_pterror[i] = (*Muons)[i].globalTrack()->ptError();
				muon_chi2[i] = (*Muons)[i].globalTrack()->chi2();
				muon_ndof[i] = (*Muons)[i].globalTrack()->ndof();
				muon_nchambers[i] = (*Muons)[i].globalTrack()->hitPattern().numberOfValidMuonHits();
			}
			else
			{	
				muon_pterror[i] = -1.;
				muon_chi2[i] = -1.;
				muon_ndof[i] = 0;
				muon_nchambers[i] = 0;
			}
			muon_isolationr3track[i]  = (*Muons)[i].isolationR03().sumPt;
			muon_isolationr3ntrack[i] = (*Muons)[i].isolationR03().nTracks;
			muon_isolationr3ecal[i]   = (*Muons)[i].isolationR03().emEt;
			muon_isolationr3hcal[i]   = (*Muons)[i].isolationR03().hadEt;

			edm::Ref<std::vector<reco::Muon> > muonref(Muons, i);
			reco::IsoDeposit muonIsoDepositChargedHadrons = (*isoDepsChargedHadrons)[muonref];
			reco::IsoDeposit muonIsoDepositChargedParticles = (*isoDepsChargedParticles)[muonref];
			reco::IsoDeposit muonIsoDepositNeutralHadrons = (*isoDepsNeutralHadrons)[muonref];
			reco::IsoDeposit muonIsoDepositPhotons = (*isoDepsPhotons)[muonref];
			reco::IsoDeposit muonIsoDepositPU = (*isoDepsPU)[muonref];

			reco::isodeposit::Direction dir = reco::isodeposit::Direction((*Muons)[i].eta(), (*Muons)[i].phi());
			reco::isodeposit::ConeVeto pf_cone_veto_charged(dir, 0.0001);
			reco::isodeposit::ThresholdVeto pf_threshold_veto_charged(0.0);
			reco::isodeposit::ConeVeto pf_cone_veto(dir, 0.01);
			reco::isodeposit::ThresholdVeto pf_threshold_veto(0.5);
		
			std::vector<reco::isodeposit::AbsVeto*> vetosPFCharged;
			vetosPFCharged.push_back(&pf_cone_veto_charged);
			vetosPFCharged.push_back(&pf_threshold_veto_charged);
			std::vector<reco::isodeposit::AbsVeto*> vetosPF;
			vetosPF.push_back(&pf_cone_veto);
			vetosPF.push_back(&pf_threshold_veto);

			muon_pfisolationr3_sumchargedhadronpt[i] = muonIsoDepositChargedHadrons.depositWithin(0.3, vetosPFCharged); //(*Muons)[i].pfIsolationR03().sumChargedHadronPt;
			muon_pfisolationr3_sumchargedparticlept[i] = muonIsoDepositChargedParticles.depositWithin(0.3, vetosPFCharged); //(*Muons)[i].pfIsolationR03().sumChargedParticlePt;
			muon_pfisolationr3_sumneutralhadronet[i] = muonIsoDepositNeutralHadrons.depositWithin(0.3, vetosPF); // (*Muons)[i].pfIsolationR03().sumNeutralHadronEt;
			muon_pfisolationr3_sumphotonet[i] = muonIsoDepositPhotons.depositWithin(0.3, vetosPF); //(*Muons)[i].pfIsolationR03().sumPhotonEt;
			muon_pfisolationr3_sumPUpt[i] = muonIsoDepositPU.depositWithin(0.3, vetosPF); // (*Muons)[i].pfIsolationR03().sumPUPt;
			muon_pfisolationr4_sumchargedhadronpt[i] = muonIsoDepositChargedHadrons.depositWithin(0.4, vetosPFCharged); //(*Muons)[i].pfIsolationR04().sumChargedHadronPt;
			muon_pfisolationr4_sumchargedparticlept[i] = muonIsoDepositChargedParticles.depositWithin(0.4, vetosPFCharged); //(*Muons)[i].pfIsolationR04().sumChargedParticlePt;
			muon_pfisolationr4_sumneutralhadronet[i] = muonIsoDepositNeutralHadrons.depositWithin(0.4, vetosPF); // (*Muons)[i].pfIsolationR04().sumNeutralHadronEt;
			muon_pfisolationr4_sumphotonet[i] = muonIsoDepositPhotons.depositWithin(0.4, vetosPF); //(*Muons)[i].pfIsolationR04().sumPhotonEt;
			muon_pfisolationr4_sumPUpt[i] = muonIsoDepositPU.depositWithin(0.4, vetosPF); // (*Muons)[i].pfIsolationR04().sumPUPt;
			muon_ecalenergy[i] = (*Muons)[i].calEnergy().em;
			muon_hcalenergy[i] = (*Muons)[i].calEnergy().had;
			muon_charge[i] = (*Muons)[i].charge();
			muon_type[i] = 0;
			muon_trackermuonquality[i] = 0;
			if((*Muons)[i].isGlobalMuon()){muon_type[i] |= 1 << 0;}
			if((*Muons)[i].isTrackerMuon()){muon_type[i] |= 1 << 1;}
			if((*Muons)[i].isStandAloneMuon()){muon_type[i] |= 1 << 2;}
			if((*Muons)[i].isCaloMuon()){muon_type[i] |= 1 << 3;}
			//if((*Muons)[i].isPFMuon()){muon_type[i] |= 1 << 8;}
			if(isPFMuon((*Muons)[i], pfCandidates)){muon_type[i] |= 1 << 8;}

			{
				using namespace muon;
				if(isGoodMuon((*Muons)[i], All)){muon_trackermuonquality[i] |= 1 << 0;}
				if(isGoodMuon((*Muons)[i], AllGlobalMuons)){muon_trackermuonquality[i] |= 1 << 1;}
				if(isGoodMuon((*Muons)[i], AllStandAloneMuons)){muon_trackermuonquality[i] |= 1 << 2;}
				if(isGoodMuon((*Muons)[i], AllTrackerMuons)){muon_trackermuonquality[i] |= 1 << 3;}
				if(isGoodMuon((*Muons)[i], TrackerMuonArbitrated)){muon_trackermuonquality[i] |= 1 << 4;}
				if(isGoodMuon((*Muons)[i], AllArbitrated)){muon_trackermuonquality[i] |= 1 << 5;}
				if(isGoodMuon((*Muons)[i], GlobalMuonPromptTight)){muon_trackermuonquality[i] |= 1 << 6;}
				if(isGoodMuon((*Muons)[i], TMLastStationLoose)){muon_trackermuonquality[i] |= 1 << 7;}
				if(isGoodMuon((*Muons)[i], TMLastStationTight)){muon_trackermuonquality[i] |= 1 << 8;}
				if(isGoodMuon((*Muons)[i], TM2DCompatibilityLoose)){muon_trackermuonquality[i] |= 1 << 9;}
				if(isGoodMuon((*Muons)[i], TM2DCompatibilityTight)){muon_trackermuonquality[i] |= 1 << 10;}
				if(isGoodMuon((*Muons)[i], TMOneStationLoose)){muon_trackermuonquality[i] |= 1 << 11;}
				if(isGoodMuon((*Muons)[i], TMOneStationTight)){muon_trackermuonquality[i] |= 1 << 12;}
				if(isGoodMuon((*Muons)[i], TMLastStationOptimizedLowPtLoose)){muon_trackermuonquality[i] |= 1 << 13;}
				if(isGoodMuon((*Muons)[i], TMLastStationOptimizedLowPtTight)){muon_trackermuonquality[i] |= 1 << 14;}
				if(isGoodMuon((*Muons)[i], GMTkChiCompatibility)){muon_trackermuonquality[i] |= 1 << 15;}
				if(isGoodMuon((*Muons)[i], GMStaChiCompatibility)){muon_trackermuonquality[i] |= 1 << 16;}
				if(isGoodMuon((*Muons)[i], GMTkKinkTight)){muon_trackermuonquality[i] |= 1 << 17;}
				if(isGoodMuon((*Muons)[i], TMLastStationAngLoose)){muon_trackermuonquality[i] |= 1 << 18;}
				if(isGoodMuon((*Muons)[i], TMLastStationAngTight)){muon_trackermuonquality[i] |= 1 << 19;}
				if(isGoodMuon((*Muons)[i], TMOneStationAngLoose)){muon_trackermuonquality[i] |= 1 << 20;}
				if(isGoodMuon((*Muons)[i], TMOneStationAngTight)){muon_trackermuonquality[i] |= 1 << 21;}
				if(isGoodMuon((*Muons)[i], TMLastStationOptimizedBarrelLowPtLoose)){muon_trackermuonquality[i] |= 1 << 22;}
				if(isGoodMuon((*Muons)[i], TMLastStationOptimizedBarrelLowPtTight)){muon_trackermuonquality[i] |= 1 << 23;}
				if((*Muons)[i].stationMask(Muon::SegmentAndTrackArbitration) & 1<<0) {muon_trackermuonquality[i] |= 1 << 24;}
				if((*Muons)[i].stationMask(Muon::SegmentAndTrackArbitration) & 1<<1) {muon_trackermuonquality[i] |= 1 << 25;}
				if((*Muons)[i].stationMask(Muon::SegmentAndTrackArbitration) & 1<<2) {muon_trackermuonquality[i] |= 1 << 26;}
				if((*Muons)[i].stationMask(Muon::SegmentAndTrackArbitration) & 1<<3) {muon_trackermuonquality[i] |= 1 << 27;}
				if((*Muons)[i].stationMask(Muon::SegmentAndTrackArbitration) & 1<<4) {muon_trackermuonquality[i] |= 1 << 28;}
				if((*Muons)[i].stationMask(Muon::SegmentAndTrackArbitration) & 1<<5) {muon_trackermuonquality[i] |= 1 << 29;}
				if((*Muons)[i].stationMask(Muon::SegmentAndTrackArbitration) & 1<<6) {muon_trackermuonquality[i] |= 1 << 30;}
				if((*Muons)[i].stationMask(Muon::SegmentAndTrackArbitration) & 1<<7) {muon_trackermuonquality[i] |= 1 << 31;}
				muon_numchambers[i] = (*Muons)[i].numberOfChambers();
				muon_numchamberswithsegments[i] = (*Muons)[i].numberOfMatches(Muon::SegmentAndTrackArbitration);
			}

			if((*Muons)[i].time().direction() == 1)
			{
				muon_type[i] |= 1<<6;
			}
			else if((*Muons)[i].time().direction() == -1)
			{
				muon_type[i] |= 1<<7;
			}


			TrackRef innertrack = (*Muons)[i].innerTrack();

			if(innertrack.isNonnull())
			{
				muon_type[i] |= 1 << 4;
				TransientTrack TTrack = TTrackBuilder->build(innertrack);
				TrajectoryStateClosestToPoint TTrackState = TTrack.trajectoryStateClosestToPoint(GlobalPoint(pv_position.x(), pv_position.y(), pv_position.z()));
				edm::Handle<edm::ValueMap<DeDxData> > dEdxharmonic2;
				iEvent.getByLabel(edm::InputTag("dedxHarmonic2"), dEdxharmonic2);

				muon_innertrack_px[muon_count] = TTrackState.momentum().x();
				muon_innertrack_py[muon_count] = TTrackState.momentum().y();
				muon_innertrack_pz[muon_count] = TTrackState.momentum().z();
				muon_innertrack_closestpointx[muon_count] = TTrackState.position().x();
				muon_innertrack_closestpointy[muon_count] = TTrackState.position().y();
				muon_innertrack_closestpointz[muon_count] = TTrackState.position().z();
				muon_innertrack_dxy[muon_count]    = TTrackState.perigeeParameters().transverseImpactParameter();
				muon_innertrack_dxyerr[muon_count] = TTrackState.perigeeError().transverseImpactParameterError();
				muon_innertrack_dz[muon_count]     = TTrackState.perigeeParameters().longitudinalImpactParameter();
				muon_innertrack_dzerr[muon_count]  = TTrackState.perigeeError().longitudinalImpactParameterError();
				muon_innertrack_chi2[muon_count]   = innertrack->chi2();
				muon_innertrack_ndof[muon_count]   = innertrack->ndof();
				muon_innertrack_charge[muon_count] = innertrack->charge();
				muon_innertrack_nhits[muon_count] = innertrack->numberOfValidHits();
				muon_innertrack_nmissinghits[muon_count] = innertrack->numberOfLostHits();
				muon_innertrack_npixelhits[muon_count] = innertrack->hitPattern().numberOfValidPixelHits();
				muon_innertrack_npixellayers[muon_count] = innertrack->hitPattern().pixelLayersWithMeasurement();
				muon_innertrack_nstriplayers[muon_count] = innertrack->hitPattern().stripLayersWithMeasurement();
				muon_innertrack_dedxharmonic2[muon_count] = (*dEdxharmonic2)[innertrack].dEdx();
				math::XYZPoint ecalPos = PositionOnECalSurface(TTrack);
				muon_innertrack_outerx[muon_count] = ecalPos.x();
				muon_innertrack_outery[muon_count] = ecalPos.y();
				muon_innertrack_outerz[muon_count] = ecalPos.z();
			}

			TrackRef muontrack = (*Muons)[i].outerTrack();

			if(muontrack.isNonnull())
			{
				muon_type[i] |= 1 << 5;
				muon_outertrack_px[i] = muontrack->px();
				muon_outertrack_py[i] = muontrack->py();
				muon_outertrack_pz[i] = muontrack->pz();
				muon_outertrack_hits[i] = muontrack->numberOfValidHits();
				muon_outertrack_missinghits[i] = muontrack->numberOfLostHits();
				muon_outertrack_chi2[i] = muontrack->chi2();
				muon_outertrack_ndof[i] = muontrack->ndof();
			}

			muon_trigger[i] = GetTriggerMatch((*Muons)[i], muontriggers);

			muon_mva_id[i] = fMuonIsoMVA.mvaValue( (*Muons)[i], primaryVtx, pfCandidates, *hRho, currentAreaTarget, IdentifiedElectrons, IdentifiedMuons);
			muon_mva_iso[i] = fMuonIDMVA.mvaValue( (*Muons)[i], primaryVtx, pfCandidates, *hRho, currentAreaTarget, IdentifiedElectrons, IdentifiedMuons);

			muon_count++;
			if(muon_count == M_muonmaxcount)
			{
				cerr << "number of muon > M_muonmaxcount. They are missing." << endl; 
				errors |= 1<<0;
				break;
			}

			if((*Muons)[i].pt() >= 10 && fabs((*Muons)[i].eta()) <= 2.5)
			{
				++nGoodMuons;
			}
		}
	}

	return nGoodMuons;
}


bool RootMaker::GetL1ExtraTriggerMatch(const l1extra::L1JetParticleCollection* l1jets,  const l1extra::L1JetParticleCollection* l1taus, const LeafCandidate& leg2) 
{
  bool matched = false;
  //check matching to l1tau 44 or l1jet 64
  if(l1taus)
    {
      //check matching with l1tau Pt>44 |eta|<2.172
      matched = false;
      for(unsigned int i=0; i<l1taus->size(); ++i)
	{
	  if( (*l1taus)[i].pt() < 44 || fabs((*l1taus)[i].eta() ) > 2.172 ) continue;
	  if( ROOT::Math::VectorUtil::DeltaR( (*l1taus)[i].p4(), leg2.p4() )  < 0.5 )
	    {
	      matched = true;
	      break;
	    }
	}// for(unsigned int i=0; i<l1taus->size(); ++i) 
      
      if(!matched)
	{ 
	  if(l1jets){//check matching with l1jet Pt>64 |eta|<2.172
	    for(unsigned int i=0; i < l1jets->size(); ++i)
	      {
		if( (*l1jets)[i].pt() < 64 || fabs((*l1jets)[i].eta() ) > 2.172 ) continue;
		if( ROOT::Math::VectorUtil::DeltaR((*l1jets)[i].p4(), leg2.p4() ) < 0.5 ) {
		  matched = true;
		  break;
		}
	      }//for(unsigned int i=0; i<l1jets->size(); ++i)
	  }
	}
    } //if(l1taus)
  return matched;
}


UInt_t RootMaker::GetTriggerMatch(const LeafCandidate& particle, const vector<std::pair<string, string> >& triggers)
{
  UInt_t result = 0;
  assert(HLTrigger.isValid() && HLTriggerEvent.isValid());
  const trigger::TriggerObjectCollection TOC(HLTriggerEvent->getObjects());
  for(unsigned n = 0 ; n < triggers.size() ; n++)
    {
      const std::string& triggerName = triggers[n].first;
      const std::string& triggerFilter = triggers[n].second;
      const unsigned int triggerIndex = HLTConfiguration.triggerIndex(triggerName);
      if(triggerIndex < HLTrigger->size()) // Trigger might not exist
	{
	  const unsigned int FilterIndex = HLTriggerEvent->filterIndex(edm::InputTag(triggerFilter, "", cTriggerProcess));
	  // If the trigger has fired the filter must have been run
	  assert(FilterIndex < HLTriggerEvent->sizeFilters() || !HLTrigger->accept(triggerIndex));
	  if(FilterIndex < HLTriggerEvent->sizeFilters())
	    {
	      const trigger::Keys KEYS(HLTriggerEvent->filterKeys(FilterIndex));
	      for(unsigned int j = 0; j < KEYS.size(); ++j)
		{
		  const double dr = ROOT::Math::VectorUtil::DeltaR(particle.p4(), TOC[KEYS[j]].particle().p4());
		  const double dpt = fabs(particle.pt() - TOC[KEYS[j]].pt())/particle.pt();
		  if(dr < 0.5/* && dpt < 0.5*/){result |= 1<<n; break;}
		}
	    }
	}
    }
  return(result);
}

bool RootMaker::AddPhotons(const edm::Event& iEvent)
{
	int NumGood = 0;
	//edm::Handle<pat::PhotonCollection> Photons;
	//iEvent.getByLabel(edm::InputTag("patPhotons"), Photons);
	edm::Handle<PhotonCollection> Photons;
	iEvent.getByLabel(edm::InputTag("photons"), Photons);

	if(Photons.isValid())
	{
		for(unsigned i = 0 ; i < Photons->size() ; i++)
		{
			photon_px[i] = (*Photons)[i].px();
			photon_py[i] = (*Photons)[i].py();
			photon_pz[i] = (*Photons)[i].pz();
			photon_e1x5[i] = (*Photons)[i].e1x5();
			photon_e2x5[i] = (*Photons)[i].e2x5();
			photon_e3x3[i] = (*Photons)[i].e3x3();
			photon_e5x5[i] = (*Photons)[i].e5x5();
			photon_maxenergyxtal[i] = (*Photons)[i].maxEnergyXtal();
			photon_sigmaetaeta[i] = (*Photons)[i].sigmaEtaEta();
			photon_sigmaietaieta[i] = (*Photons)[i].sigmaIetaIeta();
			photon_ehcaloverecal[i] = (*Photons)[i].hadronicOverEm();
			photon_ehcaloverecaldepth1[i] = (*Photons)[i].hadronicDepth1OverEm();
			photon_ehcaloverecaldepth2[i] = (*Photons)[i].hadronicDepth2OverEm();
			photon_isolationr3track[i] = (*Photons)[i].trkSumPtSolidConeDR03();
			photon_isolationr3trackhollow[i] = (*Photons)[i].trkSumPtHollowConeDR03();
			photon_isolationr3ecal[i] = (*Photons)[i].ecalRecHitSumEtConeDR03();
			photon_isolationr3hcal[i] = (*Photons)[i].hcalTowerSumEtConeDR03();
			photon_isolationr3ntrack[i] = (*Photons)[i].nTrkSolidConeDR03();
			photon_isolationr3ntrackhollow[i] = (*Photons)[i].nTrkHollowConeDR03();
			photon_isolationr4track[i] = (*Photons)[i].trkSumPtSolidConeDR04();
			photon_isolationr4trackhollow[i] = (*Photons)[i].trkSumPtHollowConeDR04();
			photon_isolationr4ecal[i] = (*Photons)[i].ecalRecHitSumEtConeDR04();
			photon_isolationr4hcal[i] = (*Photons)[i].hcalTowerSumEtConeDR04();
			photon_isolationr4ntrack[i] = (*Photons)[i].nTrkSolidConeDR04();
			photon_isolationr4ntrackhollow[i] = (*Photons)[i].nTrkHollowConeDR04();
			photon_superclusterindex[i] = getSuperCluster((*Photons)[i].superCluster()->energy(), (*Photons)[i].superCluster()->x(), (*Photons)[i].superCluster()->y(), (*Photons)[i].superCluster()->z()); 
			photon_info[i] = 0;
			photon_info[i] |= (*Photons)[i].isPhoton() << 0;
			photon_info[i] |= (*Photons)[i].hasConversionTracks() << 1;
			photon_info[i] |= (*Photons)[i].hasPixelSeed() << 2;
			photon_gapinfo[i] = 0;
			photon_gapinfo[i] |= (*Photons)[i].isEB() << 0;
			photon_gapinfo[i] |= (*Photons)[i].isEE() << 1;
			photon_gapinfo[i] |= (*Photons)[i].isEBGap() << 2;
			photon_gapinfo[i] |= (*Photons)[i].isEBEtaGap() << 3;
			photon_gapinfo[i] |= (*Photons)[i].isEBPhiGap() << 4;
			photon_gapinfo[i] |= (*Photons)[i].isEEGap() << 5;
			photon_gapinfo[i] |= (*Photons)[i].isEERingGap() << 6;
			photon_gapinfo[i] |= (*Photons)[i].isEEDeeGap() << 7;
			photon_gapinfo[i] |= (*Photons)[i].isEBEEGap() << 8;
			photon_conversionbegin[i] = conversion_count;
			photon_trigger[i] = GetTriggerMatch((*Photons)[i], photontriggers); 
		
			ConversionRefVector conversions = (*Photons)[i].conversions();
			for(unsigned j = 0 ; j < conversions.size() ; j++)
			{
				ConversionRef currconv = conversions[j];

				conversion_info[conversion_count] = 0;
				conversion_info[conversion_count] |= currconv->isConverted() << 0;
				conversion_info[conversion_count] |= currconv->conversionVertex().isValid() << 1;
				conversion_vx[conversion_count] = currconv->conversionVertex().x();
				conversion_vy[conversion_count] = currconv->conversionVertex().y();
				conversion_vz[conversion_count] = currconv->conversionVertex().z();
				conversion_chi2[conversion_count] = currconv->conversionVertex().chi2();
				conversion_ndof[conversion_count] = currconv->conversionVertex().ndof();
				conversion_cov[conversion_count][0] = currconv->conversionVertex().covariance(0,0);
				conversion_cov[conversion_count][1] = currconv->conversionVertex().covariance(0,1);
				conversion_cov[conversion_count][2] = currconv->conversionVertex().covariance(0,2);
				conversion_cov[conversion_count][3] = currconv->conversionVertex().covariance(1,1);
				conversion_cov[conversion_count][4] = currconv->conversionVertex().covariance(1,2);
				conversion_cov[conversion_count][5] = currconv->conversionVertex().covariance(2,2);
				conversion_mvaout[conversion_count] = currconv->MVAout();

				conversion_trackndof[conversion_count][0] = -1.;
				conversion_trackndof[conversion_count][1] = -1.;
				if(currconv->nTracks() == 2)
				{
					edm::RefToBase<Track> trA  = currconv->tracks()[0];
					TransientTrack SVTTrackA = TTrackBuilder->build(*trA);
					TrajectoryStateClosestToPoint TTrackStateA = SVTTrackA.trajectoryStateClosestToPoint(GlobalPoint(currconv->conversionVertex().x(), currconv->conversionVertex().y(), currconv->conversionVertex().z()));
					edm::RefToBase<Track> trB  = currconv->tracks()[1];
					TransientTrack SVTTrackB = TTrackBuilder->build(*trB);
					TrajectoryStateClosestToPoint TTrackStateB = SVTTrackB.trajectoryStateClosestToPoint(GlobalPoint(currconv->conversionVertex().x(), currconv->conversionVertex().y(), currconv->conversionVertex().z()));

					if(TTrackStateB.isValid() && TTrackStateA.isValid())
					{

						conversion_trackecalpointx[conversion_count][0] = currconv->ecalImpactPosition()[0].X();
						conversion_trackecalpointy[conversion_count][0] = currconv->ecalImpactPosition()[0].Y();
						conversion_trackecalpointz[conversion_count][0] = currconv->ecalImpactPosition()[0].Z();
						conversion_trackpx[conversion_count][0] = TTrackStateA.momentum().x();
						conversion_trackpy[conversion_count][0] = TTrackStateA.momentum().y();
						conversion_trackpz[conversion_count][0] = TTrackStateA.momentum().z();
						conversion_trackclosestpointx[conversion_count][0] =  TTrackStateA.position().x();
						conversion_trackclosestpointy[conversion_count][0] =  TTrackStateA.position().y();
						conversion_trackclosestpointz[conversion_count][0] =  TTrackStateA.position().z();
						conversion_trackchi2[conversion_count][0] = currconv->tracks()[0]->chi2();
						conversion_trackndof[conversion_count][0] = currconv->tracks()[0]->ndof();
						conversion_trackdxy[conversion_count][0] = TTrackStateA.perigeeParameters().transverseImpactParameter();
						conversion_trackdxyerr[conversion_count][0] = TTrackStateA.perigeeError().transverseImpactParameterError();
						conversion_trackdz[conversion_count][0] = TTrackStateA.perigeeParameters().longitudinalImpactParameter();
						conversion_trackdzerr[conversion_count][0] = TTrackStateA.perigeeError().longitudinalImpactParameterError();
						conversion_trackcharge[conversion_count][0] = currconv->tracks()[0]->charge();
						conversion_tracknhits[conversion_count][0] = currconv->tracks()[0]->numberOfValidHits();
						conversion_tracknmissinghits[conversion_count][0] = currconv->tracks()[0]->numberOfLostHits();
						conversion_tracknpixelhits[conversion_count][0] = currconv->tracks()[0]->hitPattern().numberOfValidPixelHits();
						conversion_tracknpixellayers[conversion_count][0] = currconv->tracks()[0]->hitPattern().pixelLayersWithMeasurement();
						conversion_tracknstriplayers[conversion_count][0] = currconv->tracks()[0]->hitPattern().stripLayersWithMeasurement();
						conversion_trackecalpointx[conversion_count][1] = currconv->ecalImpactPosition()[1].X();
						conversion_trackecalpointy[conversion_count][1] = currconv->ecalImpactPosition()[1].Y();
						conversion_trackecalpointz[conversion_count][1] = currconv->ecalImpactPosition()[1].Z();
						conversion_trackpx[conversion_count][1] = TTrackStateB.momentum().x();
						conversion_trackpy[conversion_count][1] = TTrackStateB.momentum().y();
						conversion_trackpz[conversion_count][1] = TTrackStateB.momentum().z();
						conversion_trackclosestpointx[conversion_count][1] =  TTrackStateB.position().x();
						conversion_trackclosestpointy[conversion_count][1] =  TTrackStateB.position().y();
						conversion_trackclosestpointz[conversion_count][1] =  TTrackStateB.position().z();
						conversion_trackchi2[conversion_count][1] = currconv->tracks()[1]->chi2();
						conversion_trackndof[conversion_count][1] = currconv->tracks()[1]->ndof();
						conversion_trackdxy[conversion_count][1] = TTrackStateB.perigeeParameters().transverseImpactParameter();
						conversion_trackdxyerr[conversion_count][1] = TTrackStateB.perigeeError().transverseImpactParameterError();
						conversion_trackdz[conversion_count][1] = TTrackStateB.perigeeParameters().longitudinalImpactParameter();
						conversion_trackdzerr[conversion_count][1] = TTrackStateB.perigeeError().longitudinalImpactParameterError();
						conversion_trackcharge[conversion_count][1] = currconv->tracks()[1]->charge();
						conversion_tracknhits[conversion_count][1] = currconv->tracks()[1]->numberOfValidHits();
						conversion_tracknmissinghits[conversion_count][1] = currconv->tracks()[1]->numberOfLostHits();
						conversion_tracknpixelhits[conversion_count][1] = currconv->tracks()[1]->hitPattern().numberOfValidPixelHits();
						conversion_tracknpixellayers[conversion_count][1] = currconv->tracks()[1]->hitPattern().pixelLayersWithMeasurement();
						conversion_tracknstriplayers[conversion_count][1] = currconv->tracks()[1]->hitPattern().stripLayersWithMeasurement();
					}
				}
				conversion_count++;
				if(conversion_count == M_conversionmaxcount){cerr << "number of conversions > M_conversionmaxcount. They are missing." << endl; errors |= 1<<4; break;}
			}
			photon_count++;
			if(photon_count == M_photonmaxcount || conversion_count == M_conversionmaxcount){cerr << "number of photon > M_photonmaxcount. They are missing." << endl; errors |= 1<<3; break;}
			if((*Photons)[i].pt() >= cPhotonPtMin && fabs((*Photons)[i].eta()) <= cPhotonEtaMax) NumGood++;
		}
	}

	if(NumGood >= cPhotonNum) return(true);
	return(false);
}

bool RootMaker::AddAllConversions(const edm::Event& iEvent){
  edm::Handle<ConversionCollection> Conversions;
  iEvent.getByLabel(edm::InputTag("allConversions"), Conversions);
  if(Conversions.isValid()){
    for(unsigned i = 0 ; i < Conversions->size() ; i++){
      allconversion_info[allconversion_count] = 0;
      allconversion_info[allconversion_count] |= (*Conversions)[i].isConverted() << 0;
      allconversion_info[allconversion_count] |= (*Conversions)[i].conversionVertex().isValid() << 1;
      allconversion_vx[allconversion_count] = (*Conversions)[i].conversionVertex().x();
      allconversion_vy[allconversion_count] = (*Conversions)[i].conversionVertex().y();
      allconversion_vz[allconversion_count] = (*Conversions)[i].conversionVertex().z();
      allconversion_chi2[allconversion_count] = (*Conversions)[i].conversionVertex().chi2();
      allconversion_ndof[allconversion_count] = (*Conversions)[i].conversionVertex().ndof();
      allconversion_cov[allconversion_count][0] = (*Conversions)[i].conversionVertex().covariance(0,0);
      allconversion_cov[allconversion_count][1] = (*Conversions)[i].conversionVertex().covariance(0,1);
      allconversion_cov[allconversion_count][2] = (*Conversions)[i].conversionVertex().covariance(0,2);
      allconversion_cov[allconversion_count][3] = (*Conversions)[i].conversionVertex().covariance(1,1);
      allconversion_cov[allconversion_count][4] = (*Conversions)[i].conversionVertex().covariance(1,2);
      allconversion_cov[allconversion_count][5] = (*Conversions)[i].conversionVertex().covariance(2,2);
      allconversion_mvaout[allconversion_count] = (*Conversions)[i].MVAout();
      
      allconversion_trackndof[allconversion_count][0] = -1.;
      allconversion_trackndof[allconversion_count][1] = -1.;
      if((*Conversions)[i].nTracks() == 2)
	{
	  edm::RefToBase<Track> trA  = (*Conversions)[i].tracks()[0];
	  TransientTrack SVTTrackA = TTrackBuilder->build(*trA);
	  TrajectoryStateClosestToPoint 
	    TTrackStateA = SVTTrackA.trajectoryStateClosestToPoint(GlobalPoint((*Conversions)[i].conversionVertex().x(), 
									       (*Conversions)[i].conversionVertex().y(), 
									       (*Conversions)[i].conversionVertex().z()));
	  edm::RefToBase<Track> trB  = (*Conversions)[i].tracks()[1];
	  TransientTrack SVTTrackB = TTrackBuilder->build(*trB);
	  TrajectoryStateClosestToPoint 
	    TTrackStateB = SVTTrackB.trajectoryStateClosestToPoint(GlobalPoint((*Conversions)[i].conversionVertex().x(), 
									       (*Conversions)[i].conversionVertex().y(), 
									       (*Conversions)[i].conversionVertex().z()));
	  if(TTrackStateB.isValid() && TTrackStateA.isValid()){
	    allconversion_trackecalpointx[allconversion_count][0] = (*Conversions)[i].ecalImpactPosition()[0].X();
	    allconversion_trackecalpointy[allconversion_count][0] = (*Conversions)[i].ecalImpactPosition()[0].Y();
	    allconversion_trackecalpointz[allconversion_count][0] = (*Conversions)[i].ecalImpactPosition()[0].Z();
	    allconversion_trackpx[allconversion_count][0] = TTrackStateA.momentum().x();
	    allconversion_trackpy[allconversion_count][0] = TTrackStateA.momentum().y();
	    allconversion_trackpz[allconversion_count][0] = TTrackStateA.momentum().z();
	    allconversion_trackclosestpointx[allconversion_count][0] =  TTrackStateA.position().x();
	    allconversion_trackclosestpointy[allconversion_count][0] =  TTrackStateA.position().y();
	    allconversion_trackclosestpointz[allconversion_count][0] =  TTrackStateA.position().z();
	    allconversion_trackchi2[allconversion_count][0] = (*Conversions)[i].tracks()[0]->chi2();
	    allconversion_trackndof[allconversion_count][0] = (*Conversions)[i].tracks()[0]->ndof();
	    allconversion_trackdxy[allconversion_count][0] = TTrackStateA.perigeeParameters().transverseImpactParameter();
	    allconversion_trackdxyerr[allconversion_count][0] = TTrackStateA.perigeeError().transverseImpactParameterError();
	    allconversion_trackdz[allconversion_count][0] = TTrackStateA.perigeeParameters().longitudinalImpactParameter();
	    allconversion_trackdzerr[allconversion_count][0] = TTrackStateA.perigeeError().longitudinalImpactParameterError();
	    allconversion_trackcharge[allconversion_count][0] = (*Conversions)[i].tracks()[0]->charge();
	    allconversion_tracknhits[allconversion_count][0] = (*Conversions)[i].tracks()[0]->numberOfValidHits();
	    allconversion_tracknmissinghits[allconversion_count][0] = (*Conversions)[i].tracks()[0]->numberOfLostHits();
	    allconversion_tracknpixelhits[allconversion_count][0] = (*Conversions)[i].tracks()[0]->hitPattern().numberOfValidPixelHits();
	    allconversion_tracknpixellayers[allconversion_count][0] = (*Conversions)[i].tracks()[0]->hitPattern().pixelLayersWithMeasurement();
	    allconversion_tracknstriplayers[allconversion_count][0] = (*Conversions)[i].tracks()[0]->hitPattern().stripLayersWithMeasurement();
	    allconversion_trackecalpointx[allconversion_count][1] = (*Conversions)[i].ecalImpactPosition()[1].X();
	    allconversion_trackecalpointy[allconversion_count][1] = (*Conversions)[i].ecalImpactPosition()[1].Y();
	    allconversion_trackecalpointz[allconversion_count][1] = (*Conversions)[i].ecalImpactPosition()[1].Z();
	    allconversion_trackpx[allconversion_count][1] = TTrackStateB.momentum().x();
	    allconversion_trackpy[allconversion_count][1] = TTrackStateB.momentum().y();
	    allconversion_trackpz[allconversion_count][1] = TTrackStateB.momentum().z();
	    allconversion_trackclosestpointx[allconversion_count][1] =  TTrackStateB.position().x();
	    allconversion_trackclosestpointy[allconversion_count][1] =  TTrackStateB.position().y();
	    allconversion_trackclosestpointz[allconversion_count][1] =  TTrackStateB.position().z();
	    allconversion_trackchi2[allconversion_count][1] = (*Conversions)[i].tracks()[1]->chi2();
	    allconversion_trackndof[allconversion_count][1] = (*Conversions)[i].tracks()[1]->ndof();
	    allconversion_trackdxy[allconversion_count][1] = TTrackStateB.perigeeParameters().transverseImpactParameter();
	    allconversion_trackdxyerr[allconversion_count][1] = TTrackStateB.perigeeError().transverseImpactParameterError();
	    allconversion_trackdz[allconversion_count][1] = TTrackStateB.perigeeParameters().longitudinalImpactParameter();
	    allconversion_trackdzerr[allconversion_count][1] = TTrackStateB.perigeeError().longitudinalImpactParameterError();
	    allconversion_trackcharge[allconversion_count][1] = (*Conversions)[i].tracks()[1]->charge();
	    allconversion_tracknhits[allconversion_count][1] = (*Conversions)[i].tracks()[1]->numberOfValidHits();
	    allconversion_tracknmissinghits[allconversion_count][1] = (*Conversions)[i].tracks()[1]->numberOfLostHits();
	    allconversion_tracknpixelhits[allconversion_count][1] = (*Conversions)[i].tracks()[1]->hitPattern().numberOfValidPixelHits();
	    allconversion_tracknpixellayers[allconversion_count][1] = (*Conversions)[i].tracks()[1]->hitPattern().pixelLayersWithMeasurement();
	    allconversion_tracknstriplayers[allconversion_count][1] = (*Conversions)[i].tracks()[1]->hitPattern().stripLayersWithMeasurement();
	  }
	}
      allconversion_count++;
      if(allconversion_count == M_conversionmaxcount){cerr << "number of conversions > M_conversionmaxcount. They are missing." << endl; errors |= 1<<4; break;}
    }
  }
  return(true);
}

double newEmFraction(const pat::Tau* a)
{
  double EcalE = 0;
  double E = 0;
  
  //for (reco::PFCandidateRefVector::const_iterator cand = a->signalPFCands().begin(); cand != a->signalPFCands().end(); ++cand) {
  for (std::vector<reco::PFCandidatePtr>::const_iterator cand = a->signalPFCands().begin(); cand != a->signalPFCands().end(); ++cand) {
    EcalE += (*cand)->ecalEnergy();
    E += (*cand)->ecalEnergy();
    E += (*cand)->hcalEnergy();
  }
  return EcalE/E;
}



int RootMaker::AddTaus(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(doDebug) cout<<"inside the AddTaus()"<< endl;
  edm::Handle<l1extra::L1JetParticleCollection> l1jetsHandle;
  const l1extra::L1JetParticleCollection* l1jets = 0;
  
  edm::Handle<l1extra::L1JetParticleCollection> l1tausHandle;
  const l1extra::L1JetParticleCollection* l1taus = 0;
  
  edm::Handle<BeamSpot> beamSpot;
  iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpot);
  
  // Obtain Collections 
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);

  iEvent.getByLabel(edm::InputTag("l1extraParticles","Central"), l1jetsHandle);
  if( !l1jetsHandle.isValid() )  edm::LogError("DataNotAvailable")  << "No L1CentralJets collection available \n";
  else  l1jets = l1jetsHandle.product();
  
  iEvent.getByLabel(edm::InputTag("l1extraParticles","Tau"), l1tausHandle);
  if( !l1jetsHandle.isValid() )  edm::LogError("DataNotAvailable")  << "No L1TauJets collection available \n";
  else  l1taus = l1tausHandle.product(); 
 
  // writing the piece to refit the vertex
  edm::Handle<reco::TrackCollection>  trackCollection;
  iEvent.getByLabel(TrackCollectionTag_,trackCollection);
  reco::TrackCollection nonTauTracks;
  std::vector<reco::TrackBaseRef> SignalTracks;
  
  edm::Handle<reco::VertexCollection > PV;
  iEvent.getByLabel(PVTag_,PV);
  const reco::VertexCollection* vertexes = PV.product();

  //tau collection
  edm::Handle<pat::TauCollection> Taus;
  iEvent.getByLabel(edm::InputTag("selectedPatTaus"), Taus);
  //iEvent.getByLabel(edm::InputTag("shrinkingConePFTauProducer"), Taus);

  unsigned int nGoodTaus = 0;
  if(Taus.isValid())
    {
      //vector<edm::Handle<PFTauDiscriminator> > PFTauDiscriminatiors(cTauDiscriminators.size());
      //for(unsigned n = 0 ; n < cTauDiscriminators.size() ; n++)
      //{
      //iEvent.getByLabel(cTauDiscriminators[n].c_str(), PFTauDiscriminatiors[n]);
      //}
      //edm::Handle<PFTauDiscriminator> hpsPFTauDiscriminationByLooseIsolation;
      //iEvent.getByLabel("hpsPFTauDiscriminationByLooseIsolation", hpsPFTauDiscriminationByLooseIsolation);
      tauIndexSelection.clear();

      for(unsigned i = 0 ; i < Taus->size() ; i++)
	{
	  // skiming  selection
	  if((cSkim & SKIM_ALL) == 0 && (*Taus)[i].pt() < 10) continue; 
	  if((*Taus)[i].eta() > 2.2 )                         continue;
	  if((*Taus)[i].tauID("decayModeFinding")     < 0.5 ) continue;

	  //pat::Tau aTauLeg( (*Taus)[i] );
	  //if( aTauLeg.tauID("againstMuonTight")       < 0.5 ) continue;
	  //if( aTauLeg.tauID("againstElectronTight")   < 0.5 ) continue;
	  //if( aTauLeg.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") < 0.5) continue;

	  if(doDebug) cout << "Skimmed events..."<< endl;
	  if((*Taus)[i].pt() > 20 && (*Taus)[i].tauID("decayModeFinding") > 0.5) ++nGoodTaus;
	  tauIndexSelection.push_back(i);
	  
	  // int numtrack = (*Taus)[i].signalPFChargedHadrCands().size();
	  // PFTauRef tauCandidate(Taus, i);
	  // if((*hpsPFTauDiscriminationByLooseIsolation)[tauCandidate] < 0.5 || !(numtrack == 3 || numtrack == 1)) continue;

	  tau_dishps[tau_count] = 0;
	  for(unsigned n = 0 ; n < cTauDiscriminators.size() ; n++)
	    {
	      if((*Taus)[i].tauID(cTauDiscriminators[n]) > 0.5) {
		tau_dishps[tau_count] |= (ULong64_t)1 << (ULong64_t)n;
	      }
	    }
	  
	  tau_e[tau_count]                                        = (*Taus)[i].energy();
	  tau_px[tau_count]                                       = (*Taus)[i].px();
	  tau_py[tau_count]                                       = (*Taus)[i].py();
	  tau_pz[tau_count]                                       = (*Taus)[i].pz();
	  tau_emfraction[tau_count]                               = (*Taus)[i].emFraction();
	  tau_newemfraction[tau_count]                            = newEmFraction(&(*Taus)[i]);
	  tau_hcaltotoverplead[tau_count]                         = (*Taus)[i].hcalTotOverPLead();
	  tau_hcal3x3overplead[tau_count]                         = (*Taus)[i].hcalMaxOverPLead();
	  tau_ecalstripsumeoverplead[tau_count]                   = (*Taus)[i].hcal3x3OverPLead();
	  tau_bremsrecoveryeoverplead[tau_count]                  = (*Taus)[i].ecalStripSumEOverPLead();
	  tau_calocomp[tau_count]                                 = (*Taus)[i].caloComp();
	  tau_segcomp[tau_count]                                  = (*Taus)[i].segComp();

	  tau_againstelectronmva5raw[tau_count]                   = (*Taus)[i].tauID("againstElectronMVA5raw");      //OK
	  tau_byIsolationmva3newDMwoLTraw[tau_count]              = (*Taus)[i].tauID("byIsolationMVA3newDMwoLTraw"); //OK
	  tau_byIsolationmva3newDMwLTraw[tau_count]               = (*Taus)[i].tauID("byIsolationMVA3newDMwLTraw");  //OK
	  tau_againstelectronmva5raw[tau_count]                   = (*Taus)[i].tauID("againstElectronMVA5raw");      //OK
	  tau_againstelectronVLoosemva5[tau_count]                = (*Taus)[i].tauID("againstElectronVLooseMVA5");
	  tau_againstelectronLoosemva5[tau_count]                 = (*Taus)[i].tauID("againstElectronLooseMVA5");
	  tau_againstelectronMediummva5[tau_count]                = (*Taus)[i].tauID("againstElectronMediumMVA5");
	  tau_againstelectronTightmva5[tau_count]                 = (*Taus)[i].tauID("againstElectronTightMVA5");
	  tau_againstelectronDeadECAL[tau_count]                  = (*Taus)[i].tauID("againstElectronDeadECAL");
	  tau_againstelectronmva5category[tau_count]              = (*Taus)[i].tauID("againstElectronMVA5category"); //OK
	  tau_bycombinedisolationdeltabetacorrraw3hits[tau_count] = (*Taus)[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"); //OK
	  tau_signalPFChargedHadrCands_size[tau_count]            = (*Taus)[i].signalPFChargedHadrCands().size(); 	  
	  tau_signalPFGammaCands_size[tau_count]	          = (*Taus)[i].signalPFGammaCands().size();
	  if( ((*Taus)[i].genJet())) 
	    tau_genTaudecayMode[tau_count]                        = JetMCTagUtils::genTauDecayMode(*((*Taus)[i].genJet())); //added on Jul 27
	  else tau_genTaudecayMode[tau_count]                     = -1.;
	  
	  SignedImpactParameter3D signed_ip3D;
	  // std::vector<reco::PFCandidate> leadPFChargedHadrCand_;
	  // const reco::PFCandidateRef leadPFChargedHadrCand() const;

	  if((*Taus)[i].leadPFChargedHadrCand().isNonnull())
	    {
	      tau_leadpfchargedhadrcandecalenergy[tau_count] = (*Taus)[i].leadPFChargedHadrCand()->ecalEnergy();
	      tau_leadpfchargedhadrcandhcalenergy[tau_count] = (*Taus)[i].leadPFChargedHadrCand()->hcalEnergy();
	      tau_leadpfchargedhadrcandp[tau_count]          = (*Taus)[i].leadPFChargedHadrCand()->p();
	      
	      TrackRef track = (*Taus)[i].leadPFChargedHadrCand()->trackRef();
	      if(track.isNonnull())
		{
		  const reco::TransientTrack transientTrack = TTrackBuilder->build(*track);
		  const Measurement1D meas = signed_ip3D.apply(transientTrack, GlobalVector(track->px(), track->py(), track->pz()), primvertex).second;
		  
		  tau_dxy[tau_count] = track->dxy(pv_position);
		  tau_dz[tau_count] = track->dz(pv_position);
		  tau_ip3d[tau_count] = meas.value();
		  tau_ip3dSig[tau_count] = meas.significance();
		}
	      else
		{
		  tau_dxy[tau_count] = -100.0f;
		  tau_dz[tau_count] = -100.0f;
		  tau_ip3d[tau_count] = -1.0f;
		  tau_ip3dSig[tau_count] = -1.0f;
		}
	    }
	  else
	    {
	      tau_leadpfchargedhadrcandecalenergy[tau_count] = -1.0f;
	      tau_leadpfchargedhadrcandhcalenergy[tau_count] = -1.0f;
	      tau_leadpfchargedhadrcandp[tau_count] = -1.0f;
	      tau_dxy[tau_count] = -100.0f;
	      tau_dz[tau_count] = -100.0f;
	      tau_ip3d[tau_count] = -1.0f;
	      tau_ip3dSig[tau_count] = -1.0f;
	    }
	  
	  tau_nprongs[tau_count]      = (*Taus)[i].signalPFCands().size();
	  tau_charge[tau_count]       = (*Taus)[i].charge();
	  tau_chargedbegin[tau_count] = tau_charged_count;
	  
	  // matching with the jet collection can be done offline based on these information:
	  PFJetRef thejet = (*Taus)[i].pfJetRef();
	  tau_ak5pfjet_e[tau_count] = thejet->energy();
	  tau_ak5pfjet_px[tau_count] = thejet->px();
	  tau_ak5pfjet_py[tau_count] = thejet->py();
	  tau_ak5pfjet_pz[tau_count] = thejet->pz();
	  
	  tau_isolationchargednum[tau_count] = 0;
	  tau_isolationchargedpt[tau_count] = 0.;
	  tau_isolationneutralsnum[tau_count] = 0;
	  tau_isolationneutralspt[tau_count] = 0.;
	  tau_isolationgammanum[tau_count] = (*Taus)[i].isolationPFGammaCands().size();
	  tau_isolationgammapt[tau_count] = (*Taus)[i].isolationPFGammaCandsEtSum();
	  
	  for(unsigned n = 0 ; n < (*Taus)[i].isolationPFCands().size() ; n++)
	    {
	      //PFCandidateRef isocand = (*Taus)[i].isolationPFCands()[n];
	      PFCandidatePtr isocand = (*Taus)[i].isolationPFCands()[n];
	      if(isocand->charge() == 0)
		{
		  tau_isolationneutralsnum[tau_count]++;
		  tau_isolationneutralspt[tau_count]+=isocand->pt();
		}
	      else
		{
		  tau_isolationchargednum[tau_count]++;
		  tau_isolationchargedpt[tau_count]+=isocand->pt();
		}
	    }
	  tau_trigger[tau_count] = GetTriggerMatch((*Taus)[i], tautriggers);
	  tau_L1trigger_match[tau_count] = GetL1ExtraTriggerMatch(l1jets, l1taus,  (*Taus)[i] );
	  
	  //tau track
	  TrackRef track;
	  for(unsigned n = 0 ; n < (*Taus)[i].signalPFCands().size() ; n++)
	    {
	      //PFCandidateRef cand = (*Taus)[i].signalPFCands()[n];
	      PFCandidatePtr cand = (*Taus)[i].signalPFCands()[n];
	      track = cand->trackRef();
	      if(!track.isNull())
		{
		  edm::Handle<edm::ValueMap<DeDxData> > dEdxharmonic2;
		  iEvent.getByLabel(edm::InputTag("dedxHarmonic2"), dEdxharmonic2);
		  TransientTrack TTrack = TTrackBuilder->build(track);
		  TrajectoryStateClosestToPoint TTrackState = TTrack.trajectoryStateClosestToPoint(GlobalPoint(pv_position.x(), pv_position.y(), pv_position.z()));
		  tau_charged_px[tau_charged_count] = TTrackState.momentum().x();
		  tau_charged_py[tau_charged_count] = TTrackState.momentum().y();
		  tau_charged_pz[tau_charged_count] = TTrackState.momentum().z();
		  tau_charged_closestpointx[tau_charged_count] = TTrackState.position().x();
		  tau_charged_closestpointy[tau_charged_count] = TTrackState.position().y();
		  tau_charged_closestpointz[tau_charged_count] = TTrackState.position().z();
		  tau_charged_dxy[tau_charged_count]    = TTrackState.perigeeParameters().transverseImpactParameter();
		  tau_charged_dxyerr[tau_charged_count] = TTrackState.perigeeError().transverseImpactParameterError();
		  tau_charged_dz[tau_charged_count]     = TTrackState.perigeeParameters().longitudinalImpactParameter();
		  tau_charged_dzerr[tau_charged_count]  = TTrackState.perigeeError().longitudinalImpactParameterError();
		  tau_charged_chi2[tau_charged_count]   = track->chi2();
		  tau_charged_ndof[tau_charged_count]   = track->ndof();
		  tau_charged_charge[tau_charged_count] = track->charge();
		  tau_charged_nhits[tau_charged_count] = track->numberOfValidHits();
		  tau_charged_nmissinghits[tau_charged_count] = track->numberOfLostHits();
		  tau_charged_npixelhits[tau_charged_count] = track->hitPattern().numberOfValidPixelHits();
		  tau_charged_npixellayers[tau_charged_count] = track->hitPattern().pixelLayersWithMeasurement();
		  tau_charged_nstriplayers[tau_charged_count] = track->hitPattern().stripLayersWithMeasurement();
		  tau_charged_dedxharmonic2[tau_charged_count] = (*dEdxharmonic2)[track].dEdx();
		  math::XYZPoint ecalPos = PositionOnECalSurface(TTrack);
		  tau_charged_outerx[tau_charged_count] = ecalPos.x();
		  tau_charged_outery[tau_charged_count] = ecalPos.y();
		  tau_charged_outerz[tau_charged_count] = ecalPos.z();
		  tau_charged_count++;
		  if(tau_charged_count == M_taumaxcount*10) break;
		} //if(!track.isNull())
	    } // for(unsigned n = 0 ; n < (*Taus)[i].signalPFCands().size() ; n++)
	  tau_count++;

	  if(tau_count == M_taumaxcount || tau_charged_count == M_taumaxcount*10) {
	    cerr << "number of taus > M_taumaxcount. They are missing." << endl; errors |= 1<<10; 
	    break;
	  }
	} // for(unsigned i = 0 ; i < Taus->size() ; i++)
 
      if(doDebug) cout << "Starting the diTau Study"<< endl;     
      // diTau Vertex Study
      reco::Vertex thePV;
      bool tau_vtxFound = false;

      //to make sure if vertex is found, else break
      if(vertexes->size() > 0) 
	{
	  for(size_t ivtx = 0; ivtx < vertexes->size(); ivtx++){
	    reco::Vertex theVtx = vertexes->at(ivtx);
	    if(!theVtx.isValid())                      continue;
	    if(theVtx.ndof() < 4)                      continue;
	    if(theVtx.z()    < -24 || theVtx.z() > 24) continue;
	    if(theVtx.position().Rho() > 2.)           continue;
	    thePV=vertexes->at(ivtx);
	    tau_vtxFound = true;
	    break;
	  }
	}
      
      if(tau_vtxFound) { 
	ditau_VtxX    = thePV.x();
	ditau_VtxY    = thePV.y();
	ditau_VtxZ    = thePV.z();
	ditau_VtxXErr = thePV.xError();
	ditau_VtxYErr = thePV.yError();
	ditau_VtxZErr = thePV.zError();
	
	// loop over ditau pair
	for (unsigned indxTau1 =0; indxTau1 < tauIndexSelection.size() ; indxTau1++) 
	  {
	    pat::Tau aTauLeg1((*Taus)[indxTau1]);
	    
	    for (unsigned indxTau2 = indxTau1+1; indxTau2 < tauIndexSelection.size(); indxTau2++) 
	      {
		if(indxTau2 <= indxTau1)  continue;
		ditau_leg1_index[ditau_Index] = -99;	
		ditau_leg2_index[ditau_Index] = -99;	
		
		pat::Tau aTauLeg2( (*Taus)[indxTau2]);
		ditau_leg1_index[ditau_Index] = indxTau1;
		ditau_leg2_index[ditau_Index] = indxTau2;
		
		ditau_reFitVtxX[ditau_Index]       = -99;  
		ditau_reFitVtxY[ditau_Index]       = -99;
		ditau_reFitVtxZ[ditau_Index]       = -99;
		ditau_reFitVtxXErr[ditau_Index]    = -99;
		ditau_reFitVtxYErr[ditau_Index]    = -99;
		ditau_reFitVtxZErr[ditau_Index]    = -99;
		ditau_reFitVtxRho[ditau_Index]     = -99;
		ditau_reFitVtxNdof[ditau_Index]    = -99;
		ditau_leg1_dxy[ditau_Index]        = -99;
		ditau_leg1_dz[ditau_Index]         = -99;
		ditau_leg1_dxyErr[ditau_Index]     = -99;
		ditau_leg1_dzErr[ditau_Index]      = -99;
		ditau_leg1_dxy_OPV[ditau_Index]    = -99;
		ditau_leg1_dz_OPV[ditau_Index]     = -99;
		ditau_leg1_dxyErr_OPV[ditau_Index] = -99;
		ditau_leg1_dzErr_OPV[ditau_Index]  = -99;
		ditau_leg1_X[ditau_Index]          = -99;
		ditau_leg1_Y[ditau_Index]          = -99;
		ditau_leg1_Z[ditau_Index]          = -99;
		ditau_leg1_X_OPV[ditau_Index]      = -99;
		ditau_leg1_Y_OPV[ditau_Index]      = -99;
		ditau_leg1_Z_OPV[ditau_Index]      = -99;
		ditau_leg2_dxy[ditau_Index]        = -99;
		ditau_leg2_dz[ditau_Index]         = -99;
		ditau_leg2_dxyErr[ditau_Index]     = -99;
		ditau_leg2_dzErr[ditau_Index]      = -99;
		ditau_leg2_dxy_OPV[ditau_Index]    = -99;
		ditau_leg2_dz_OPV[ditau_Index]     = -99;
		ditau_leg2_dxyErr_OPV[ditau_Index] = -99;
		ditau_leg2_dzErr_OPV[ditau_Index]  = -99;
		ditau_leg2_X_OPV[ditau_Index]      = -99;
		ditau_leg2_Y_OPV[ditau_Index]      = -99;
		ditau_leg2_Z_OPV[ditau_Index]      = -99;
		
		// Get signal tracks from DiTaus
		//const reco::PFCandidateRefVector & candsLeg1 = aTauLeg1.signalPFChargedHadrCands();
		const std::vector<reco::PFCandidatePtr> & candsLeg1 = aTauLeg1.signalPFChargedHadrCands();
		//for (reco::PFCandidateRefVector::const_iterator iter = candsLeg1.begin(); iter!=candsLeg1.end(); iter++)
		for (std::vector<reco::PFCandidatePtr>::const_iterator iter = candsLeg1.begin(); iter!=candsLeg1.end(); iter++)
		  {
		    if(iter->get()->trackRef().isNonnull()) 
		      SignalTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
		    else if(iter->get()->gsfTrackRef().isNonnull())
		      {
			SignalTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));
		      }
		  }
		
		// Get tracks form PFTau daugthers of tauLeg2 
		//const reco::PFCandidateRefVector & candsLeg2 = aTauLeg2.signalPFChargedHadrCands();
		const std::vector<reco::PFCandidatePtr> & candsLeg2 = aTauLeg2.signalPFChargedHadrCands();
		//for (reco::PFCandidateRefVector::const_iterator iter = candsLeg2.begin(); iter!=candsLeg2.end(); iter++)
		for (std::vector<reco::PFCandidatePtr>::const_iterator iter = candsLeg2.begin(); iter!=candsLeg2.end(); iter++)
		  {
		    if(iter->get()->trackRef().isNonnull()) 
		      SignalTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
		    else if(iter->get()->gsfTrackRef().isNonnull())
		      {
			SignalTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));
		      }
		  }
		
		// Get Non-Tau tracks 
		reco::TrackCollection nonTauTracks;
		if (trackCollection.isValid()) {
		  
		  // remove tau tracks and only tracks associated with the vertex                                                                    
		  unsigned int idx = 0;
		  for (reco::TrackCollection::const_iterator iTrk = trackCollection->begin(); iTrk != trackCollection->end(); ++iTrk, idx++) {
		    reco::TrackRef tmpRef(trackCollection, idx);
		    reco::TrackRef tmpRefForBase=tmpRef;
		    bool isSigTrk = false;
		    bool fromVertex=false;
		    
		    for (unsigned int sigTrk = 0; sigTrk < SignalTracks.size(); sigTrk++) {
		      if (reco::TrackBaseRef(tmpRefForBase) == SignalTracks.at(sigTrk)) 
			{
			  isSigTrk = true;
			  break;
			}
		    }
		    
		    for(std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef=thePV.tracks_begin(); vtxTrkRef < thePV.tracks_end(); vtxTrkRef++ )
		      {
			if ( thePV.trackWeight(*vtxTrkRef) > 0 ) {
			  if( (*vtxTrkRef) == reco::TrackBaseRef(tmpRefForBase) )
			    {
			      fromVertex = true;
			      break;
			    }
			}
		      }
		    if (!isSigTrk && fromVertex) {
		      nonTauTracks.push_back(*iTrk);
		    }
		  }
		} // if (trackCollection.isValid())
		
		
		if(nonTauTracks.size() >= 2) {
		  
		  // Refit the vertex
		  TransientVertex transVtx;
		  std::vector<reco::TransientTrack> transTracks;
		  for (reco::TrackCollection::iterator iter=nonTauTracks.begin(); iter!=nonTauTracks.end(); ++iter){
		    transTracks.push_back(TTrackBuilder->build(*iter));
		  }
		  
		  bool FitOk(true);
		  AdaptiveVertexFitter avf;	    
		  avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception                                                       
		  try   {
		    if(!cbeamspot) transVtx = avf.vertex(transTracks);
		    else transVtx = avf.vertex(transTracks,*beamSpot);
		  } 
		  catch(...)  {
		    FitOk = false;
		  }
		  
		  reco::Vertex primaryVertexReFit;
		  
		  if(FitOk) {
		    primaryVertexReFit = transVtx;
		    
		    ditau_reFitVtxX[ditau_Index] = primaryVertexReFit.x();
		    ditau_reFitVtxY[ditau_Index] = primaryVertexReFit.y();
		    ditau_reFitVtxZ[ditau_Index] = primaryVertexReFit.z();
		    ditau_reFitVtxXErr[ditau_Index] = primaryVertexReFit.xError();
		    ditau_reFitVtxYErr[ditau_Index] = primaryVertexReFit.yError();
		    ditau_reFitVtxZErr[ditau_Index] = primaryVertexReFit.zError();
		    ditau_reFitVtxRho[ditau_Index]  = primaryVertexReFit.position().Rho();
		    ditau_reFitVtxNdof[ditau_Index] = primaryVertexReFit.ndof();
		    
		    if(aTauLeg1.leadPFChargedHadrCand().isNonnull()){
		      if(aTauLeg1.leadPFChargedHadrCand()->trackRef().isNonnull()){
			
			reco::TransientTrack transTrk=transTrackBuilder->build(aTauLeg1.leadPFChargedHadrCand()->trackRef());
			
			GlobalPoint pv(primaryVertexReFit.position().x(), primaryVertexReFit.position().y(), primaryVertexReFit.position().z());
			ditau_leg1_dxy[ditau_Index]       = transTrk.trajectoryStateClosestToPoint(pv).perigeeParameters().transverseImpactParameter();
			ditau_leg1_dz[ditau_Index]        = transTrk.trajectoryStateClosestToPoint(pv).perigeeParameters().longitudinalImpactParameter();
			ditau_leg1_dxyErr[ditau_Index]    = transTrk.trajectoryStateClosestToPoint(pv).perigeeError().transverseImpactParameterError();
			ditau_leg1_dzErr[ditau_Index]     = transTrk.trajectoryStateClosestToPoint(pv).perigeeError().longitudinalImpactParameterError();

			GlobalPoint pos                   = transTrk.trajectoryStateClosestToPoint(pv).position();
			ditau_leg1_X[ditau_Index]         = pos.x();
			ditau_leg1_Y[ditau_Index]         = pos.y();
			ditau_leg1_Z[ditau_Index]         = pos.z();

			//wrt original vertex 
			GlobalPoint pvOrig(thePV.position().x(),thePV.position().y(),thePV.position().z());
			ditau_leg1_dxy_OPV[ditau_Index]       = transTrk.trajectoryStateClosestToPoint(pvOrig).perigeeParameters().transverseImpactParameter();
			ditau_leg1_dz_OPV[ditau_Index]        = transTrk.trajectoryStateClosestToPoint(pvOrig).perigeeParameters().longitudinalImpactParameter();
			ditau_leg1_dxyErr_OPV[ditau_Index]    = transTrk.trajectoryStateClosestToPoint(pvOrig).perigeeError().transverseImpactParameterError();
			ditau_leg1_dzErr_OPV[ditau_Index]     = transTrk.trajectoryStateClosestToPoint(pvOrig).perigeeError().longitudinalImpactParameterError();

			GlobalPoint posOPV                    = transTrk.trajectoryStateClosestToPoint(pvOrig).position();
			ditau_leg1_X[ditau_Index]             = posOPV.x();
			ditau_leg1_Y[ditau_Index]             = posOPV.y();
			ditau_leg1_Z[ditau_Index]             = posOPV.z();
		      }
		      
		      if(aTauLeg2.leadPFChargedHadrCand().isNonnull())
			{
			  if(aTauLeg2.leadPFChargedHadrCand()->trackRef().isNonnull())
			    {
			      reco::TransientTrack transTrk=transTrackBuilder->build(aTauLeg2.leadPFChargedHadrCand()->trackRef());
			      GlobalPoint pv(primaryVertexReFit.position().x(),primaryVertexReFit.position().y(),primaryVertexReFit.position().z());
			      ditau_leg2_dxy[ditau_Index]        = transTrk.trajectoryStateClosestToPoint(pv).perigeeParameters().transverseImpactParameter();
			      ditau_leg2_dz[ditau_Index]         = transTrk.trajectoryStateClosestToPoint(pv).perigeeParameters().longitudinalImpactParameter();
			      ditau_leg2_dxyErr[ditau_Index]     = transTrk.trajectoryStateClosestToPoint(pv).perigeeError().transverseImpactParameterError();
			      ditau_leg2_dzErr[ditau_Index]      = transTrk.trajectoryStateClosestToPoint(pv).perigeeError().longitudinalImpactParameterError();
			      GlobalPoint pos                  = transTrk.trajectoryStateClosestToPoint(pv).position();
			      ditau_leg2_X[ditau_Index]        = pos.x();
			      ditau_leg2_Y[ditau_Index]        = pos.y();
			      ditau_leg2_Z[ditau_Index]        = pos.z();
			      
			      //wrt original vertex
			      GlobalPoint pvOrig(thePV.position().x(),thePV.position().y(),thePV.position().z());
			      ditau_leg2_dxy_OPV[ditau_Index]         = transTrk.trajectoryStateClosestToPoint(pvOrig).perigeeParameters().transverseImpactParameter();
			      ditau_leg2_dz_OPV[ditau_Index]          = transTrk.trajectoryStateClosestToPoint(pvOrig).perigeeParameters().longitudinalImpactParameter();
			      ditau_leg2_dxyErr_OPV[ditau_Index]      = transTrk.trajectoryStateClosestToPoint(pvOrig).perigeeError().transverseImpactParameterError();
			      ditau_leg2_dzErr_OPV[ditau_Index]       = transTrk.trajectoryStateClosestToPoint(pvOrig).perigeeError().longitudinalImpactParameterError();
			      GlobalPoint posOPV                      = transTrk.trajectoryStateClosestToPoint(pvOrig).position();
			      ditau_leg2_X_OPV[ditau_Index]           = posOPV.x();
			      ditau_leg2_Y_OPV[ditau_Index]           = posOPV.y();
			      ditau_leg2_Z_OPV[ditau_Index]           = posOPV.z();
			    } // if(aTauLeg2.leadPFChargedHadrCand()->trackRef().isNonnull())
			} //  if(aTauLeg2.leadPFChargedHadrCand().isNonnull())
		    } //  if(aTauLeg1.leadPFChargedHadrCand().isNonnull())
		  } // if(FitOk)
		} // if(nonTauTracks.size() >= 2) 
		
		ditau_Index++; 
	      } // for (unsigned indxTau2 = indxTau1+1; indxTau2 < tauIndexSelection.size(); indxTau2++)
	  } // for (unsigned indxTau1 =0; indxTau1 < tauIndexSelection.size() ; indxTau1++) 
      } // if(tau_vtxFound)
    } // if(Taus.isValid())

  if(doDebug) cout<<"Leaving the AddTaus()"<< endl;
  return nGoodTaus;
}
 

#if 0
template<typename TCollection>
const NSVfitEventHypothesisBase* matchSVfitHypothesis(const NSVfitEventHypothesisBase& hypo, const TCollection& coll) //const NSVfitEventHypothesisBaseCollection& coll)
{
	for(unsigned int i = 0; i < coll.size(); ++i)
	{
		const NSVfitEventHypothesisBase& cand = coll[i];
		assert(hypo.numResonances() == cand.numResonances());

		bool particle_match = true;
		for(unsigned int j = 0; j < hypo.numResonances(); ++j)
		{
			const NSVfitResonanceHypothesisBase* hypo_res = hypo.resonance(j);
			const NSVfitResonanceHypothesisBase* cand_res = cand.resonance(j);
			assert(hypo_res && cand_res);

			assert(hypo_res->numDaughters() == cand_res->numDaughters());
			for(unsigned int k = 0; k < hypo_res->numDaughters(); ++k)
			{
				if(hypo_res->daughter(k)->particle() != cand_res->daughter(k)->particle())
					particle_match = false;
			}
		}

		if(particle_match) return &cand;
	}

	assert(false);
	return NULL;
}
#endif

RootMaker::DCA RootMaker::calculateDCA(const pat::Tau& tau1, const pat::Tau& tau2)
{
        // TODO: Use the reconstructed decay mode, and then only the mass of the decay products?
	float tauMass = 1.777;
	// TODO: What is this sigma?
	float tauSigma = tauMass*1e-6;

	DCA dca = { -1.0, -1.0f, -1.0f, -1.0f };
	reco::TrackRef track1 = tau1.leadPFChargedHadrCand()->trackRef();
	reco::TrackRef track2 = tau2.leadPFChargedHadrCand()->trackRef();
	if(track1.isNonnull() && track2.isNonnull())
	{
		reco::TransientTrack transientTrack1 = TTrackBuilder->build(*track1);
		reco::TransientTrack transientTrack2 = TTrackBuilder->build(*track2);
		if(transientTrack1.impactPointTSCP().isValid() && transientTrack2.impactPointTSCP().isValid())
		{
			FreeTrajectoryState state1 = transientTrack1.impactPointTSCP().theState();
			FreeTrajectoryState state2 = transientTrack2.impactPointTSCP().theState();
			TwoTrackMinimumDistance minDist;
			minDist.calculate(state1, state2);
			if(minDist.status())
			{
				const float dist3D = minDist.distance();
				std::pair<GlobalPoint,GlobalPoint> pcas = minDist.points();
				GlobalPoint pca1 = pcas.first;
				GlobalPoint pca2 = pcas.second;

				ROOT::Math::SVector<double, 3> distanceVector(pca1.x()-pca2.x(), pca1.y()-pca2.y(), pca1.z()-pca2.z());
				const float twoTauDist3D = ROOT::Math::Mag(distanceVector);
				assert(fabs(dist3D - twoTauDist3D) < 1e-3);

				float chi2 = 0.0f, ndf = 0.0f;
				KinematicParticleFactoryFromTransientTrack pFactory;
				RefCountedKinematicParticle tau1Particle = pFactory.particle(transientTrack1, tauMass, chi2, ndf, tauSigma);
				RefCountedKinematicParticle tau2Particle = pFactory.particle(transientTrack2, tauMass, chi2, ndf, tauSigma);

				float sig[6];
				sig[0] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,0);
				sig[1] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,1);
				sig[2] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,1);
				sig[3] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,2);
				sig[4] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,2);
				sig[5] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(2,2);
				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca1Cov(sig, sig+6);

				sig[0] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,0);
				sig[1] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,1);
				sig[2] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,1);
				sig[3] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,2);
				sig[4] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,2);
				sig[5] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(2,2);
				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca2Cov(sig, sig+6);

				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > totCov = pca1Cov + pca2Cov;
				const float twoTauDist3DErr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTauDist3D;

				distanceVector(2) = 0.0;
				const float twoTauDist2D = ROOT::Math::Mag(distanceVector);
				const float twoTauDist2DErr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTauDist2D;

				dca.dca2d = twoTauDist2D;
				dca.dca2dErr = twoTauDist2DErr;
				dca.dca3d = twoTauDist3D;
				dca.dca3dErr = twoTauDist3DErr;
			}
		}
	}

	return dca;
}


bool RootMaker::AddMuTauTauPairs(const edm::Event& iEvent)
{
	edm::Handle<NSVfitEventHypothesisByIntegrationCollection> TauPairsInt;
	iEvent.getByLabel(edm::InputTag("MuTauTauPairsWHInt"), TauPairsInt);

	mutautaupair_count = TauPairsInt->size();
	for(unsigned i = 0 ; i < TauPairsInt->size() ; i++)
	{
		const NSVfitEventHypothesisBase& TauPair = (*TauPairsInt)[i];

		const NSVfitResonanceHypothesisBase* svFitResonanceHypothesis = TauPair.resonance("Higgs");
		assert(svFitResonanceHypothesis);

		assert(svFitResonanceHypothesis->numDaughters() == 2);
		const NSVfitSingleParticleHypothesisBase* svFitDaughter1 = svFitResonanceHypothesis->daughter("leg1");
		const NSVfitSingleParticleHypothesisBase* svFitDaughter2 = svFitResonanceHypothesis->daughter("leg2");
		assert(svFitDaughter1 && svFitDaughter2);

		const NSVfitResonanceHypothesisBase* svFitResonanceHypothesisW = TauPair.NSVfitEventHypothesisBase::resonance("W");
		assert(svFitResonanceHypothesisW);

		assert(svFitResonanceHypothesisW->numDaughters() >= 1);
		const NSVfitSingleParticleHypothesisBase* svFitDaughterMu = svFitResonanceHypothesisW->daughter("chargedLepton");
		assert(svFitDaughterMu);

		const pat::Tau* tau1 = dynamic_cast<const pat::Tau*>(&*svFitDaughter1->particle());
		const pat::Tau* tau2 = dynamic_cast<const pat::Tau*>(&*svFitDaughter2->particle());
		if(!tau1 || !tau2) throw cms::Exception("RootMaker::AddMuTauTauPairs") << "Candidates are not PAT taus";
		DCA dca = calculateDCA(*tau1, *tau2);

		mutautaupair_leg1_px[i] = svFitDaughter1->particle()->px();
		mutautaupair_leg1_py[i] = svFitDaughter1->particle()->py();
		mutautaupair_leg1_pz[i] = svFitDaughter1->particle()->pz();
		mutautaupair_leg1_energy[i] = svFitDaughter1->particle()->energy();
		mutautaupair_leg2_px[i] = svFitDaughter2->particle()->px();
		mutautaupair_leg2_py[i] = svFitDaughter2->particle()->py();
		mutautaupair_leg2_pz[i] = svFitDaughter2->particle()->pz();
		mutautaupair_leg2_energy[i] = svFitDaughter2->particle()->energy();
		mutautaupair_mu_px[i] = svFitDaughterMu->particle()->px();
		mutautaupair_mu_py[i] = svFitDaughterMu->particle()->py();
		mutautaupair_mu_pz[i] = svFitDaughterMu->particle()->pz();
		mutautaupair_mu_energy[i] = svFitDaughterMu->particle()->energy();
		mutautaupair_svfit_int_valid[i] = svFitResonanceHypothesis->isValidSolution();
		mutautaupair_svfit_mass_int[i] = svFitResonanceHypothesis->mass();
		mutautaupair_svfit_mass_int_err_up[i] = svFitResonanceHypothesis->massErrUp();
		mutautaupair_svfit_mass_int_err_down[i] = svFitResonanceHypothesis->massErrDown();

		mutautaupair_dca2d[i] = dca.dca2d;
		mutautaupair_dca2d_err[i] = dca.dca2dErr;
		mutautaupair_dca3d[i] = dca.dca3d;
		mutautaupair_dca3d_err[i] = dca.dca3dErr;
	}

	return true;
}

bool RootMaker::AddElTauTauPairs(const edm::Event& iEvent)
{
	edm::Handle<NSVfitEventHypothesisByIntegrationCollection> TauPairsInt;
	iEvent.getByLabel(edm::InputTag("ElTauTauPairsWHInt"), TauPairsInt);

	eltautaupair_count = TauPairsInt->size();
	for(unsigned i = 0 ; i < TauPairsInt->size() ; i++)
	{
		const NSVfitEventHypothesisBase& TauPair = (*TauPairsInt)[i];

		const NSVfitResonanceHypothesisBase* svFitResonanceHypothesis = TauPair.resonance("Higgs");
		assert(svFitResonanceHypothesis);

		assert(svFitResonanceHypothesis->numDaughters() == 2);
		const NSVfitSingleParticleHypothesisBase* svFitDaughter1 = svFitResonanceHypothesis->daughter("leg1");
		const NSVfitSingleParticleHypothesisBase* svFitDaughter2 = svFitResonanceHypothesis->daughter("leg2");
		assert(svFitDaughter1 && svFitDaughter2);

		const NSVfitResonanceHypothesisBase* svFitResonanceHypothesisW = TauPair.NSVfitEventHypothesisBase::resonance("W");
		assert(svFitResonanceHypothesisW);

		assert(svFitResonanceHypothesisW->numDaughters() >= 1);
		const NSVfitSingleParticleHypothesisBase* svFitDaughterEl = svFitResonanceHypothesisW->daughter("chargedLepton");
		assert(svFitDaughterEl);

		const pat::Tau* tau1 = dynamic_cast<const pat::Tau*>(&*svFitDaughter1->particle());
		const pat::Tau* tau2 = dynamic_cast<const pat::Tau*>(&*svFitDaughter2->particle());
		if(!tau1 || !tau2) throw cms::Exception("RootMaker::AddElTauTauPairs") << "Candidates are not PAT taus";
		DCA dca = calculateDCA(*tau1, *tau2);

		eltautaupair_leg1_px[i] = svFitDaughter1->particle()->px();
		eltautaupair_leg1_py[i] = svFitDaughter1->particle()->py();
		eltautaupair_leg1_pz[i] = svFitDaughter1->particle()->pz();
		eltautaupair_leg1_energy[i] = svFitDaughter1->particle()->energy();
		eltautaupair_leg2_px[i] = svFitDaughter2->particle()->px();
		eltautaupair_leg2_py[i] = svFitDaughter2->particle()->py();
		eltautaupair_leg2_pz[i] = svFitDaughter2->particle()->pz();
		eltautaupair_leg2_energy[i] = svFitDaughter2->particle()->energy();
		eltautaupair_el_px[i] = svFitDaughterEl->particle()->px();
		eltautaupair_el_py[i] = svFitDaughterEl->particle()->py();
		eltautaupair_el_pz[i] = svFitDaughterEl->particle()->pz();
		eltautaupair_el_energy[i] = svFitDaughterEl->particle()->energy();
		eltautaupair_svfit_int_valid[i] = svFitResonanceHypothesis->isValidSolution();
		eltautaupair_svfit_mass_int[i] = svFitResonanceHypothesis->mass();
		eltautaupair_svfit_mass_int_err_up[i] = svFitResonanceHypothesis->massErrUp();
		eltautaupair_svfit_mass_int_err_down[i] = svFitResonanceHypothesis->massErrDown();

		eltautaupair_dca2d[i] = dca.dca2d;
		eltautaupair_dca2d_err[i] = dca.dca2dErr;
		eltautaupair_dca3d[i] = dca.dca3d;
		eltautaupair_dca3d_err[i] = dca.dca3dErr;
	}

	return true;
}

bool RootMaker::AddTracks(const edm::Event& iEvent)
{
  int NumGood = 0;

  edm::Handle<TrackCollection> Tracks;
  iEvent.getByLabel(edm::InputTag("generalTracks"), Tracks);

  edm::Handle<edm::ValueMap<DeDxData> > dEdxharmonic2;
  iEvent.getByLabel(edm::InputTag("dedxHarmonic2"), dEdxharmonic2);

  if(Tracks.isValid()){
    for(unsigned i = 0 ; i < Tracks->size() ; i++)
      {
	TransientTrack TTrack = TTrackBuilder->build((*Tracks)[i]);
	TrajectoryStateClosestToPoint TTrackState = TTrack.trajectoryStateClosestToPoint(GlobalPoint(pv_position.x(), pv_position.y(), pv_position.z()));

	if(TTrackState.pt() >= cTrackFilterPtMin)
	  {
	    track_px[track_count] = TTrackState.momentum().x();
	    track_py[track_count] = TTrackState.momentum().y();
	    track_pz[track_count] = TTrackState.momentum().z();
	    track_closestpointx[track_count] = TTrackState.position().x();
	    track_closestpointy[track_count] = TTrackState.position().y();
	    track_closestpointz[track_count] = TTrackState.position().z();
	    track_dxy[track_count]    = TTrackState.perigeeParameters().transverseImpactParameter();
	    track_dxyerr[track_count] = TTrackState.perigeeError().transverseImpactParameterError();
	    track_dz[track_count]     = TTrackState.perigeeParameters().longitudinalImpactParameter();
	    track_dzerr[track_count]  = TTrackState.perigeeError().longitudinalImpactParameterError();
	    track_chi2[track_count]   = (*Tracks)[i].chi2();
	    track_ndof[track_count]   = (*Tracks)[i].ndof();
	    track_charge[track_count] = (*Tracks)[i].charge();
	    track_nhits[track_count] = (*Tracks)[i].numberOfValidHits();
	    track_nmissinghits[track_count] = (*Tracks)[i].numberOfLostHits();
	    track_npixelhits[track_count] = (*Tracks)[i].hitPattern().numberOfValidPixelHits();
	    track_npixellayers[track_count] = (*Tracks)[i].hitPattern().pixelLayersWithMeasurement();
	    track_nstriplayers[track_count] = (*Tracks)[i].hitPattern().stripLayersWithMeasurement();

	    TrackRef tr  = TrackRef(Tracks, i);
	    track_dedxharmonic2[track_count] = (*dEdxharmonic2)[tr].dEdx();
	    math::XYZPoint ecalPos = PositionOnECalSurface(TTrack);
	    track_outerx[track_count] = ecalPos.x();
	    track_outery[track_count] = ecalPos.y();
	    track_outerz[track_count] = ecalPos.z();
	    track_count++;

	    if(track_count == M_trackmaxcount)
	      {
		cerr << "number of tracks > M_trackmaxcount. They are missing." << endl;
		errors |= 1<<2;
		break;
	      }
	    
	    if(TTrackState.pt() >= cTrackPtMin && fabs((*Tracks)[i].eta()) <= cTrackEtaMax)
	      {
		NumGood++;
		//double energy = sqrt(pow(track_px[i],2) + pow(track_py[i],2) + pow(track_pz[i],2));
		//TrackVector.push_back(TLorentzVector(track_px[i], track_py[i], track_pz[i], energy));
	      }
	  }
      }
    
  }
  if(NumGood >= cTrackNum) return(true);
  return(false);
}

bool RootMaker::AddAK5CaloJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	int NumGood = 0;

	edm::Handle<pat::JetCollection> ak5caloJets;
	iEvent.getByLabel(edm::InputTag("patJetsAK5Calo"), ak5caloJets);
	if(ak5caloJets.isValid())
	{
		for(unsigned i = 0 ; i < ak5caloJets->size() ; i++)
		{
			if((*ak5caloJets)[i].pt() >= cAK5CaloFilterPtMin)
			{
				ak5calojet_e[ak5calojet_count] = (*ak5caloJets)[i].energy();
				ak5calojet_px[ak5calojet_count] = (*ak5caloJets)[i].px();
				ak5calojet_py[ak5calojet_count] = (*ak5caloJets)[i].py();
				ak5calojet_pz[ak5calojet_count] = (*ak5caloJets)[i].pz();
				ak5calojet_hadronicenergy[ak5calojet_count] = (*ak5caloJets)[i].energyFractionHadronic() * (*ak5caloJets)[i].energy();
				ak5calojet_emenergy[ak5calojet_count] = (*ak5caloJets)[i].emEnergyFraction() * (*ak5caloJets)[i].energy();
				if((*ak5caloJets)[i].jecSetsAvailable())
				{
					ak5calojet_energycorr[ak5calojet_count] = (*ak5caloJets)[i].jecFactor("Uncorrected");
				}
				else
				{
					ak5calojet_energycorr[ak5calojet_count] = -1.;
				}
				for(unsigned n = 0 ; n < bdisclabel.size() ; n++)
				{
					if(bdisclabel[n] != "F"){
						ak5calojet_btag[ak5calojet_count][n] = (*ak5caloJets)[i].bDiscriminator(bdisclabel[n]);
					}
					else
					{
						ak5calojet_btag[ak5calojet_count][n] = -1000;
					}
				}

				ak5calojet_fhpd[ak5calojet_count] = (*ak5caloJets)[i].jetID().fHPD;
				ak5calojet_restrictedemf[ak5calojet_count] = (*ak5caloJets)[i].jetID().restrictedEMF;
				/*for(unsigned n = 0 ; n < min((*ak5caloJets)[i].getPairDiscri().size(), unsigned(M_btagmax)) ; n++)
				  {
				  cout << (*ak5caloJets)[i].getPairDiscri()[n].first << " " << (*ak5caloJets)[i].getPairDiscri()[n].second << " " << ak5calojet_btag[i][n] << endl;
				  }*/
				ak5calojet_n90[ak5calojet_count] = (*ak5caloJets)[i].n90();
				ak5calojet_n60[ak5calojet_count] = (*ak5caloJets)[i].n60();
				ak5calojet_count++;
				if(ak5calojet_count == M_jetmaxcount){cerr << "number of ak5calojet > M_jetmaxcount. They are missing." << endl; errors |= 1<<8; break;}
				if((*ak5caloJets)[i].pt() >= cAK5CaloPtMin && fabs((*ak5caloJets)[i].eta()) < cAK5CaloEtaMax) NumGood++;
			}
		}
	}
	if(NumGood >= cAK5CaloNum) return(true);
	return(false);
}

bool RootMaker::AddAK5JPTJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	int NumGood = 0;

	edm::Handle<pat::JetCollection> ak5jptJets;
	iEvent.getByLabel(edm::InputTag("patJetsAK5JPT"), ak5jptJets);
	if(ak5jptJets.isValid())
	{
		for(unsigned i = 0 ; i < ak5jptJets->size() ; i++)
		{
			if((*ak5jptJets)[i].pt() >= cAK5JPTFilterPtMin)
			{
				ak5jptjet_e[ak5jptjet_count] = (*ak5jptJets)[i].energy();
				ak5jptjet_px[ak5jptjet_count] = (*ak5jptJets)[i].px();
				ak5jptjet_py[ak5jptjet_count] = (*ak5jptJets)[i].py();
				ak5jptjet_pz[ak5jptjet_count] = (*ak5jptJets)[i].pz();
				ak5jptjet_hadronicenergy[ak5jptjet_count] = (*ak5jptJets)[i].chargedHadronEnergy() + (*ak5jptJets)[i].neutralHadronEnergy();
				ak5jptjet_chargedhadronicenergy[ak5jptjet_count] = (*ak5jptJets)[i].chargedHadronEnergy();
				ak5jptjet_emenergy[ak5jptjet_count] = (*ak5jptJets)[i].chargedEmEnergy() + (*ak5jptJets)[i].neutralEmEnergy();
				ak5jptjet_chargedemenergy[ak5jptjet_count] = (*ak5jptJets)[i].chargedEmEnergy();
				ak5jptjet_chargedmulti[ak5jptjet_count] = (*ak5jptJets)[i].chargedMultiplicity();
				if((*ak5jptJets)[i].jecSetsAvailable())
				{
					ak5jptjet_energycorr[ak5jptjet_count] = (*ak5jptJets)[i].jecFactor("Uncorrected");
				}
				else
				{
					ak5jptjet_energycorr[ak5jptjet_count] = -1.;
				}
				for(unsigned n = 0 ; n < bdisclabel.size() ; n++)
				{
					if(bdisclabel[n] != "F"){
						ak5jptjet_btag[ak5jptjet_count][n] = (*ak5jptJets)[i].bDiscriminator(bdisclabel[n]);
					}
					else
					{
						ak5jptjet_btag[ak5jptjet_count][n] = -1000;
					}
				}

				ak5jptjet_fhpd[ak5jptjet_count] = (*ak5jptJets)[i].jetID().fHPD;
				ak5jptjet_restrictedemf[ak5jptjet_count] = (*ak5jptJets)[i].jetID().restrictedEMF;
				ak5jptjet_n90[ak5jptjet_count] = (*ak5jptJets)[i].jetID().n90Hits;
				ak5jptjet_count++;
				if(ak5jptjet_count == M_jetmaxcount){cerr << "number of ak5jptjet > M_jetmaxcount. They are missing." << endl; errors |= 1<<9; break;}
				if((*ak5jptJets)[i].pt() >= cAK5JPTPtMin && fabs((*ak5jptJets)[i].eta()) < cAK5JPTEtaMax) NumGood++;
			}
		}
	}

	if(NumGood >= cAK5JPTNum) return(true);
	return(false);
}

bool RootMaker::AddAK5PFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	int NumGood = 0;

	edm::Handle<pat::JetCollection> ak5pfJets;
	iEvent.getByLabel(edm::InputTag("patJets"), ak5pfJets);

	edm::Handle<std::vector<reco::SecondaryVertexTagInfo> > svInfos;
	iEvent.getByLabel(edm::InputTag("secondaryVertexTagInfosAOD"), svInfos);
	assert(svInfos.isValid());

	edm::Handle<edm::ValueMap<float> > puJetIdMVAFull;
	iEvent.getByLabel(edm::InputTag("pileupJetIdProducer","fullDiscriminant"), puJetIdMVAFull);

	edm::Handle<edm::ValueMap<int> > puJetIdFlagFull;
	iEvent.getByLabel(edm::InputTag("pileupJetIdProducer","fullId"), puJetIdFlagFull);

	if(ak5pfJets.isValid())
	  {
	    for(unsigned i = 0 ; i < ak5pfJets->size() ; i++)
	      {
		if((*ak5pfJets)[i].pt() >= cAK5PFFilterPtMin)
		  {
		    //std::cout << "Jet, pT=" <<  (*ak5pfJets)[i].pt() << std::endl;
		    
		    ak5pfjet_e[ak5pfjet_count] = (*ak5pfJets)[i].energy();
		    ak5pfjet_px[ak5pfjet_count] = (*ak5pfJets)[i].px();
		    ak5pfjet_py[ak5pfjet_count] = (*ak5pfJets)[i].py();
		    ak5pfjet_pz[ak5pfjet_count] = (*ak5pfJets)[i].pz();
		    ak5pfjet_hadronicenergy[ak5pfjet_count] = (*ak5pfJets)[i].chargedHadronEnergy() + (*ak5pfJets)[i].neutralHadronEnergy();
		    ak5pfjet_chargedhadronicenergy[ak5pfjet_count] = (*ak5pfJets)[i].chargedHadronEnergy();
		    ak5pfjet_emenergy[ak5pfjet_count] = (*ak5pfJets)[i].chargedEmEnergy() + (*ak5pfJets)[i].neutralEmEnergy();
		    ak5pfjet_chargedemenergy[ak5pfjet_count] = (*ak5pfJets)[i].chargedEmEnergy();
		    ak5pfjet_chargedmulti[ak5pfjet_count] = (*ak5pfJets)[i].chargedMultiplicity();
		    ak5pfjet_neutralmulti[ak5pfjet_count] = (*ak5pfJets)[i].neutralMultiplicity();
		    //float pt = (*ak5pfJets)[i].pt();
		    //cout << (*ak5pfJets)[i].availableJECSets().size() << endl;
		    //cout << (*ak5pfJets)[i].jecFactor("Uncorrected")*pt << " " << (*ak5pfJets)[i].jecFactor("L1FastJet")*pt << " " << (*ak5pfJets)[i].jecFactor("L2Relative")*pt << " " << (*ak5pfJets)[i].jecFactor("L3Absolute")*pt << " " << (*ak5pfJets)[i].jecFactor("L2L3Residual")*pt << endl;
		    if((*ak5pfJets)[i].jecSetsAvailable())
		      {
			ak5pfjet_energycorr[ak5pfjet_count] = (*ak5pfJets)[i].jecFactor("Uncorrected");
		      }
		    else
		      {
			ak5pfjet_energycorr[ak5pfjet_count] = -1.;
		      }
		    
		    // Secondary vertex information
		    ak5pfjet_sv_count[ak5pfjet_count] = 0;
		    ak5pfjet_sv_x[ak5pfjet_count] = 0.0f;
		    ak5pfjet_sv_y[ak5pfjet_count] = 0.0f;
		    ak5pfjet_sv_z[ak5pfjet_count] = 0.0f;
		    ak5pfjet_sv_dx[ak5pfjet_count] = 0.0f;
		    ak5pfjet_sv_dy[ak5pfjet_count] = 0.0f;
		    ak5pfjet_sv_dz[ak5pfjet_count] = 0.0f;
		    for(unsigned n = 0; n < svInfos->size(); ++n)
		      {
			//if(ROOT::Math::VectorUtil::DeltaR((*svInfos)[n].jet()->p4(), (*ak5pfJets)[i].p4()) < 0.0001)
			if((*svInfos)[n].jet().get() == (*ak5pfJets)[i].originalObject())
			  {
			    ak5pfjet_sv_count[ak5pfjet_count] = (*svInfos)[n].nVertices();
			    //std::cout << "Matched " << i << " and " << n << std::endl;
			    //for(int k = 0; k < min((*svInfos)[n].nVertices(), M_svmax); ++k)
			    if((*svInfos)[n].nVertices() > 0)
			      {
				ak5pfjet_sv_x[ak5pfjet_count] = (*svInfos)[n].secondaryVertex(0).x();
				ak5pfjet_sv_y[ak5pfjet_count] = (*svInfos)[n].secondaryVertex(0).y();
				ak5pfjet_sv_z[ak5pfjet_count] = (*svInfos)[n].secondaryVertex(0).z();
				ak5pfjet_sv_dx[ak5pfjet_count] = (*svInfos)[n].secondaryVertex(0).xError();
				ak5pfjet_sv_dy[ak5pfjet_count] = (*svInfos)[n].secondaryVertex(0).yError();
				ak5pfjet_sv_dz[ak5pfjet_count] = (*svInfos)[n].secondaryVertex(0).zError();
			      }
			  }
		      }
		    
		    //edm::Ref<pat::JetCollection> jetRef(ak5pfJets, i);
		    
		    ak5pfjet_pu_jet_simple_loose[ak5pfjet_count] = false;
		    ak5pfjet_pu_jet_simple_medium[ak5pfjet_count] = false;
		    ak5pfjet_pu_jet_simple_tight[ak5pfjet_count] = false;
		    ak5pfjet_pu_jet_simple_mva[ak5pfjet_count] = -1.0f;

		    ak5pfjet_pu_jet_full_loose[ak5pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*ak5pfJets)[i].originalObjectRef()], PileupJetIdentifier::kLoose);
		    ak5pfjet_pu_jet_full_medium[ak5pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*ak5pfJets)[i].originalObjectRef()], PileupJetIdentifier::kMedium);
		    ak5pfjet_pu_jet_full_tight[ak5pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*ak5pfJets)[i].originalObjectRef()], PileupJetIdentifier::kTight);
		    ak5pfjet_pu_jet_full_mva[ak5pfjet_count] = (*puJetIdMVAFull)[(*ak5pfJets)[i].originalObjectRef()];
		    
		    for(unsigned n = 0 ; n < bdisclabel.size() ; n++)
		      {
			if(bdisclabel[n] != "F"){
			  ak5pfjet_btag[ak5pfjet_count][n] = (*ak5pfJets)[i].bDiscriminator(bdisclabel[n]);
			}
			else
			  {
			    ak5pfjet_btag[ak5pfjet_count][n] = -1000;
			  }
		      }
		    ak5pfjet_trigger[ak5pfjet_count] = GetTriggerMatch((*ak5pfJets)[i], jettriggers); 
		    ak5pfjet_count++;
		    if(ak5pfjet_count == M_jetmaxcount){cerr << "number of ak5pfjet > M_jetmaxcount. They are missing." << endl; errors |= 1<<10; break;}
		    if((*ak5pfJets)[i].pt() >= cAK5PFPtMin && fabs((*ak5pfJets)[i].eta()) < cAK5PFEtaMax) NumGood++;
		  }
	      }
	  }
	
	if(NumGood >= cAK5PFNum) return(true);
	return(false);
}

/*bool RootMaker::AddAK5PFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	int NumGood = 0;
	edm::Handle<PFJetCollection> ak5pfJets;
	iEvent.getByLabel(edm::InputTag("ak5PFJets"), ak5pfJets);
	if(ak5pfJets.isValid())
	{
		vector<edm::Handle<JetFloatAssociation::Container> > btag(bdisclabel.size());
		for(unsigned n = 0 ; n < bdisclabel.size() ; n++)
		{
			iEvent.getByLabel(edm::InputTag(bdisclabel[n].c_str()), btag[n]);
		}

		const JetCorrector* corrector;
		if(cdata)
		{
			corrector = JetCorrector::getJetCorrector("ak5PFL1L2L3Residual", iSetup);
		}
		else
		{
			corrector = JetCorrector::getJetCorrector("ak5PFL1L2L3", iSetup);
		}
		
		for(unsigned i = 0 ; i < ak5pfJets->size() ; i++)
		{
			if((*ak5pfJets)[i].pt() >= cAK5PFFilterPtMin)
			{
				ak5pfjet_e[ak5pfjet_count] = (*ak5pfJets)[i].energy();
				ak5pfjet_px[ak5pfjet_count] = (*ak5pfJets)[i].px();
				ak5pfjet_py[ak5pfjet_count] = (*ak5pfJets)[i].py();
				ak5pfjet_pz[ak5pfjet_count] = (*ak5pfJets)[i].pz();
				ak5pfjet_hadronicenergy[ak5pfjet_count] = (*ak5pfJets)[i].chargedHadronEnergy() + (*ak5pfJets)[i].neutralHadronEnergy();
				ak5pfjet_chargedhadronicenergy[ak5pfjet_count] = (*ak5pfJets)[i].chargedHadronEnergy();
				ak5pfjet_emenergy[ak5pfjet_count] = (*ak5pfJets)[i].chargedEmEnergy() + (*ak5pfJets)[i].neutralEmEnergy();
				ak5pfjet_chargedemenergy[ak5pfjet_count] = (*ak5pfJets)[i].chargedEmEnergy();
				ak5pfjet_chargedmulti[ak5pfjet_count] = (*ak5pfJets)[i].chargedMultiplicity();
				ak5pfjet_neutralmulti[ak5pfjet_count] = (*ak5pfJets)[i].neutralMultiplicity();
				JetBaseRef jetRef(edm::Ref<PFJetCollection>(ak5pfJets, i));
				ak5pfjet_energycorr[ak5pfjet_count] = corrector->correction((*ak5pfJets)[i], jetRef, iEvent, iSetup);

				for(unsigned n = 0 ; n < bdisclabel.size() ; n++)
				{
					if(btag[n].isValid())
					{
						//ak5pfjet_btag[ak5pfjet_count][n] = (*(btag[n]))[jetRef];
						//ak5pfjet_btag[ak5pfjet_count][n] = JetFloatAssociation::getValue(*(btag[n]), (*ak5pfJets)[i]);
						cout << "btag1: " << n << " " << ak5pfjet_btag[ak5pfjet_count][n] << endl;
					}
					else
					{
						ak5pfjet_btag[ak5pfjet_count][n] = -1000;
					}
				}

				if(crecjettrigger)
				{
					ak5pfjet_trigger[i] = GetTriggerMatch((*ak5pfJets)[i], jettriggers); 
				}
				else
				{
					ak5pfjet_trigger[i] = 0;
				}
				ak5pfjet_count++;
				if(ak5pfjet_count == M_jetmaxcount){cerr << "number of ak5pfjet > M_jetmaxcount. They are missing." << endl; errors |= 1<<10; break;}
				if((*ak5pfJets)[i].pt() >= cAK5PFPtMin && fabs((*ak5pfJets)[i].eta()) < cAK5PFEtaMax) NumGood++;
			}
		}
	}

	if(NumGood >= cAK5PFNum) return(true);
	return(false);
}*/

int RootMaker::AddElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	edm::Handle<double> hRho;
	edm::InputTag tag("kt6PFJets","rho");
	iEvent.getByLabel(tag,hRho);

	edm::Handle<reco::BeamSpot> hBeamSpot;
	iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), hBeamSpot);

	edm::Handle<reco::ConversionCollection> hConversions;
	iEvent.getByLabel(edm::InputTag("allConversions"), hConversions);

	edm::Handle<std::vector<reco::Vertex> > hVertices;
	iEvent.getByLabel(edm::InputTag("offlinePrimaryVertices"), hVertices);
	const reco::Vertex& primaryVtx = (*hVertices)[0];

	edm::Handle<std::vector<reco::PFCandidate> > hPfCandidates;
	iEvent.getByLabel(edm::InputTag("particleFlow"), hPfCandidates);
	//const std::vector<reco::PFCandidate>& pfCandidates = *hPfCandidates;

#if 0
	std::vector<reco::Muon> IdentifiedMuons;
	std::vector<reco::GsfElectron> IdentifiedElectrons;
#endif

	edm::Handle<GsfElectronCollection> Electrons;
	iEvent.getByLabel(edm::InputTag("gsfElectrons"), Electrons);

	// ecal cluster shape variables
	EcalClusterLazyTools lazyTools(iEvent, iSetup,
		edm::InputTag("reducedEcalRecHitsEB"),
		edm::InputTag("reducedEcalRecHitsEE"));

	//edm::Handle<pat::MuonCollection> Muons;
	//iEvent.getByLabel(edm::InputTag("patMuons"), Muons);
	edm::Handle<MuonCollection> inMuons;
	iEvent.getByLabel(edm::InputTag("muons"), inMuons);

#if 0
	for (std::vector<reco::GsfElectron>::const_iterator iE = Electrons->begin(); iE != Electrons->end(); ++iE)
	{
	        double electronTrackZ = 0;
	        if (iE->gsfTrack().isNonnull())
	        {
	                electronTrackZ = iE->gsfTrack()->dz(primaryVtx.position());
	        }
	        else if (iE->closestCtfTrackRef().isNonnull())
	        {
	                electronTrackZ = iE->closestCtfTrackRef()->dz(primaryVtx.position());
	        }
	        if(fabs(electronTrackZ) > 0.2)
	                continue;
	        if(fabs(iE->superCluster()->eta())<1.479)
	        {
	                if(iE->pt() > 20)
	                {
	                        if(iE->sigmaIetaIeta() > 0.01) continue;
	                        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
	                        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
	                        if(iE->hadronicOverEm() > 0.15) continue;
	                }
	                else
	                {
	                        if(iE->sigmaIetaIeta() > 0.012) continue;
	                        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
	                        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
	                        if(iE->hadronicOverEm() > 0.15) continue;
	                }
	        }
	        else
	        {
	                if(iE->pt() > 20)
	                {
	                        if(iE->sigmaIetaIeta() > 0.03) continue;
	                        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
	                        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8) continue;
	                }
	                else
	                {
	                        if(iE->sigmaIetaIeta() > 0.032) continue;
	                        if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
	                        if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8) continue;
	                }
	        }
	        IdentifiedElectrons.push_back(*iE);
	}

	for (MuonCollection::const_iterator iM = inMuons->begin(); iM != inMuons->end(); ++iM)
	{
	        if(!(iM->innerTrack().isNonnull()))
	                continue;
	        if(!(iM->isGlobalMuon() || iM->isTrackerMuon()))
	                continue;
	        if(iM->innerTrack()->numberOfValidHits() < 11 )
	                continue;
	        IdentifiedMuons.push_back(*iM);
	}

	ElectronEffectiveArea::ElectronEffectiveAreaTarget currentAreaTarget = ElectronEffectiveArea::kEleEANoCorr;
	if (cYear == 2011)
	{
	        if (iEvent.isRealData())
	                currentAreaTarget = ElectronEffectiveArea::kEleEAData2011;
	        else
	        {
	                if (cPeriod == "Summer11")
	                        currentAreaTarget = ElectronEffectiveArea::kEleEASummer11MC;
	                if (cPeriod == "Fall11")
	                        currentAreaTarget = ElectronEffectiveArea::kEleEAFall11MC;
	        }
	}
	if (cYear == 2012)
	{
	        if (iEvent.isRealData())
	                currentAreaTarget = ElectronEffectiveArea::kEleEAData2012;
	}
#endif

	edm::Handle<edm::ValueMap<reco::IsoDeposit> > isoDepsChargedHadrons, isoDepsChargedParticles, isoDepsNeutralHadrons, isoDepsPhotons, isoDepsPU;
	iEvent.getByLabel(edm::InputTag("pfElIsoDepositChargedHadrons"), isoDepsChargedHadrons);
	iEvent.getByLabel(edm::InputTag("pfElIsoDepositChargedParticles"), isoDepsChargedParticles);
	iEvent.getByLabel(edm::InputTag("pfElIsoDepositNeutralHadrons"), isoDepsNeutralHadrons);
	iEvent.getByLabel(edm::InputTag("pfElIsoDepositPhotons"), isoDepsPhotons);
	iEvent.getByLabel(edm::InputTag("pfElIsoDepositPU"), isoDepsPU);

	edm::Handle<edm::ValueMap<float> > mvaIdTrig;
	edm::Handle<edm::ValueMap<float> > mvaIdNonTrig;
	iEvent.getByLabel(edm::InputTag("mvaTrigV0"), mvaIdTrig);
	iEvent.getByLabel(edm::InputTag("mvaNonTrigV0"), mvaIdNonTrig);

	int nGoodElectrons = 0;
	/*if(crecelectrontrigger)
	{
		iEvent.getByLabel(edm::InputTag("l1extraParticles", "NonIsolated"), L1Electrons);
		iEvent.getByLabel(edm::InputTag("l1extraParticles", "Isolated"), L1ElectronsIso);
	}*/

	assert(Electrons.isValid());

	for(unsigned i = 0 ; i < Electrons->size() ; i++){
	  electron_px[i] = (*Electrons)[i].px();
	  electron_py[i] = (*Electrons)[i].py();
	  electron_pz[i] = (*Electrons)[i].pz();
	  electron_charge[i] = (*Electrons)[i].charge();
	  electron_esuperclusterovertrack[i] = (*Electrons)[i].eSuperClusterOverP();
	  electron_eseedclusterovertrack[i] = (*Electrons)[i].eSeedClusterOverP();
	  electron_deltaetasuperclustertrack[i] = (*Electrons)[i].deltaEtaSuperClusterTrackAtVtx();
	  electron_deltaphisuperclustertrack[i] = (*Electrons)[i].deltaEtaSuperClusterTrackAtVtx();
	  electron_e1x5[i] = (*Electrons)[i].e1x5();
	  electron_e2x5[i] = (*Electrons)[i].e2x5Max();
	  electron_e5x5[i] = (*Electrons)[i].e5x5();
	  electron_sigmaetaeta[i] = (*Electrons)[i].sigmaEtaEta();
	  electron_sigmaietaieta[i] = (*Electrons)[i].sigmaIetaIeta();
	  electron_ehcaloverecal[i] = (*Electrons)[i].hcalOverEcal();
	  electron_ehcaloverecaldepth1[i] = (*Electrons)[i].hcalDepth1OverEcal();
	  electron_ehcaloverecaldepth2[i] = (*Electrons)[i].hcalDepth2OverEcal();
	  electron_isolationr3track[i] = (*Electrons)[i].dr03TkSumPt();
	  electron_isolationr3ecal[i] = (*Electrons)[i].dr03EcalRecHitSumEt();
	  electron_isolationr3hcal[i] = (*Electrons)[i].dr03HcalTowerSumEt();
	  electron_isolationr4track[i] = (*Electrons)[i].dr04TkSumPt();
	  electron_isolationr4ecal[i] = (*Electrons)[i].dr04EcalRecHitSumEt();
	  electron_isolationr4hcal[i] = (*Electrons)[i].dr04HcalTowerSumEt();
	  electron_info[i] = 0;
	  electron_info[i] |= (*Electrons)[i].isElectron(); 
	  
	  edm::Ref<std::vector<reco::GsfElectron> > electronref(Electrons, i);
	  reco::IsoDeposit electronIsoDepositChargedHadrons = (*isoDepsChargedHadrons)[electronref];
	  reco::IsoDeposit electronIsoDepositChargedParticles = (*isoDepsChargedParticles)[electronref];
	  reco::IsoDeposit electronIsoDepositNeutralHadrons = (*isoDepsNeutralHadrons)[electronref];
	  reco::IsoDeposit electronIsoDepositPhotons = (*isoDepsPhotons)[electronref];
	  reco::IsoDeposit electronIsoDepositPU = (*isoDepsPU)[electronref];
	  
	  reco::isodeposit::Direction dir = reco::isodeposit::Direction((*Electrons)[i].eta(), (*Electrons)[i].phi());
	  reco::isodeposit::ConeVeto pf_cone_veto_charged(dir, (*Electrons)[i].isEB() ? 0.01 : 0.015);
	  reco::isodeposit::ConeVeto pf_cone_veto_photons(dir, (*Electrons)[i].isEB() ? 0.08 : 0.08);
	  reco::isodeposit::ConeVeto pf_cone_veto_other(dir, 0.0);
	  reco::isodeposit::ThresholdVeto pf_threshold_veto(0.0);
	  
	  std::vector<reco::isodeposit::AbsVeto*> vetosPFCharged, vetosPFPhotons, vetosPFOther;
	  vetosPFCharged.push_back(&pf_cone_veto_charged);
	  vetosPFCharged.push_back(&pf_threshold_veto);
	  vetosPFPhotons.push_back(&pf_cone_veto_photons);
	  vetosPFPhotons.push_back(&pf_threshold_veto);
	  vetosPFOther.push_back(&pf_threshold_veto);
	  
	  electron_pfisolationr3_sumchargedhadronpt[i] = electronIsoDepositChargedHadrons.depositWithin(0.3, vetosPFCharged);
	  electron_pfisolationr3_sumchargedparticlept[i] = electronIsoDepositChargedParticles.depositWithin(0.3, vetosPFCharged);
	  electron_pfisolationr3_sumneutralhadronet[i] = electronIsoDepositNeutralHadrons.depositWithin(0.3, vetosPFOther);
	  electron_pfisolationr3_sumphotonet[i] = electronIsoDepositPhotons.depositWithin(0.3, vetosPFPhotons);
	  electron_pfisolationr3_sumPUpt[i] = electronIsoDepositPU.depositWithin(0.3, vetosPFOther);
	  electron_pfisolationr4_sumchargedhadronpt[i] = electronIsoDepositChargedHadrons.depositWithin(0.4, vetosPFCharged);
	  electron_pfisolationr4_sumchargedparticlept[i] = electronIsoDepositChargedParticles.depositWithin(0.4, vetosPFCharged);
	  electron_pfisolationr4_sumneutralhadronet[i] = electronIsoDepositNeutralHadrons.depositWithin(0.4, vetosPFOther);
	  electron_pfisolationr4_sumphotonet[i] = electronIsoDepositPhotons.depositWithin(0.4, vetosPFPhotons);
	  electron_pfisolationr4_sumPUpt[i] = electronIsoDepositPU.depositWithin(0.4, vetosPFOther);
	  
	  electron_gapinfo[i] = 0;
	  electron_gapinfo[i] |= (*Electrons)[i].isEB() << 0;
	  electron_gapinfo[i] |= (*Electrons)[i].isEE() << 1;
	  electron_gapinfo[i] |= (*Electrons)[i].isEBGap() << 2;
	  electron_gapinfo[i] |= (*Electrons)[i].isEBEtaGap() << 3;
	  electron_gapinfo[i] |= (*Electrons)[i].isEBPhiGap() << 4;
	  electron_gapinfo[i] |= (*Electrons)[i].isEEGap() << 5;
	  electron_gapinfo[i] |= (*Electrons)[i].isEERingGap() << 6;
	  electron_gapinfo[i] |= (*Electrons)[i].isEEDeeGap() << 7;
	  electron_gapinfo[i] |= (*Electrons)[i].isEBEEGap() << 8;
	  
	  electron_chargeinfo[i] = 0;
	  if((*Electrons)[i].isGsfCtfChargeConsistent()) electron_chargeinfo[i] |= (1 << 0);
	  if((*Electrons)[i].isGsfCtfScPixChargeConsistent()) electron_chargeinfo[i] |= (1 << 1);
	  if((*Electrons)[i].isGsfScPixChargeConsistent()) electron_chargeinfo[i] |= (1 << 2);
	  
	  edm::Handle<reco::TrackCollection> ctfTracks;
	  iEvent.getByLabel("generalTracks", ctfTracks);
	  ConversionFinder convFinder;
	  ConversionInfo convInfo = convFinder.getConversionInfo((*Electrons)[i], ctfTracks, magneticField->inTesla(GlobalPoint(0.,0.,0.)).z());
	  
	  electron_convdist[i] = convInfo.dist();
	  electron_convdcot[i] = convInfo.dcot();
	  electron_convradius[i] = convInfo.radiusOfConversion();
	  electron_fbrems[i] = (*Electrons)[i].fbrem();
	  electron_numbrems[i] = (*Electrons)[i].numberOfBrems();
	  
	  GsfTrackRef gsfTr_e = (*Electrons)[i].gsfTrack();
	  TransientTrack TTrack = TTrackBuilder->build(gsfTr_e);
	  math::XYZPoint ecalPos = PositionOnECalSurface(TTrack);
	  electron_outerx[i] = ecalPos.x();
	  electron_outery[i] = ecalPos.y();
	  electron_outerz[i] = ecalPos.z();
	  //TrajectoryStateClosestToPoint TTrackState = TTrack.trajectoryStateClosestToPoint(GlobalPoint(pv_position.x(), pv_position.y(), pv_position.z()));
	  
	  electron_trackchi2[i] = gsfTr_e->chi2();
	  electron_trackndof[i] = gsfTr_e->ndof();
	  electron_closestpointx[i] = gsfTr_e->vx();
	  electron_closestpointy[i] = gsfTr_e->vy();
	  electron_closestpointz[i] = gsfTr_e->vz();
	  
	  electron_nhits[i]        = gsfTr_e->numberOfValidHits();
	  electron_nmissinghits[i] = gsfTr_e->numberOfLostHits();
	  electron_npixelhits[i]   = (gsfTr_e->hitPattern()).numberOfValidPixelHits();
	  electron_npixellayers[i]   = (gsfTr_e->hitPattern()).pixelLayersWithMeasurement();
	  electron_nstriplayers[i]   = (gsfTr_e->hitPattern()).stripLayersWithMeasurement();
	  
	  electron_dxy[i]          = gsfTr_e->dxy(pv_position);
	  electron_dxyerr[i]       = gsfTr_e->dxyError();
	  electron_dz[i]           = gsfTr_e->dz(pv_position);
	  electron_dzerr[i]        = gsfTr_e->dzError();
	  
	  electron_superclusterindex[i] = getSuperCluster((*Electrons)[i].superCluster()->energy(), (*Electrons)[i].superCluster()->x(), (*Electrons)[i].superCluster()->y(), (*Electrons)[i].superCluster()->z()); 
	  
	  electron_trigger[i] = GetTriggerMatch((*Electrons)[i], electrontriggers); 
	  
	  
	  electron_mva_id_trig[i] = (*mvaIdTrig)[electronref];
	  electron_mva_id_nontrig[i] = (*mvaIdNonTrig)[electronref];
	  electron_mva_iso[i] = 0.0f;
	  
	  electron_has_conversion[i] = ConversionTools::hasMatchedConversion((*Electrons)[i], hConversions, hBeamSpot->position(), true, 2.0, 1e-6, 0);
	  
	  electron_count++;
	  if(electron_count == M_electronmaxcount)
	    {
	      cerr << "number of electron > M_electronmaxcount. They are missing." << endl;
	      errors |= 1<<1;
	      break;
	    }
	  if((*Electrons)[i].pt() >= 10 && fabs((*Electrons)[i].eta()) <= 2.5)
	    {
	      ++nGoodElectrons;
	    }
	}
	
	return nGoodElectrons;
}

Int_t RootMaker::getSuperCluster(float e, float x, float y, float z)
{
  Int_t result = -1;
  for(UInt_t i = 0 ; i < supercluster_count ; i++)
    {
      if(supercluster_e[i] == e && supercluster_x[i] == x && supercluster_y[i] == y && supercluster_z[i] == z) {result = i; break;}
    }
  return(result);
}



bool RootMaker::AddVertices(const edm::Event& iEvent)
{
	double mpi = 0.13957018;
	double mp = 0.93827203;
	edm::Handle<TrackCollection> Tracks;
	iEvent.getByLabel(edm::InputTag("generalTracks"), Tracks);
	edm::Handle<edm::ValueMap<DeDxData> > dEdxharmonic2;
	iEvent.getByLabel(edm::InputTag("dedxHarmonic2"), dEdxharmonic2);
	if(Tracks.isValid() && primvertex_count > 0)
	{
		KalmanVertexFitter fitter(true);

		for(unsigned a = 0 ; a < Tracks->size() ; a++)
		{
			if(!((*Tracks)[a].normalizedChi2() < cVertexTRKChi2 && (*Tracks)[a].numberOfValidHits() >= cVertexTRKHitsMin)){ continue;}
			for(unsigned b = 0 ; b < a ; b++)
			{
				if(!((*Tracks)[a].charge() != (*Tracks)[b].charge() && (*Tracks)[b].normalizedChi2() < cVertexTRKChi2 && (*Tracks)[b].numberOfValidHits() >= cVertexTRKHitsMin)){ continue;}
				vector<TransientTrack> testvec;
				testvec.resize(2);
				testvec[0] = TTrackBuilder->build((*Tracks)[a]);
				testvec[1] = TTrackBuilder->build((*Tracks)[b]);

				TransientVertex testvertex = fitter.vertex(testvec);
				if(testvertex.isValid() && testvertex.normalisedChiSquared() < cVertexChi2)
				{
					const Vertex& secvertex = (Vertex)testvertex;
					SMatrixSym3D covmatrix = primvertex.covariance() + secvertex.covariance();
					SVector3 dist(primvertex.x() - secvertex.x(), primvertex.y() - secvertex.y(), 0.);
					double vertexsig2d =  ROOT::Math::Mag2(dist) / sqrt(ROOT::Math::Similarity(covmatrix, dist));

					if(vertexsig2d > cVertexSig2D)
					{
						using namespace TMath;
						TLorentzVector A((*Tracks)[a].px(), (*Tracks)[a].py(), (*Tracks)[a].pz(), Sqrt(Power((*Tracks)[a].px(), 2) + Power((*Tracks)[a].py(),2) + Power((*Tracks)[a].pz(), 2) + mpi*mpi));
						TLorentzVector B((*Tracks)[b].px(), (*Tracks)[b].py(), (*Tracks)[b].pz(), Sqrt(Power((*Tracks)[b].px(), 2) + Power((*Tracks)[b].py(),2) + Power((*Tracks)[b].pz(), 2) + mpi*mpi));
						double mkaon = (A+B).M();
						A.SetE(Sqrt(A.E()*A.E() - mpi*mpi + mp*mp));
						double mlambdaa = (A+B).M();
						A.SetE(Sqrt(A.E()*A.E() + mpi*mpi - mp*mp));
						B.SetE(Sqrt(B.E()*B.E() - mpi*mpi + mp*mp));
						double mlambdab = (A+B).M();
						if(Abs(mkaon - 0.49765) < cKaonMassWindow || Abs(mlambdaa - 1.115683) < cLambdaMassWindow || Abs(mlambdab - 1.115683) < cLambdaMassWindow)
						{

							secvertices_vx[secvertices_count] = secvertex.x();
							secvertices_vy[secvertices_count] = secvertex.y();
							secvertices_vz[secvertices_count] = secvertex.z();
							secvertices_chi2[secvertices_count] = secvertex.chi2();
							secvertices_ndof[secvertices_count] = secvertex.ndof();
							secvertices_cov[secvertices_count][0] = secvertex.covariance(0,0);
							secvertices_cov[secvertices_count][1] = secvertex.covariance(0,1);
							secvertices_cov[secvertices_count][2] = secvertex.covariance(0,2);
							secvertices_cov[secvertices_count][3] = secvertex.covariance(1,1);
							secvertices_cov[secvertices_count][4] = secvertex.covariance(1,2);
							secvertices_cov[secvertices_count][5] = secvertex.covariance(2,2);
							for(Int_t j = 0 ; j < 2 ; j++)
							{
								TrajectoryStateClosestToPoint TTrackState = testvec[j].trajectoryStateClosestToPoint(GlobalPoint(secvertex.x(), secvertex.y(), secvertex.z()));									
								TrackRef tr;
								if(j==0) tr = TrackRef(Tracks, a);
								if(j==1) tr = TrackRef(Tracks, b);
								secvertices_track_px[secvertices_count][j] = TTrackState.momentum().x();
								secvertices_track_py[secvertices_count][j] = TTrackState.momentum().y();
								secvertices_track_pz[secvertices_count][j] = TTrackState.momentum().z();
								secvertices_track_closestpointx[secvertices_count][j] = TTrackState.position().x();
								secvertices_track_closestpointy[secvertices_count][j] = TTrackState.position().y();
								secvertices_track_closestpointz[secvertices_count][j] = TTrackState.position().z();
								secvertices_track_chi2[secvertices_count][j]   = tr->chi2();
								secvertices_track_ndof[secvertices_count][j]   = tr->ndof();
								secvertices_track_dxy[secvertices_count][j]    = TTrackState.perigeeParameters().transverseImpactParameter();
								secvertices_track_dxyerr[secvertices_count][j] = TTrackState.perigeeError().transverseImpactParameterError();
								secvertices_track_dz[secvertices_count][j]     = TTrackState.perigeeParameters().longitudinalImpactParameter();
								secvertices_track_dzerr[secvertices_count][j]  = TTrackState.perigeeError().longitudinalImpactParameterError();
								secvertices_track_dedxharmonic2[secvertices_count][j] = (*dEdxharmonic2)[tr].dEdx();
								secvertices_track_charge[secvertices_count][j] = tr->charge();
								secvertices_track_nhits[secvertices_count][j] = tr->numberOfValidHits();
								secvertices_track_nmissinghits[secvertices_count][j] = tr->numberOfLostHits();
								secvertices_track_npixelhits[secvertices_count][j] = tr->hitPattern().numberOfValidPixelHits();
								secvertices_track_npixellayers[secvertices_count][j] = tr->hitPattern().pixelLayersWithMeasurement();
								secvertices_track_nstriplayers[secvertices_count][j] = tr->hitPattern().stripLayersWithMeasurement();
							}
							secvertices_count++;
							if(secvertices_count == M_secverticesmaxcount){cerr << "number of secvertices > M_secverticesmaxcount. They are missing." << endl; break;}
						}

					}
				}
			}
		}
	}
	return(true);
}

bool RootMaker::AddMuVertices(const edm::Event& iEvent)
{
	edm::Handle<edm::ValueMap<DeDxData> > dEdxharmonic2;
	iEvent.getByLabel(edm::InputTag("dedxHarmonic2"), dEdxharmonic2);
	edm::Handle<MuonCollection> Muons;
	iEvent.getByLabel(edm::InputTag("muons"), Muons);
	int count = 0;
	if(Muons.isValid() && primvertex_count > 0)
	{
		vector<TransientTrack> testvec(2);
		KalmanVertexFitter fitter(true);
		// 1st loop over muons
		for(unsigned int i = 0; i < Muons->size(); ++i)
		{
			if(!((*Muons)[i].isGlobalMuon() && (*Muons)[i].isTrackerMuon() && (*Muons)[i].innerTrack().isNonnull())) continue;
			TrackRef trA = (*Muons)[i].innerTrack();
			// 2nd loop over muons
			for(unsigned int j = i+1; j < Muons->size(); ++j)
			{
				if(!((*Muons)[j].isGlobalMuon() && (*Muons)[j].isTrackerMuon() && (*Muons)[j].innerTrack().isNonnull())) continue;
				if(!((*Muons)[i].charge()*(*Muons)[j].charge() < 0)) continue;
				TrackRef trB = (*Muons)[j].innerTrack();
				testvec[0] = TTrackBuilder->build(trA);
				testvec[1] = TTrackBuilder->build(trB);

				TransientVertex testvertex = fitter.vertex(testvec);
				if( testvertex.isValid())
				{
					const Vertex& secvertex = (Vertex)testvertex;
					using namespace TMath;
					TLorentzVector A(trA->px(), trA->py(), trA->pz(), Sqrt(Power(trA->px(), 2) + Power(trA->py(), 2) + Power(trA->pz(), 2)));
					TLorentzVector B(trB->px(), trB->py(), trB->pz(), Sqrt(Power(trB->px(), 2) + Power(trB->py(), 2) + Power(trB->pz(), 2)));
					double mass = (A+B).M();
					if(Abs(mass - 3.097) < 1.5 || Abs(mass - 9.46) < 1.5) //J/Psi and Y
					{
						count++;
						musecvertices_vx[musecvertices_count] = secvertex.x();
						musecvertices_vy[musecvertices_count] = secvertex.y();
						musecvertices_vz[musecvertices_count] = secvertex.z();
						musecvertices_chi2[musecvertices_count] = secvertex.chi2();
						musecvertices_ndof[musecvertices_count] = secvertex.ndof();
						musecvertices_cov[musecvertices_count][0] = secvertex.covariance(0,0);
						musecvertices_cov[musecvertices_count][1] = secvertex.covariance(0,1);
						musecvertices_cov[musecvertices_count][2] = secvertex.covariance(0,2);
						musecvertices_cov[musecvertices_count][3] = secvertex.covariance(1,1);
						musecvertices_cov[musecvertices_count][4] = secvertex.covariance(1,2);
						musecvertices_cov[musecvertices_count][5] = secvertex.covariance(2,2);

						for( int j = 0; j < 2; ++j)
						{
							TrajectoryStateClosestToPoint TTrackState = testvec[j].trajectoryStateClosestToPoint(GlobalPoint(secvertex.x(), secvertex.y(), secvertex.z()));
							TrackRef tr;
							if(j==0) tr = (*Muons)[i].innerTrack();
							if(j==1) tr = (*Muons)[j].innerTrack();

							musecvertices_track_px[musecvertices_count][j]            = TTrackState.momentum().x();
							musecvertices_track_py[musecvertices_count][j]            = TTrackState.momentum().y();
							musecvertices_track_pz[musecvertices_count][j]            = TTrackState.momentum().z();
							musecvertices_track_closestpointx[musecvertices_count][j] = TTrackState.position().x();
							musecvertices_track_closestpointy[musecvertices_count][j] = TTrackState.position().y();
							musecvertices_track_closestpointz[musecvertices_count][j] = TTrackState.position().z();
							musecvertices_track_chi2[musecvertices_count][j]          = tr->chi2();
							musecvertices_track_ndof[musecvertices_count][j]          = tr->ndof();
							musecvertices_track_dxy[musecvertices_count][j]           = TTrackState.perigeeParameters().transverseImpactParameter();
							musecvertices_track_dxyerr[musecvertices_count][j]        = TTrackState.perigeeError().transverseImpactParameterError();
							musecvertices_track_dz[musecvertices_count][j]            = TTrackState.perigeeParameters().longitudinalImpactParameter();
							musecvertices_track_dzerr[musecvertices_count][j]         = TTrackState.perigeeError().longitudinalImpactParameterError();
							musecvertices_track_dedxharmonic2[musecvertices_count][j] = (*dEdxharmonic2)[tr].dEdx();
							musecvertices_track_charge[musecvertices_count][j]        = tr->charge();
							musecvertices_track_nhits[musecvertices_count][j]         = tr->numberOfValidHits();
							musecvertices_track_nmissinghits[musecvertices_count][j]  = tr->numberOfLostHits();
							musecvertices_track_npixelhits[musecvertices_count][j]    = tr->hitPattern().numberOfValidPixelHits();
							musecvertices_track_npixellayers[musecvertices_count][j]  = tr->hitPattern().pixelLayersWithMeasurement();
							musecvertices_track_nstriplayers[musecvertices_count][j]  = tr->hitPattern().stripLayersWithMeasurement();
						}
						musecvertices_count++;
						if(musecvertices_count == M_musecverticesmaxcount)
						{
							cerr << "number of secvertices > M_musecverticesmaxcount. They are missing." << endl; 
							break;
						}
					}
				}
			}// end of 2nd loop over muons
		}// end of 1st loop over muons
	}

	if(crecsecvertices && count > 0)
	{
		AddVertices(iEvent);
	}
	return(true);
}
